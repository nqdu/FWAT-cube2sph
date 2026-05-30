#!/usr/bin/env python3
"""
Job monitor for the FWAT shell pipeline.

Responsibilities
----------------
1. Track every job the pipeline submits (SEM arrays, POST, WOLFE) in a JSON
   state file so we survive restarts of the monitor itself.
2. Poll the scheduler (Slurm or PBS) and detect failed jobs / failed
   array-task subsets.
3. Resubmit failed work:
       * For arrays: resubmit only the failed task IDs.
       * For single jobs: resubmit the whole job.
4. Rebuild dependency chains: when a parent SEM is resubmitted, any pending
   downstream POST/WOLFE that depended on the old job ID is cancelled and
   resubmitted with a dependency on the new ID(s).
5. Enforce a hard retry cap (MAX_RETRIES) per logical job.

Designed to be called from `run.sh` in two modes:
    * Registration mode (one-shot):   monitor.py --cmd register ...
    * Wait/monitor loop:              monitor.py --cmd wait

The script intentionally has no third-party dependencies: stdlib only.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional


# -------- Configuration (override via env vars) --------
STATE_FILE   = Path(os.environ.get("MONITOR_STATE",   "LOG/monitor_state.json"))
MAX_RETRIES  = int(os.environ.get("MONITOR_MAX_RETRIES", "3"))
POLL_SECS    = int(os.environ.get("MONITOR_POLL_SECS",   "60"))
PLATFORM     = os.environ.get("PLATFORM", "slurm").lower()   # slurm | pbs
LOG_DIR      = Path("LOG")


# -------- Logging --------
def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[monitor {ts}] {msg}"
    print(line, flush=True)
    with open(LOG_DIR / "monitor.log", "a") as f:
        f.write(line + "\n")


# ==========================================================================
#  Scheduler backend abstraction
# ==========================================================================
class Backend:
    """Abstract scheduler interface."""

    def submit(self, script: str, args: list[str],
               dependency: Optional[list[str]] = None,
               array_spec: Optional[str] = None) -> str:
        """Submit a job and return its ID as a string."""
        raise NotImplementedError

    def cancel(self, jobid: str) -> None:
        raise NotImplementedError

    def state(self, jobid: str) -> dict:
        """
        Return:
          {
            'status'      : 'PENDING'|'RUNNING'|'COMPLETED'|'FAILED'|'CANCELLED'|'UNKNOWN',
            'failed_tasks': [int, ...]      # only for arrays; empty otherwise
          }
        """
        raise NotImplementedError


# -------- Slurm --------
class SlurmBackend(Backend):
    # States that mean the job (or array task) finished badly
    BAD = {"FAILED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY",
           "CANCELLED", "BOOT_FAIL", "PREEMPTED", "DEADLINE"}
    OK  = {"COMPLETED"}
    ACTIVE = {"PENDING", "RUNNING", "REQUEUED", "RESIZING",
              "SUSPENDED", "CONFIGURING"}

    def submit(self, script, args, dependency=None, array_spec=None):
        cmd = ["sbatch", "--parsable"]
        if dependency:
            cmd.append(f"--dependency=afterok:{':'.join(dependency)}")
        if array_spec:
            # override the #SBATCH --array line baked into the script
            cmd.append(f"--array={array_spec}")
        cmd.append(script)
        cmd.extend(args)
        out = subprocess.check_output(cmd, text=True).strip()
        # --parsable returns "JOBID" or "JOBID;CLUSTER"
        return out.split(";")[0]

    def cancel(self, jobid):
        subprocess.run(["scancel", jobid], check=False)

    def state(self, jobid):
        # -P pipe-delimited, -n no header, include array task rows (no -X)
        try:
            out = subprocess.check_output(
                ["sacct", "-j", jobid, "-n", "-P",
                 "-o", "JobID,State,ExitCode"],
                text=True, stderr=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError:
            return {"status": "UNKNOWN", "failed_tasks": []}

        rows = [r for r in out.splitlines() if r.strip()]
        if not rows:
            # Try squeue: job may still be queued and not yet in sacct
            sq = subprocess.run(["squeue", "-h", "-j", jobid],
                                capture_output=True, text=True)
            if sq.returncode == 0 and sq.stdout.strip():
                return {"status": "PENDING", "failed_tasks": []}
            return {"status": "UNKNOWN", "failed_tasks": []}

        # Look at array-task rows only (e.g. "12345_3"); skip ".batch" / ".extern"
        # If no task rows exist (non-array job), use the master row.
        task_rows, master_state = [], None
        for r in rows:
            jid, raw_state, _ = (r.split("|") + ["", ""])[:3]
            user_cancelled = raw_state.startswith("CANCELLED by ")
            st = raw_state.split()[0]
            if "." in jid:                              # .batch / .extern
                continue
            if "_" in jid:
                m = re.match(rf"{re.escape(jobid)}_(\d+)$", jid)
                if m:
                    task_rows.append((int(m.group(1)), st, user_cancelled))
            else:
                master_state = (st, user_cancelled)

        if task_rows:
            cancelled = sorted({i for i, _, cancelled in task_rows if cancelled})
            failed = sorted({i for i, s, _ in task_rows if s in self.BAD})
            active = [s for _, s, _ in task_rows if s in self.ACTIVE]
            if active:
                return {"status": "RUNNING", "failed_tasks": []}
            if cancelled:
                return {"status": "CANCELLED", "failed_tasks": cancelled}
            if failed:
                return {"status": "FAILED", "failed_tasks": failed}
            if all(s in self.OK for _, s, _ in task_rows):
                return {"status": "COMPLETED", "failed_tasks": []}
            return {"status": "UNKNOWN", "failed_tasks": failed}

        # Single (non-array) job
        if master_state is None:
            return {"status": "UNKNOWN", "failed_tasks": []}

        state_name, user_cancelled = master_state
        if state_name in self.ACTIVE:
            return {"status": "RUNNING", "failed_tasks": []}
        if user_cancelled:
            return {"status": "CANCELLED", "failed_tasks": []}
        if state_name in self.BAD:
            return {"status": "FAILED", "failed_tasks": []}
        if state_name in self.OK:
            return {"status": "COMPLETED", "failed_tasks": []}
        return {"status": "UNKNOWN", "failed_tasks": []}


# -------- PBS --------
def normalize_pbs_array_spec(array_spec: str) -> str:
    spec = array_spec.strip()

    single_index = re.fullmatch(r"(\d+)", spec)
    if single_index:
        start = int(single_index.group(1))
        return f"{start}-{start + 1}:2"

    single_range = re.fullmatch(r"(\d+)-(\d+)(?::(\d+))?", spec)
    if single_range:
        start = int(single_range.group(1))
        end = int(single_range.group(2))
        if start == end:
            return f"{start}-{start + 1}:2"

    return spec


class PBSBackend(Backend):
    def submit(self, script, args, dependency=None, array_spec=None):
        cmd = ["qsub"]
        if dependency:
            cmd += ["-W", f"depend=afterok:{':'.join(dependency)}"]
        if array_spec:
            # PBS rejects singleton ranges like 1-1, so normalize them.
            cmd += ["-J", normalize_pbs_array_spec(array_spec)]
        # Positional args have to come via -F on PBS
        if args:
            cmd += ["-F", " ".join(args)]
        cmd.append(script)
        out = subprocess.check_output(cmd, text=True).strip()
        # qsub returns e.g. "12345.server" or "12345[].server"
        return out.split(".")[0].rstrip("[]")

    def cancel(self, jobid):
        subprocess.run(["qdel", jobid], check=False)

    @staticmethod
    def was_user_cancelled(info: dict) -> bool:
        text_fields = [
            info.get("comment"),
            info.get("Comment"),
            info.get("obit_comment"),
            info.get("Obit_comment"),
            info.get("Submit_arguments"),
        ]
        combined = " ".join(str(value).lower() for value in text_fields if value)
        return any(
            phrase in combined
            for phrase in (
                "deleted as requested",
                "job deleted",
                "requestor=user",
                "cancelled",
                "canceled",
            )
        )

    def state(self, jobid):
        # -x includes completed jobs; -f -F json for structured output if avail.
        try:
            out = subprocess.check_output(
                ["qstat", "-x", "-f", "-F", "json", jobid],
                text=True, stderr=subprocess.DEVNULL,
            )
            data = json.loads(out)
        except (subprocess.CalledProcessError, json.JSONDecodeError):
            return {"status": "UNKNOWN", "failed_tasks": []}

        jobs = data.get("Jobs", {})
        if not jobs:
            return {"status": "UNKNOWN", "failed_tasks": []}

        # An array master has key "12345[].server"; subjobs appear as
        # "12345[1].server", "12345[2].server", ...
        sub_pattern = re.compile(rf"^{re.escape(jobid)}\[(\d+)\]")
        master_pattern = re.compile(rf"^{re.escape(jobid)}\[\]")
        sub_states, master_state = {}, None

        for key, info in jobs.items():
            st = info.get("job_state", "?")          # Q,R,H,E,F,X
            exit_status = info.get("Exit_status")
            user_cancelled = self.was_user_cancelled(info)
            ms = sub_pattern.match(key.split(".")[0] + "." + key.split(".",1)[1]
                                   if "." in key else key)
            if ms:
                idx = int(ms.group(1))
                sub_states[idx] = (st, exit_status, user_cancelled)
            elif master_pattern.match(key) or key.startswith(jobid):
                master_state = (st, exit_status, user_cancelled)

        def classify(st, exit_status, user_cancelled):
            if st in ("Q", "H", "W"):
                return "PENDING"
            if st in ("R", "E", "B"):                # B = array begun
                return "RUNNING"
            if st in ("F", "X"):                      # finished / exited
                if user_cancelled:
                    return "CANCELLED"
                if exit_status in (0, "0", None):
                    return "COMPLETED" if exit_status == 0 else "UNKNOWN"
                return "FAILED"
            return "UNKNOWN"

        if sub_states:
            classified = {
                i: classify(s, x, cancelled)
                for i, (s, x, cancelled) in sub_states.items()
            }
            cancelled = sorted(i for i, c in classified.items() if c == "CANCELLED")
            failed = sorted(i for i, c in classified.items() if c == "FAILED")
            if any(c in ("PENDING", "RUNNING") for c in classified.values()):
                return {"status": "RUNNING", "failed_tasks": []}
            if cancelled:
                return {"status": "CANCELLED", "failed_tasks": cancelled}
            if failed:
                return {"status": "FAILED", "failed_tasks": failed}
            if all(c == "COMPLETED" for c in classified.values()):
                return {"status": "COMPLETED", "failed_tasks": []}
            return {"status": "UNKNOWN", "failed_tasks": failed}

        if master_state:
            st, x, cancelled = master_state
            return {"status": classify(st, x, cancelled), "failed_tasks": []}
        return {"status": "UNKNOWN", "failed_tasks": []}


def get_backend() -> Backend:
    if PLATFORM == "slurm":
        return SlurmBackend()
    if PLATFORM == "pbs":
        return PBSBackend()
    raise SystemExit(f"monitor: unsupported PLATFORM={PLATFORM!r}")


# ==========================================================================
#  State model
# ==========================================================================
@dataclass
class JobRecord:
    name: str                     # logical name, e.g. "sem_noise", "post", "wolfe"
    jobid: str                    # current scheduler job ID
    script: str                   # script path that was submitted
    args: list[str] = field(default_factory=list)
    array_spec: Optional[str] = None     # current array spec (None for non-array)
    full_array_spec: Optional[str] = None # original full spec for reference
    parents: list[str] = field(default_factory=list)   # names of parent jobs
    retries: int = 0
    status: str = "PENDING"
    failed_tasks: list[int] = field(default_factory=list)


@dataclass
class State:
    jobs: dict[str, JobRecord] = field(default_factory=dict)

    @classmethod
    def load(cls) -> "State":
        if STATE_FILE.exists():
            data = json.loads(STATE_FILE.read_text())
            return cls(jobs={k: JobRecord(**v) for k, v in data.get("jobs", {}).items()})
        return cls()

    def save(self) -> None:
        STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
        tmp = STATE_FILE.with_suffix(".tmp")
        tmp.write_text(json.dumps(
            {"jobs": {k: asdict(v) for k, v in self.jobs.items()}},
            indent=2,
        ))
        tmp.replace(STATE_FILE)

    def parent_ids(self, names: list[str]) -> list[str]:
        return [self.jobs[n].jobid for n in names if n in self.jobs]

    def children_of(self, name: str) -> list[JobRecord]:
        return [j for j in self.jobs.values() if name in j.parents]


# ==========================================================================
#  Core operations
# ==========================================================================
def register(args) -> None:
    """Record a job that the shell script already submitted."""
    st = State.load()
    rec = JobRecord(
        name=args.name,
        jobid=args.jobid,
        script=args.script,
        args=args.script_args or [],
        array_spec=args.array_spec,
        full_array_spec=args.array_spec,
        parents=args.parents or [],
    )
    st.jobs[args.name] = rec
    st.save()
    log(f"registered {args.name} jobid={args.jobid} "
        f"array={args.array_spec or '-'} parents={rec.parents}")


def resubmit_failed_array(backend: Backend, rec: JobRecord) -> None:
    """Resubmit only failed array indices."""
    if not rec.failed_tasks:
        return
    spec = ",".join(str(i) for i in rec.failed_tasks)
    log(f"resubmitting {rec.name}: failed tasks [{spec}] "
        f"(retry {rec.retries + 1}/{MAX_RETRIES})")
    new_id = backend.submit(rec.script, rec.args, array_spec=spec)
    rec.jobid       = new_id
    rec.array_spec  = spec
    rec.retries    += 1
    rec.status      = "PENDING"
    rec.failed_tasks = []
    log(f"  -> new jobid={new_id}")


def resubmit_single(backend: Backend, rec: JobRecord,
                    state: State) -> None:
    """Resubmit a non-array job (POST/WOLFE)."""
    parent_ids = state.parent_ids(rec.parents)
    log(f"resubmitting {rec.name} "
        f"(retry {rec.retries + 1}/{MAX_RETRIES}) "
        f"deps={parent_ids or '-'}")
    new_id = backend.submit(rec.script, rec.args,
                            dependency=parent_ids or None)
    rec.jobid    = new_id
    rec.retries += 1
    rec.status   = "PENDING"
    log(f"  -> new jobid={new_id}")


def rewire_children(backend: Backend, parent: JobRecord,
                    state: State,
                    allow_cancel: bool = True) -> None:
    """
    When `parent` was just resubmitted (new jobid), any child whose dependency
    referred to the old jobid is now broken. Cancel and resubmit each pending
    child with a dependency on the NEW parent jobid.
    """
    if not allow_cancel:
        log(f"  skipping child rewiring for {parent.name}: parent retry in progress")
        return

    for child in state.children_of(parent.name):
        if child.status in ("COMPLETED",):
            continue
        if child.retries >= MAX_RETRIES:
            log(f"  child {child.name} already at retry cap; not rewiring")
            continue
        log(f"  rewiring child {child.name} (cancel {child.jobid}, resubmit)")
        backend.cancel(child.jobid)
        resubmit_single(backend, child, state)


def poll_once(backend: Backend, state: State) -> bool:
    """
    Update statuses and trigger resubmissions.
    Returns True if all jobs reached a terminal good state.
    """
    all_done = True
    changed  = False

    # Snapshot the keys: we mutate state.jobs during iteration only by
    # replacing fields on existing records, not by adding/removing.
    for name in list(state.jobs.keys()):
        rec = state.jobs[name]
        if rec.status == "COMPLETED":
            continue

        info = backend.state(rec.jobid)
        new_status = info["status"]
        if new_status != rec.status:
            log(f"{name} ({rec.jobid}): {rec.status} -> {new_status}")
            rec.status = new_status
            changed = True

        if new_status == "FAILED":
            rec.failed_tasks = info["failed_tasks"]
            if rec.retries >= MAX_RETRIES:
                log(f"ABORT: {name} exceeded MAX_RETRIES={MAX_RETRIES}. "
                    f"jobid={rec.jobid}")
                state.save()
                sys.exit(2)

            if rec.failed_tasks:                         # array job
                resubmit_failed_array(backend, rec)
            else:                                        # single job
                resubmit_single(backend, rec, state)
            rewire_children(backend, rec, state, allow_cancel=False)
            all_done = False
            changed = True
        elif new_status == "CANCELLED":
            rec.failed_tasks = info["failed_tasks"]
            if rec.parents:
                parents = [state.jobs[p] for p in rec.parents if p in state.jobs]

                if any(p.status == "CANCELLED" for p in parents):
                    cancelled_parents = [p.name for p in parents if p.status == "CANCELLED"]
                    log(f"ABORT: {name} cancelled and parent(s) cancelled: "
                        f"{','.join(cancelled_parents)}")
                    state.save()
                    sys.exit(3)

                if any(p.status == "FAILED" and p.retries < MAX_RETRIES for p in parents):
                    log(f"{name} ({rec.jobid}) cancelled while parent retrying; "
                        "waiting for parent resubmission")
                    all_done = False
                    continue

                if parents and any(p.status != "COMPLETED" for p in parents):
                    if rec.retries >= MAX_RETRIES:
                        log(f"ABORT: {name} exceeded MAX_RETRIES={MAX_RETRIES}. "
                            f"jobid={rec.jobid}")
                        state.save()
                        sys.exit(2)
                    log(f"{name} ({rec.jobid}) cancelled while parent chain is "
                        f"still resolving; resubmitting")
                    resubmit_single(backend, rec, state)
                    all_done = False
                    changed = True
                    continue

            log(f"ABORT: {name} was cancelled. jobid={rec.jobid}")
            state.save()
            sys.exit(3)
        elif new_status != "COMPLETED":
            all_done = False

    if changed:
        state.save()
    return all_done


def wait_loop() -> None:
    backend = get_backend()
    log(f"wait loop start (platform={PLATFORM}, poll={POLL_SECS}s, "
        f"max_retries={MAX_RETRIES})")
    while True:
        state = State.load()
        if not state.jobs:
            log("no jobs registered; exiting")
            return
        if poll_once(backend, state):
            log("all jobs COMPLETED")
            return
        time.sleep(POLL_SECS)


def show() -> None:
    st = State.load()
    if not st.jobs:
        print("(no jobs registered)")
        return
    print(f"{'NAME':<20} {'JOBID':<14} {'STATUS':<10} "
          f"{'RETRIES':<8} {'PARENTS'}")
    for j in st.jobs.values():
        print(f"{j.name:<20} {j.jobid:<14} {j.status:<10} "
              f"{j.retries}/{MAX_RETRIES:<5} {','.join(j.parents) or '-'}")


def reset() -> None:
    if STATE_FILE.exists():
        STATE_FILE.unlink()
        log("state cleared")


# ==========================================================================
#  CLI
# ==========================================================================
def main(argv: Optional[list[str]] = None) -> None:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--cmd",
        required=True,
        choices=["register", "wait", "show", "reset"],
        help="monitor action to run",
    )
    p.add_argument("--name")
    p.add_argument("--jobid")
    p.add_argument("--script")
    p.add_argument("--array-spec", default=None,
                   help="e.g. '1-10' or '1,3,7'")
    p.add_argument("--parents", nargs="*", default=[],
                   help="logical names of parent jobs")
    p.add_argument("--script-args", nargs="*", default=[])

    args = p.parse_args(argv)

    # create log dir if missing
    LOG_DIR.mkdir(exist_ok=True)

    if args.cmd == "register":
        missing = [
            flag for flag, value in (
                ("--name", args.name),
                ("--jobid", args.jobid),
                ("--script", args.script),
            )
            if not value
        ]
        if missing:
            p.error(f"{' '.join(missing)} required when --cmd=register")

    if   args.cmd == "register": register(args)
    elif args.cmd == "wait":     wait_loop()
    elif args.cmd == "show":     show()
    elif args.cmd == "reset":    reset()


if __name__ == "__main__":
    main()