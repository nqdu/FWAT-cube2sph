"""
Auto-generate the base <-> user model/kernel conversions for the ``dtti`` model
type used by ``FwatModel`` (``_cijkl2dtti`` and ``_cijkl_kl2dtti``).

Everything is driven by a single symbolic stiffness matrix ``C_ij(params)`` per
kltype (see ``_case``), so the model conversion and the kernel conversion can
never drift apart:

* backward model  (user params -> cijkl) : evaluate ``C_ij(params)``
* forward  model  (cijkl -> user params) : the analytic inverse/projection
                                           (hand-specified in ``_forward_exprs``)
* kernels (base cijkl kernels -> optimization space) : chain rule
                                           ``K_p = sum_ij (dC_ij/dp) * md_kl[k]``

We differentiate with sympy, run common-subexpression elimination (CSE), and
write plain numpy code.

Usage
-----
    python auto_kergen.py --write [path]   # (re)generate the module below
    python auto_kergen.py dtti <kltype>    # print one kltype kernel body (inspect)
"""
import os
import re
import sys

import sympy as sp

# module written by --write and imported by model.py
GENERATED_NAME = "dtti_generated.py"


def _index(m, n):
    """Upper-triangular (m<=n) index of c_{m,n} in the length-22 cijkl vector."""
    return m * 6 + n - (m * (m + 1)) // 2


def _c66mat(A, C, L, N, F, gc, gs):
    C0 = sp.Array([
        [A, A - 2 * N, F, 0, 0, 0],
        [A - 2 * N, A, F, 0, 0, 0],
        [F, F, C, 0, 0, 0],
        [0, 0, 0, L - gc, -gs, 0],
        [0, 0, 0, -gs, L + gc, 0],
        [0, 0, 0, 0, 0, N],
    ])

    return C0


def _case(kltype):
    """
    Return (params, C0, rho, n_log) for a given kltype.

    params : user-defined parameters in output order; the same order as the
             ``md_usr`` layout, so md_usr[i] corresponds to params[i].
    n_log  : number of leading parameters carried in log space (velocities and
             rho); their kernels are multiplied by md_usr to become d/dln(m).
    """
    if kltype == 1:  # vp, vs, rho, gcp, gsp
        vp, vs, rho, gcp, gsp = sp.symbols("vp vs rho gcp gsp")
        A = rho * vp**2
        C = A
        N = rho * vs**2
        L = N
        gc = gcp * L
        gs = gsp * L
        F = A - 2 * L                       # eta = 1
        params, n_log = [vp, vs, rho, gcp, gsp], 3

    elif kltype == 2:  # vph, vpv, vsh, vsv, rho, eta, gcp, gsp
        vph, vpv, vsh, vsv, rho, eta, gcp, gsp = sp.symbols("vph vpv vsh vsv rho eta gcp gsp")
        A = rho * vph**2
        L = rho * vsv**2
        C = rho * vpv**2
        N = rho * vsh**2
        gc = gcp * L
        gs = gsp * L
        F = eta * (A - 2 * L)
        params, n_log = [vph, vpv, vsh, vsv, rho, eta, gcp, gsp], 5

    elif kltype == 3:  # vp, vs, rho, (vph-vpv)/vpv, (vsh-vsv)/vsv, eta, gcp, gsp
        vp, vs, rho, kappaa, kappab, eta, gcp, gsp = sp.symbols("vp vs rho kappaa kappab eta gcp gsp")
        # vp = sqrt((2*vph**2 + vpv**2)/3),  kappaa = (vph - vpv)/vpv
        vpv_sq = 3 * vp**2 / (2 * (1 + kappaa)**2 + 1)
        vph_sq = (1 + kappaa)**2 * vpv_sq
        vsv_sq = 3 * vs**2 / (2 * (1 + kappab)**2 + 1)
        vsh_sq = (1 + kappab)**2 * vsv_sq
        A = rho * vph_sq
        L = rho * vsv_sq
        C = rho * vpv_sq
        N = rho * vsh_sq
        gc = gcp * L
        gs = gsp * L
        F = eta * (A - 2 * L)
        params, n_log = [vp, vs, rho, kappaa, kappab, eta, gcp, gsp], 3

    else:
        raise NotImplementedError(f"kltype = {kltype} is not implemented")

    return params, _c66mat(A, C, L, N, F, gc, gs), rho, n_log


def _forward_exprs(kltype):
    """
    Analytic inverse/projection: recover the user parameters from a base cijkl
    vector ``model`` (length 22, rho at index 21). Returns the output-order list
    of expressions. These encode a modelling choice (e.g. L = (c44+c55)/2) and
    are therefore hand-specified rather than derived.
    """
    model = sp.IndexedBase('model')
    rho = model[21]

    def c(i, j):
        return model[_index(i, j)]

    if kltype == 1:
        vp = sp.sqrt(c(0, 0) / rho)
        vs = sp.sqrt(c(5, 5) / rho)
        gc = sp.Rational(1, 2) * (c(4, 4) - c(3, 3))
        gs = -c(3, 4)
        out = [vp, vs, rho, gc / (rho * vs**2), gs / (rho * vs**2)]

    elif kltype == 2:
        A = c(0, 0)
        C = c(2, 2)
        L = sp.Rational(1, 2) * (c(3, 3) + c(4, 4))
        N = c(5, 5)
        F = c(0, 2)
        eta = F / (A - 2 * L)
        gcp = (c(4, 4) - L) / L
        gsp = -c(3, 4) / L
        out = [sp.sqrt(A / rho), sp.sqrt(C / rho), sp.sqrt(N / rho), sp.sqrt(L / rho),
               rho, eta, gcp, gsp]

    elif kltype == 3:
        A = c(0, 0)
        C = c(2, 2)
        L = sp.Rational(1, 2) * (c(3, 3) + c(4, 4))
        N = c(5, 5)
        F = c(0, 2)
        eta = F / (A - 2 * L)
        gcp = (c(4, 4) - L) / L
        gsp = -c(3, 4) / L
        vph = sp.sqrt(A / rho)
        vpv = sp.sqrt(C / rho)
        vsh = sp.sqrt(N / rho)
        vsv = sp.sqrt(L / rho)
        vp = sp.sqrt((2 * vph**2 + vpv**2) / 3)
        vs = sp.sqrt((2 * vsh**2 + vsv**2) / 3)
        kappaa = (vph - vpv) / vpv
        kappab = (vsh - vsv) / vsv
        out = [vp, vs, rho, kappaa, kappab, eta, gcp, gsp]

    else:
        raise NotImplementedError(f"kltype = {kltype} is not implemented")

    return out


_PRINTER = sp.printing.numpy.NumPyPrinter()


def _fmt(e):
    code = _PRINTER.doprint(e)
    code = code.replace('numpy.', 'np.')
    # md_kl[3] / model[5] -> md_kl[3,...] / model[5,...]  (broadcast over gridpoints)
    return re.sub(r'\b(md_kl|model)\[(-?\d+)\]', r'\1[\2,...]', code)


def _cse_block(exprs, pad):
    """CSE a list of expressions; return (temp_lines, reduced_exprs)."""
    replacements, reduced = sp.cse(exprs)
    lines = [f"{pad}{sym} = {_fmt(sub)}" for sym, sub in replacements]
    return lines, reduced


def build_kernel_body(kltype, indent=4):
    """numpy source lines for the base-kernel -> optimization-space conversion."""
    params, C0, rho, n_log = _case(kltype)
    md_kl = sp.IndexedBase('md_kl')

    # chain rule: one expression per output parameter
    exprs = []
    for p in params:
        expr = sp.diff(rho, p) * md_kl[21]   # rho kernel is the last entry
        k = 0
        for i in range(6):
            for j in range(i, 6):
                expr += sp.diff(C0[i, j], p) * md_kl[k]
                k += 1
        exprs.append(expr)

    pad = ' ' * indent
    lines = [f"{pad}{p.name} = md_usr[{i},...]" for i, p in enumerate(params)]
    lines.append(f"{pad}kl_opt = md_usr * 0")

    temps, reduced = _cse_block(exprs, pad)
    lines += temps
    for i, r in enumerate(reduced):
        lines.append(f"{pad}kl_opt[{i},...] = {_fmt(r)}")

    # chain rule to kernels w.r.t. log parameters (velocities and rho)
    lines.append(f"{pad}kl_opt[0:{n_log},...] *= md_usr[0:{n_log},...]")
    lines.append(f"{pad}return kl_opt")
    return lines


def build_model_body(kltype, indent=4):
    """numpy source lines for the cijkl <-> user model conversion (both directions)."""
    params, C0, rho, n_log = _case(kltype)
    n = len(params)
    pad = ' ' * indent
    pad2 = ' ' * (indent + 4)
    lines = []

    # ---- forward: cijkl -> user params ----
    lines.append(f"{pad}if not backward:")
    lines.append(f"{pad2}model_new = np.zeros(({n},) + model.shape[1:])")
    temps, reduced = _cse_block(_forward_exprs(kltype), pad2)
    lines += temps
    for i, r in enumerate(reduced):
        lines.append(f"{pad2}model_new[{i},...] = {_fmt(r)}")

    # ---- backward: user params -> cijkl (length 22) ----
    lines.append(f"{pad}else:")
    for i, p in enumerate(params):
        lines.append(f"{pad2}{p.name} = model[{i},...]")
    lines.append(f"{pad2}model_new = np.zeros((22,) + model.shape[1:])")
    # only the nonzero upper-triangle stiffnesses, plus rho at index 21
    entries = [(_index(i, j), C0[i, j])
               for i in range(6) for j in range(i, 6) if C0[i, j] != 0]
    entries.append((21, rho))
    temps, reduced = _cse_block([e for _, e in entries], pad2)
    lines += temps
    for (idx, _), r in zip(entries, reduced):
        lines.append(f"{pad2}model_new[{idx},...] = {_fmt(r)}")

    lines.append(f"{pad}return model_new")
    return lines


def write_module(path):
    """Write the full importable module with model + kernel converters."""
    out = [
        '"""',
        "AUTO-GENERATED by auto_kergen.py -- do not edit by hand.",
        "Regenerate with:  python auto_kergen.py --write",
        "",
        "Base <-> user model and kernel conversions for the dtti model type.",
        '"""',
        "import numpy as np",
        "",
        "",
    ]

    for kltype in (1, 2, 3):
        out.append(f"def cijkl2dtti_kltype{kltype}(model, backward=False):")
        out.extend(build_model_body(kltype))
        out += ["", ""]
        out.append(f"def kl2dtti_kltype{kltype}(md_usr, md_kl):")
        out.extend(build_kernel_body(kltype))
        out += ["", ""]

    out.append("MODEL_CONVERTERS = {")
    for kltype in (1, 2, 3):
        out.append(f"    {kltype}: cijkl2dtti_kltype{kltype},")
    out.append("}")
    out.append("")
    out.append("KERNEL_CONVERTERS = {")
    for kltype in (1, 2, 3):
        out.append(f"    {kltype}: kl2dtti_kltype{kltype},")
    out.append("}")
    out.append("")

    with open(path, "w") as f:
        f.write("\n".join(out))
    print(f"wrote {path}")


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "--write":
        path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
            os.path.dirname(os.path.abspath(__file__)), GENERATED_NAME)
        write_module(path)
        return

    if len(sys.argv) != 3:
        print("Usage: python auto_kergen.py --write [path]")
        print("       python auto_kergen.py dtti <kltype>")
        exit(1)

    mdtype = sys.argv[1]
    kltype = int(sys.argv[2])
    if (mdtype not in ['dtti']) or (kltype not in [1, 2, 3]):
        print(f"mdtype = {mdtype} and kltype = {kltype} is not implemented")
        return

    print("\n".join(build_kernel_body(kltype)))


if __name__ == "__main__":
    main()
