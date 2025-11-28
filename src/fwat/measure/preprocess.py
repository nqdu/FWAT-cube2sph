"""Preprocessing module for different measurement types.

This module provides a unified interface for preprocessing various seismic
measurement types including teleseismic, SKS, noise, and receiver functions.
"""

import sys
from typing import List, Protocol
from mpi4py import MPI


class PreprocessOperator(Protocol):
    """Protocol defining the interface for preprocessing operators."""
    
    def execute(self) -> None:
        """Execute the preprocessing operation."""
        ...


class PreprocessError(Exception):
    """Exception raised for preprocessing errors."""
    pass


def create_operator(
    mtype: str, 
    iter0: int, 
    evtid: str, 
    run_opt: int
) -> PreprocessOperator:
    """
    Factory function to create appropriate preprocessing operator.
    
    Args:
        mtype: Measurement type ('tele', 'sks', 'noise', 'rf')
        iter0: Iteration number
        evtid: Event ID
        run_opt: Run option (1=forward, 2=line search, 3=adjoint)
        
    Returns:
        Preprocessing operator instance
        
    Raises:
        PreprocessError: If measurement type is not supported
    """
    operators = {
        'tele': lambda: _create_tele_operator(mtype, iter0, evtid, run_opt),
        'sks': lambda: _create_sks_operator(mtype, iter0, evtid, run_opt),
        'noise': lambda: _create_noise_operator(mtype, iter0, evtid, run_opt),
        'rf': lambda: _create_rf_operator(mtype, iter0, evtid, run_opt),
    }
    
    if mtype not in operators:
        supported = ', '.join(operators.keys())
        raise PreprocessError(
            f"Unsupported measurement type '{mtype}'. "
            f"Supported types: {supported}"
        )
    
    return operators[mtype]()


def _create_tele_operator(mtype: str, iter0: int, evtid: str, run_opt: int):
    """Create teleseismic preprocessing operator."""
    from .tele_preproc import Tele_PreOP
    return Tele_PreOP(mtype, iter0, evtid, run_opt)


def _create_sks_operator(mtype: str, iter0: int, evtid: str, run_opt: int):
    """Create SKS preprocessing operator."""
    from .sks_preproc import SKS_PreOP
    return SKS_PreOP(mtype, iter0, evtid, run_opt)


def _create_noise_operator(mtype: str, iter0: int, evtid: str, run_opt: int):
    """Create noise preprocessing operator."""
    from .noise_mc_preproc import NoiseMC_PreOP
    return NoiseMC_PreOP(mtype, iter0, evtid, run_opt)


def _create_rf_operator(mtype: str, iter0: int, evtid: str, run_opt: int):
    """Create receiver function preprocessing operator."""
    from .rf_preproc import RF_PreOP
    return RF_PreOP(mtype, iter0, evtid, run_opt)


def run(argv: List[str]) -> None:
    """
    Run preprocessing with given arguments.
    
    Args:
        argv: Command line arguments [mtype, iter, evtid, run_opt]
        
    Raises:
        SystemExit: If arguments are invalid or execution fails
    """
    # Validate arguments
    if len(argv) != 4:
        print("Usage: fwat measure <measure_type> <iter> <evtid> <run_opt>")
        print()
        print("Arguments:")
        print("  measure_type: Type of measurement (tele, sks, noise, rf)")
        print("  iter:         Iteration number (integer)")
        print("  evtid:        Event ID (string)")
        print("  run_opt:      Run option (1=forward, 2=line_search, 3=adjoint)")
        print()
        print("Example:")
        print("  mpirun -np 31 fwat measure tele 0 XZ.FAF 3")
        sys.exit(1)
    
    # Parse arguments
    mtype = argv[0]
    
    try:
        iter0 = int(argv[1])
    except ValueError:
        print(f"Error: iter must be an integer, got '{argv[1]}'")
        sys.exit(1)
    
    evtid = argv[2]
    
    try:
        run_opt = int(argv[3])
    except ValueError:
        print(f"Error: run_opt must be an integer, got '{argv[3]}'")
        sys.exit(1)
    
    # Validate run_opt
    if run_opt not in {1, 2, 3}:
        print(f"Error: run_opt must be 1, 2, or 3, got {run_opt}")
        sys.exit(1)
    
    # Create and execute operator
    try:
        operator = create_operator(mtype, iter0, evtid, run_opt)
        operator.execute()
    except PreprocessError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error during preprocessing: {e}")
        raise
    finally:
        # Ensure MPI is finalized
        if MPI.Is_initialized() and not MPI.Is_finalized():
            MPI.Finalize()


def main() -> None:
    """Main entry point for the preprocessing module."""
    if len(sys.argv) != 5:
        print("Usage: python run_preprocess.py <measure_type> <iter> <evtid> <run_opt>")
        sys.exit(1)
    
    run(sys.argv[1:])


if __name__ == "__main__":
    main()