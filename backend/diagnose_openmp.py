#!/usr/bin/env python3
"""
Diagnose OpenMP library conflicts.

This script checks which packages are loading libiomp5md.dll and identifies
the source of the conflict.
"""

import os
import sys
from pathlib import Path


def find_openmp_libraries():
    """Find all OpenMP libraries in the Python environment."""
    print('=' * 70)
    print('DIAGNOSING OPENMP LIBRARY CONFLICT')
    print('=' * 70)
    print()

    # Check Python environment
    print(f'Python executable: {sys.executable}')
    print(f'Python version: {sys.version}')
    print()

    # Import packages that might load OpenMP
    print('Importing packages...')
    packages_to_check = []

    try:
        import numpy
        packages_to_check.append(('numpy', numpy.__version__, numpy.__file__))
        print(f'  NumPy {numpy.__version__}: OK')
    except ImportError as e:
        print(f'  NumPy: FAILED ({e})')

    try:
        import scipy
        packages_to_check.append(('scipy', scipy.__version__, scipy.__file__))
        print(f'  SciPy {scipy.__version__}: OK')
    except ImportError as e:
        print(f'  SciPy: FAILED ({e})')

    try:
        import numba
        packages_to_check.append(('numba', numba.__version__, numba.__file__))
        print(f'  Numba {numba.__version__}: OK')
    except ImportError as e:
        print(f'  Numba: FAILED ({e})')

    print()
    print('=' * 70)
    print('SEARCHING FOR OPENMP LIBRARIES')
    print('=' * 70)
    print()

    # Search for OpenMP libraries in each package
    openmp_libs = []

    for pkg_name, pkg_version, pkg_file in packages_to_check:
        pkg_dir = Path(pkg_file).parent

        # Search for libiomp5md.dll (Windows) or libomp.so (Linux)
        if sys.platform == 'win32':
            lib_patterns = ['**/libiomp5md.dll', '**/libomp*.dll']
        else:
            lib_patterns = ['**/libomp*.so', '**/libgomp*.so']

        for pattern in lib_patterns:
            for lib_file in pkg_dir.glob(pattern):
                openmp_libs.append((pkg_name, str(lib_file)))
                print(f'[FOUND] {pkg_name}:')
                print(f'        {lib_file}')
                print(f'        Size: {lib_file.stat().st_size:,} bytes')
                print()

    if not openmp_libs:
        print('[INFO] No OpenMP libraries found in package directories')
        print('       (Libraries may be loaded from system paths)')
        print()

    print('=' * 70)
    print('LOADED MODULES')
    print('=' * 70)
    print()

    # Check which modules are loaded
    import_order = []

    for module_name in sorted(sys.modules.keys()):
        if any(pkg in module_name for pkg in ['numpy', 'scipy', 'numba', 'mkl', 'omp']):
            import_order.append(module_name)

    print('Relevant modules loaded:')
    for i, module_name in enumerate(import_order[:20], 1):
        print(f'  {i:2d}. {module_name}')

    if len(import_order) > 20:
        print(f'  ... and {len(import_order) - 20} more')

    print()
    print('=' * 70)
    print('DIAGNOSIS')
    print('=' * 70)
    print()

    if len(openmp_libs) > 1:
        print('[ISSUE] Multiple OpenMP libraries found!')
        print()
        print('This happens when multiple packages (NumPy, SciPy, Numba) each')
        print('bundle their own copy of Intel MKL, which includes OpenMP.')
        print()
        print('Solutions:')
        print()
        print('1. Set environment variable (quick workaround):')
        print('   export KMP_DUPLICATE_LIB_OK=TRUE  # Linux/Mac')
        print('   set KMP_DUPLICATE_LIB_OK=TRUE     # Windows CMD')
        print('   $env:KMP_DUPLICATE_LIB_OK="TRUE"  # Windows PowerShell')
        print()
        print('2. Use conda-forge packages (recommended):')
        print('   conda install -c conda-forge numpy scipy numba')
        print('   (conda-forge ensures single MKL installation)')
        print()
        print('3. Use Intel Distribution (most robust):')
        print('   conda install -c intel numpy scipy')
        print('   pip install numba')
        print()
    else:
        print('[INFO] No obvious conflict detected in package directories.')
        print('       The conflict may be coming from system-wide libraries.')
        print()
        print('To prevent this warning, set:')
        print('   KMP_DUPLICATE_LIB_OK=TRUE')

    print()


if __name__ == '__main__':
    find_openmp_libraries()
