"""
Comprehensive test runner for rRNA-Phylo package.

Runs all tests and provides a summary of results.
"""

import subprocess
import sys
from pathlib import Path

def run_test(test_file):
    """Run a single test file and return results."""
    print(f"\n{'='*80}")
    print(f"Running: {test_file.name}")
    print('='*80)

    try:
        result = subprocess.run(
            [sys.executable, str(test_file)],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=test_file.parent.parent
        )

        # Check if test passed (exit code 0)
        if result.returncode == 0:
            print(f"PASSED: {test_file.name}")
            return True, None
        else:
            print(f"FAILED: {test_file.name}")
            print("\nError output:")
            print(result.stderr[-500:] if len(result.stderr) > 500 else result.stderr)
            return False, result.stderr
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT: {test_file.name} (exceeded 120 seconds)")
        return False, "Timeout exceeded"
    except Exception as e:
        print(f"ERROR: {test_file.name} - {str(e)}")
        return False, str(e)

def main():
    """Run all tests and generate summary."""
    print("="*80)
    print("rRNA-Phylo Test Suite")
    print("="*80)

    # Find all test files
    tests_dir = Path(__file__).parent / "tests"
    test_files = sorted(tests_dir.glob("test_*.py"))

    if not test_files:
        print("ERROR: No test files found!")
        return 1

    print(f"\nFound {len(test_files)} test files:")
    for test_file in test_files:
        print(f"  - {test_file.name}")

    # Run all tests
    results = {}
    for test_file in test_files:
        passed, error = run_test(test_file)
        results[test_file.name] = (passed, error)

    # Print summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)

    passed_tests = []
    failed_tests = []

    for test_name, (passed, error) in results.items():
        if passed:
            passed_tests.append(test_name)
            print(f"PASS: {test_name}")
        else:
            failed_tests.append(test_name)
            print(f"FAIL: {test_name}")

    print("\n" + "="*80)
    print(f"Total: {len(test_files)} tests")
    print(f"Passed: {len(passed_tests)} tests")
    print(f"Failed: {len(failed_tests)} tests")
    print("="*80)

    if failed_tests:
        print("\nFailed tests:")
        for test_name in failed_tests:
            print(f"  - {test_name}")

    # Return exit code
    return 0 if len(failed_tests) == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
