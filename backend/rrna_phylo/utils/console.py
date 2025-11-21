"""
Console output utilities with proper encoding handling.

Handles Unicode characters gracefully on Windows consoles.
"""

import sys
import io

# Safe print function that handles Unicode
def safe_print(*args, **kwargs):
    """
    Print function that safely handles Unicode characters on Windows.

    Replaces problematic Unicode characters with ASCII equivalents if needed.
    """
    try:
        print(*args, **kwargs)
    except UnicodeEncodeError:
        # Convert all args to strings and replace Unicode characters
        safe_args = []
        for arg in args:
            s = str(arg)
            # Replace common Unicode characters with ASCII equivalents
            s = s.replace('\u2713', '[OK]')  # [OK]
            s = s.replace('\u2717', '[X]')   # ✗
            s = s.replace('\u2192', '->')    # ->
            s = s.replace('\u2190', '<-')    # <-
            s = s.replace('\u2514', '`-')    # └
            s = s.replace('\u2500', '-')     # ─
            s = s.replace('\u251c', '|-')    # ├
            s = s.replace('\u2502', '|')     # │
            safe_args.append(s)
        print(*safe_args, **kwargs)

# Configure stdout to handle UTF-8 on Windows
def configure_console():
    """
    Configure console for UTF-8 output on Windows.

    Call this at the start of scripts to enable Unicode output.
    """
    if sys.platform == 'win32':
        try:
            # Try to reconfigure stdout to use UTF-8
            if hasattr(sys.stdout, 'reconfigure'):
                sys.stdout.reconfigure(encoding='utf-8')
            else:
                # Fallback for older Python versions
                sys.stdout = io.TextIOWrapper(
                    sys.stdout.buffer,
                    encoding='utf-8',
                    errors='replace'
                )
        except (AttributeError, io.UnsupportedOperation):
            # If reconfigure fails, we'll use safe_print instead
            pass

# ASCII-safe versions of tree drawing characters
TREE_CHARS_ASCII = {
    'branch': '|-',
    'last_branch': '`-',
    'vertical': '| ',
    'space': '  ',
    'horizontal': '-'
}

TREE_CHARS_UNICODE = {
    'branch': '├─',
    'last_branch': '└─',
    'vertical': '│ ',
    'space': '  ',
    'horizontal': '─'
}

def get_tree_chars():
    """Get tree drawing characters appropriate for the console."""
    # Try to detect if Unicode is supported
    try:
        test = '├─└│'
        test.encode(sys.stdout.encoding or 'ascii')
        return TREE_CHARS_UNICODE
    except (UnicodeEncodeError, AttributeError):
        return TREE_CHARS_ASCII
