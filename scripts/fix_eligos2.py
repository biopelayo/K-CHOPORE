#!/usr/bin/env python3
"""Fix ELIGOS2 compatibility issues.

1. Wrap all pandas2ri.activate() calls in try/except for rpy2 >= 3.6
2. Fix BED merge losing strand column when processing BED12 input
"""
import re

# Fix 1: pandas2ri.activate() DeprecationWarning
for f in ['/home/eligos2/_rna_mod.py', '/home/eligos2/_eligos_func.py']:
    try:
        with open(f) as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        continue
    new_lines = []
    for line in lines:
        stripped = line.rstrip('\n')
        if 'pandas2ri.activate()' in stripped and 'try:' not in stripped and 'except' not in stripped:
            indent = len(stripped) - len(stripped.lstrip())
            spaces = ' ' * indent
            new_lines.append(spaces + 'try:\n')
            new_lines.append(spaces + '    pandas2ri.activate()\n')
            new_lines.append(spaces + 'except (DeprecationWarning, Exception):\n')
            new_lines.append(spaces + '    pass\n')
        else:
            new_lines.append(line)
    with open(f, 'w') as fh:
        fh.writelines(new_lines)
    print(f'Fixed pandas2ri in {f}')

# Fix 2: BED merge strand column
misc_file = '/home/eligos2/_misc.py'
try:
    with open(misc_file) as fh:
        content = fh.read()
    content = content.replace(
        "mergedBed = beds.sort().merge(s=True,c='4',o='distinct')",
        "mergedBed = beds.sort().merge(s=True,c='6,4',o='distinct,distinct')"
    )
    with open(misc_file, 'w') as fh:
        fh.write(content)
    print(f'Fixed BED merge in {misc_file}')
except FileNotFoundError:
    pass
