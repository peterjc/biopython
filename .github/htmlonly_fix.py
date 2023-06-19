#!/usr/bin/env python
r"""Python script to discard LaTeX htmlonly blocks, keeping latexonly.

Input arguments are filenames (in situ), or - for stdin to stdout.

Will discard these blocks, used by hevea for HTML output::

    \begin{htmlonly}
    ...
    \end{htmlonly}

Will retain the contents of these blocks, used to tell hevea to ignore
bits:

    \begin{latexonly}
    ...
    \end{htmlonly}

These are typically used in the Biopython Tutorial to have two versions
of an image - for the PDF output figures may not be placed as is in the
text but delayed to a later page, requiring a label and reference, while
for the HTML output we can assume the image appears in line and there is
no need for a label and reference.

NOTE: Assumes these LaTeX commands are alone on their own line
NOTE: Assumes these LaTeX environments are not nested
NOTE: Does not remove any comments about the syntax removed
"""

import sys


def fix_lines(lines):
    """Loop over lines applying simple state tracking."""
    html = False
    for line in lines:
        if line.strip() == r"\begin{htmlonly}":
            html = True
            continue
        elif line.strip() == r"\end{htmlonly}":
            html = False
            continue
        elif line.strip() in (r"\begin{latexonly}", r"\end{latexonly}"):
            pass
        elif html:
            pass
        else:
            yield line


def fix_file(filename):
    """Edit a TEX file in-situ."""
    with open(filename) as handle:
        lines = list(handle)
    with open(filename, "w") as handle:
        for line in fix_lines(lines):
            handle.write(line)


for f in sys.argv[1:]:
    if f == "-":
        sys.stderr.write("Removing hevea syntax on stdin\n")
        for line in fix_lines(sys.stdin):
            sys.stdout.write(line)
    else:
        sys.stderr.write("Fixing %s\n" % f)
        fix_file(f)
