# Copyright 2012 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Biopython sequence objects (Seq) and related code (EXPERIMENTAL).

This is Bio.seq, a new module intended in the long term to replace
Bio.Seq (which provides the core object Seq, and related classes
like MutableSeq and UnknownSeq), and provide a namespace for other
related objects like the SeqRecord and SeqFeature and modules like
Bio.SeqIO.

For now Bio.seq and anything under it is considered EXPERIMENTAL.
This means you will see the BiopythonExperimentalWarning warning
when importing it, and the code is still expected to change before
being declared stable -- hopefully in one or two releases time.
"""
from Bio.Seq import * #Expose old module under new namespace

import warnings
from Bio import BiopythonExperimentalWarning
warnings.warn("Bio.seq is experimental code and likely to change! "
              "(Please read the release notes, NEWS and README files)",
              BiopythonExperimentalWarning)
