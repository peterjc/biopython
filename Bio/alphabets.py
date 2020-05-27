# Copyright 2020 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Provide a labelling system for biological sequence alphabets/molecule-tyes.

The alphabet system has been greatly simplified in Biopython 1.78, and is now
just a Python enumeration offering only:

 - DNA (including both single and double stranded)
 - RNA
 - Nucleotide (for use when DNA versus RNA is unknown)
 - Protein (including any modified residues)
 - Other/Unknown (default)

The two primary uses of this are to track the intended meaning of sequences
from file formats where the sequence type is explicit (e.g. GenBank and EMBL,
but not FASTA or FASTQ), and to prevent improper use of sequence methods (e.g.
you cannot transcribe a protein sequence).

This is essential a restricted molecule type, making no attempt to record or
restrict the allowed letters in the associated sequence (e.g. upper/lower or
mixed case, special characters like alignment gaps or stop codons). The now
deprecated ``Bio.Alphabet`` system could record but did not enforce this.
"""

import enum

class Alphabets(enum.Enum):
    DNA = enum.auto()
    RNA = enum.auto()
    Protein = enum.auto()
    Nucleotide = enum.auto()
    Other = enum.auto()

    def __repr__(self):
        """String representation for debugging."""
        return f"{self.__class__.__name__}.{self.name}"

    def __or__(self, other):
        """Implement the alphabet consensus rules.

        >>> Alphabets.DNA | Alphabets.RNA
        Alphabets.Nucleotide
        >>> Alphabets.Nucleotide | Alphabets.Protein
        Alphabets.Other

        NOTE: Follows the rules from Bio.Alphabet._consensus_base_alphabet
        """
        if not isinstance(other, self.__class__):
            return NotImplemented
        if self is other:
            return self
        if self is Alphabets.Other or other is Alphabets.Other:
            # Other | Any = Any | Other = Other
            return Alphabets.Other
        elif self is Alphabets.Protein or other is Alphabets.Protein:
            # Any-nucl | Protein = Protein | any-nucl = Other
            return Alphabets.Other
        else:
            # Any-nucl | Any-nucl = Nucleotide
            # e.g. DNA | RNA = Nucleotide
            return Alphabets.Nucleotide

    def __add__(self, other):
        """Implement the alphabet based sequence addition rules.

        This will raise TypeError on DNA+RNA, or any nucleotide+protein,
        or other+anything-else.

        >>> Alphabets.DNA + Alphabets.Nucleotide
        Alphabets.Nucleotide
        >>> Alphabets.Nucleotide + Alphabets.RNA
        Alphabets.Nucleotide
        >>> Alphabets.DNA + Alphabets.Other
        Alphabets.Other
        >>> Alphabets.RNA + Alphabets.DNA
        Traceback (most recent call last):
        ...
        TypeError: Incompatible alphabets RNA and DNA

        NOTE: Follows the rules from Bio.Alphabet._check_type_compatible
        """
        if not isinstance(other, self.__class__):
            return NotImplemented
        try:
            other = other._as_enum
            # Should only work on legacy Alphabet objects, issue warning?
        except AttributeError:
            pass
        if self is other:
            return self
        elif self is Alphabets.RNA and other in (Alphabets.DNA, Alphabets.Protein):
            raise TypeError(f"Incompatible alphabets {self.name} and {other.name}")
        elif self is Alphabets.DNA and other in (Alphabets.RNA, Alphabets.Protein):
            raise TypeError(f"Incompatible alphabets {self.name} and {other.name}")
        elif self is Alphabets.Protein and other is not Alphabets.Other:
            # i.e. other is RNA/DNA/Nucleotide
            raise TypeError(f"Incompatible alphabets {self.name} and {other.name}")
        elif self is Alphabets.Nucleotide and other is Alphabets.Protein:
            raise TypeError(
                f"Incompatible alphabets {self.name} and {other.name}"
            )
        return self | other
