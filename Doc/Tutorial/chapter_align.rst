.. _`chapter:align`:

Multiple Sequence Alignment objects
===================================

This chapter is about Multiple Sequence Alignments, by which we mean a
collection of multiple sequences which have been aligned together –
usually with the insertion of gap characters, and addition of leading or
trailing gaps – such that all the sequence strings are the same length.
Such an alignment can be regarded as a matrix of letters, where each row
is held as a ``SeqRecord`` object internally.

We will introduce the ``MultipleSeqAlignment`` object which holds this
kind of data, and the ``Bio.AlignIO`` module for reading and writing
them as various file formats (following the design of the ``Bio.SeqIO``
module from the previous chapter). Note that both ``Bio.SeqIO`` and
``Bio.AlignIO`` can read and write sequence alignment files. The
appropriate choice will depend largely on what you want to do with the
data.

The final part of this chapter is about our command line wrappers for
common multiple sequence alignment tools like ClustalW and MUSCLE.

Parsing or Reading Sequence Alignments
--------------------------------------

We have two functions for reading in sequence alignments,
``Bio.AlignIO.read()`` and ``Bio.AlignIO.parse()`` which following the
convention introduced in ``Bio.SeqIO`` are for files containing one or
multiple alignments respectively.

Using ``Bio.AlignIO.parse()`` will return an *iterator* which gives
``MultipleSeqAlignment`` objects. Iterators are typically used in a for
loop. Examples of situations where you will have multiple different
alignments include resampled alignments from the PHYLIP tool
``seqboot``, or multiple pairwise alignments from the EMBOSS tools
``water`` or ``needle``, or Bill Pearson’s FASTA tools.

However, in many situations you will be dealing with files which contain
only a single alignment. In this case, you should use the
``Bio.AlignIO.read()`` function which returns a single
``MultipleSeqAlignment`` object.

Both functions expect two mandatory arguments:

#. The first argument is a *handle* to read the data from, typically an
   open file (see
   Section :ref:`sec:appendix-handles`), or a
   filename.

#. The second argument is a lower case string specifying the alignment
   format. As in ``Bio.SeqIO`` we don’t try and guess the file format
   for you! See http://biopython.org/wiki/AlignIO for a full listing of
   supported formats.

There is also an optional ``seq_count`` argument which is discussed in
Section :ref:`sec:AlignIO-count-argument` below for dealing with
ambiguous file formats which may contain more than one alignment.

Single Alignments
~~~~~~~~~~~~~~~~~

As an example, consider the following annotation rich protein alignment
in the PFAM or Stockholm file format:

.. code:: text

   # STOCKHOLM 1.0
   #=GS COATB_BPIKE/30-81  AC P03620.1
   #=GS COATB_BPIKE/30-81  DR PDB; 1ifl ; 1-52;
   #=GS Q9T0Q8_BPIKE/1-52  AC Q9T0Q8.1
   #=GS COATB_BPI22/32-83  AC P15416.1
   #=GS COATB_BPM13/24-72  AC P69541.1
   #=GS COATB_BPM13/24-72  DR PDB; 2cpb ; 1-49;
   #=GS COATB_BPM13/24-72  DR PDB; 2cps ; 1-49;
   #=GS COATB_BPZJ2/1-49   AC P03618.1
   #=GS Q9T0Q9_BPFD/1-49   AC Q9T0Q9.1
   #=GS Q9T0Q9_BPFD/1-49   DR PDB; 1nh4 A; 1-49;
   #=GS COATB_BPIF1/22-73  AC P03619.2
   #=GS COATB_BPIF1/22-73  DR PDB; 1ifk ; 1-50;
   COATB_BPIKE/30-81             AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA
   #=GR COATB_BPIKE/30-81  SS    -HHHHHHHHHHHHHH--HHHHHHHH--HHHHHHHHHHHHHHHHHHHHH----
   Q9T0Q8_BPIKE/1-52             AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA
   COATB_BPI22/32-83             DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA
   COATB_BPM13/24-72             AEGDDP...AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
   #=GR COATB_BPM13/24-72  SS    ---S-T...CHCHHHHCCCCTCCCTTCHHHHHHHHHHHHHHHHHHHHCTT--
   COATB_BPZJ2/1-49              AEGDDP...AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA
   Q9T0Q9_BPFD/1-49              AEGDDP...AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
   #=GR Q9T0Q9_BPFD/1-49   SS    ------...-HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH--
   COATB_BPIF1/22-73             FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA
   #=GR COATB_BPIF1/22-73  SS    XX-HHHH--HHHHHH--HHHHHHH--HHHHHHHHHHHHHHHHHHHHHHH---
   #=GC SS_cons                  XHHHHHHHHHHHHHHHCHHHHHHHHCHHHHHHHHHHHHHHHHHHHHHHHC--
   #=GC seq_cons                 AEssss...AptAhDSLpspAT-hIu.sWshVsslVsAsluIKLFKKFsSKA
   //

This is the seed alignment for the Phage_Coat_Gp8 (PF05371) PFAM entry,
downloaded from a now out of date release of PFAM from
https://pfam.xfam.org/. We can load this file as follows (assuming it
has been saved to disk as “PF05371_seed.sth” in the current working
directory):

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")

This code will print out a summary of the alignment:

.. code:: pycon

   >>> print(alignment)
   Alignment with 7 rows and 52 columns
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73

You’ll notice in the above output the sequences have been truncated. We
could instead write our own code to format this as we please by
iterating over the rows as ``SeqRecord`` objects:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> print("Alignment length %i" % alignment.get_alignment_length())
   Alignment length 52
   >>> for record in alignment:
   ...     print("%s - %s" % (record.seq, record.id))
   ...
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73

You could also call Python’s built-in ``format`` function on the
alignment object to show it in a particular file format – see
Section :ref:`sec:alignment-format` for details.

Did you notice in the raw file above that several of the sequences
include database cross-references to the PDB and the associated known
secondary structure? Try this:

.. code:: pycon

   >>> for record in alignment:
   ...     if record.dbxrefs:
   ...         print("%s %s" % (record.id, record.dbxrefs))
   ...
   COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
   COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
   Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
   COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']

To have a look at all the sequence annotation, try this:

.. code:: pycon

   >>> for record in alignment:
   ...     print(record)
   ...

PFAM provide a nice web interface at http://pfam.xfam.org/family/PF05371
which will actually let you download this alignment in several other
formats. This is what the file looks like in the FASTA file format:

.. code:: text

   >COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA
   >Q9T0Q8_BPIKE/1-52
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA
   >COATB_BPI22/32-83
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA
   >COATB_BPM13/24-72
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
   >COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA
   >Q9T0Q9_BPFD/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
   >COATB_BPIF1/22-73
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA

Note the website should have an option about showing gaps as periods
(dots) or dashes, we’ve shown dashes above. Assuming you download and
save this as file “PF05371_seed.faa” then you can load it with almost
exactly the same code:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.faa", "fasta")
   >>> print(alignment)

All that has changed in this code is the filename and the format string.
You’ll get the same output as before, the sequences and record
identifiers are the same. However, as you should expect, if you check
each ``SeqRecord`` there is no annotation nor database cross-references
because these are not included in the FASTA file format.

Note that rather than using the Sanger website, you could have used
``Bio.AlignIO`` to convert the original Stockholm format file into a
FASTA file yourself (see below).

With any supported file format, you can load an alignment in exactly the
same way just by changing the format string. For example, use “phylip”
for PHYLIP files, “nexus” for NEXUS files or “emboss” for the alignments
output by the EMBOSS tools. There is a full listing on the wiki page
(http://biopython.org/wiki/AlignIO) and in the built in documentation
(also
`online <http://biopython.org/docs/\bpversion/api/Bio.AlignIO.html>`__):

.. code:: pycon

   >>> from Bio import AlignIO
   >>> help(AlignIO)

Multiple Alignments
~~~~~~~~~~~~~~~~~~~

The previous section focused on reading files containing a single
alignment. In general however, files can contain more than one
alignment, and to read these files we must use the
``Bio.AlignIO.parse()`` function.

Suppose you have a small alignment in PHYLIP format:

.. code:: text

       5    6
   Alpha     AACAAC
   Beta      AACCCC
   Gamma     ACCAAC
   Delta     CCACCA
   Epsilon   CCAAAC

If you wanted to bootstrap a phylogenetic tree using the PHYLIP tools,
one of the steps would be to create a set of many resampled alignments
using the tool ``bootseq``. This would give output something like this,
which has been abbreviated for conciseness:

.. code:: text

       5     6
   Alpha     AAACCA
   Beta      AAACCC
   Gamma     ACCCCA
   Delta     CCCAAC
   Epsilon   CCCAAA
       5     6
   Alpha     AAACAA
   Beta      AAACCC
   Gamma     ACCCAA
   Delta     CCCACC
   Epsilon   CCCAAA
       5     6
   Alpha     AAAAAC
   Beta      AAACCC
   Gamma     AACAAC
   Delta     CCCCCA
   Epsilon   CCCAAC
   ...
       5     6
   Alpha     AAAACC
   Beta      ACCCCC
   Gamma     AAAACC
   Delta     CCCCAA
   Epsilon   CAAACC

If you wanted to read this in using ``Bio.AlignIO`` you could use:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignments = AlignIO.parse("resampled.phy", "phylip")
   >>> for alignment in alignments:
   ...     print(alignment)
   ...     print()
   ...

This would give the following output, again abbreviated for display:

.. code:: text

   Alignment with 5 rows and 6 columns
   AAACCA Alpha
   AAACCC Beta
   ACCCCA Gamma
   CCCAAC Delta
   CCCAAA Epsilon

   Alignment with 5 rows and 6 columns
   AAACAA Alpha
   AAACCC Beta
   ACCCAA Gamma
   CCCACC Delta
   CCCAAA Epsilon

   Alignment with 5 rows and 6 columns
   AAAAAC Alpha
   AAACCC Beta
   AACAAC Gamma
   CCCCCA Delta
   CCCAAC Epsilon

   ...

   Alignment with 5 rows and 6 columns
   AAAACC Alpha
   ACCCCC Beta
   AAAACC Gamma
   CCCCAA Delta
   CAAACC Epsilon

As with the function ``Bio.SeqIO.parse()``, using
``Bio.AlignIO.parse()`` returns an iterator. If you want to keep all the
alignments in memory at once, which will allow you to access them in any
order, then turn the iterator into a list:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignments = list(AlignIO.parse("resampled.phy", "phylip"))
   >>> last_align = alignments[-1]
   >>> first_align = alignments[0]

.. _`sec:AlignIO-count-argument`:

Ambiguous Alignments
~~~~~~~~~~~~~~~~~~~~

Many alignment file formats can explicitly store more than one
alignment, and the division between each alignment is clear. However,
when a general sequence file format has been used there is no such block
structure. The most common such situation is when alignments have been
saved in the FASTA file format. For example consider the following:

.. code:: text

   >Alpha
   ACTACGACTAGCTCAG--G
   >Beta
   ACTACCGCTAGCTCAGAAG
   >Gamma
   ACTACGGCTAGCACAGAAG
   >Alpha
   ACTACGACTAGCTCAGG--
   >Beta
   ACTACCGCTAGCTCAGAAG
   >Gamma
   ACTACGGCTAGCACAGAAG

This could be a single alignment containing six sequences (with repeated
identifiers). Or, judging from the identifiers, this is probably two
different alignments each with three sequences, which happen to all have
the same length.

What about this next example?

.. code:: text

   >Alpha
   ACTACGACTAGCTCAG--G
   >Beta
   ACTACCGCTAGCTCAGAAG
   >Alpha
   ACTACGACTAGCTCAGG--
   >Gamma
   ACTACGGCTAGCACAGAAG
   >Alpha
   ACTACGACTAGCTCAGG--
   >Delta
   ACTACGGCTAGCACAGAAG

Again, this could be a single alignment with six sequences. However this
time based on the identifiers we might guess this is three pairwise
alignments which by chance have all got the same lengths.

This final example is similar:

.. code:: text

   >Alpha
   ACTACGACTAGCTCAG--G
   >XXX
   ACTACCGCTAGCTCAGAAG
   >Alpha
   ACTACGACTAGCTCAGG
   >YYY
   ACTACGGCAAGCACAGG
   >Alpha
   --ACTACGAC--TAGCTCAGG
   >ZZZ
   GGACTACGACAATAGCTCAGG

In this third example, because of the differing lengths, this cannot be
treated as a single alignment containing all six records. However, it
could be three pairwise alignments.

Clearly trying to store more than one alignment in a FASTA file is not
ideal. However, if you are forced to deal with these as input files
``Bio.AlignIO`` can cope with the most common situation where all the
alignments have the same number of records. One example of this is a
collection of pairwise alignments, which can be produced by the EMBOSS
tools ``needle`` and ``water`` – although in this situation,
``Bio.AlignIO`` should be able to understand their native output using
“emboss” as the format string.

To interpret these FASTA examples as several separate alignments, we can
use ``Bio.AlignIO.parse()`` with the optional ``seq_count`` argument
which specifies how many sequences are expected in each alignment (in
these examples, 3, 2 and 2 respectively). For example, using the third
example as the input data:

.. code:: pycon

   >>> for alignment in AlignIO.parse(handle, "fasta", seq_count=2):
   ...     print("Alignment length %i" % alignment.get_alignment_length())
   ...     for record in alignment:
   ...         print("%s - %s" % (record.seq, record.id))
   ...     print()
   ...

giving:

.. code:: text

   Alignment length 19
   ACTACGACTAGCTCAG--G - Alpha
   ACTACCGCTAGCTCAGAAG - XXX

   Alignment length 17
   ACTACGACTAGCTCAGG - Alpha
   ACTACGGCAAGCACAGG - YYY

   Alignment length 21
   --ACTACGAC--TAGCTCAGG - Alpha
   GGACTACGACAATAGCTCAGG - ZZZ

Using ``Bio.AlignIO.read()`` or ``Bio.AlignIO.parse()`` without the
``seq_count`` argument would give a single alignment containing all six
records for the first two examples. For the third example, an exception
would be raised because the lengths differ preventing them being turned
into a single alignment.

If the file format itself has a block structure allowing ``Bio.AlignIO``
to determine the number of sequences in each alignment directly, then
the ``seq_count`` argument is not needed. If it is supplied, and doesn’t
agree with the file contents, an error is raised.

Note that this optional ``seq_count`` argument assumes each alignment in
the file has the same number of sequences. Hypothetically you may come
across stranger situations, for example a FASTA file containing several
alignments each with a different number of sequences – although I would
love to hear of a real world example of this. Assuming you cannot get
the data in a nicer file format, there is no straight forward way to
deal with this using ``Bio.AlignIO``. In this case, you could consider
reading in the sequences themselves using ``Bio.SeqIO`` and batching
them together to create the alignments as appropriate.

Writing Alignments
------------------

We’ve talked about using ``Bio.AlignIO.read()`` and
``Bio.AlignIO.parse()`` for alignment input (reading files), and now
we’ll look at ``Bio.AlignIO.write()`` which is for alignment output
(writing files). This is a function taking three arguments: some
``MultipleSeqAlignment`` objects (or for backwards compatibility the
obsolete ``Alignment`` objects), a handle or filename to write to, and a
sequence format.

Here is an example, where we start by creating a few
``MultipleSeqAlignment`` objects the hard way (by hand, rather than by
loading them from a file). Note we create some ``SeqRecord`` objects to
construct the alignment from.

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord
   >>> from Bio.Align import MultipleSeqAlignment
   >>> align1 = MultipleSeqAlignment(
   ...     [
   ...         SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
   ...         SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
   ...         SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
   ...     ]
   ... )
   >>> align2 = MultipleSeqAlignment(
   ...     [
   ...         SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
   ...         SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
   ...         SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
   ...     ]
   ... )
   >>> align3 = MultipleSeqAlignment(
   ...     [
   ...         SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
   ...         SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
   ...         SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
   ...     ]
   ... )
   >>> my_alignments = [align1, align2, align3]

Now we have a list of ``Alignment`` objects, we’ll write them to a
PHYLIP format file:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> AlignIO.write(my_alignments, "my_example.phy", "phylip")

And if you open this file in your favorite text editor it should look
like this:

.. code:: text

    3 12
   Alpha      ACTGCTAGCT AG
   Beta       ACT-CTAGCT AG
   Gamma      ACTGCTAGDT AG
    3 9
   Delta      GTCAGC-AG
   Epislon    GACAGCTAG
   Zeta       GTCAGCTAG
    3 13
   Eta        ACTAGTACAG CTG
   Theta      ACTAGTACAG CT-
   Iota       -CTACTACAG GTG

Its more common to want to load an existing alignment, and save that,
perhaps after some simple manipulation like removing certain rows or
columns.

Suppose you wanted to know how many alignments the
``Bio.AlignIO.write()`` function wrote to the handle? If your alignments
were in a list like the example above, you could just use
``len(my_alignments)``, however you can’t do that when your records come
from a generator/iterator. Therefore the ``Bio.AlignIO.write()``
function returns the number of alignments written to the file.

*Note* - If you tell the ``Bio.AlignIO.write()`` function to write to a
file that already exists, the old file will be overwritten without any
warning.

.. _`sec:converting-alignments`:

Converting between sequence alignment file formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Converting between sequence alignment file formats with ``Bio.AlignIO``
works in the same way as converting between sequence file formats with
``Bio.SeqIO``
(Section :ref:`sec:SeqIO-conversion`). We load
generally the alignment(s) using ``Bio.AlignIO.parse()`` and then save
them using the ``Bio.AlignIO.write()`` – or just use the
``Bio.AlignIO.convert()`` helper function.

For this example, we’ll load the PFAM/Stockholm format file used earlier
and save it as a Clustal W format file:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
   >>> print("Converted %i alignments" % count)
   Converted 1 alignments

Or, using ``Bio.AlignIO.parse()`` and ``Bio.AlignIO.write()``:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
   >>> count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
   >>> print("Converted %i alignments" % count)
   Converted 1 alignments

The ``Bio.AlignIO.write()`` function expects to be given multiple
alignment objects. In the example above we gave it the alignment
iterator returned by ``Bio.AlignIO.parse()``.

In this case, we know there is only one alignment in the file so we
could have used ``Bio.AlignIO.read()`` instead, but notice we have to
pass this alignment to ``Bio.AlignIO.write()`` as a single element list:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> AlignIO.write([alignment], "PF05371_seed.aln", "clustal")

Either way, you should end up with the same new Clustal W format file
“PF05371_seed.aln” with the following content:

.. code:: text

   CLUSTAL X (1.81) multiple sequence alignment


   COATB_BPIKE/30-81                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSS
   Q9T0Q8_BPIKE/1-52                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVS
   COATB_BPI22/32-83                   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSS
   COATB_BPM13/24-72                   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS
   COATB_BPZJ2/1-49                    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFAS
   Q9T0Q9_BPFD/1-49                    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS
   COATB_BPIF1/22-73                   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVS

   COATB_BPIKE/30-81                   KA
   Q9T0Q8_BPIKE/1-52                   RA
   COATB_BPI22/32-83                   KA
   COATB_BPM13/24-72                   KA
   COATB_BPZJ2/1-49                    KA
   Q9T0Q9_BPFD/1-49                    KA
   COATB_BPIF1/22-73                   RA

Alternatively, you could make a PHYLIP format file which we’ll name
“PF05371_seed.phy”:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")

This time the output looks like this:

.. code:: text

    7 52
   COATB_BPIK AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIRLFKKFSS
   Q9T0Q8_BPI AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIKLFKKFVS
   COATB_BPI2 DGTSTATSYA TEAMNSLKTQ ATDLIDQTWP VVTSVAVAGL AIRLFKKFSS
   COATB_BPM1 AEGDDP---A KAAFNSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   COATB_BPZJ AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFAS
   Q9T0Q9_BPF AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   COATB_BPIF FAADDATSQA KAAFDSLTAQ ATEMSGYAWA LVVLVVGATV GIKLFKKFVS

              KA
              RA
              KA
              KA
              KA
              KA
              RA

One of the big handicaps of the original PHYLIP alignment file format is
that the sequence identifiers are strictly truncated at ten characters.
In this example, as you can see the resulting names are still unique -
but they are not very readable. As a result, a more relaxed variant of
the original PHYLIP format is now quite widely used:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")

This time the output looks like this, using a longer indentation to
allow all the identifiers to be given in full:

.. code:: text

    7 52
   COATB_BPIKE/30-81  AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIRLFKKFSS
   Q9T0Q8_BPIKE/1-52  AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIKLFKKFVS
   COATB_BPI22/32-83  DGTSTATSYA TEAMNSLKTQ ATDLIDQTWP VVTSVAVAGL AIRLFKKFSS
   COATB_BPM13/24-72  AEGDDP---A KAAFNSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   COATB_BPZJ2/1-49   AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFAS
   Q9T0Q9_BPFD/1-49   AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   COATB_BPIF1/22-73  FAADDATSQA KAAFDSLTAQ ATEMSGYAWA LVVLVVGATV GIKLFKKFVS

                      KA
                      RA
                      KA
                      KA
                      KA
                      KA
                      RA

If you have to work with the original strict PHYLIP format, then you may
need to compress the identifiers somehow – or assign your own names or
numbering system. This following bit of code manipulates the record
identifiers before saving the output:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> name_mapping = {}
   >>> for i, record in enumerate(alignment):
   ...     name_mapping[i] = record.id
   ...     record.id = "seq%i" % i
   ...
   >>> print(name_mapping)
   {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}
   >>> AlignIO.write([alignment], "PF05371_seed.phy", "phylip")

This code used a Python dictionary to record a simple mapping from the
new sequence system to the original identifier:

.. code:: python

   {
       0: "COATB_BPIKE/30-81",
       1: "Q9T0Q8_BPIKE/1-52",
       2: "COATB_BPI22/32-83",
       # ...
   }

Here is the new (strict) PHYLIP format output:

.. code:: text

    7 52
   seq0       AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIRLFKKFSS
   seq1       AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIKLFKKFVS
   seq2       DGTSTATSYA TEAMNSLKTQ ATDLIDQTWP VVTSVAVAGL AIRLFKKFSS
   seq3       AEGDDP---A KAAFNSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   seq4       AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFAS
   seq5       AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   seq6       FAADDATSQA KAAFDSLTAQ ATEMSGYAWA LVVLVVGATV GIKLFKKFVS

              KA
              RA
              KA
              KA
              KA
              KA
              RA

In general, because of the identifier limitation, working with *strict*
PHYLIP file formats shouldn’t be your first choice. Using the
PFAM/Stockholm format on the other hand allows you to record a lot of
additional annotation too.

.. _`sec:alignment-format`:

Getting your alignment objects as formatted strings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Bio.AlignIO`` interface is based on handles, which means if you
want to get your alignment(s) into a string in a particular file format
you need to do a little bit more work (see below). However, you will
probably prefer to call Python’s built-in ``format`` function on the
alignment object. This takes an output format specification as a single
argument, a lower case string which is supported by ``Bio.AlignIO`` as
an output format. For example:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> print(format(alignment, "clustal"))
   CLUSTAL X (1.81) multiple sequence alignment


   COATB_BPIKE/30-81                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSS
   Q9T0Q8_BPIKE/1-52                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVS
   COATB_BPI22/32-83                   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSS
   ...

Without an output format specification, ``format`` returns the same
output as ``str``.

As described in
Section :ref:`sec:SeqRecord-format`, the
``SeqRecord`` object has a similar method using output formats supported
by ``Bio.SeqIO``.

Internally ``format`` is calling ``Bio.AlignIO.write()`` with a
``StringIO`` handle. You can do this in your own code if for example you
are using an older version of Biopython:

.. code:: pycon

   >>> from io import StringIO
   >>> from Bio import AlignIO
   >>> alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
   >>> out_handle = StringIO()
   >>> AlignIO.write(alignments, out_handle, "clustal")
   1
   >>> clustal_data = out_handle.getvalue()
   >>> print(clustal_data)
   CLUSTAL X (1.81) multiple sequence alignment


   COATB_BPIKE/30-81                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSS
   Q9T0Q8_BPIKE/1-52                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVS
   COATB_BPI22/32-83                   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSS
   COATB_BPM13/24-72                   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS
   ...

.. _`sec:manipulating-alignments`:

Manipulating Alignments
-----------------------

Now that we’ve covered loading and saving alignments, we’ll look at what
else you can do with them.

Slicing alignments
~~~~~~~~~~~~~~~~~~

First of all, in some senses the alignment objects act like a Python
``list`` of ``SeqRecord`` objects (the rows). With this model in mind
hopefully the actions of ``len()`` (the number of rows) and iteration
(each row as a ``SeqRecord``) make sense:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> print("Number of rows: %i" % len(alignment))
   Number of rows: 7
   >>> for record in alignment:
   ...     print("%s - %s" % (record.seq, record.id))
   ...
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73

You can also use the list-like ``append`` and ``extend`` methods to add
more rows to the alignment (as ``SeqRecord`` objects). Keeping the list
metaphor in mind, simple slicing of the alignment should also make sense
- it selects some of the rows giving back another alignment object:

.. code:: pycon

   >>> print(alignment)
   Alignment with 7 rows and 52 columns
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73
   >>> print(alignment[3:7])
   Alignment with 4 rows and 52 columns
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73

What if you wanted to select by column? Those of you who have used the
NumPy matrix or array objects won’t be surprised at this - you use a
double index.

.. code:: pycon

   >>> print(alignment[2, 6])
   T

Using two integer indices pulls out a single letter, short hand for
this:

.. code:: pycon

   >>> print(alignment[2].seq[6])
   T

You can pull out a single column as a string like this:

.. code:: pycon

   >>> print(alignment[:, 6])
   TTT---T

You can also select a range of columns. For example, to pick out those
same three rows we extracted earlier, but take just their first six
columns:

.. code:: pycon

   >>> print(alignment[3:6, :6])
   Alignment with 3 rows and 6 columns
   AEGDDP COATB_BPM13/24-72
   AEGDDP COATB_BPZJ2/1-49
   AEGDDP Q9T0Q9_BPFD/1-49

Leaving the first index as ``:`` means take all the rows:

.. code:: pycon

   >>> print(alignment[:, :6])
   Alignment with 7 rows and 6 columns
   AEPNAA COATB_BPIKE/30-81
   AEPNAA Q9T0Q8_BPIKE/1-52
   DGTSTA COATB_BPI22/32-83
   AEGDDP COATB_BPM13/24-72
   AEGDDP COATB_BPZJ2/1-49
   AEGDDP Q9T0Q9_BPFD/1-49
   FAADDA COATB_BPIF1/22-73

This brings us to a neat way to remove a section. Notice columns 7, 8
and 9 which are gaps in three of the seven sequences:

.. code:: pycon

   >>> print(alignment[:, 6:9])
   Alignment with 7 rows and 3 columns
   TNY COATB_BPIKE/30-81
   TNY Q9T0Q8_BPIKE/1-52
   TSY COATB_BPI22/32-83
   --- COATB_BPM13/24-72
   --- COATB_BPZJ2/1-49
   --- Q9T0Q9_BPFD/1-49
   TSQ COATB_BPIF1/22-73

Again, you can slice to get everything after the ninth column:

.. code:: pycon

   >>> print(alignment[:, 9:])
   Alignment with 7 rows and 43 columns
   ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
   ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
   ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
   AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
   AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
   AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
   AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73

Now, the interesting thing is that addition of alignment objects works
by column. This lets you do this as a way to remove a block of columns:

.. code:: pycon

   >>> edited = alignment[:, :6] + alignment[:, 9:]
   >>> print(edited)
   Alignment with 7 rows and 49 columns
   AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
   AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
   DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
   AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
   AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
   AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
   FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73

Another common use of alignment addition would be to combine alignments
for several different genes into a meta-alignment. Watch out though -
the identifiers need to match up (see
Section :ref:`sec:SeqRecord-addition` for how
adding ``SeqRecord`` objects works). You may find it helpful to first
sort the alignment rows alphabetically by id:

.. code:: pycon

   >>> edited.sort()
   >>> print(edited)
   Alignment with 7 rows and 49 columns
   DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
   FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
   AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
   AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
   AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
   AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
   AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49

Note that you can only add two alignments together if they have the same
number of rows.

Alignments as arrays
~~~~~~~~~~~~~~~~~~~~

Depending on what you are doing, it can be more useful to turn the
alignment object into an array of letters – and you can do this with
NumPy:

.. code:: pycon

   >>> import numpy as np
   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> align_array = np.array(alignment)
   >>> print("Array shape %i by %i" % align_array.shape)
   Array shape 7 by 52
   >>> align_array[:, :10]  # doctest:+ELLIPSIS
   array([['A', 'E', 'P', 'N', 'A', 'A', 'T', 'N', 'Y', 'A'],
          ['A', 'E', 'P', 'N', 'A', 'A', 'T', 'N', 'Y', 'A'],
          ['D', 'G', 'T', 'S', 'T', 'A', 'T', 'S', 'Y', 'A'],
          ['A', 'E', 'G', 'D', 'D', 'P', '-', '-', '-', 'A'],
          ['A', 'E', 'G', 'D', 'D', 'P', '-', '-', '-', 'A'],
          ['A', 'E', 'G', 'D', 'D', 'P', '-', '-', '-', 'A'],
          ['F', 'A', 'A', 'D', 'D', 'A', 'T', 'S', 'Q', 'A']],...

Note that this leaves the original Biopython alignment object and the
NumPy array in memory as separate objects - editing one will not update
the other!

Getting information on the alignment
------------------------------------

Substitutions
~~~~~~~~~~~~~

The ``substitutions`` property of an alignment reports how often letters
in the alignment are substituted for each other. This is calculated by
taking all pairs of rows in the alignment, counting the number of times
two letters are aligned to each other, and summing this over all pairs.
For example,

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord
   >>> from Bio.Align import MultipleSeqAlignment
   >>> alignment = MultipleSeqAlignment(
   ...     [
   ...         SeqRecord(Seq("ACTCCTA"), id="seq1"),
   ...         SeqRecord(Seq("AAT-CTA"), id="seq2"),
   ...         SeqRecord(Seq("CCTACT-"), id="seq3"),
   ...         SeqRecord(Seq("TCTCCTC"), id="seq4"),
   ...     ]
   ... )
   >>> print(alignment)
   Alignment with 4 rows and 7 columns
   ACTCCTA seq1
   AAT-CTA seq2
   CCTACT- seq3
   TCTCCTC seq4
   >>> substitutions = alignment.substitutions
   >>> print(substitutions)
       A    C    T
   A 2.0  4.5  1.0
   C 4.5 10.0  0.5
   T 1.0  0.5 12.0
   <BLANKLINE>

As the ordering of pairs is arbitrary, counts are divided equally above
and below the diagonal. For example, the 9 alignments of ``A`` to ``C``
are stored as 4.5 at position ``['A', 'C']`` and 4.5 at position
``['C', 'A']``. This arrangement helps to make the math easier when
calculating a substitution matrix from these counts, as described in
Section :ref:`sec:subs_mat_ex`.

Note that ``alignment.substitutions`` contains entries for the letters
appearing in the alignment only. You can use the ``select`` method to
add entries for missing letters, for example

.. code:: pycon

   >>> m = substitutions.select("ATCG")
   >>> print(m)
       A    T    C   G
   A 2.0  1.0  4.5 0.0
   T 1.0 12.0  0.5 0.0
   C 4.5  0.5 10.0 0.0
   G 0.0  0.0  0.0 0.0
   <BLANKLINE>

This also allows you to change the order of letters in the alphabet.

.. _`sec:alignment-tools`:

Alignment Tools
---------------

There are *lots* of algorithms out there for aligning sequences, both
pairwise alignments and multiple sequence alignments. These calculations
are relatively slow, and you generally wouldn’t want to write such an
algorithm in Python. For pairwise alignments Biopython contains
``PairwiseAligner`` (see Section :ref:`sec:pairwise`). In addition,
you can use Biopython to invoke a command line tool on your behalf.
Normally you would:

#. Prepare an input file of your unaligned sequences, typically this
   will be a FASTA file which you might create using ``Bio.SeqIO`` (see
   Chapter :ref:`chapter:seqio`).

#. Call the command line tool to process this input file, typically via
   one of Biopython’s command line wrappers (which we’ll discuss here).

#. Read the output from the tool, i.e. your aligned sequences, typically
   using ``Bio.AlignIO`` (see earlier in this chapter).

All the command line wrappers we’re going to talk about in this chapter
follow the same style. You create a command line object specifying the
options (e.g. the input filename and the output filename), then invoke
this command line via a Python operating system call (e.g. using the
``subprocess`` module).

*WARNING:* We have decided to drop these command line wrappers in a
future Biopython release. We will be updating this documentation to
instead build the command line directly, and invoke it with the
``subprocess`` module.

Most of these wrappers are defined in the ``Bio.Align.Applications``
module:

.. code:: pycon

   >>> import Bio.Align.Applications
   >>> dir(Bio.Align.Applications)  # doctest:+ELLIPSIS
   ['ClustalOmegaCommandline', 'ClustalwCommandline', 'DialignCommandline', 'MSAProbsCommandline', 'MafftCommandline', 'MuscleCommandline', 'PrankCommandline', 'ProbconsCommandline', 'TCoffeeCommandline', ...]

(Ignore the entries starting with an underscore – these have special
meaning in Python.) The module ``Bio.Emboss.Applications`` has wrappers
for some of the `EMBOSS suite <http://emboss.sourceforge.net/>`__,
including ``needle`` and ``water``, which are described below in
Section :ref:`sec:emboss-needle-water`, and wrappers for the
EMBOSS packaged versions of the PHYLIP tools (which EMBOSS refer to as
one of their EMBASSY packages - third party tools with an EMBOSS style
interface). We won’t explore all these alignment tools here in the
section, just a sample, but the same principles apply.

.. _`sec:align_clustal`:

ClustalW
~~~~~~~~

ClustalW is a popular command line tool for multiple sequence alignment
(there is also a graphical interface called ClustalX). Biopython’s
``Bio.Align.Applications`` module has a wrapper for this alignment tool
(and several others).

Before trying to use ClustalW from within Python, you should first try
running the ClustalW tool yourself by hand at the command line, to
familiarize yourself the other options. You’ll find the Biopython
wrapper is very faithful to the actual command line API:

.. code:: pycon

   >>> from Bio.Align.Applications import ClustalwCommandline
   >>> help(ClustalwCommandline)

For the most basic usage, all you need is to have a FASTA input file,
such as
`opuntia.fasta <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.fasta>`__
(available online or in the Doc/examples subdirectory of the Biopython
source code). This is a small FASTA file containing seven prickly-pear
DNA sequences (from the cactus family *Opuntia*).

By default ClustalW will generate an alignment and guide tree file with
names based on the input FASTA file, in this case ``opuntia.aln`` and
``opuntia.dnd``, but you can override this or make it explicit:

.. code:: pycon

   >>> from Bio.Align.Applications import ClustalwCommandline
   >>> cline = ClustalwCommandline("clustalw2", infile="opuntia.fasta")
   >>> print(cline)
   clustalw2 -infile=opuntia.fasta

Notice here we have given the executable name as ``clustalw2``,
indicating we have version two installed, which has a different filename
to version one (``clustalw``, the default). Fortunately both versions
support the same set of arguments at the command line (and indeed,
should be functionally identical).

You may find that even though you have ClustalW installed, the above
command doesn’t work – you may get a message about “command not found”
(especially on Windows). This indicated that the ClustalW executable is
not on your PATH (an environment variable, a list of directories to be
searched). You can either update your PATH setting to include the
location of your copy of ClustalW tools (how you do this will depend on
your OS), or simply type in the full path of the tool. For example:

.. code:: pycon

   >>> import os
   >>> from Bio.Align.Applications import ClustalwCommandline
   >>> clustalw_exe = r"C:\Program Files\new clustal\clustalw2.exe"
   >>> clustalw_cline = ClustalwCommandline(clustalw_exe, infile="opuntia.fasta")

.. code:: pycon

   >>> assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
   >>> stdout, stderr = clustalw_cline()

Remember, in Python strings ``\n`` and ``\t`` are by default interpreted
as a new line and a tab – which is why we’re put a letter “r” at the
start for a raw string that isn’t translated in this way. This is
generally good practice when specifying a Windows style file name.

Internally this uses the ``subprocess`` module which is now the
recommended way to run another program in Python. This replaces older
options like the ``os.system()`` and the ``os.popen*`` functions.

Now, at this point it helps to know about how command line tools “work”.
When you run a tool at the command line, it will often print text output
directly to screen. This text can be captured or redirected, via two
“pipes”, called standard output (the normal results) and standard error
(for error messages and debug messages). There is also standard input,
which is any text fed into the tool. These names get shortened to stdin,
stdout and stderr. When the tool finishes, it has a return code (an
integer), which by convention is zero for success.

When you run the command line tool like this via the Biopython wrapper,
it will wait for it to finish, and check the return code. If this is non
zero (indicating an error), an exception is raised. The wrapper then
returns two strings, stdout and stderr.

In the case of ClustalW, when run at the command line all the important
output is written directly to the output files. Everything normally
printed to screen while you wait (via stdout or stderr) is boring and
can be ignored (assuming it worked).

What we care about are the two output files, the alignment and the guide
tree. We didn’t tell ClustalW what filenames to use, but it defaults to
picking names based on the input file. In this case the output should be
in the file ``opuntia.aln``. You should be able to work out how to read
in the alignment using ``Bio.AlignIO`` by now:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> align = AlignIO.read("opuntia.aln", "clustal")
   >>> print(align)
   Alignment with 7 rows and 906 columns
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
   TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191

In case you are interested (and this is an aside from the main thrust of
this chapter), the ``opuntia.dnd`` file ClustalW creates is just a
standard Newick tree file, and ``Bio.Phylo`` can parse these:

.. code:: pycon

   >>> from Bio import Phylo
   >>> tree = Phylo.read("opuntia.dnd", "newick")
   >>> Phylo.draw_ascii(tree)
                                _______________ gi|6273291|gb|AF191665.1|AF191665
     __________________________|
    |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
    |                          |__|
    |                             |_____ gi|6273289|gb|AF191663.1|AF191663
    |
   _|_________________ gi|6273287|gb|AF191661.1|AF191661
    |
    |__________ gi|6273286|gb|AF191660.1|AF191660
    |
    |    __ gi|6273285|gb|AF191659.1|AF191659
    |___|
        | gi|6273284|gb|AF191658.1|AF191658
   <BLANKLINE>

Chapter :ref:`chapter:phylo` covers Biopython’s support
for phylogenetic trees in more depth.

MUSCLE
~~~~~~

MUSCLE is a more recent multiple sequence alignment tool than ClustalW,
and Biopython also has a wrapper for it under the
``Bio.Align.Applications`` module. As before, we recommend you try using
MUSCLE from the command line before trying it from within Python, as the
Biopython wrapper is very faithful to the actual command line API:

.. code:: pycon

   >>> from Bio.Align.Applications import MuscleCommandline
   >>> help(MuscleCommandline)

For the most basic usage, all you need is to have a FASTA input file,
such as
`opuntia.fasta <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.fasta>`__
(available online or in the Doc/examples subdirectory of the Biopython
source code). You can then tell MUSCLE to read in this FASTA file, and
write the alignment to an output file:

.. code:: pycon

   >>> from Bio.Align.Applications import MuscleCommandline
   >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.txt")
   >>> print(cline)
   muscle -in opuntia.fasta -out opuntia.txt

Note that MUSCLE uses “-in” and “-out” but in Biopython we have to use
“input” and “out” as the keyword arguments or property names. This is
because “in” is a reserved word in Python.

By default MUSCLE will output the alignment as a FASTA file (using
gapped sequences). The ``Bio.AlignIO`` module should be able to read
this alignment using ``format="fasta"``. You can also ask for
ClustalW-like output:

.. code:: pycon

   >>> from Bio.Align.Applications import MuscleCommandline
   >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.aln", clw=True)
   >>> print(cline)
   muscle -in opuntia.fasta -out opuntia.aln -clw

Or, strict ClustalW output where the original ClustalW header line is
used for maximum compatibility:

.. code:: pycon

   >>> from Bio.Align.Applications import MuscleCommandline
   >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.aln", clwstrict=True)
   >>> print(cline)
   muscle -in opuntia.fasta -out opuntia.aln -clwstrict

The ``Bio.AlignIO`` module should be able to read these alignments using
``format="clustal"``.

MUSCLE can also output in GCG MSF format (using the ``msf`` argument),
but Biopython can’t currently parse that, or using HTML which would give
a human readable web page (not suitable for parsing).

You can also set the other optional parameters, for example the maximum
number of iterations. See the built in help for details.

You would then run MUSCLE command line string as described above for
ClustalW, and parse the output using ``Bio.AlignIO`` to get an alignment
object.

MUSCLE using stdout
~~~~~~~~~~~~~~~~~~~

Using a MUSCLE command line as in the examples above will write the
alignment to a file. This means there will be no important information
written to the standard out (stdout) or standard error (stderr) handles.
However, by default MUSCLE will write the alignment to standard output
(stdout). We can take advantage of this to avoid having a temporary
output file! For example:

.. code:: pycon

   >>> from Bio.Align.Applications import MuscleCommandline
   >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
   >>> print(muscle_cline)
   muscle -in opuntia.fasta

If we run this via the wrapper, we get back the output as a string. In
order to parse this we can use ``StringIO`` to turn it into a handle.
Remember that MUSCLE defaults to using FASTA as the output format:

.. code:: pycon

   >>> from Bio.Align.Applications import MuscleCommandline
   >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
   >>> stdout, stderr = muscle_cline()
   >>> from io import StringIO
   >>> from Bio import AlignIO
   >>> align = AlignIO.read(StringIO(stdout), "fasta")
   >>> print(align)
   Alignment with 7 rows and 906 columns
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191663
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191665
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191664
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191661
   TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191660
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191659
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191658

The above approach is fairly simple, but if you are dealing with very
large output text the fact that all of stdout and stderr is loaded into
memory as a string can be a potential drawback. Using the ``subprocess``
module we can work directly with handles instead:

.. code:: pycon

   >>> import subprocess
   >>> from Bio.Align.Applications import MuscleCommandline
   >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
   >>> child = subprocess.Popen(
   ...     str(muscle_cline),
   ...     stdout=subprocess.PIPE,
   ...     stderr=subprocess.PIPE,
   ...     text=True,
   ...     shell=(sys.platform != "win32"),
   ... )
   >>> from Bio import AlignIO
   >>> align = AlignIO.read(child.stdout, "fasta")
   >>> print(align)
   Alignment with 7 rows and 906 columns
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191663
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191665
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191664
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191661
   TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191660
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191659
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191658

MUSCLE using stdin and stdout
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We don’t actually *need* to have our FASTA input sequences prepared in a
file, because by default MUSCLE will read in the input sequence from
standard input! Note this is a bit more advanced and fiddly, so don’t
bother with this technique unless you need to.

First, we’ll need some unaligned sequences in memory as ``SeqRecord``
objects. For this demonstration I’m going to use a filtered version of
the original FASTA file (using a generator expression), taking just six
of the seven sequences:

.. code:: pycon

   >>> from Bio import SeqIO
   >>> records = (r for r in SeqIO.parse("opuntia.fasta", "fasta") if len(r) < 900)

Then we create the MUSCLE command line, leaving the input and output to
their defaults (stdin and stdout). I’m also going to ask for strict
ClustalW format as for the output.

.. code:: pycon

   >>> from Bio.Align.Applications import MuscleCommandline
   >>> muscle_cline = MuscleCommandline(clwstrict=True)
   >>> print(muscle_cline)
   muscle -clwstrict

Now for the fiddly bits using the ``subprocess`` module, stdin and
stdout:

.. code:: pycon

   >>> import subprocess
   >>> import sys
   >>> child = subprocess.Popen(
   ...     str(cline),
   ...     stdin=subprocess.PIPE,
   ...     stdout=subprocess.PIPE,
   ...     stderr=subprocess.PIPE,
   ...     text=True,
   ...     shell=(sys.platform != "win32"),
   ... )

That should start MUSCLE, but it will be sitting waiting for its FASTA
input sequences, which we must supply via its stdin handle:

.. code:: pycon

   >>> SeqIO.write(records, child.stdin, "fasta")
   6
   >>> child.stdin.close()

After writing the six sequences to the handle, MUSCLE will still be
waiting to see if that is all the FASTA sequences or not – so we must
signal that this is all the input data by closing the handle. At that
point MUSCLE should start to run, and we can ask for the output:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> align = AlignIO.read(child.stdout, "clustal")
   >>> print(align)
   Alignment with 6 rows and 900 columns
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF19166
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF19166
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF19166
   TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF19166
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF19165
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF19165

Wow! There we are with a new alignment of just the six records, without
having created a temporary FASTA input file, or a temporary alignment
output file. However, a word of caution: Dealing with errors with this
style of calling external programs is much more complicated. It also
becomes far harder to diagnose problems, because you can’t try running
MUSCLE manually outside of Biopython (because you don’t have the input
file to supply). There can also be subtle cross platform issues (e.g.
Windows versus Linux), and how you run your script can have an impact
(e.g. at the command line, from IDLE or an IDE, or as a GUI script).
These are all generic Python issues though, and not specific to
Biopython.

If you find working directly with ``subprocess`` like this scary, there
is an alternative. If you execute the tool with ``muscle_cline()`` you
can supply any standard input as a big string,
``muscle_cline(stdin=...)``. So, provided your data isn’t very big, you
can prepare the FASTA input in memory as a string using ``StringIO``
(see Section :ref:`sec:appendix-handles`):

.. code:: pycon

   >>> from Bio import SeqIO
   >>> records = (r for r in SeqIO.parse("opuntia.fasta", "fasta") if len(r) < 900)
   >>> from io import StringIO
   >>> handle = StringIO()
   >>> SeqIO.write(records, handle, "fasta")
   6
   >>> data = handle.getvalue()

You can then run the tool and parse the alignment as follows:

.. code:: pycon

   >>> stdout, stderr = muscle_cline(stdin=data)
   >>> from Bio import AlignIO
   >>> align = AlignIO.read(StringIO(stdout), "clustal")
   >>> print(align)
   Alignment with 6 rows and 900 columns
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF19166
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF19166
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF19166
   TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF19166
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF19165
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF19165

You might find this easier, but it does require more memory (RAM) for
the strings used for the input FASTA and output Clustal formatted data.

.. _`sec:emboss-needle-water`:

EMBOSS needle and water
~~~~~~~~~~~~~~~~~~~~~~~

The `EMBOSS <http://emboss.sourceforge.net/>`__ suite includes the
``water`` and ``needle`` tools for Smith-Waterman algorithm local
alignment, and Needleman-Wunsch global alignment. The tools share the
same style interface, so switching between the two is trivial – we’ll
just use ``needle`` here.

Suppose you want to do a global pairwise alignment between two
sequences, prepared in FASTA format as follows:

.. code:: text

   >HBA_HUMAN
   MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
   KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
   AVHASLDKFLASVSTVLTSKYR

in a file ``alpha.faa``, and secondly in a file ``beta.faa``:

.. code:: text

   >HBB_HUMAN
   MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
   VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
   KEFTPPVQAAYQKVVAGVANALAHKYH

You can find copies of these example files with the Biopython source
code under the ``Doc/examples/`` directory.

Let’s start by creating a complete ``needle`` command line object in one
go:

.. code:: pycon

   >>> from Bio.Emboss.Applications import NeedleCommandline
   >>> needle_cline = NeedleCommandline(
   ...     asequence="alpha.faa",
   ...     bsequence="beta.faa",
   ...     gapopen=10,
   ...     gapextend=0.5,
   ...     outfile="needle.txt",
   ... )
   >>> print(needle_cline)
   needle -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5

Why not try running this by hand at the command prompt? You should see
it does a pairwise comparison and records the output in the file
``needle.txt`` (in the default EMBOSS alignment file format).

Even if you have EMBOSS installed, running this command may not work –
you might get a message about “command not found” (especially on
Windows). This probably means that the EMBOSS tools are not on your PATH
environment variable. You can either update your PATH setting, or simply
tell Biopython the full path to the tool, for example:

.. code:: pycon

   >>> from Bio.Emboss.Applications import NeedleCommandline
   >>> needle_cline = NeedleCommandline(
   ...     r"C:\EMBOSS\needle.exe",
   ...     asequence="alpha.faa",
   ...     bsequence="beta.faa",
   ...     gapopen=10,
   ...     gapextend=0.5,
   ...     outfile="needle.txt",
   ... )

Remember in Python that for a default string ``\n`` or ``\t`` means a
new line or a tab – which is why we’re put a letter “r” at the start for
a raw string.

At this point it might help to try running the EMBOSS tools yourself by
hand at the command line, to familiarize yourself the other options and
compare them to the Biopython help text:

.. code:: pycon

   >>> from Bio.Emboss.Applications import NeedleCommandline
   >>> help(NeedleCommandline)

Note that you can also specify (or change or look at) the settings like
this:

.. code:: pycon

   >>> from Bio.Emboss.Applications import NeedleCommandline
   >>> needle_cline = NeedleCommandline()
   >>> needle_cline.asequence = "alpha.faa"
   >>> needle_cline.bsequence = "beta.faa"
   >>> needle_cline.gapopen = 10
   >>> needle_cline.gapextend = 0.5
   >>> needle_cline.outfile = "needle.txt"
   >>> print(needle_cline)
   needle -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5
   >>> print(needle_cline.outfile)
   needle.txt

Next we want to use Python to run this command for us. As explained
above, for full control, we recommend you use the built in Python
``subprocess`` module, but for simple usage the wrapper object usually
suffices:

.. code:: pycon

   >>> stdout, stderr = needle_cline()
   >>> print(stdout + stderr)
   Needleman-Wunsch global alignment of two sequences

Next we can load the output file with ``Bio.AlignIO`` as discussed
earlier in this chapter, as the ``emboss`` format:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> align = AlignIO.read("needle.txt", "emboss")
   >>> print(align)
   Alignment with 2 rows and 149 columns
   MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY...KYR HBA_HUMAN
   MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRF...KYH HBB_HUMAN

In this example, we told EMBOSS to write the output to a file, but you
*can* tell it to write the output to stdout instead (useful if you don’t
want a temporary output file to get rid of – use ``stdout=True`` rather
than the ``outfile`` argument), and also to read *one* of the one of the
inputs from stdin (e.g. ``asequence="stdin"``, much like in the MUSCLE
example in the section above).

This has only scratched the surface of what you can do with ``needle``
and ``water``. One useful trick is that the second file can contain
multiple sequences (say five), and then EMBOSS will do five pairwise
alignments.

.. _`sec:pairwise`:

Pairwise sequence alignment
---------------------------

Pairwise sequence alignment is the process of aligning two sequences to
each other by optimizing the similarity score between them. The
``Bio.Align`` module contains the ``PairwiseAligner`` class for global
and local alignments using the Needleman-Wunsch, Smith-Waterman, Gotoh
(three-state), and Waterman-Smith-Beyer global and local pairwise
alignment algorithms, with numerous options to change the alignment
parameters. We refer to Durbin *et al.* :raw-latex:`\cite{durbin1998}`
for in-depth information on sequence alignment algorithms.

.. _`sec:pairwise-basic`:

Basic usage
~~~~~~~~~~~

To generate pairwise alignments, first create a ``PairwiseAligner``
object:

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()

The ``PairwiseAligner`` object ``aligner`` (see
Section :ref:`sec:pairwise-aligner`) stores the alignment
parameters to be used for the pairwise alignments.

These attributes can be set in the constructor of the object or after
the object is made.

.. code:: pycon

   >>> aligner = Align.PairwiseAligner(match_score=1.0)

Or, equivalently:

.. code:: pycon

   >>> aligner.match_score = 1.0

Use the ``aligner.score`` method to calculate the alignment score
between two sequences:

.. code:: pycon

   >>> target = "GAACT"
   >>> query = "GAT"
   >>> score = aligner.score(target, query)
   >>> score
   3.0

To see the actual alignments, use the ``aligner.align`` method and
iterate over the ``Alignment`` objects returned:

.. code:: pycon

   >>> alignments = aligner.align(target, query)
   >>> for alignment in alignments:
   ...     print(alignment)
   ...
   target            0 GAACT 5
                     0 ||--| 5
   query             0 GA--T 3
   <BLANKLINE>
   target            0 GAACT 5
                     0 |-|-| 5
   query             0 G-A-T 3
   <BLANKLINE>

By default, a global pairwise alignment is performed, which finds the
optimal alignment over the whole length of ``target`` and ``query``.
Instead, a local alignment will find the subsequence of ``target`` and
``query`` with the highest alignment score. Local alignments can be
generated by setting ``aligner.mode`` to ``"local"``:

.. code:: pycon

   >>> aligner.mode = "local"
   >>> target = "AGAACTC"
   >>> query = "GAACT"
   >>> score = aligner.score(target, query)
   >>> score
   5.0
   >>> alignments = aligner.align(target, query)
   >>> for alignment in alignments:
   ...     print(alignment)
   ...
   target            1 GAACT 6
                     0 ||||| 5
   query             0 GAACT 5
   <BLANKLINE>

Note that there is some ambiguity in the definition of the best local
alignments if segments with a score 0 can be added to the alignment. We
follow the suggestion by Waterman & Eggert
:raw-latex:`\cite{waterman1987}` and disallow such extensions.

.. _`sec:pairwise-aligner`:

The pairwise aligner object
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``PairwiseAligner`` object stores all alignment parameters to be
used for the pairwise alignments. To see an overview of the values for
all parameters, use

.. code:: pycon

   >>> print(aligner)
   Pairwise sequence aligner with parameters
     wildcard: None
     match_score: 1.000000
     mismatch_score: 0.000000
     target_internal_open_gap_score: 0.000000
     target_internal_extend_gap_score: 0.000000
     target_left_open_gap_score: 0.000000
     target_left_extend_gap_score: 0.000000
     target_right_open_gap_score: 0.000000
     target_right_extend_gap_score: 0.000000
     query_internal_open_gap_score: 0.000000
     query_internal_extend_gap_score: 0.000000
     query_left_open_gap_score: 0.000000
     query_left_extend_gap_score: 0.000000
     query_right_open_gap_score: 0.000000
     query_right_extend_gap_score: 0.000000
     mode: local
   <BLANKLINE>

See Sections :ref:`sec:pairwise-substitution-scores`,
:ref:`sec:pairwise-affine-gapscores`, and
:ref:`sec:pairwise-general-gapscores` below for the definition of
these parameters. The attribute ``mode`` (described above in
Section :ref:`sec:pairwise-basic`) can be set equal to
``"global"`` or ``"local"`` to specify global or local pairwise
alignment, respectively.

Depending on the gap scoring parameters (see
Sections :ref:`sec:pairwise-affine-gapscores` and
:ref:`sec:pairwise-general-gapscores`) and mode, a
``PairwiseAligner`` object automatically chooses the appropriate
algorithm to use for pairwise sequence alignment. To verify the selected
algorithm, use

.. code:: pycon

   >>> aligner.algorithm
   'Smith-Waterman'

This attribute is read-only.

A ``PairwiseAligner`` object also stores the precision :math:`\epsilon`
to be used during alignment. The value of :math:`\epsilon` is stored in
the attribute ``aligner.epsilon``, and by default is equal to
:math:`10^{-6}`:

.. code:: pycon

   >>> aligner.epsilon
   1e-06

Two scores will be considered equal to each other for the purpose of the
alignment if the absolute difference between them is less than
:math:`\epsilon`.

.. _`sec:pairwise-substitution-scores`:

Substitution scores
~~~~~~~~~~~~~~~~~~~

Substitution scores define the value to be added to the total score when
two letters (nucleotides or amino acids) are aligned to each other. The
substitution scores to be used by the ``PairwiseAligner`` can be
specified in two ways:

-  By specifying a match score for identical letters, and a mismatch
   scores for mismatched letters. Nucleotide sequence alignments are
   typically based on match and mismatch scores. For example, by default
   BLAST :raw-latex:`\cite{altschul1990}` uses a match score of
   :math:`+1` and a mismatch score of :math:`-2` for nucleotide
   alignments by ``megablast``, with a gap penalty of 2.5 (see section
   :ref:`sec:pairwise-affine-gapscores` for more information on
   gap scores). Match and mismatch scores can be specified by setting
   the ``match`` and ``mismatch`` attributes of the ``PairwiseAligner``
   object:

   .. code:: pycon

      >>> from Bio import Align
      >>> aligner = Align.PairwiseAligner()
      >>> aligner.match_score
      1.0
      >>> aligner.mismatch_score
      0.0
      >>> score = aligner.score("ACGT", "ACAT")
      >>> print(score)
      3.0
      >>> aligner.match_score = 1.0
      >>> aligner.mismatch_score = -2.0
      >>> aligner.gap_score = -2.5
      >>> score = aligner.score("ACGT", "ACAT")
      >>> print(score)
      1.0

   When using match and mismatch scores, you can specify a wildcard
   character (``None`` by default) for unknown letters. These will get a
   zero score in alignments, irrespective of the value of the match or
   mismatch score:

   .. code:: pycon

      >>> aligner.wildcard = "?"
      >>> score = aligner.score("ACGT", "AC?T")
      >>> print(score)
      3.0

-  Alternatively, you can use the ``substitution_matrix`` attribute of
   the ``PairwiseAligner`` object to specify a substitution matrix. This
   allows you to apply different scores for different pairs of matched
   and mismatched letters. This is typically used for amino acid
   sequence alignments. For example, by default BLAST
   :raw-latex:`\cite{altschul1990}` uses the BLOSUM62 substitution
   matrix for protein alignments by ``blastp``. This substitution matrix
   is available from Biopython:

   .. code:: pycon

      >>> from Bio.Align import substitution_matrices
      >>> substitution_matrices.load()  # doctest: +ELLIPSIS
      ['BENNER22', 'BENNER6', 'BENNER74', 'BLASTN', 'BLASTP', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62', ..., 'TRANS']
      >>> matrix = substitution_matrices.load("BLOSUM62")
      >>> print(matrix)  # doctest: +ELLIPSIS
      #  Matrix made by matblas from blosum62.iij
      ...
           A    R    N    D    C    Q ...
      A  4.0 -1.0 -2.0 -2.0  0.0 -1.0 ...
      R -1.0  5.0  0.0 -2.0 -3.0  1.0 ...
      N -2.0  0.0  6.0  1.0 -3.0  0.0 ...
      D -2.0 -2.0  1.0  6.0 -3.0  0.0 ...
      C  0.0 -3.0 -3.0 -3.0  9.0 -3.0 ...
      Q -1.0  1.0  0.0  0.0 -3.0  5.0 ...
      ...
      >>> aligner.substitution_matrix = matrix
      >>> score = aligner.score("ACDQ", "ACDQ")
      >>> score
      24.0
      >>> score = aligner.score("ACDQ", "ACNQ")
      >>> score
      19.0

   When using a substitution matrix, ``X`` is *not* interpreted as an
   unknown character. Instead, the score provided by the substitution
   matrix will be used:

   .. code:: pycon

      >>> matrix["D", "X"]
      -1.0
      >>> score = aligner.score("ACDQ", "ACXQ")
      >>> score
      17.0

By default, ``aligner.substitution_matrix`` is ``None``. The attributes
``aligner.match_score`` and ``aligner.mismatch_score`` are ignored if
``aligner.substitution_matrix`` is not ``None``. Setting
``aligner.match_score`` or ``aligner.mismatch_score`` to valid values
will reset ``aligner.substitution_matrix`` to ``None``.

.. _`sec:pairwise-affine-gapscores`:

Affine gap scores
~~~~~~~~~~~~~~~~~

Affine gap scores are defined by a score to open a gap, and a score to
extend an existing gap:

:math:`\textrm{gap score} = \textrm{open gap score} + (n-1) \times \textrm{extend gap score}`,

where :math:`n` is the length of the gap. Biopython’s pairwise sequence
aligner allows fine-grained control over the gap scoring scheme by
specifying the following twelve attributes of a ``PairwiseAligner``
object:

================================== ====================================
**Opening scores**                 **Extending scores**
================================== ====================================
``query_left_open_gap_score``      ``query_left_extend_gap_score``
``query_internal_open_gap_score``  ``query_internal_extend_gap_score``
``query_right_open_gap_score``     ``query_right_extend_gap_score``
``target_left_open_gap_score``     ``target_left_extend_gap_score``
``target_internal_open_gap_score`` ``target_internal_extend_gap_score``
``target_right_open_gap_score``    ``target_right_extend_gap_score``
================================== ====================================

These attributes allow for different gap scores for internal gaps and on
either end of the sequence, as shown in this example:

========== ========= ================================
**target** **query** **score**
========== ========= ================================
A          -         query left open gap score
C          -         query left extend gap score
C          -         query left extend gap score
G          G         match score
G          T         mismatch score
G          -         query internal open gap score
A          -         query internal extend gap score
A          -         query internal extend gap score
T          T         match score
A          A         match score
G          -         query internal open gap score
C          C         match score
-          C         target internal open gap score
-          C         target internal extend gap score
C          C         match score
T          G         mismatch score
C          C         match score
-          C         target internal open gap score
A          A         match score
-          T         target right open gap score
-          A         target right extend gap score
-          A         target right extend gap score
========== ========= ================================

For convenience, ``PairwiseAligner`` objects have additional attributes
that refer to a number of these values collectively, as shown
(hierarchically) in Table :ref:`table:align-meta-attributes`.

.. container::
   :name: table:align-meta-attributes

   .. table:: Meta-attributes of the pairwise aligner objects.

      +---------------------------------+-----------------------------------+
      | **Meta-attribute**              | **Attributes it maps to**         |
      +=================================+===================================+
      | ``gap_score``                   | ``target_gap_score``,             |
      |                                 | ``query_gap_score``               |
      +---------------------------------+-----------------------------------+
      | ``open_gap_score``              | ``target_open_gap_score``,        |
      |                                 | ``query_open_gap_score``          |
      +---------------------------------+-----------------------------------+
      | ``extend_gap_score``            | ``target_extend_gap_score``,      |
      |                                 | ``query_extend_gap_score``        |
      +---------------------------------+-----------------------------------+
      | ``internal_gap_score``          | ``target_internal_gap_score``,    |
      |                                 | ``query_internal_gap_score``      |
      +---------------------------------+-----------------------------------+
      | ``internal_open_gap_score``     | ``                                |
      |                                 | target_internal_open_gap_score``, |
      |                                 | ``query_internal_open_gap_score`` |
      +---------------------------------+-----------------------------------+
      | ``internal_extend_gap_score``   | ``ta                              |
      |                                 | rget_internal_extend_gap_score``, |
      |                                 | ``                                |
      |                                 | query_internal_extend_gap_score`` |
      +---------------------------------+-----------------------------------+
      | ``end_gap_score``               | ``target_end_gap_score``,         |
      |                                 | ``query_end_gap_score``           |
      +---------------------------------+-----------------------------------+
      | ``end_open_gap_score``          | ``target_end_open_gap_score``,    |
      |                                 | ``query_end_open_gap_score``      |
      +---------------------------------+-----------------------------------+
      | ``end_extend_gap_score``        | ``target_end_extend_gap_score``,  |
      |                                 | ``query_end_extend_gap_score``    |
      +---------------------------------+-----------------------------------+
      | ``left_gap_score``              | ``target_left_gap_score``,        |
      |                                 | ``query_left_gap_score``          |
      +---------------------------------+-----------------------------------+
      | ``right_gap_score``             | ``target_right_gap_score``,       |
      |                                 | ``query_right_gap_score``         |
      +---------------------------------+-----------------------------------+
      | ``left_open_gap_score``         | ``target_left_open_gap_score``,   |
      |                                 | ``query_left_open_gap_score``     |
      +---------------------------------+-----------------------------------+
      | ``left_extend_gap_score``       | ``target_left_extend_gap_score``, |
      |                                 | ``query_left_extend_gap_score``   |
      +---------------------------------+-----------------------------------+
      | ``right_open_gap_score``        | ``target_right_open_gap_score``,  |
      |                                 | ``query_right_open_gap_score``    |
      +---------------------------------+-----------------------------------+
      | ``right_extend_gap_score``      | `                                 |
      |                                 | `target_right_extend_gap_score``, |
      |                                 | ``query_right_extend_gap_score``  |
      +---------------------------------+-----------------------------------+
      | ``target_open_gap_score``       | ``                                |
      |                                 | target_internal_open_gap_score``, |
      |                                 | ``target_left_open_gap_score``,   |
      +---------------------------------+-----------------------------------+
      |                                 | ``target_right_open_gap_score``   |
      +---------------------------------+-----------------------------------+
      | ``target_extend_gap_score``     | ``ta                              |
      |                                 | rget_internal_extend_gap_score``, |
      |                                 | ``target_left_extend_gap_score``, |
      +---------------------------------+-----------------------------------+
      |                                 | ``target_right_extend_gap_score`` |
      +---------------------------------+-----------------------------------+
      | ``target_gap_score``            | ``target_open_gap_score``,        |
      |                                 | ``target_extend_gap_score``       |
      +---------------------------------+-----------------------------------+
      | ``query_open_gap_score``        | `                                 |
      |                                 | `query_internal_open_gap_score``, |
      |                                 | ``query_left_open_gap_score``,    |
      +---------------------------------+-----------------------------------+
      |                                 | ``query_right_open_gap_score``    |
      +---------------------------------+-----------------------------------+
      | ``query_extend_gap_score``      | ``q                               |
      |                                 | uery_internal_extend_gap_score``, |
      |                                 | ``query_left_extend_gap_score``,  |
      +---------------------------------+-----------------------------------+
      |                                 | ``query_right_extend_gap_score``  |
      +---------------------------------+-----------------------------------+
      | ``query_gap_score``             | ``query_open_gap_score``,         |
      |                                 | ``query_extend_gap_score``        |
      +---------------------------------+-----------------------------------+
      | ``target_internal_gap_score``   | ``                                |
      |                                 | target_internal_open_gap_score``, |
      |                                 | ``t                               |
      |                                 | arget_internal_extend_gap_score`` |
      +---------------------------------+-----------------------------------+
      | ``target_end_gap_score``        | ``target_end_open_gap_score``,    |
      |                                 | ``target_end_extend_gap_score``   |
      +---------------------------------+-----------------------------------+
      | ``target_end_open_gap_score``   | ``target_left_open_gap_score``,   |
      |                                 | ``target_right_open_gap_score``   |
      +---------------------------------+-----------------------------------+
      | ``target_end_extend_gap_score`` | ``target_left_extend_gap_score``, |
      |                                 | ``target_right_extend_gap_score`` |
      +---------------------------------+-----------------------------------+
      | ``target_left_gap_score``       | ``target_left_open_gap_score``,   |
      |                                 | ``target_left_extend_gap_score``  |
      +---------------------------------+-----------------------------------+
      | ``target_right_gap_score``      | ``target_right_open_gap_score``,  |
      |                                 | ``target_right_extend_gap_score`` |
      +---------------------------------+-----------------------------------+
      | ``query_end_gap_score``         | ``query_end_open_gap_score``,     |
      |                                 | ``query_end_extend_gap_score``    |
      +---------------------------------+-----------------------------------+
      | ``query_end_open_gap_score``    | ``query_left_open_gap_score``,    |
      |                                 | ``query_right_open_gap_score``    |
      +---------------------------------+-----------------------------------+
      | ``query_end_extend_gap_score``  | ``query_left_extend_gap_score``,  |
      |                                 | ``query_right_extend_gap_score``  |
      +---------------------------------+-----------------------------------+
      | ``query_internal_gap_score``    | `                                 |
      |                                 | `query_internal_open_gap_score``, |
      |                                 | ``                                |
      |                                 | query_internal_extend_gap_score`` |
      +---------------------------------+-----------------------------------+
      | ``query_left_gap_score``        | ``query_left_open_gap_score``,    |
      |                                 | ``query_left_extend_gap_score``   |
      +---------------------------------+-----------------------------------+
      | ``query_right_gap_score``       | ``query_right_open_gap_score``,   |
      |                                 | ``query_right_extend_gap_score``  |
      +---------------------------------+-----------------------------------+

.. _`sec:pairwise-general-gapscores`:

General gap scores
~~~~~~~~~~~~~~~~~~

For even more fine-grained control over the gap scores, you can specify
a gap scoring function. For example, the gap scoring function below
disallows a gap after two nucleotides in the query sequence:

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()
   >>> def my_gap_score_function(start, length):
   ...     if start == 2:
   ...         return -1000
   ...     else:
   ...         return -1 * length
   ...
   >>> aligner.query_gap_score = my_gap_score_function
   >>> alignments = aligner.align("AACTT", "AATT")
   >>> for alignment in alignments:
   ...     print(alignment)
   ...
   target            0 AACTT 5
                     0 -|.|| 5
   query             0 -AATT 4
   <BLANKLINE>
   target            0 AACTT 5
                     0 |-.|| 5
   query             0 A-ATT 4
   <BLANKLINE>
   target            0 AACTT 5
                     0 ||.-| 5
   query             0 AAT-T 4
   <BLANKLINE>
   target            0 AACTT 5
                     0 ||.|- 5
   query             0 AATT- 4
   <BLANKLINE>

.. _`sec:pairwise-predefined-scoring`:

Using a pre-defined substitution matrix and gap scores
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, a ``PairwiseAligner`` object is initialized with a match
score of +1.0, a mismatch score of 0.0, and all gap scores equal to 0.0,
While this has the benefit of being a simple scoring scheme, in general
it does not give the best performance. Instead, you can use the argument
``scoring`` to select a predefined scoring scheme when initializing a
``PairwiseAligner`` object. Currently, the provided scoring schemes are
``blastn`` and ``megablast``, which are suitable for nucleotide
alignments, and ``blastp``, which is suitable for protein alignments.
Selecting these scoring schemes will initialize the ``PairwiseAligner``
object to the default scoring parameters used by BLASTN, MegaBLAST, and
BLASTP, respectively.

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner(scoring="blastn")
   >>> print(aligner)  # doctest:+ELLIPSIS
   Pairwise sequence aligner with parameters
     substitution_matrix: <Array object at ...>
     target_internal_open_gap_score: -7.000000
     target_internal_extend_gap_score: -2.000000
     target_left_open_gap_score: -7.000000
     target_left_extend_gap_score: -2.000000
     target_right_open_gap_score: -7.000000
     target_right_extend_gap_score: -2.000000
     query_internal_open_gap_score: -7.000000
     query_internal_extend_gap_score: -2.000000
     query_left_open_gap_score: -7.000000
     query_left_extend_gap_score: -2.000000
     query_right_open_gap_score: -7.000000
     query_right_extend_gap_score: -2.000000
     mode: global
   <BLANKLINE>
   >>> print(aligner.substitution_matrix[:, :])
        A    T    G    C    S    W    R    Y    K    M    B    V    H    D    N
   A  2.0 -3.0 -3.0 -3.0 -3.0 -1.0 -1.0 -3.0 -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -2.0
   T -3.0  2.0 -3.0 -3.0 -3.0 -1.0 -3.0 -1.0 -1.0 -3.0 -1.0 -3.0 -1.0 -1.0 -2.0
   G -3.0 -3.0  2.0 -3.0 -1.0 -3.0 -1.0 -3.0 -1.0 -3.0 -1.0 -1.0 -3.0 -1.0 -2.0
   C -3.0 -3.0 -3.0  2.0 -1.0 -3.0 -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -3.0 -2.0
   S -3.0 -3.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   W -1.0 -1.0 -3.0 -3.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   R -1.0 -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   Y -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   K -3.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -2.0
   M -1.0 -3.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   B -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   V -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   H -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   D -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   N -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0
   <BLANKLINE>

Iterating over alignments
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``alignments`` returned by ``aligner.align`` are a kind of immutable
iterable objects (similar to ``range``). While they appear similar to a
``tuple`` or ``list`` of ``Alignment`` objects, they are different in
the sense that each ``Alignment`` object is created dynamically when it
is needed. This approach was chosen because the number of alignments can
be extremely large, in particular for poor alignments (see
Section :ref:`sec:pairwise-examples` for an example).

You can perform the following operations on ``alignments``:

-  ``len(alignments)`` returns the number of alignments stored. This
   function returns quickly, even if the number of alignments is huge.
   If the number of alignments is extremely large (typically, larger
   than 9,223,372,036,854,775,807, which is the largest integer that can
   be stored as a ``long int`` on 64 bit machines), ``len(alignments)``
   will raise an ``OverflowError``. A large number of alignments
   suggests that the alignment quality is low.

   .. code:: pycon

      >>> from Bio import Align
      >>> aligner = Align.PairwiseAligner()
      >>> alignments = aligner.align("AAA", "AA")
      >>> len(alignments)
      3

-  You can extract a specific alignment by index:

   .. code:: pycon

      >>> from Bio import Align
      >>> aligner = Align.PairwiseAligner()
      >>> alignments = aligner.align("AAA", "AA")
      >>> print(alignments[2])
      target            0 AAA 3
                        0 -|| 3
      query             0 -AA 2
      <BLANKLINE>
      >>> print(alignments[0])
      target            0 AAA 3
                        0 ||- 3
      query             0 AA- 2
      <BLANKLINE>

-  You can iterate over alignments, for example as in

   .. code:: pycon

      >>> for alignment in alignments:
      ...     print(alignment)
      ...

   By calling ``alignments.rewind``, you can rewind the ``alignments``
   iterator to the first alignment and iterate over the alignments from
   the beginning:

   .. code:: pycon

      >>> from Bio import Align
      >>> aligner = Align.PairwiseAligner()
      >>> alignments = aligner.align("AAA", "AA")
      >>> for alignment in alignments:
      ...     print(alignment)
      ...
      target            0 AAA 3
                        0 ||- 3
      query             0 AA- 2
      <BLANKLINE>
      target            0 AAA 3
                        0 |-| 3
      query             0 A-A 2
      <BLANKLINE>
      target            0 AAA 3
                        0 -|| 3
      query             0 -AA 2
      <BLANKLINE>
      >>> alignments.rewind()
      >>> for alignment in alignments:
      ...     print(alignment)
      ...
      target            0 AAA 3
                        0 ||- 3
      query             0 AA- 2
      <BLANKLINE>
      target            0 AAA 3
                        0 |-| 3
      query             0 A-A 2
      <BLANKLINE>
      target            0 AAA 3
                        0 -|| 3
      query             0 -AA 2
      <BLANKLINE>

   You can also convert the ``alignments`` iterator into a ``list`` or
   ``tuple``:

   .. code:: pycon

      >>> alignments = list(alignments)

   It is wise to check the number of alignments by calling
   ``len(alignments)`` before attempting to call ``list(alignments)`` to
   save all alignments as a list.

-  The alignment score (which has the same value for each alignment in
   ``alignments``) is stored as an attribute. This allows you to check
   the alignment score before proceeding to extract individual
   alignments:

   .. code:: pycon

      >>> print(alignments.score)
      2.0

Alignment objects
~~~~~~~~~~~~~~~~~

The ``aligner.align`` method returns ``Alignment`` objects, each
representing one alignment between the two sequences.

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()
   >>> target = "GAACT"
   >>> query = "GAT"
   >>> alignments = aligner.align(target, query)
   >>> alignment = alignments[0]
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (2 rows x 5 columns) at ...>

Each alignment stores the alignment score:

.. code:: pycon

   >>> alignment.score
   3.0

as well as pointers to the sequences that were aligned:

.. code:: pycon

   >>> alignment.target
   'GAACT'
   >>> alignment.query
   'GAT'

Print the ``Alignment`` object to show the alignment explicitly:

.. code:: pycon

   >>> print(alignment)
   target            0 GAACT 5
                     0 ||--| 5
   query             0 GA--T 3
   <BLANKLINE>

Internally, the alignment is stored in terms of the sequence
coordinates:

.. code:: pycon

   >>> alignment.coordinates
   array([[0, 2, 4, 5],
          [0, 2, 2, 3]])

Here, the two rows refer to the target and query sequence. These
coordinates show that the alignment consists of the following three
blocks:

-  ``target[0:2]`` aligned to ``query[0:2]``;

-  ``target[2:4]`` aligned to a gap, since ``query[2:2]`` is an empty
   string;

-  ``target[4:5]`` aligned to ``query[2:3]``.

The length of the alignment is defined as the number of aligned
sequences, which is always 2 for a pairwise alignment:

.. code:: pycon

   >>> len(alignment)
   2

The ``shape`` property returns a tuple consisting of the length of the
alignment and the number of columns in the alignment as printed:

.. code:: pycon

   >>> alignment.shape
   (2, 5)

For local alignments, sections that are not aligned are not included in
the number of columns:

.. code:: pycon

   >>> aligner.mode = "local"
   >>> local_alignments = aligner.align("TGAACT", "GAC")
   >>> local_alignment = local_alignments[0]
   >>> print(local_alignment)
   target            1 GAAC 5
                     0 ||-| 4
   query             0 GA-C 3
   <BLANKLINE>
   >>> local_alignment.shape
   (2, 4)

Use the ``aligned`` property to find the start and end indices of
subsequences in the target and query sequence that were aligned to each
other. Generally, if the alignment between target (t) and query (q)
consists of :math:`N` chunks, you get a numpy array with dimensions
:math:`2 \times N \times 2`:

.. code:: python

   (
       ((t_start1, t_end1), (t_start2, t_end2), ..., (t_startN, t_endN)),
       ((q_start1, q_end1), (q_start2, q_end2), ..., (q_startN, q_endN)),
   )

In the current example, ``alignment.aligned`` returns two tuples of
length 2:

.. code:: pycon

   >>> alignment.aligned
   array([[[0, 2],
           [4, 5]],
   <BLANKLINE>
          [[0, 2],
           [2, 3]]])

while for the alternative alignment, two tuples of length 3 are
returned:

.. code:: pycon

   >>> alignment = alignments[1]
   >>> print(alignment)
   target            0 GAACT 5
                     0 |-|-| 5
   query             0 G-A-T 3
   <BLANKLINE>
   >>> alignment.aligned
   array([[[0, 1],
           [2, 3],
           [4, 5]],
   <BLANKLINE>
          [[0, 1],
           [1, 2],
           [2, 3]]])

Note that different alignments may have the same subsequences aligned to
each other. In particular, this may occur if alignments differ from each
other in terms of their gap placement only:

.. code:: pycon

   >>> aligner.mode = "global"
   >>> aligner.mismatch_score = -10
   >>> alignments = aligner.align("AAACAAA", "AAAGAAA")
   >>> len(alignments)
   2
   >>> print(alignments[0])
   target            0 AAAC-AAA 7
                     0 |||--||| 8
   query             0 AAA-GAAA 7
   <BLANKLINE>
   >>> alignments[0].aligned
   array([[[0, 3],
           [4, 7]],
   <BLANKLINE>
          [[0, 3],
           [4, 7]]])
   >>> print(alignments[1])
   target            0 AAA-CAAA 7
                     0 |||--||| 8
   query             0 AAAG-AAA 7
   <BLANKLINE>
   >>> alignments[1].aligned
   array([[[0, 3],
           [4, 7]],
   <BLANKLINE>
          [[0, 3],
           [4, 7]]])

The ``aligned`` property can be used to identify alignments that are
identical to each other in terms of their aligned sequences.

The ``sort`` method sorts the alignment sequences. By default, sorting
is done based on the ``id`` attribute of each sequence if available, or
the sequence contents otherwise.

.. code:: pycon

   >>> print(local_alignment)
   target            1 GAAC 5
                     0 ||-| 4
   query             0 GA-C 3
   <BLANKLINE>
   >>> local_alignment.sort()
   >>> print(local_alignment)
   target            0 GA-C 3
                     0 ||-| 4
   query             1 GAAC 5
   <BLANKLINE>

Alternatively, you can supply a ``key`` function to determine the sort
order. For example, you can sort the sequences by increasing GC content:

.. code:: pycon

   >>> from Bio.SeqUtils import gc_fraction
   >>> local_alignment.sort(key=gc_fraction)
   >>> print(local_alignment)
   target            1 GAAC 5
                     0 ||-| 4
   query             0 GA-C 3
   <BLANKLINE>

The ``reverse`` argument lets you reverse the sort order to obtain the
sequences in decreasing GC content:

.. code:: pycon

   >>> local_alignment.sort(key=gc_fraction, reverse=True)
   >>> print(local_alignment)
   target            0 GA-C 3
                     0 ||-| 4
   query             1 GAAC 5
   <BLANKLINE>

Use the ``substitutions`` method to find the number of substitutions
between each pair of nucleotides:

.. code:: pycon

   >>> target = "AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT"
   >>> query = "AAAAAAACCCTCCCCGGCCGGGGTTTAGTTT"
   >>> alignments = aligner.align(target, query)
   >>> aligner.mismatch_score = -1
   >>> aligner.gap_score = -1
   >>> alignments = aligner.align(target, query)
   >>> len(alignments)
   8
   >>> print(alignments[0])
   target            0 AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT 32
                     0 |||||||-|||.||||||..|||||||..||| 32
   query             0 AAAAAAA-CCCTCCCCGGCCGGGGTTTAGTTT 31
   <BLANKLINE>
   >>> m = alignments[0].substitutions
   >>> print(m)
       A   C   G   T
   A 7.0 0.0 0.0 0.0
   C 0.0 7.0 0.0 1.0
   G 0.0 2.0 6.0 0.0
   T 1.0 0.0 1.0 6.0
   <BLANKLINE>

Note that the matrix is not symmetric: rows correspond to the target
sequence, and columns to the query sequence. For example, the number of
G’s in the target sequence that are aligned to a C in the query sequence
is

.. code:: pycon

   >>> m["G", "C"]
   2.0

and the number of C’s in the query sequence tat are aligned to a T in
the query sequence is

.. code:: pycon

   >>> m["C", "G"]
   0.0

To get a symmetric matrix, use

.. code:: pycon

   >>> m += m.transpose()
   >>> m /= 2.0
   >>> print(m)
       A   C   G   T
   A 7.0 0.0 0.0 0.5
   C 0.0 7.0 1.0 0.5
   G 0.0 1.0 6.0 0.5
   T 0.5 0.5 0.5 6.0
   <BLANKLINE>
   >>> m["G", "C"]
   1.0
   >>> m["C", "G"]
   1.0

The total number of substitutions between C’s and G’s in the alignment
is 1.0 + 1.0 = 2.

The ``map`` method can be applied on a pairwise alignment ``alignment1``
to find the pairwise alignment of the query of ``alignment2`` to the
target of ``alignment1``, where the target of ``alignment2`` and the
query of ``alignment1`` are identical. A typical example is where
``alignment1`` is the pairwise alignment between a chromosome and a
transcript, ``alignment2`` is the pairwise alignment between the
transcript and a sequence (e.g., an RNA-seq read), and we want to find
the alignment of the sequence to the chromosome:

.. code:: pycon

   >>> aligner.mode = "local"
   >>> aligner.open_gap_score = -1
   >>> aligner.extend_gap_score = 0
   >>> chromosome = "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA"
   >>> transcript = "CCCCCCCGGGGGG"
   >>> alignments1 = aligner.align(chromosome, transcript)
   >>> len(alignments1)
   1
   >>> alignment1 = alignments1[0]
   >>> print(alignment1)
   target            8 CCCCCCCAAAAAAAAAAAGGGGGG 32
                     0 |||||||-----------|||||| 24
   query             0 CCCCCCC-----------GGGGGG 13
   <BLANKLINE>
   >>> sequence = "CCCCGGGG"
   >>> alignments2 = aligner.align(transcript, sequence)
   >>> len(alignments2)
   1
   >>> alignment2 = alignments2[0]
   >>> print(alignment2)
   target            3 CCCCGGGG 11
                     0 ||||||||  8
   query             0 CCCCGGGG  8
   <BLANKLINE>
   >>> mapped_alignment = alignment1.map(alignment2)
   >>> print(mapped_alignment)
   target           11 CCCCAAAAAAAAAAAGGGG 30
                     0 ||||-----------|||| 19
   query             0 CCCC-----------GGGG  8
   <BLANKLINE>
   >>> format(mapped_alignment, "psl")
   '8\t0\t0\t0\t0\t0\t1\t11\t+\tquery\t8\t0\t8\ttarget\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n'

Mapping the alignment does not depend on the sequence contents. If we
delete the sequence contents, the same alignment is found in PSL format
(though we obviously lose the ability to print the sequence alignment):

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> alignment1.target = Seq(None, len(alignment1.target))
   >>> alignment1.query = Seq(None, len(alignment1.query))
   >>> alignment2.target = Seq(None, len(alignment2.target))
   >>> alignment2.query = Seq(None, len(alignment2.query))
   >>> mapped_alignment = alignment1.map(alignment2)
   >>> format(mapped_alignment, "psl")
   '8\t0\t0\t0\t0\t0\t1\t11\t+\tquery\t8\t0\t8\ttarget\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n'

Slicing and indexing a pairwise alignment
'''''''''''''''''''''''''''''''''''''''''

Slices of the form ``alignment[k, i:j]``, where ``k`` is an integer and
``i`` and ``j`` are integers or are absent, return a string showing the
aligned sequence (including gaps) for the target (if ``k=0``) or the
query (if ``k=1``) that includes only the columns ``i`` through ``j`` in
the printed alignment.

To illustrate this, in the following example the printed alignment has 5
columns:

.. code:: pycon

   >>> print(alignment)
   target            0 GAACT 5
                     0 |-|-| 5
   query             0 G-A-T 3
   <BLANKLINE>

To get the aligned sequence strings individually, use

.. code:: pycon

   >>> alignment[0]
   'GAACT'
   >>> alignment[1]
   'G-A-T'
   >>> alignment[0, :]
   'GAACT'
   >>> alignment[1, :]
   'G-A-T'
   >>> alignment[0, 1:-1]
   'AAC'
   >>> alignment[1, 1:-1]
   '-A-'

Columns to be included can also be selected using an iterable over
integers:

.. code:: pycon

   >>> alignment[0, (1, 3, 4)]
   'ACT'
   >>> alignment[1, range(0, 5, 2)]
   'GAT'

To get specific columns in the alignment, use

.. code:: pycon

   >>> alignment[:, 0]
   'GG'
   >>> alignment[:, 1]
   'A-'
   >>> alignment[:, 2]
   'AA'

Slices of the form ``alignment[:, i:j]``, where ``i`` and ``j`` are
integers or are absent, return a new ``Alignment`` object that includes
only the columns ``i`` through ``j`` in the printed alignment.

Extracting the first 4 columns for the example alignment above gives:

.. code:: pycon

   >>> alignment[:, :4]  # doctest:+ELLIPSIS
   <Alignment object (2 rows x 4 columns) at ...>
   >>> print(alignment[:, :4])
   target            0 GAAC 4
                     0 |-|- 4
   query             0 G-A- 2
   <BLANKLINE>

Here, the final ``T`` nucleotides are still shown, but they are not
aligned to each other. Note that ``alignment`` is a global alignment,
but ``alignment[:, :4]`` is a local alignment.

Similarly, extracting the last 3 columns gives:

.. code:: pycon

   >>> alignment[:, -3:]  # doctest:+ELLIPSIS
   <Alignment object (2 rows x 3 columns) at ...>
   >>> print(alignment[:, -3:])
   target            2 ACT 5
                     0 |-| 3
   query             1 A-T 3
   <BLANKLINE>

This is also now a local alignment, with the initial ``GA`` nucleotides
in the target and ``G`` nucleotide in the query not aligned to each
other.

The column index can also be an iterable of integers:

.. code:: pycon

   >>> alignment[:, -3:]  # doctest:+ELLIPSIS
   <Alignment object (2 rows x 3 columns) at ...>
   >>> print(alignment[:, (1, 3, 0)])
   target            0 ACG 3
                     0 --| 3
   query             0 --G 1
   <BLANKLINE>

Reverse-complementing the alignment
'''''''''''''''''''''''''''''''''''

Reverse_complementing an alignment will take the reverse complement of
each sequence, and recalculate the coordinates:

.. code:: pycon

   >>> print(alignment.sequences)
   ['GAACT', 'GAT']
   >>> rc_alignment = alignment.reverse_complement()
   >>> print(rc_alignment.sequences)
   ['AGTTC', 'ATC']
   >>> print(rc_alignment)
   target            0 AGTTC 5
                     0 |-|-| 5
   query             0 A-T-C 3
   <BLANKLINE>
   >>> print(alignment[:, :4].sequences)
   ['GAACT', 'GAT']
   >>> print(alignment[:, :4])
   target            0 GAAC 4
                     0 |-|- 4
   query             0 G-A- 2
   <BLANKLINE>
   >>> rc_alignment = alignment[:, :4].reverse_complement()
   >>> print(rc_alignment[:, :4].sequences)
   ['AGTTC', 'ATC']
   >>> print(rc_alignment[:, :4])
   target            1 GTTC 5
                     0 -|-| 4
   query             1 -T-C 3
   <BLANKLINE>

Reverse-complementing an alignment preserves its column annotations (in
reverse order), but discards all other annotations.

Exporting alignments
''''''''''''''''''''

Use the ``format`` method to create a string representation of the
alignment in various file formats. This method takes an argument ``fmt``
specifying the file format, and may take additional keyword arguments
depending on file type. The following values for ``fmt`` are supported:

-  ``""`` (empty string; default): Create a human-readable
   representation of the alignment (same as when you ``print`` the
   alignment).

-  ``"SAM"``: Create a line representing the alignment in the Sequence
   Alignment/Map (SAM) format:

   .. code:: pycon

      >>> alignment.format("sam")
      'query\t0\ttarget\t1\t255\t1M1D1M1D1M\t*\t0\t0\tGAT\t*\tAS:i:3\n'

-  ``"BED"``: Create a line representing the alignment in the Browser
   Extensible Data (BED) file format:

   .. code:: pycon

      >>> alignment.format("bed")
      'target\t0\t5\tquery\t3.0\t+\t0\t5\t0\t3\t1,1,1,\t0,2,4,\n'

-  ``"PSL"``: Create a line representing the alignment in the Pattern
   Space Layout (PSL) file format as generated by BLAT
   :raw-latex:`\cite{kent2002}`).

   .. code:: pycon

      >>> alignment.format("psl")
      '3\t0\t0\t0\t0\t0\t2\t2\t+\tquery\t3\t0\t3\ttarget\t5\t0\t5\t3\t1,1,1,\t0,1,2,\t0,2,4,\n'

   The first four columns in the PSL output contain the number of
   matched and mismatched characters, the number of matches to repeat
   regions, and the number of matches to unknown nucleotides. Repeat
   regions in the target sequence are indicated by masking the sequence
   as lower-case or upper-case characters, as defined by the following
   values for the ``mask`` keyword argument:

   -  ``False`` (default): Do not count matches to masked sequences
      separately;

   -  ``"lower"``: Count and report matches to lower-case characters as
      matches to repeat regions;

   -  ``"upper"``: Count and report matches to upper-case characters as
      matches to repeat regions;

   The character used for unknown nucleotides is defined by the
   ``wildcard`` argument. For consistency with BLAT, the wildcard
   character is ``"N"`` by default. Use ``wildcard=None`` if you don’t
   want to count matches to any unknown nucleotides separately.

   .. code:: pycon

      >>> from Bio import Align
      >>> aligner = Align.PairwiseAligner()
      >>> aligner.mismatch_score = -1
      >>> aligner.internal_gap_score = -5
      >>> aligner.wildcard = "N"
      >>> target = "AAAAAAAggggGGNGAAAAA"
      >>> query = "GGTGGGGG"
      >>> alignments = aligner.align(target.upper(), query)
      >>> print(len(alignments))
      1
      >>> alignment = alignments[0]
      >>> print(alignment)
      target            0 AAAAAAAGGGGGGNGAAAAA 20
                        0 -------||.|||.|----- 20
      query             0 -------GGTGGGGG-----  8
      <BLANKLINE>
      >>> alignment.score
      5.0
      >>> alignment.target
      'AAAAAAAGGGGGGNGAAAAA'
      >>> alignment.target = target
      >>> alignment.target
      'AAAAAAAggggGGNGAAAAA'
      >>> print(alignment)
      target            0 AAAAAAAggggGGNGAAAAA 20
                        0 -------....||.|----- 20
      query             0 -------GGTGGGGG-----  8
      <BLANKLINE>
      >>> print(alignment.format("psl"))  # doctest: +NORMALIZE_WHITESPACE
      6   1   0   1   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,
      >>> print(alignment.format("psl", mask="lower"))  # doctest: +NORMALIZE_WHITESPACE
      3   1   3   1   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,
      >>> print(
      ...     alignment.format("psl", mask="lower", wildcard=None)
      ... )  # doctest: +NORMALIZE_WHITESPACE
      3   2   3   0   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,

In addition to the ``format`` method, you can use Python’s built-in
``format`` function:

.. code:: pycon

   >>> print(format(alignment, "psl"))  # doctest: +NORMALIZE_WHITESPACE
   6   1   0   1   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,

allowing ``Alignment`` objects to be used in formatted (f-) strings in
Python:

.. code:: pycon

   >>> print(
   ...     f"The alignment in PSL format is '{alignment:psl}'."
   ... )  # doctest: +NORMALIZE_WHITESPACE
   The alignment in PSL format is '6   1   0   1   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,
   '

Note that optional keyword arguments cannot be used with the ``format``
function or with formatted strings.

Aligning to the reverse strand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the pairwise aligner aligns the forward strand of the query
to the forward strand of the target. To calculate the alignment score
for ``query`` to the reverse strand of ``target``, use ``strand="-"``:

.. code:: pycon

   >>> from Bio import Align
   >>> from Bio.Seq import reverse_complement
   >>> target = "AAAACCC"
   >>> query = "AACC"
   >>> aligner = Align.PairwiseAligner()
   >>> aligner.mismatch_score = -1
   >>> aligner.internal_gap_score = -1
   >>> aligner.score(target, query)  # strand is "+" by default
   4.0
   >>> aligner.score(target, reverse_complement(query), strand="-")
   4.0
   >>> aligner.score(target, query, strand="-")
   0.0
   >>> aligner.score(target, reverse_complement(query))
   0.0

The alignments against the reverse strand can be obtained by specifying
``strand="-"`` when calling ``aligner.align``:

.. code:: pycon

   >>> alignments = aligner.align(target, query)
   >>> len(alignments)
   1
   >>> print(alignments[0])
   target            0 AAAACCC 7
                     0 --||||- 7
   query             0 --AACC- 4
   <BLANKLINE>
   >>> print(alignments[0].format("bed"))  # doctest: +NORMALIZE_WHITESPACE
   target   2   6   query   4.0   +   2   6   0   1   4,   0,
   <BLANKLINE>
   >>> alignments = aligner.align(target, reverse_complement(query), strand="-")
   >>> len(alignments)
   1
   >>> print(alignments[0])
   target            0 AAAACCC 7
                     0 --||||- 7
   query             4 --AACC- 0
   <BLANKLINE>
   >>> print(alignments[0].format("bed"))  # doctest: +NORMALIZE_WHITESPACE
   target   2   6   query   4.0   -   2   6   0   1   4,   0,
   <BLANKLINE>
   >>> alignments = aligner.align(target, query, strand="-")
   >>> len(alignments)
   2
   >>> print(alignments[0])
   target            0 AAAACCC----  7
                     0 ----------- 11
   query             4 -------GGTT  0
   <BLANKLINE>
   >>> print(alignments[1])
   target            0 ----AAAACCC  7
                     0 ----------- 11
   query             4 GGTT-------  0
   <BLANKLINE>

Note that the score for aligning ``query`` to the reverse strand of
``target`` may be different from the score for aligning the reverse
complement of ``query`` to the forward strand of ``target`` if the left
and right gap scores are different:

.. code:: pycon

   >>> aligner.left_gap_score = -0.5
   >>> aligner.right_gap_score = -0.2
   >>> aligner.score(target, query)
   2.8
   >>> alignments = aligner.align(target, query)
   >>> len(alignments)
   1
   >>> print(alignments[0])
   target            0 AAAACCC 7
                     0 --||||- 7
   query             0 --AACC- 4
   <BLANKLINE>
   >>> aligner.score(target, reverse_complement(query), strand="-")
   3.1
   >>> alignments = aligner.align(target, reverse_complement(query), strand="-")
   >>> len(alignments)
   1
   >>> print(alignments[0])
   target            0 AAAACCC 7
                     0 --||||- 7
   query             4 --AACC- 0
   <BLANKLINE>

.. _`sec:pairwise-examples`:

Examples
~~~~~~~~

Suppose you want to do a global pairwise alignment between the same two
hemoglobin sequences from above (``HBA_HUMAN``, ``HBB_HUMAN``) stored in
``alpha.faa`` and ``beta.faa``:

.. code:: pycon

   >>> from Bio import Align
   >>> from Bio import SeqIO
   >>> seq1 = SeqIO.read("alpha.faa", "fasta")
   >>> seq2 = SeqIO.read("beta.faa", "fasta")
   >>> aligner = Align.PairwiseAligner()
   >>> score = aligner.score(seq1.seq, seq2.seq)
   >>> print(score)
   72.0

showing an alignment score of 72.0. To see the individual alignments, do

.. code:: pycon

   >>> alignments = aligner.align(seq1.seq, seq2.seq)

In this example, the total number of optimal alignments is huge (more
than :math:`4 \times 10^{37}`), and calling ``len(alignments)`` will
raise an ``OverflowError``:

.. code:: pycon

   >>> len(alignments)
   Traceback (most recent call last):
   ...
   OverflowError: number of optimal alignments is larger than 9223372036854775807

Let’s have a look at the first alignment:

.. code:: pycon

   >>> alignment = alignments[0]

The alignment object stores the alignment score, as well as the
alignment itself:

.. code:: pycon

   >>> print(alignment.score)
   72.0
   >>> print(alignment)
   target            0 MV-LS-PAD--KTN--VK-AA-WGKV-----GAHAGEYGAEALE-RMFLSF----P-TTK
                     0 ||-|--|----|----|--|--||||-----|---||--|--|--|--|------|-|--
   query             0 MVHL-TP--EEK--SAV-TA-LWGKVNVDEVG---GE--A--L-GR--L--LVVYPWT--
   <BLANKLINE>
   target           41 TY--FPHF----DLSHGS---AQVK-G------HGKKV--A--DA-LTNAVAHV-DDMPN
                    60 ----|--|----|||------|-|--|------|||||--|--|--|--|--|--|---|
   query            39 --QRF--FESFGDLS---TPDA-V-MGNPKVKAHGKKVLGAFSD-GL--A--H-LD---N
   <BLANKLINE>
   target           79 ALS----A-LSD-LHAH--KLR-VDPV-NFK-LLSHC---LLVT--LAAHLPA----EFT
                   120 -|-----|-||--||----||--|||--||--||------|-|---||-|-------|||
   query            81 -L-KGTFATLS-ELH--CDKL-HVDP-ENF-RLL---GNVL-V-CVLA-H---HFGKEFT
   <BLANKLINE>
   target          119 PA-VH-ASLDKFLAS---VSTV------LTS--KYR- 142
                   180 |--|--|------|----|--|------|----||-- 217
   query           124 P-PV-QA------A-YQKV--VAGVANAL--AHKY-H 147
   <BLANKLINE>

Better alignments are usually obtained by penalizing gaps: higher costs
for opening a gap and lower costs for extending an existing gap. For
amino acid sequences match scores are usually encoded in matrices like
``PAM`` or ``BLOSUM``. Thus, a more meaningful alignment for our example
can be obtained by using the BLOSUM62 matrix, together with a gap open
penalty of 10 and a gap extension penalty of 0.5:

.. code:: pycon

   >>> from Bio import Align
   >>> from Bio import SeqIO
   >>> from Bio.Align import substitution_matrices
   >>> seq1 = SeqIO.read("alpha.faa", "fasta")
   >>> seq2 = SeqIO.read("beta.faa", "fasta")
   >>> aligner = Align.PairwiseAligner()
   >>> aligner.open_gap_score = -10
   >>> aligner.extend_gap_score = -0.5
   >>> aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
   >>> score = aligner.score(seq1.seq, seq2.seq)
   >>> print(score)
   292.5
   >>> alignments = aligner.align(seq1.seq, seq2.seq)
   >>> len(alignments)
   2
   >>> print(alignments[0].score)
   292.5
   >>> print(alignments[0])
   target            0 MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGS
                     0 ||-|.|..|..|.|.||||--...|.|.|||.|.....|.|...|..|-|||-----.|.
   query             0 MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGN
   <BLANKLINE>
   target           53 AQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAH
                    60 ..||.|||||..|.....||.|........||.||..||.|||.||.||...|...||.|
   query            58 PKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHH
   <BLANKLINE>
   target          113 LPAEFTPAVHASLDKFLASVSTVLTSKYR 142
                   120 ...||||.|.|...|..|.|...|..||. 149
   query           118 FGKEFTPPVQAAYQKVVAGVANALAHKYH 147
   <BLANKLINE>

This alignment has the same score that we obtained earlier with EMBOSS
needle using the same sequences and the same parameters.

To perform a local alignment, set ``aligner.mode`` to ``'local'``:

.. code:: pycon

   >>> aligner.mode = "local"
   >>> aligner.open_gap_score = -10
   >>> aligner.extend_gap_score = -1
   >>> alignments = aligner.align("LSPADKTNVKAA", "PEEKSAV")
   >>> print(len(alignments))
   1
   >>> alignment = alignments[0]
   >>> print(alignment)
   target            2 PADKTNV 9
                     0 |..|..| 7
   query             0 PEEKSAV 7
   <BLANKLINE>
   >>> print(alignment.score)
   16.0

.. _`sec:generalized-pairwise`:

Generalized pairwise alignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In most cases, ``PairwiseAligner`` is used to perform alignments of
sequences (strings or ``Seq`` objects) consisting of single-letter
nucleotides or amino acids. More generally, ``PairwiseAligner`` can also
be applied to lists or tuples of arbitrary objects. This section will
describe some examples of such generalized pairwise alignments.

Generalized pairwise alignments using a substitution matrix and alphabet
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Schneider *et al.* :raw-latex:`\cite{schneider2005}` created a
substitution matrix for aligning three-nucleotide codons (see
`below <#codonmatrix>`__ in section :ref:`sec:substitution_matrices`
for more information). This substitution matrix is associated with an
alphabet consisting of all three-letter codons:

.. code:: pycon

   >>> from Bio.Align import substitution_matrices
   >>> m = substitution_matrices.load("SCHNEIDER")
   >>> m.alphabet  # doctest: +ELLIPSIS
   ('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', ..., 'TTG', 'TTT')

We can use this matrix to align codon sequences to each other:

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()
   >>> aligner.substitution_matrix = m
   >>> aligner.gap_score = -1.0
   >>> s1 = ("AAT", "CTG", "TTT", "TTT")
   >>> s2 = ("AAT", "TTA", "TTT")
   >>> alignments = aligner.align(s1, s2)
   >>> len(alignments)
   2
   >>> print(alignments[0])
   AAT CTG TTT TTT
   ||| ... ||| ---
   AAT TTA TTT ---
   <BLANKLINE>
   >>> print(alignments[1])
   AAT CTG TTT TTT
   ||| ... --- |||
   AAT TTA --- TTT
   <BLANKLINE>

Note that aligning ``TTT`` to ``TTA``, as in this example:

.. code:: pycon

   AAT CTG TTT TTT
   ||| --- ... |||
   AAT --- TTA TTT

would get a much lower score:

.. code:: pycon

   >>> print(m["CTG", "TTA"])
   7.6
   >>> print(m["TTT", "TTA"])
   -0.3

presumably because ``CTG`` and ``TTA`` both code for leucine, while
``TTT`` codes for phenylalanine. The three-letter codon substitution
matrix also reveals a preference among codons representing the same
amino acid. For example, ``TTA`` has a preference for ``CTG`` preferred
compared to ``CTC``, though all three code for leucine:

.. code:: pycon

   >>> s1 = ("AAT", "CTG", "CTC", "TTT")
   >>> s2 = ("AAT", "TTA", "TTT")
   >>> alignments = aligner.align(s1, s2)
   >>> len(alignments)
   1
   >>> print(alignments[0])
   AAT CTG CTC TTT
   ||| ... --- |||
   AAT TTA --- TTT
   <BLANKLINE>
   >>> print(m["CTC", "TTA"])
   6.5

Generalized pairwise alignments using match/mismatch scores and an alphabet
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Using the three-letter amino acid symbols, the sequences above translate
to

.. code:: pycon

   >>> s1 = ("Asn", "Leu", "Leu", "Phe")
   >>> s2 = ("Asn", "Leu", "Phe")

We can align these sequences directly to each other by using a
three-letter amino acid alphabet:

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()
   >>> aligner.alphabet = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys',
   ...                     'Gln', 'Glu', 'Gly', 'His', 'Ile',
   ...                     'Leu', 'Lys', 'Met', 'Phe', 'Pro',
   ...                     'Ser', 'Thr', 'Trp', 'Tyr', 'Val']  # fmt: skip
   ...

We use +6/-1 match and mismatch scores as an approximation of the
BLOSUM62 matrix, and align these sequences to each other:

.. code:: pycon

   >>> aligner.match = +6
   >>> aligner.mismatch = -1
   >>> alignments = aligner.align(s1, s2)
   >>> print(len(alignments))
   2
   >>> print(alignments[0])
   Asn Leu Leu Phe
   ||| ||| --- |||
   Asn Leu --- Phe
   <BLANKLINE>
   >>> print(alignments[1])
   Asn Leu Leu Phe
   ||| --- ||| |||
   Asn --- Leu Phe
   <BLANKLINE>
   >>> print(alignments.score)
   18.0

Generalized pairwise alignments using match/mismatch scores and integer sequences
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Internally, the first step when performing an alignment is to replace
the two sequences by integer arrays consisting of the indices of each
letter in each sequence in the alphabet associated with the aligner.
This step can be bypassed by passing integer arrays directly:

.. code:: pycon

   >>> import numpy as np
   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()
   >>> s1 = np.array([2, 10, 10, 13], np.int32)
   >>> s2 = np.array([2, 10, 13], np.int32)
   >>> aligner.match = +6
   >>> aligner.mismatch = -1
   >>> alignments = aligner.align(s1, s2)
   >>> print(len(alignments))
   2
   >>> print(alignments[0])
   2 10 10 13
   | || -- ||
   2 10 -- 13
   <BLANKLINE>
   >>> print(alignments[1])
   2 10 10 13
   | -- || ||
   2 -- 10 13
   <BLANKLINE>
   >>> print(alignments.score)
   18.0

Note that the indices should consist of 32-bit integers, as specified in
this example by ``numpy.int32``.

Unknown letters can again be included by defining a wildcard character,
and using the corresponding Unicode code point number as the index:

.. code:: pycon

   >>> aligner.wildcard = "?"
   >>> ord(aligner.wildcard)
   63
   >>> s2 = np.array([2, 63, 13], np.int32)
   >>> aligner.gap_score = -3
   >>> alignments = aligner.align(s1, s2)
   >>> print(len(alignments))
   2
   >>> print(alignments[0])
   2 10 10 13
   | .. -- ||
   2 63 -- 13
   <BLANKLINE>
   >>> print(alignments[1])
   2 10 10 13
   | -- .. ||
   2 -- 63 13
   <BLANKLINE>
   >>> print(alignments.score)
   9.0

Generalized pairwise alignments using a substitution matrix and integer sequences
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Integer sequences can also be aligned using a substitution matrix, in
this case a numpy square array without an alphabet associated with it.
In this case, all index values must be non-negative, and smaller than
the size of the substitution matrix:

.. code:: pycon

   >>> from Bio import Align
   >>> import numpy as np
   >>> aligner = Align.PairwiseAligner()
   >>> m = np.eye(5)
   >>> m[0, 1:] = m[1:, 0] = -2
   >>> m[2, 2] = 3
   >>> print(m)
   [[ 1. -2. -2. -2. -2.]
    [-2.  1.  0.  0.  0.]
    [-2.  0.  3.  0.  0.]
    [-2.  0.  0.  1.  0.]
    [-2.  0.  0.  0.  1.]]
   >>> aligner.substitution_matrix = m
   >>> aligner.gap_score = -1
   >>> s1 = np.array([0, 2, 3, 4], np.int32)
   >>> s2 = np.array([0, 3, 2, 1], np.int32)
   >>> alignments = aligner.align(s1, s2)
   >>> print(len(alignments))
   2
   >>> print(alignments[0])
   0 - 2 3 4
   | - | . -
   0 3 2 1 -
   <BLANKLINE>
   >>> print(alignments[1])
   0 - 2 3 4
   | - | - .
   0 3 2 - 1
   <BLANKLINE>
   >>> print(alignments.score)
   2.0

.. _`sec:substitution_matrices`:

Substitution matrices
---------------------

The ``Array`` class in ``Bio.Align.substitution_matrices`` is a subclass
of numpy arrays that supports indexing both by integers and by specific
strings. An ``Array`` instance can either be a one-dimensional array or
a square two-dimensional arrays. A one-dimensional ``Array`` object can
for example be used to store the nucleotide frequency of a DNA sequence,
while a two-dimensional ``Array`` object can be used to represent a
scoring matrix for sequence alignments.

Creating an Array object
~~~~~~~~~~~~~~~~~~~~~~~~

To create a one-dimensional ``Array``, only the alphabet of allowed
letters needs to be specified:

.. code:: pycon

   >>> from Bio.Align.substitution_matrices import Array
   >>> counts = Array("ACGT")
   >>> print(counts)
   A 0.0
   C 0.0
   G 0.0
   T 0.0
   <BLANKLINE>

The allowed letters are stored in the ``alphabet`` property:

.. code:: pycon

   >>> counts.alphabet
   'ACGT'

This property is read-only; modifying the underlying ``_alphabet``
attribute may lead to unexpected results. Elements can be accessed both
by letter and by integer index:

.. code:: pycon

   >>> counts["C"] = -3
   >>> counts[2] = 7
   >>> print(counts)
   A  0.0
   C -3.0
   G  7.0
   T  0.0
   <BLANKLINE>
   >>> counts[1]
   -3.0

Using a letter that is not in the alphabet, or an index that is out of
bounds, will cause a ``IndexError``:

.. code:: pycon

   >>> counts["U"]
   Traceback (most recent call last):
       ...
   IndexError: 'U'
   >>> counts["X"] = 6
   Traceback (most recent call last):
       ...
   IndexError: 'X'
   >>> counts[7]
   Traceback (most recent call last):
       ...
   IndexError: index 7 is out of bounds for axis 0 with size 4

A two-dimensional ``Array`` can be created by specifying ``dims=2``:

.. code:: pycon

   >>> from Bio.Align.substitution_matrices import Array
   >>> counts = Array("ACGT", dims=2)
   >>> print(counts)
       A   C   G   T
   A 0.0 0.0 0.0 0.0
   C 0.0 0.0 0.0 0.0
   G 0.0 0.0 0.0 0.0
   T 0.0 0.0 0.0 0.0
   <BLANKLINE>

Again, both letters and integers can be used for indexing, and
specifying a letter that is not in the alphabet will cause an
``IndexError``:

.. code:: pycon

   >>> counts["A", "C"] = 12.0
   >>> counts[2, 1] = 5.0
   >>> counts[3, "T"] = -2
   >>> print(counts)
       A    C   G    T
   A 0.0 12.0 0.0  0.0
   C 0.0  0.0 0.0  0.0
   G 0.0  5.0 0.0  0.0
   T 0.0  0.0 0.0 -2.0
   <BLANKLINE>
   >>> counts["X", 1]
   Traceback (most recent call last):
       ...
   IndexError: 'X'
   >>> counts["A", 5]
   Traceback (most recent call last):
       ...
   IndexError: index 5 is out of bounds for axis 1 with size 4

Selecting a row or column from the two-dimensional array will return a
one-dimensional ``Array``:

.. code:: pycon

   >>> counts = Array("ACGT", dims=2)
   >>> counts["A", "C"] = 12.0
   >>> counts[2, 1] = 5.0
   >>> counts[3, "T"] = -2

.. code:: pycon

   >>> counts["G"]
   Array([0., 5., 0., 0.],
         alphabet='ACGT')
   >>> counts[:, "C"]
   Array([12.,  0.,  5.,  0.],
         alphabet='ACGT')

``Array`` objects can thus be used as an array and as a dictionary. They
can be converted to plain numpy arrays or plain dictionary objects:

.. code:: pycon

   >>> import numpy as np
   >>> x = Array("ACGT")
   >>> x["C"] = 5

.. code:: pycon

   >>> x
   Array([0., 5., 0., 0.],
         alphabet='ACGT')
   >>> a = np.array(x)  # create a plain numpy array
   >>> a
   array([0., 5., 0., 0.])
   >>> d = dict(x)  # create a plain dictionary
   >>> d
   {'A': 0.0, 'C': 5.0, 'G': 0.0, 'T': 0.0}

While the alphabet of an ``Array`` is usually a string, you may also use
a tuple of (immutable) objects. This is used for example for a `codon
substitution matrix <#codonmatrix>`__, where the keys are not individual
nucleotides or amino acids but instead three-nucleotide codons.

While the ``alphabet`` property of an ``Array`` is immutable, you can
create a new ``Array`` object by selecting the letters you are
interested in from the alphabet. For example,

.. code:: pycon

   >>> a = Array("ABCD", dims=2, data=np.arange(16).reshape(4, 4))
   >>> print(a)
        A    B    C    D
   A  0.0  1.0  2.0  3.0
   B  4.0  5.0  6.0  7.0
   C  8.0  9.0 10.0 11.0
   D 12.0 13.0 14.0 15.0
   <BLANKLINE>
   >>> b = a.select("CAD")
   >>> print(b)
        C    A    D
   C 10.0  8.0 11.0
   A  2.0  0.0  3.0
   D 14.0 12.0 15.0
   <BLANKLINE>

Note that this also allows you to reorder the alphabet.

Data for letters that are not found in the alphabet are set to zero:

.. code:: pycon

   >>> c = a.select("DEC")
   >>> print(c)
        D   E    C
   D 15.0 0.0 14.0
   E  0.0 0.0  0.0
   C 11.0 0.0 10.0
   <BLANKLINE>

Calculating a substitution matrix from a pairwise sequence alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As ``Array`` is a subclass of a numpy array, you can apply mathematical
operations on an ``Array`` object in much the same way. Here, we
illustrate this by calculating a scoring matrix from the alignment of
the 16S ribosomal RNA gene sequences of *Escherichia coli* and *Bacillus
subtilis*. First, we create a ``PairwiseAligner`` and initialize it with
the default scores used by ``blastn``:

.. code:: pycon

   >>> from Bio.Align import PairwiseAligner
   >>> aligner = PairwiseAligner()
   >>> aligner.mode = "local"
   >>> aligner.match_score = 2
   >>> aligner.mismatch_score = -3
   >>> aligner.open_gap_score = -7
   >>> aligner.extend_gap_score = -2

Next, we read in the 16S ribosomal RNA gene sequence of *Escherichia
coli* and *Bacillus subtilis* (provided in
``Tests/scoring_matrices/ecoli.fa`` and
``Tests/scoring_matrices/bsubtilis.fa``), and align them to each other:

.. code:: pycon

   >>> from Bio import SeqIO
   >>> sequence1 = SeqIO.read("ecoli.fa", "fasta")
   >>> sequence2 = SeqIO.read("bsubtilis.fa", "fasta")
   >>> alignments = aligner.align(sequence1.seq, sequence2.seq)

The number of alignments generated is very large:

.. code:: pycon

   >>> len(alignments)
   1990656

However, as they only differ trivially from each other, we arbitrarily
choose the first alignment, and count the number of each substitution:

.. code:: pycon

   >>> alignment = alignments[0]
   >>> from Bio.Align.substitution_matrices import Array
   >>> frequency = Array("ACGT", dims=2)
   >>> for (start1, end1), (start2, end2) in zip(*alignment.aligned):
   ...     seq1 = sequence1[start1:end1]
   ...     seq2 = sequence2[start2:end2]
   ...     for c1, c2 in zip(seq1, seq2):
   ...         frequency[c1, c2] += 1
   ...
   >>> print(frequency)
         A     C     G     T
   A 307.0  19.0  34.0  19.0
   C  15.0 280.0  25.0  29.0
   G  34.0  24.0 401.0  20.0
   T  24.0  36.0  20.0 228.0
   <BLANKLINE>

We normalize against the total number to find the probability of each
substitution, and create a symmetric matrix:

.. code:: pycon

   >>> import numpy as np
   >>> probabilities = frequency / np.sum(frequency)
   >>> probabilities = (probabilities + probabilities.transpose()) / 2.0
   >>> print(probabilities.format("%.4f"))
          A      C      G      T
   A 0.2026 0.0112 0.0224 0.0142
   C 0.0112 0.1848 0.0162 0.0215
   G 0.0224 0.0162 0.2647 0.0132
   T 0.0142 0.0215 0.0132 0.1505
   <BLANKLINE>

The background probability is the probability of finding an A, C, G, or
T nucleotide in each sequence separately. This can be calculated as the
sum of each row or column:

.. code:: pycon

   >>> background = np.sum(probabilities, 0)
   >>> print(background.format("%.4f"))
   A 0.2505
   C 0.2337
   G 0.3165
   T 0.1993
   <BLANKLINE>

The number of substitutions expected at random is simply the product of
the background distribution with itself:

.. code:: pycon

   >>> expected = np.dot(background[:, None], background[None, :])
   >>> print(expected.format("%.4f"))
          A      C      G      T
   A 0.0627 0.0585 0.0793 0.0499
   C 0.0585 0.0546 0.0740 0.0466
   G 0.0793 0.0740 0.1002 0.0631
   T 0.0499 0.0466 0.0631 0.0397
   <BLANKLINE>

The scoring matrix can then be calculated as the logarithm of the
odds-ratio of the observed and the expected probabilities:

.. code:: pycon

   >>> oddsratios = probabilities / expected
   >>> scoring_matrix = np.log2(oddsratios)
   >>> print(scoring_matrix)
        A    C    G    T
   A  1.7 -2.4 -1.8 -1.8
   C -2.4  1.8 -2.2 -1.1
   G -1.8 -2.2  1.4 -2.3
   T -1.8 -1.1 -2.3  1.9
   <BLANKLINE>

The matrix can be used to set the substitution matrix for the pairwise
aligner:

.. code:: pycon

   >>> aligner.substitution_matrix = scoring_matrix

A ``ValueError`` is triggered if the ``Array`` objects appearing in a
mathematical operation have different alphabets:

.. code:: pycon

   >>> from Bio.Align.substitution_matrices import Array
   >>> d = Array("ACGT")
   >>> r = Array("ACGU")
   >>> d + r
   Traceback (most recent call last):
       ...
   ValueError: alphabets are inconsistent

Reading ``Array`` objects from file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Bio.Align.substitution_matrices`` includes a parser to read one- and
two-dimensional ``Array`` objects from file. One-dimensional arrays are
represented by a simple two-column format, with the first column
containing the key and the second column the corresponding value. For
example, the file ``hg38.chrom.sizes`` (obtained from UCSC), available
in the ``Tests/Align`` subdirectory of the Biopython distribution,
contains the size in nucleotides of each chromosome in human genome
assembly hg38:

.. code:: text

   chr1    248956422
   chr2    242193529
   chr3    198295559
   chr4    190214555
   ...
   chrUn_KI270385v1    990
   chrUn_KI270423v1    981
   chrUn_KI270392v1    971
   chrUn_KI270394v1    970

To parse this file, use

.. code:: pycon

   >>> from Bio.Align import substitution_matrices
   >>> with open("hg38.chrom.sizes") as handle:
   ...     table = substitution_matrices.read(handle)
   ...
   >>> print(table)  # doctest: +ELLIPSIS
   chr1 248956422.0
   chr2 242193529.0
   chr3 198295559.0
   chr4 190214555.0
   ...
   chrUn_KI270423v1       981.0
   chrUn_KI270392v1       971.0
   chrUn_KI270394v1       970.0
   <BLANKLINE>

Use ``dtype=int`` to read the values as integers:

.. code:: pycon

   >>> with open("hg38.chrom.sizes") as handle:
   ...     table = substitution_matrices.read(handle, int)
   ...
   >>> print(table)  # doctest: +ELLIPSIS
   chr1 248956422
   chr2 242193529
   chr3 198295559
   chr4 190214555
   ...
   chrUn_KI270423v1       981
   chrUn_KI270392v1       971
   chrUn_KI270394v1       970
   <BLANKLINE>

For two-dimensional arrays, we follow the file format of substitution
matrices provided by NCBI. For example, the BLOSUM62 matrix, which is
the default substitution matrix for NCBI’s protein-protein BLAST
:raw-latex:`\cite{altschul1990}` program ``blastp``, is stored as
follows:

.. code:: text

   #  Matrix made by matblas from blosum62.iij
   #  * column uses minimum score
   #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
   #  Blocks Database = /data/blocks_5.0/blocks.dat
   #  Cluster Percentage: >= 62
   #  Entropy =   0.6979, Expected =  -0.5209
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
   ...

This file is included in the Biopython distribution under
``Bio/Align/substitution_matrices/data``. To parse this file, use

.. code:: pycon

   >>> from Bio.Align import substitution_matrices
   >>> with open("BLOSUM62") as handle:
   ...     matrix = substitution_matrices.read(handle)
   ...
   >>> print(matrix.alphabet)
   ARNDCQEGHILKMFPSTWYVBZX*
   >>> print(matrix["A", "D"])
   -2.0

The header lines starting with ``#`` are stored in the attribute
``header``:

.. code:: pycon

   >>> matrix.header[0]
   'Matrix made by matblas from blosum62.iij'

We can now use this matrix as the substitution matrix on an aligner
object:

.. code:: pycon

   >>> from Bio.Align import PairwiseAligner
   >>> aligner = PairwiseAligner()
   >>> aligner.substitution_matrix = matrix

To save an Array object, create a string first:

.. code:: pycon

   >>> text = str(matrix)
   >>> print(text)  # doctest: +ELLIPSIS
   #  Matrix made by matblas from blosum62.iij
   #  * column uses minimum score
   #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
   #  Blocks Database = /data/blocks_5.0/blocks.dat
   #  Cluster Percentage: >= 62
   #  Entropy =   0.6979, Expected =  -0.5209
        A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S ...
   A  4.0 -1.0 -2.0 -2.0  0.0 -1.0 -1.0  0.0 -2.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0  1.0 ...
   R -1.0  5.0  0.0 -2.0 -3.0  1.0  0.0 -2.0  0.0 -3.0 -2.0  2.0 -1.0 -3.0 -2.0 -1.0 ...
   N -2.0  0.0  6.0  1.0 -3.0  0.0  0.0  0.0  1.0 -3.0 -3.0  0.0 -2.0 -3.0 -2.0  1.0 ...
   D -2.0 -2.0  1.0  6.0 -3.0  0.0  2.0 -1.0 -1.0 -3.0 -4.0 -1.0 -3.0 -3.0 -1.0  0.0 ...
   C  0.0 -3.0 -3.0 -3.0  9.0 -3.0 -4.0 -3.0 -3.0 -1.0 -1.0 -3.0 -1.0 -2.0 -3.0 -1.0 ...
   ...

and write the ``text`` to a file.

Loading predefined substitution matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Biopython contains a large set of substitution matrices defined in the
literature, including BLOSUM (Blocks Substitution Matrix)
:raw-latex:`\cite{henikoff1992}` and PAM (Point Accepted Mutation)
matrices :raw-latex:`\cite{dayhoff1978}`. These matrices are available
as flat files in the ``Bio/Align/scoring_matrices/data`` directory, and
can be loaded into Python using the ``load`` function in the
``scoring_matrices`` submodule. For example, the BLOSUM62 matrix can be
loaded by running

.. code:: pycon

   >>> from Bio.Align import substitution_matrices
   >>> m = substitution_matrices.load("BLOSUM62")

This substitution matrix has an alphabet consisting of the 20 amino
acids used in the genetic code, the three ambiguous amino acids B
(asparagine or aspartic acid), Z (glutamine or glutamic acid), and X
(representing any amino acid), and the stop codon represented by an
asterisk:

.. code:: pycon

   >>> m.alphabet
   'ARNDCQEGHILKMFPSTWYVBZX*'

To get a full list of available substitution matrices, use ``load``
without an argument:

.. code:: pycon

   >>> substitution_matrices.load()  # doctest: +ELLIPSIS
   ['BENNER22', 'BENNER6', 'BENNER74', 'BLASTN', 'BLASTP', 'BLOSUM45', 'BLOSUM50', ..., 'TRANS']

.. container::
   :name: codonmatrix

   Note that the substitution matrix provided by Schneider *et al.*
   :raw-latex:`\cite{schneider2005}` uses an alphabet consisting of
   three-nucleotide codons:

.. code:: pycon

   >>> m = substitution_matrices.load("SCHNEIDER")
   >>> m.alphabet  # doctest: +ELLIPSIS
   ('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', ..., 'TTG', 'TTT')

.. _`sec:pairwise2`:

Pairwise alignments using pairwise2
-----------------------------------

**Please note that Bio.pairwise2 was deprecated in Release 1.80.** As an
alternative, please consider using ``Bio.Align.PairwiseAligner``
(described in section :ref:`sec:pairwise`).

``Bio.pairwise2`` contains essentially the same algorithms as ``water``
(local) and ``needle`` (global) from the
`EMBOSS <http://emboss.sourceforge.net/>`__ suite (see above) and should
return the same results. The ``pairwise2`` module has undergone some
optimization regarding speed and memory consumption recently (Biopython
versions >1.67) so that for short sequences (global alignments: ~2000
residues, local alignments ~600 residues) it’s faster (or equally fast)
to use ``pairwise2`` than calling EMBOSS’ ``water`` or ``needle`` via
the command line tools.

Suppose you want to do a global pairwise alignment between the same two
hemoglobin sequences from above (``HBA_HUMAN``, ``HBB_HUMAN``) stored in
``alpha.faa`` and ``beta.faa``:

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio import SeqIO
   >>> seq1 = SeqIO.read("alpha.faa", "fasta")
   >>> seq2 = SeqIO.read("beta.faa", "fasta")
   >>> alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)

As you see, we call the alignment function with ``align.globalxx``. The
tricky part are the last two letters of the function name (here:
``xx``), which are used for decoding the scores and penalties for
matches (and mismatches) and gaps. The first letter decodes the match
score, e.g. ``x`` means that a match counts 1 while mismatches have no
costs. With ``m`` general values for either matches or mismatches can be
defined (for more options see `Biopython’s
API <http://biopython.org/docs/1.77/api/Bio.pairwise2.html>`__). The
second letter decodes the cost for gaps; ``x`` means no gap costs at
all, with ``s`` different penalties for opening and extending a gap can
be assigned. So, ``globalxx`` means that only matches between both
sequences are counted.

Our variable ``alignments`` now contains a list of alignments (at least
one) which have the same optimal score for the given conditions. In our
example this are 80 different alignments with the score 72
(``Bio.pairwise2`` will return up to 1000 alignments). Have a look at
one of these alignments:

.. code:: pycon

   >>> len(alignments)
   80
   >>> print(alignments[0])  # doctest:+ELLIPSIS
   Alignment(seqA='MV-LSPADKTNV---K-A--A-WGKVGAHAG...YR-', seqB='MVHL-----T--PEEKSAVTALWGKV----...Y-H', score=72.0, start=0, end=217)

Each alignment is a named tuple consisting of the two aligned sequences,
the score, the start and the end positions of the alignment (in global
alignments the start is always 0 and the end the length of the
alignment). ``Bio.pairwise2`` has a function ``format_alignment`` for a
nicer printout:

.. code:: pycon

   >>> print(pairwise2.format_alignment(*alignments[0]))  # doctest:+ELLIPSIS
   MV-LSPADKTNV---K-A--A-WGKVGAHAG---EY-GA-EALE-RMFLSF----PTTK-TY--F...YR-
   || |     |     | |  | ||||        |  |  |||  |  |      |    |   |...|  
   MVHL-----T--PEEKSAVTALWGKV-----NVDE-VG-GEAL-GR--L--LVVYP---WT-QRF...Y-H
     Score=72
   <BLANKLINE>

Since Biopython 1.77 the required parameters can be supplied with
keywords. The last example can now also be written as:

.. code:: pycon

   >>> alignments = pairwise2.align.globalxx(sequenceA=seq1.seq, sequenceB=seq2.seq)

Better alignments are usually obtained by penalizing gaps: higher costs
for opening a gap and lower costs for extending an existing gap. For
amino acid sequences match scores are usually encoded in matrices like
``PAM`` or ``BLOSUM``. Thus, a more meaningful alignment for our example
can be obtained by using the BLOSUM62 matrix, together with a gap open
penalty of 10 and a gap extension penalty of 0.5 (using ``globalds``):

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio import SeqIO
   >>> from Bio.Align import substitution_matrices
   >>> blosum62 = substitution_matrices.load("BLOSUM62")
   >>> seq1 = SeqIO.read("alpha.faa", "fasta")
   >>> seq2 = SeqIO.read("beta.faa", "fasta")
   >>> alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
   >>> len(alignments)
   2
   >>> print(pairwise2.format_alignment(*alignments[0]))
   MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY...KYR
   || |.|..|..|.|.|||| ......|............|.......||.
   MVHLTPEEKSAVTALWGKV-NVDEVGGEALGRLLVVYPWTQRFF...KYH
     Score=292.5

This alignment has the same score that we obtained earlier with EMBOSS
needle using the same sequences and the same parameters.

Local alignments are called similarly with the function
``align.localXX``, where again XX stands for a two letter code for the
match and gap functions:

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio.Align import substitution_matrices
   >>> blosum62 = substitution_matrices.load("BLOSUM62")
   >>> alignments = pairwise2.align.localds("LSPADKTNVKAA", "PEEKSAV", blosum62, -10, -1)
   >>> print(pairwise2.format_alignment(*alignments[0]))
   3 PADKTNV
     |..|..|
   1 PEEKSAV
     Score=16
   <BLANKLINE>

In recent Biopython versions, ``format_alignment`` will only print the
aligned part of a local alignment (together with the start positions in
1-based notation, as shown in the above example). If you are also
interested in the non- aligned parts of the sequences, use the
keyword-parameter ``full_sequences=True``:

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio.Align import substitution_matrices
   >>> blosum62 = substitution_matrices.load("BLOSUM62")
   >>> alignments = pairwise2.align.localds("LSPADKTNVKAA", "PEEKSAV", blosum62, -10, -1)
   >>> print(pairwise2.format_alignment(*alignments[0], full_sequences=True))
   LSPADKTNVKAA
     |..|..|   
   --PEEKSAV---
     Score=16
   <BLANKLINE>

Note that local alignments must, as defined by Smith & Waterman, have a
positive score (>0). Thus, ``pairwise2`` may return no alignments if no
score >0 has been obtained. Also, ``pairwise2`` will not report
alignments which are the result of the addition of zero-scoring
extensions on either site. In the next example, the pairs
serine/aspartic acid (S/D) and lysine/asparagine (K/N) both have a match
score of 0. As you see, the aligned part has not been extended:

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio.Align import substitution_matrices
   >>> blosum62 = substitution_matrices.load("BLOSUM62")
   >>> alignments = pairwise2.align.localds("LSSPADKTNVKKAA", "DDPEEKSAVNN", blosum62, -10, -1)
   >>> print(pairwise2.format_alignment(*alignments[0]))
   4 PADKTNV
     |..|..|
   3 PEEKSAV
     Score=16
   <BLANKLINE>

Instead of supplying a complete match/mismatch matrix, the match code
``m`` allows for easy defining general match/mismatch values. The next
example uses match/mismatch scores of 5/-4 and gap penalties
(open/extend) of 2/0.5 using ``localms``:

.. code:: pycon

   >>> alignments = pairwise2.align.localms("AGAACT", "GAC", 5, -4, -2, -0.5)
   >>> print(pairwise2.format_alignment(*alignments[0]))
   2 GAAC
     | ||
   1 G-AC
     Score=13
   <BLANKLINE>

One useful keyword argument of the ``Bio.pairwise2.align`` functions is
``score_only``. When set to ``True`` it will only return the score of
the best alignment(s), but in a significantly shorter time. It will also
allow the alignment of longer sequences before a memory error is raised.
Another useful keyword argument is ``one_alignment_only=True`` which
will also result in some speed gain.

Unfortunately, ``Bio.pairwise2`` does not work with Biopython’s multiple
sequence alignment objects (yet). However, the module has some
interesting advanced features: you can define your own match and gap
functions (interested in testing affine logarithmic gap costs?), gap
penalties and end gaps penalties can be different for both sequences,
sequences can be supplied as lists (useful if you have residues that are
encoded by more than one character), etc. These features are hard (if at
all) to realize with other alignment tools. For more details see the
modules documentation in `Biopython’s
API <http://biopython.org/docs/\bpversion/api/Bio.pairwise2.html>`__.
