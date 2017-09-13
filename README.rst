clip
----

Clips overlapping regions in read mates of SAM/BAM files.

To install: ``pip install clipoverlap``
Or::

    git clone git@github.com:innovate-invent/clip.git
    cd clip
    python3 setup.py install

To run::

    $ clip -h
    Clip Overlap v1.0
    Clip overlapping reads from SAM/BAM/CRAM file
    Use: clip [-tmabcosv] [input file path | < infile > outfile] [output file path]
    If no paths are given stdin and stdout are used.
    -t # Threads to use for processing (Default=1)
    -m # Maximum template length guaranteeing no read overlap (Default=1000)
    -a Alternate strand being clipped to avoid strand bias (RAM intensive)
    -b Trim tails of reads that extend past end of mate. Used to trim barcode remnants.
    -c Clip only, do not merge clipped region into mate.
    -o [sbuc] Output format: s=SAM (Default), b=BAM compressed, bu=BAM uncompressed, c=CRAM
    -s Maintain input order (High depth regions may fill RAM), if not set will output in arbitrary order (Minimal RAM)
    -v Verbose status output

You may notice if you just run ``clip`` with no parameters it will just sit there doing nothing.
That is because the default is to listen to stdin for input.

Notes:
------
``clip`` uses a minimum of two subprocesses regardless of the ``-t`` option.

``-a`` will alternate between clipping the tail of the left most strand and clipping the head of the right most strand.
This is to avoid possible strand bias later in a processing pipeline.

If you are processing reads that had barcodes ligated and removed the 5' barcode in a previous step (See `ProDuSe:trim <https://github.com/morinlab/ProDuSe>`_)
then use the ``-b`` option to remove any possible 3' barcode sequence that would be appended if sequencing ran to the end of the molecule.

Using ``-s`` and ``-a`` together will force ``clip`` to try and sort by start reference coordinate.
If unsorted data is the input then this could potentially run out of RAM.

Merge Algorithm
---------------
The mate read cigars are assumed to align 1-1 with an offset determined by the difference in the reference start positions.

* If ``-c`` is unset then ``clip`` will retain the highest quality base at a given position in the overlapping region of the mate pairs.
* If the base qualities are equal then it will keep the base that does not match the reference.
* If base qualities are equal and both bases are different variants, then the quality score is set to 3 (3 = Phred 50% probability of either base).
* If the operations between the aligned cigars do not match then the operations from the mate with the lowest alignment cost are retained.

The alignment cost is calculated for the overlapping region only.
The cost is summed using these values:

===========  ==========================
Operation    Value
===========  ==========================
M, X, =, N   -1
I            6 to start, +1 to lengthen
D            3 to start, +1 to lengthen
===========  ==========================

TODO:
-----
Significant speed and memory optimisations are planned.
Need to eliminate the pysam dependency first.
