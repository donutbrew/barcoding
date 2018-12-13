barcoding
=========

This is a demultiplexing script for Nanopore-generated fastq files. It was forked from ONT, and I added parallel processing and several other options

*Instructions:*
```
split_barcodes.pl --barcodes ONT_NativeBarcodes.fasta --threads 12 --stringency 6  allreads.fastq

    --barcodes            FASTA file of Barcodes
    --threads             Threads to use
    --stringency          Edit distance to tolerate (default:6)
    --verbose             Extra information on the reads
    --require_both        Require both ends to have a barcode
    --enforce_orientation Require each barcode to be the expected orientation
    --check_hybrid        Check if there are additional barcodes in middle of read
                          (throws out if this is a hybrid read) 
```

**Currently, the script expects the fastq file (or symlink) to be in the working directory.**

This repository is part of the ONT Barcoding Protocol for amplicons.

Data
----
The dataset (in t/data/) was prepared with two barcodes for illustrative
purposes.

Scripts
-------
Scripts can be found in the bin/ subdirectory.
Note: Run the script (split_barcodes.pl) as-is for now. Nothing else has been modified. 

Prerequisites
-------------
The main prerequisites are listed below, with a complete list available in Build.PL:
 - Bio::Perl
 - Text::LevenshteinXS
 - Readonly

All should be available from CPAN (UNIX/Strawberry systems: "cpan <modulename>") or PPD (for ActivePerl)

