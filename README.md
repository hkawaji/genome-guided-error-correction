# Genome-guided error correctoion 
for error prone sequencing of cDNA amplicon

Usage
---

To perform all steps of error correction, please use 'correct' subcommand as:

```
./ggct.sh correct -i INFILE.fq  -g GENOME [-p PASS(default:2, or 1)] [-t THREADS(default:10)] [options used in individual steps]
```

Each of the individual steps can be run by other subcommands,
'ggc' or 'ecc', as:

```
# genome guided clustering (via BLAT)
./ggct.sh ggc -i INFILE -g GENOME [-t THREADS(default:10)] [-l SEQLEN_MIN(default:1000, should not be changed for Canu error correction)]  [-e error_correction_seqnum_min(default:5)]

# error correction within each of the clusters (with Canu)
./ggct.sh ecc -i INFILE -c CLUSTER_FILE('ggc' output) [-t THREADS(default:10)]
```


Note that sequence IDs should not include comma (',').


How it works
------------
Reads are grouped by isoform based on their alignments with
the specified genome sequence (by minimap2), and the orphan
reads (that does not grouped with others) are included
to the grouop that has the maximum length of exons in total.
Sequence error correction is performed within each of the
groups.

For two pass mode, the corrected sequences (in the first pass)
are used for grouping, and sequence error correction itself is
performed on the uncorrected reads.

Note that Canu does not perform error correction on the sequences shorter
than 1000bp (in default), and they are discarded at the beginning.



Requirements
------------
Software below is required to run "./ggct.sh". Their versions used for
development are indicated in '( )'

* minimap2 (version 2.12)
* setqk (version 1.3)
* bedtools [bed12ToBed6, mergeBed, groupBy] (version 2.27.1)
* samtools [bed12ToBed6, mergeBed, groupBy] (version 2.27.1)
* Canu (version 1.7.1)
* R (with package 'stringr')
* GNU awk (ver.4.1.3)



Author
------
This software is written by Hideya Kawaji.

Copyright (c) 2018-2020 RIKEN & Tokyo Metropolitan Institute of Medical Science

Distributed under Apache License 2.0.
