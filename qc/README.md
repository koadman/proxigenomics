Proxigenomics QC scripts
=============

This is a collection of very basic QC scripts for metagenomic 3C data

* count_cutsites.py -- counts the number of times the given cut site sequence appears in a read set. Potentially useful for evaluating whether the restriction digest was effective by comparing 3C library sequences to shotgun library sequences

* check_cuts.py -- counts how frequently reads that are split-mapped have a restriction recognition site near the split-point of the mapping. In 3C-seq data the reads containing proximity ligations should have a high fraction of cut sites at the split points.

