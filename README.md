# `nanalogue`

Package to primarily process single-molecule base modification data
stored in the modBAM format.

# Commands

### `nanalogue read-table-no-mods`

Given a BAM file and a sequencing summary file, or just the BAM file,
this command outputs a table where each row contains the basecalled
length and alignment length per molecule.

### `nanalogue read-stats`

Given a BAM file, prints read and alignment length statistics such as N50,
and how many alignments are primary, secondary etc. For more statistics,
please consider using `samtools stats`.

### Other commands

Please run `nanalogue -h` to see what other commands are available.

# Acknowledgments

This software was developed at the Earlham Institute in the UK.
This work was supported by the Biotechnology and Biological Sciences
Research Council (BBSRC), part of UK Research and Innovation,
through the Core Capability Grant BB/CCG2220/1 at the Earlham Institute
and the Earlham Institute Strategic Programme Grant Cellular Genomics
BBX011070/1 and its constituent work packages BBS/E/ER/230001B 
(CellGen WP2 Consequences of somatic genome variation on traits).
The work was also supported by the following response-mode project grants:
BB/W006014/1 (Single molecule detection of DNA replication errors) and
BB/Y00549X/1 (Single molecule analysis of Human DNA replication).
This research was supported in part by NBI Research Computing
through use of the High-Performance Computing system and Isilon storage.
