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
