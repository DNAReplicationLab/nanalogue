# Zig htslib exploration

## Preamble

This is not part of the nanalogue rust library and executable.
These files should be ignored by Cargo.
I am exploring how to use htslib directly with zig.
We may use these scripts as a benchmark to compare `nanalogue` against in the future.
These scripts may break or produce incorrect outputs as we are just exploring!

## Installation

### Download zig

Go to https://ziglang.org/download/ and download and untar zig.
I am using `zig-x86_64-linux-0.16.0.tar.xz` and I untared it at
`~/zig/zig-x86_64-linux-0.16.0/` in my computer. You have to
substitute this path with wherever you've put zig on your computer.

## Installing htslib

* Get our nanalogue repository, and do `git submodule update --init --recursive` to get the
htslib source files from their git repository.
* Then, from `zig_htslib_exploration/htslib`, run `git submodule update --init --recursive` 
to get submodules from within htslib.
* Then follow these steps from `zig_htslib_exploration/htslib`

```
mkdir ../htslib_compiled

# if you cloned htslib from git, generate ./configure first
autoreconf -i

# change the CC argument below to wherever you have untared the zig archive
./configure CC="$(realpath ~/zig/zig-x86_64-linux-0.16.0/zig) cc" --prefix=$(realpath ../htslib_compiled)

# not sure why I have to do the command below
sed -i '/HAVE_BUILTIN_CPU_SUPPORT_SSSE3/d;/HAVE_ATTRIBUTE_TARGET_SSSE3/d' config.h

# Make and make install. The `-j4` I think instructs make to use 4 threads,
# you can use as many threads as you want.
make CC="$(realpath ~/zig/zig-x86_64-linux-0.16.0/zig) cc" -j4
make install

```

## Building our program

Cd to `zig_htslib_exploration` and run the following command.

```
~/zig/zig-x86_64-linux-0.16.0/zig  build
```

## Running our programs

From `zig_htslib_exploration`, do

```
./zig-out/bin/reads_table <some_bam_file>
```

This should produce a tabular output of read ids and sequence lengths.

To convert a coordinate-sorted BAM to CRAM, do

```
./zig-out/bin/bam_to_cram <input.bam> <output.cram> <reference.fa> [threads]
```

For example, with the example files in `../examples/`:

```
./zig-out/bin/bam_to_cram ../examples/example_3.bam /tmp/example_3.cram ../examples/contigs.fa
```
