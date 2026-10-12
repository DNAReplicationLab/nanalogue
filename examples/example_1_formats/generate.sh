#!/usr/bin/env bash
# Regenerates the committed `example_1` format variants in this directory.
#
# The outputs are goldens: tests read them but never regenerate them. Rerun
# this only when deliberately changing the variant set, from the repository
# root, with samtools on PATH (they were made with samtools 1.19.2 / HTSlib
# 1.19.1, the HTSlib version bundled with the pinned `hts-sys`):
#
#     bash examples/example_1_formats/generate.sh
#
# Every variant holds the same four records as `examples/example_1.bam`,
# coordinate-sorted so that each one can be indexed (the original file is
# unsorted). `cram_v3.1_unsorted.cram` keeps the original record order.
#
# External-reference CRAMs record the relative `UR:examples/contigs.fa`, so
# they decode only when run from the repository root or when the reference is
# found through `REF_PATH`/`REF_CACHE` using the header's M5 checksums.
#
# `tests/example_1_format_variants.rs` checks every file here and pins two
# HTSlib 1.19 behaviours: CRAM 2.x reference-free records do not decode absent
# base qualities back to `*`, and multi-reference slices cannot detect a wrong
# external reference because their reference MD5 is zero by specification.
set -euo pipefail

if [[ ! -f examples/example_1.bam ]]; then
    echo "run this script from the repository root" >&2
    exit 1
fi

out=examples/example_1_formats
ref=examples/contigs.fa
# Never consult the EBI reference server; resolve references only through
# the `UR` tags in the input headers below.
export REF_PATH=/nonexistent REF_CACHE=/nonexistent

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

rm -f "$out"/*.bam "$out"/*.sam "$out"/*.gz "$out"/*.cram \
    "$out"/*.bai "$out"/*.csi "$out"/*.crai
samtools faidx "$ref"

samtools sort --no-PG -o "$tmp/sorted.bam" examples/example_1.bam
samtools view -h --no-PG -o "$tmp/sorted.sam" "$tmp/sorted.bam"

# Annotates @SQ lines with the reference checksums and the relative reference
# path, as an aligner would. Without -T, HTSlib keeps this relative `UR`.
add_ref_tags() {
    awk -v OFS='\t' '
        NR == FNR {
            if ($1 == "@SQ") {
                split($2, name, ":")
                for (i = 3; i <= NF; i++) if ($i ~ /^M5:/) m5[name[2]] = $i
            }
            next
        }
        $1 == "@SQ" {
            split($2, name, ":")
            print $0, m5[name[2]], "UR:examples/contigs.fa"
            next
        }
        { print }
    ' <(samtools dict "$ref") "$1"
}
add_ref_tags "$tmp/sorted.sam" >"$tmp/with_ref.sam"
samtools view -h --no-PG -o "$tmp/unsorted.sam" examples/example_1.bam
add_ref_tags "$tmp/unsorted.sam" >"$tmp/unsorted_with_ref.sam"

# Converts to CRAM. Usage: cram <name> <input SAM> [output-fmt-option ...]
cram_no_index() {
    local name=$1 input=$2
    shift 2
    local options=()
    for option in "$@"; do options+=(--output-fmt-option "$option"); done
    samtools view --no-PG -C "${options[@]}" -o "$out/$name.cram" "$input"
}
cram() {
    cram_no_index "$@"
    samtools index "$out/$1.cram"
}

# BAM and SAM containers.
samtools view --no-PG -b -o "$out/bam_csi.bam" "$tmp/sorted.sam"
samtools index -c "$out/bam_csi.bam"
samtools view --no-PG -u -o "$out/bam_uncompressed.bam" "$tmp/sorted.sam"
samtools index -b "$out/bam_uncompressed.bam"
cp "$tmp/sorted.sam" "$out/sam_plain.sam"
samtools view --no-PG -h --output-fmt sam,level=6 -o "$out/sam_bgzf.sam.gz" "$tmp/sorted.sam"
samtools index -c "$out/sam_bgzf.sam.gz"

# CRAM versions crossed with reference modes. CRAM 1.0 cannot be written by
# HTSlib 1.19 and 4.0 is a draft, so both are omitted. HTSlib writes CRAM 2.0
# exactly as 2.1 apart from the minor version byte, so the 2.0 files only check
# that the version byte is accepted; real legacy 2.0 files may differ.
for version in 2.0 2.1 3.0 3.1; do
    cram "cram_v${version}_external_ref" "$tmp/with_ref.sam" \
        "version=$version" embed_ref=0 no_ref=0
    cram "cram_v${version}_embedded_ref" "$tmp/with_ref.sam" \
        "version=$version" embed_ref=1 no_ref=0
    # Like an ordinary reference-less conversion, this header has no M5/UR.
    cram "cram_v${version}_no_ref" "$tmp/sorted.sam" \
        "version=$version" embed_ref=0 no_ref=1
done

# Consensus embedding (embed_ref=2): HTSlib embeds a reference built from the
# reads instead of one read from a FASTA. The mapped reads match the reference,
# so the result resembles `cram_v3.1_embedded_ref`. HTSlib warns "Failed to
# build reference" for the container holding only the unmapped read, which
# needs no reference.
cram cram_v3.1_consensus_ref "$tmp/sorted.sam" version=3.1 embed_ref=2 \
    2> >(grep -v 'Failed to build reference' >&2)

# CRAM 3.1 compression profiles and codec choices. With only four records,
# HTSlib keeps many blocks raw; these settings change the codecs that remain
# (rANS-Nx16 by default, rANS 4x16 for fast, arithmetic coding for archive,
# gzip only without rANS, and raw blocks only at level 0).
cram cram_v3.1_profile_fast "$tmp/with_ref.sam" version=3.1 fast
cram cram_v3.1_profile_archive "$tmp/with_ref.sam" version=3.1 archive
cram cram_v3.1_gzip_only "$tmp/with_ref.sam" version=3.1 use_rans=0
cram cram_v3.1_uncompressed "$tmp/with_ref.sam" version=3.1 level=0

# CRAM 3.1 container and slice layouts.
cram cram_v3.1_multi_ref_slice "$tmp/with_ref.sam" version=3.1 \
    multi_seq_per_slice=1
cram cram_v3.1_slice_per_record "$tmp/with_ref.sam" version=3.1 \
    seqs_per_slice=1 slices_per_container=4 multi_seq_per_slice=1
cram cram_v3.1_container_per_record "$tmp/with_ref.sam" version=3.1 \
    seqs_per_slice=1 slices_per_container=1

# Original (unsorted) record order. Not indexed, like any unsorted file.
cram_no_index cram_v3.1_unsorted "$tmp/unsorted_with_ref.sam" version=3.1

# Interleaved files: each record is encoded separately with different
# settings, then `samtools cat` joins the per-record containers (CRAM) or
# BGZF blocks (BAM) without re-encoding them.
# Usage: split_records <input SAM> <output prefix>
split_records() {
    for i in 0 1 2 3; do
        { grep '^@' "$1"; grep -v '^@' "$1" | sed -n "$((i + 1))p"; } >"$2_$i.sam"
    done
}
split_records "$tmp/with_ref.sam" "$tmp/record"
split_records "$tmp/sorted.sam" "$tmp/plain_record"

# Records 0-2 are dummyI, dummyII and dummyIII; record 3 is unmapped.
for version in 2.1 3.0 3.1; do
    if [[ $version == 3.1 ]]; then second_codec=archive; else second_codec=level=9; fi
    cram_no_index "r0" "$tmp/record_0.sam" "version=$version" embed_ref=0 no_ref=0
    cram_no_index "r1" "$tmp/record_1.sam" "version=$version" embed_ref=1 no_ref=0 \
        "$second_codec"
    cram_no_index "r2" "$tmp/record_2.sam" "version=$version" embed_ref=0 no_ref=1 \
        use_rans=0
    cram_no_index "r3" "$tmp/record_3.sam" "version=$version" level=0
    samtools cat --no-PG -o "$out/cram_v${version}_interleaved.cram" \
        "$out"/r{0,1,2,3}.cram
    rm "$out"/r{0,1,2,3}.cram
    samtools index "$out/cram_v${version}_interleaved.cram"
done

# BAM with one BGZF block per record, alternating stored (level 0) and
# maximally compressed (level 9) deflate blocks.
for i in 0 1 2 3; do
    samtools view --no-PG --output-fmt "bam,level=$((i % 2 == 0 ? 0 : 9))" \
        -o "$tmp/record_$i.bam" "$tmp/plain_record_$i.sam"
done
samtools cat --no-PG -o "$out/bam_interleaved_levels.bam" "$tmp"/record_{0,1,2,3}.bam
samtools index -b "$out/bam_interleaved_levels.bam"
