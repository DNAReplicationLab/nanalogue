#!/usr/bin/env bash
set -euo pipefail

# Compare the vendored HTSlib-related crates in this repository against
# upstream source checkouts.
#
# Usage:
#   ./vendor/compare.sh [--patch] [--output-dir DIR] [rust-htslib-ref] [htslib-ref]
#
# Examples:
#   ./vendor/compare.sh
#   ./vendor/compare.sh v1.0.0 v1.19
#   ./vendor/compare.sh --patch --output-dir /tmp/vendor-diffs v1.0.0 v1.19
#
# By default, the script prints a compact summary to stdout.
# Use --patch to write full diff files into an output directory.

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

write_patches=0
output_dir=""
args=()

while (($# > 0)); do
    case "$1" in
        --patch)
            write_patches=1
            shift
            ;;
        --output-dir)
            output_dir="${2:-}"
            if [[ -z "$output_dir" ]]; then
                printf 'error: --output-dir requires a directory argument\n' >&2
                exit 2
            fi
            shift 2
            ;;
        --help|-h)
            sed -n '1,28p' "$0"
            exit 0
            ;;
        --)
            shift
            while (($# > 0)); do
                args+=("$1")
                shift
            done
            ;;
        *)
            args+=("$1")
            shift
            ;;
    esac
done

set -- "${args[@]}"

if ((write_patches)) && [[ -z "$output_dir" ]]; then
    output_dir="$repo_root/vendor-diffs"
fi

clone_and_checkout() {
    local repo_url="$1"
    local ref="$2"
    local dest="$3"

    git clone --quiet "$repo_url" "$dest"
    git -C "$dest" checkout --quiet "$ref"
}

print_summary() {
    local upstream_dir="$1"
    local vendored_dir="$2"
    local label="$3"

    printf '\n=== %s ===\n' "$label"
    git diff --no-index --stat -- "$upstream_dir" "$vendored_dir" || true
}

write_patch() {
    local upstream_dir="$1"
    local vendored_dir="$2"
    local label="$3"
    local outfile="$4"

    printf '\n=== %s ===\n' "$label"
    git diff --no-index --binary -- "$upstream_dir" "$vendored_dir" > "$outfile" || true
    printf 'wrote %s\n' "$outfile"
}

rust_htslib_ref="${1:-v1.0.0}"
htslib_ref="${2:-}"

rust_htslib_upstream="$tmpdir/rust-htslib"
clone_and_checkout "https://github.com/rust-bio/rust-htslib.git" "$rust_htslib_ref" "$rust_htslib_upstream"

if ((write_patches)); then
    mkdir -p "$output_dir"
    write_patch "$rust_htslib_upstream" "$repo_root/vendor/rust-htslib-1.0.0" "rust-htslib vs $rust_htslib_ref" "$output_dir/rust-htslib.patch"
else
    print_summary "$rust_htslib_upstream" "$repo_root/vendor/rust-htslib-1.0.0" "rust-htslib vs $rust_htslib_ref"
fi

if [[ -n "$htslib_ref" ]]; then
    htslib_upstream="$tmpdir/htslib"
    clone_and_checkout "https://github.com/samtools/htslib.git" "$htslib_ref" "$htslib_upstream"

    if ((write_patches)); then
        write_patch "$htslib_upstream" "$repo_root/vendor/hts-sys-2.2.0/htslib" "htslib vs $htslib_ref" "$output_dir/htslib.patch"
    else
        print_summary "$htslib_upstream" "$repo_root/vendor/hts-sys-2.2.0/htslib" "htslib vs $htslib_ref"
    fi
else
    printf '\nSkipping HTSlib comparison; pass a second argument to compare the vendored htslib tree.\n'
fi
