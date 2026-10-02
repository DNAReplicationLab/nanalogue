#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [[ "$(uname -s):$(uname -m)" != Linux:x86_64 ]]; then
    echo "The baseline BAM viewer build requires Linux x86_64." >&2
    exit 1
fi

for command in cargo sha256sum; do
    if ! command -v "$command" >/dev/null; then
        echo "$command is required to build the baseline BAM viewer." >&2
        exit 1
    fi
done

rust_target="x86_64-unknown-linux-gnu"
if ! command -v rustup >/dev/null; then
    echo "rustup is required to verify the $rust_target Rust target." >&2
    exit 1
fi
if ! rustup target list --installed | grep -qx "$rust_target"; then
    echo "The $rust_target Rust target is required." >&2
    echo "Install it with: rustup target add $rust_target" >&2
    exit 1
fi

target_dir="$repo_root/target/portable-bam-viewer"
release_dir="$target_dir/$rust_target/release"
output_dir="$repo_root/target/portable-dist"
artifact="$output_dir/nanalogue_bam_viewer-x86_64-baseline"

rm -rf "$target_dir"
mkdir -p "$output_dir"
rm -f "$artifact"
(
    unset CARGO_ENCODED_RUSTFLAGS
    export CARGO_TARGET_DIR="$target_dir"
    export RUSTFLAGS="-C target-cpu=x86-64"
    cargo build --locked --release --target "$rust_target" \
        --features bam-viewer --bin nanalogue_bam_viewer
)

built_binary="$release_dir/nanalogue_bam_viewer"
cp "$built_binary" "$artifact"
chmod 755 "$artifact"

echo "Artifact: $artifact"
echo "Size: $(stat -c %s "$artifact") bytes"
sha256sum "$artifact"
