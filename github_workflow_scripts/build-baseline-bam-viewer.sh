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

zig_path="$(command -v zig || true)"
if [[ -z "$zig_path" && -x "$HOME/.local/zig-0.15.2/zig" ]]; then
    zig_path="$HOME/.local/zig-0.15.2/zig"
fi
if [[ -z "$zig_path" || "$("$zig_path" version)" != 0.15.2 ]]; then
    echo "Zig 0.15.2 is required." >&2
    echo "Install it with github_workflow_scripts/install-zig-for-libghostty.sh." >&2
    exit 1
fi

target_dir="$repo_root/target/portable-bam-viewer"
rust_target="x86_64-unknown-linux-gnu"
release_dir="$target_dir/$rust_target/release"
output_dir="$repo_root/target/portable-dist"
artifact="$output_dir/nanalogue_bam_viewer-x86_64-baseline"
wrapper_dir="$(mktemp -d)"
trap 'rm -rf "$wrapper_dir"' EXIT

# libghostty-vt-sys 0.2.1 does not expose CPU selection to downstream builds.
# Its native build invokes `zig build` without -Dtarget or -Dcpu, so Zig may
# optimize Ghostty for the build machine. This temporary executable named
# `zig` intercepts that invocation and adds the exact Zig option
# `-Dcpu=baseline` before forwarding it to the real Zig 0.15.2 executable.
#
# Upstream commit f720ad74a66e333181986fd5008f34c0a788a119 added the environment
# variable LIBGHOSTTY_VT_SYS_CPU and made `baseline` its default:
# https://github.com/uzaaft/libghostty-rs/commit/f720ad74a66e333181986fd5008f34c0a788a119
# Once nanalogue upgrades to a libghostty-vt release containing that change,
# this wrapper will be unnecessary. The build can pass
# LIBGHOSTTY_VT_SYS_CPU=baseline directly (or rely on its baseline default),
# and this script can become a direct Cargo command or be replaced by README
# instructions.
cat >"$wrapper_dir/zig" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == build ]]; then
    exec "$NANALOGUE_REAL_ZIG" build -Dcpu=baseline "${@:2}"
fi
exec "$NANALOGUE_REAL_ZIG" "$@"
EOF
chmod +x "$wrapper_dir/zig"

rm -rf "$target_dir"
mkdir -p "$output_dir"
rm -f "$artifact"
(
    unset CARGO_ENCODED_RUSTFLAGS
    export CARGO_TARGET_DIR="$target_dir"
    export NANALOGUE_REAL_ZIG="$zig_path"
    export PATH="$wrapper_dir:$PATH"
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
