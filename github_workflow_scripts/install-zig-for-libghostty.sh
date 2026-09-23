#!/usr/bin/env bash
set -euo pipefail

case "$(uname -s):$(uname -m)" in
    Linux:x86_64)
        archive_name="zig-x86_64-linux-0.15.2.tar.xz"
        archive_sha256="02aa270f183da276e5b5920b1dac44a63f1a49e55050ebde3aecc9eb82f93239"
        ;;
    Darwin:x86_64)
        archive_name="zig-x86_64-macos-0.15.2.tar.xz"
        archive_sha256="375b6909fc1495d16fc2c7db9538f707456bfc3373b14ee83fdd3e22b3d43f7f"
        ;;
    Darwin:arm64)
        archive_name="zig-aarch64-macos-0.15.2.tar.xz"
        archive_sha256="3cc2bab367e185cdfb27501c4b30b1b0653c28d9f73df8dc91488e66ece5fa6b"
        ;;
    *)
        echo "Zig 0.15.2 installation requires Linux x86_64 or macOS." >&2
        exit 1
        ;;
esac

zsf_minisign_key="RWSGOq2NVecA2UPNdBUZykf1CCb147pkmdtYxgb3Ti+JO/wCYvhbAb/U"
zig_dir="$HOME/.local/zig-0.15.2"
checksum_marker="$zig_dir/.archive.sha256"
if [[ -x "$zig_dir/zig" ]] \
    && [[ "$(cat "$checksum_marker" 2>/dev/null)" == "$archive_sha256" ]]; then
    if installed_version="$("$zig_dir/zig" version)" \
        && [[ "$installed_version" == 0.15.2 ]]; then
        exit 0
    fi
fi

download_dir="$(mktemp -d)"
trap 'rm -rf "$download_dir"' EXIT
archive="$download_dir/$archive_name"
signature="$archive.minisig"
mirror_list="$download_dir/community-mirrors.txt"
staged_dir="$download_dir/zig"

if ! command -v minisign >/dev/null; then
    echo "minisign is required to authenticate community mirror downloads." >&2
    exit 1
fi

curl --proto '=https' --tlsv1.2 -fsSLo "$mirror_list" \
    https://ziglang.org/download/community-mirrors.txt

verified=false
mirror_count=0
while IFS= read -r mirror; do
    ((mirror_count += 1))
    if [[ "$mirror" != https://* || "$mirror" == *[[:space:]]* ]]; then
        echo "Invalid Zig community mirror URL: $mirror" >&2
        exit 1
    fi

    rm -f "$archive" "$signature"
    archive_url="${mirror%/}/$archive_name"
    if ! curl --proto '=https' --tlsv1.2 --connect-timeout 15 --max-time 180 \
        -fsSLo "$archive" "$archive_url?source=github-DNAReplicationLab-nanalogue"; then
        continue
    fi
    if ! curl --proto '=https' --tlsv1.2 --connect-timeout 15 --max-time 30 \
        -fsSLo "$signature" "$archive_url.minisig?source=github-DNAReplicationLab-nanalogue"; then
        continue
    fi
    if ! printf '%s  %s\n' "$archive_sha256" "$archive" \
        | shasum -a 256 -c -; then
        continue
    fi
    if ! minisign -Vm "$archive" -x "$signature" -P "$zsf_minisign_key"; then
        continue
    fi

    trusted_comment="$(sed -n '3s/^trusted comment: //p' "$signature")"
    trusted_filename=""
    IFS=$'\t' read -ra comment_fields <<<"$trusted_comment"
    for field in "${comment_fields[@]:0:10}"; do
        if [[ "$field" == file:* ]]; then
            trusted_filename="${field#file:}"
        fi
    done
    if [[ "$trusted_filename" != "$archive_name" ]]; then
        continue
    fi

    verified=true
    break
done < <(head -n 10 "$mirror_list")

if ((mirror_count == 0)); then
    echo "No Zig community mirrors were returned." >&2
    exit 1
fi

if [[ "$verified" != true ]]; then
    echo "Zig download is not possible after trying up to 10 community mirrors." >&2
    exit 1
fi

mkdir "$staged_dir"
tar -xJ --strip-components=1 -C "$staged_dir" -f "$archive"
printf '%s\n' "$archive_sha256" >"$staged_dir/.archive.sha256"
mkdir -p "$(dirname "$zig_dir")"
rm -rf "$zig_dir"
mv "$staged_dir" "$zig_dir"
