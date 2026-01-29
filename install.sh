#!/bin/sh
# Nanalogue install script - downloads and installs nanalogue binaries
# Usage: curl -fsSL https://raw.githubusercontent.com/DNAReplicationLab/nanalogue/main/install.sh | sh
#
# Options:
#   -y          Non-interactive mode: use default install directory without prompting
#
# Install directory:
#   - Press Enter at the prompt to accept the default (/usr/local/bin)
#   - Type a custom path at the prompt
#   - Use -y flag to skip the prompt and use the default:
#     curl -fsSL ... | sh -s -- -y
#   - Or set the environment variable:
#     NANALOGUE_INSTALL_DIR=/your/path curl -fsSL ... | sh
#
# Dependencies: curl or wget, unzip, jq, sha256sum or shasum

set -e

REPO="DNAReplicationLab/nanalogue"
DEFAULT_INSTALL_DIR="/usr/local/bin"
INSTALL_DIR="${NANALOGUE_INSTALL_DIR:-}"
NON_INTERACTIVE="${1:-}"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

info()    { printf "${BLUE}==>${NC} %s\n" "$1"; }
success() { printf "${GREEN}==>${NC} %s\n" "$1"; }
warn()    { printf "${YELLOW}Warning:${NC} %s\n" "$1"; }
error()   { printf "${RED}Error:${NC} %s\n" "$1" >&2; exit 1; }

has_cmd() { command -v "$1" >/dev/null 2>&1; }

check_dependencies() {
    has_cmd curl || has_cmd wget || error "Either curl or wget is required but neither is installed."
    has_cmd unzip || error "unzip is required but not installed."
    has_cmd jq || error "jq is required but not installed."
    has_cmd sha256sum || has_cmd shasum || error "sha256sum or shasum is required but not installed."
}

download() {
    if has_cmd curl; then
        curl -fsSL -H "User-Agent: nanalogue-installer" --retry 3 --retry-delay 2 --connect-timeout 10 --max-time 120 "$1" -o "$2"
    else
        wget -q --header="User-Agent: nanalogue-installer" --tries=3 --waitretry=2 --timeout=20 "$1" -O "$2"
    fi
}

fetch() {
    url="$1"
    auth_header=""

    # Add authentication for GitHub API requests when token is available
    case "$url" in
        https://api.github.com/*)
            [ -n "${GITHUB_TOKEN:-}" ] && auth_header="Authorization: Bearer $GITHUB_TOKEN"
            ;;
    esac

    if has_cmd curl; then
        if [ -n "$auth_header" ]; then
            curl -fsSL -H "User-Agent: nanalogue-installer" -H "$auth_header" --retry 3 --retry-delay 2 --connect-timeout 10 --max-time 120 "$url"
        else
            curl -fsSL -H "User-Agent: nanalogue-installer" --retry 3 --retry-delay 2 --connect-timeout 10 --max-time 120 "$url"
        fi
    else
        if [ -n "$auth_header" ]; then
            wget -qO- --header="User-Agent: nanalogue-installer" --header="$auth_header" --tries=3 --waitretry=2 --timeout=20 "$url"
        else
            wget -qO- --header="User-Agent: nanalogue-installer" --tries=3 --waitretry=2 --timeout=20 "$url"
        fi
    fi
}

fetch_release_json() {
    release_json_file="$1"
    fetch "https://api.github.com/repos/${REPO}/releases/latest" > "$release_json_file"
    test -s "$release_json_file" || error "Failed to fetch release info from GitHub API."
}

validate_json() {
    json_file="$1"
    jq -e . < "$json_file" >/dev/null 2>&1 || error "Invalid JSON received from GitHub API."
}

check_github_error() {
    json_file="$1"
    message=$(jq -r '.message // .error // .errors[0].message // empty' < "$json_file" 2>/dev/null || true)
    [ -z "$message" ] && return 0

    docs=$(jq -r '.documentation_url // empty' < "$json_file" 2>/dev/null || true)
    error "GitHub API error: $message${docs:+ (see: $docs)}"
}

parse_json_field() {
    json_file="$1"
    field="${2#.}"  # Strip leading dot if present
    jq -r ".$field // empty" < "$json_file" 2>/dev/null || true
}

ensure_asset_digest() {
    json_file="$1"
    asset_name="$2"

    asset_found=$(jq -r ".assets[] | select(.name == \"$asset_name\") | .name" < "$json_file" 2>/dev/null | head -1)
    test -n "$asset_found" || error "Release asset not found: $asset_name"

    digest=$(get_asset_digest "$json_file" "$asset_name")
    test -n "$digest" || error "Release asset is missing a sha256 digest: $asset_name"
}

get_asset_digest() {
    json_file="$1"
    asset_name="$2"
    jq -r ".assets[] | select(.name == \"$asset_name\") | .digest" < "$json_file" 2>/dev/null || true
}

compute_sha256() {
    file="$1"

    if has_cmd sha256sum; then
        sha256sum "$file" | cut -d' ' -f1
    elif has_cmd shasum; then
        shasum -a 256 "$file" | cut -d' ' -f1
    else
        echo ""
    fi
}

verify_checksum() {
    file="$1"
    expected="$2"

    if [ -z "$expected" ]; then
        error "No checksum available for verification. Release assets must include a sha256 digest."
    fi

    actual_hash=$(compute_sha256 "$file")
    if [ -z "$actual_hash" ]; then
        error "No sha256sum or shasum available, cannot verify checksum."
    fi

    expected_hash="${expected#sha256:}"  # Strip sha256: prefix if present

    if [ "$actual_hash" != "$expected_hash" ]; then
        error "Checksum verification failed!
  Expected: $expected_hash
  Got:      $actual_hash
  The downloaded file may be corrupted or tampered with."
    fi

    success "Checksum verified."
}

detect_os() {
    case "$(uname -s)" in
        Linux*)                        echo "linux" ;;
        Darwin*)                       echo "darwin" ;;
        MINGW*|MSYS*|CYGWIN*|Windows_NT) echo "windows" ;;
        *)                             error "Unsupported operating system: $(uname -s)" ;;
    esac
}

detect_arch() {
    case "$(uname -m)" in
        x86_64|amd64)   echo "x86_64" ;;
        aarch64|arm64)  echo "aarch64" ;;
        *)              error "Unsupported architecture: $(uname -m)" ;;
    esac
}

is_musl() {
    if ls /lib/libc.musl-*.so.1 >/dev/null 2>&1; then
        return 0
    fi
    if has_cmd ldd && ldd --version 2>&1 | grep -qi musl; then
        return 0
    fi
    return 1
}

get_glibc_version() {
    if ! has_cmd ldd; then
        echo ""
        return
    fi
    ldd --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+' | head -1
}

version_gte() {
    v1_major=${1%%.*}
    v1_minor=${1#*.}
    v1_minor=${v1_minor%%.*}
    v1_minor=${v1_minor:-0}

    v2_major=${2%%.*}
    v2_minor=${2#*.}
    v2_minor=${v2_minor%%.*}
    v2_minor=${v2_minor:-0}

    if [ "$v1_major" -gt "$v2_major" ]; then
        return 0
    fi
    if [ "$v1_major" -eq "$v2_major" ] && [ "$v1_minor" -ge "$v2_minor" ]; then
        return 0
    fi
    return 1
}

get_asset_name() {
    os="$1"
    arch="$2"

    case "$os" in
        darwin)
            if [ "$arch" = "aarch64" ]; then
                echo "binaries-macos-latest.zip"
            else
                echo "binaries-macos-15-intel.zip"
            fi
            ;;
        linux)
            if is_musl; then
                echo "binaries-musllinux_1_2_${arch}.zip"
                return
            fi
            glibc_version=$(get_glibc_version)
            if [ -n "$glibc_version" ] && version_gte "$glibc_version" "2.34"; then
                echo "binaries-manylinux_2_34_${arch}.zip"
            elif [ -n "$glibc_version" ] && version_gte "$glibc_version" "2.28"; then
                echo "binaries-manylinux_2_28_${arch}.zip"
            else
                echo "binaries-manylinux2014_${arch}.zip"
            fi
            ;;
        *)
            error "Cannot determine asset for OS: $os"
            ;;
    esac
}

show_windows_message() {
    printf "\n${YELLOW}========================================${NC}\n"
    printf "${YELLOW}  Windows is not directly supported${NC}\n"
    printf "${YELLOW}========================================${NC}\n\n"
    printf "You can use nanalogue via Docker instead:\n\n"
    printf "  ${GREEN}docker pull dockerofsat/nanalogue:latest${NC}\n\n"
    printf "Then run commands like:\n\n"
    printf "  ${BLUE}docker run --rm -v \$(pwd):/data dockerofsat/nanalogue:latest nanalogue --help${NC}\n\n"
    exit 0
}

prompt_install_dir() {
    # Use environment variable if set
    if [ -n "$INSTALL_DIR" ]; then
        return
    fi

    # Use default if non-interactive mode or no TTY available
    if [ "$NON_INTERACTIVE" = "-y" ] || { [ ! -t 0 ] && [ ! -e /dev/tty ]; }; then
        INSTALL_DIR="$DEFAULT_INSTALL_DIR"
        return
    fi

    # Prompt user for install directory
    printf "\n  Where would you like to install nanalogue?\n"
    printf "  Press Enter for default [${BLUE}$DEFAULT_INSTALL_DIR${NC}]: "
    read -r user_input </dev/tty || user_input=""
    INSTALL_DIR="${user_input:-$DEFAULT_INSTALL_DIR}"
}

check_existing_install() {
    existing_path=$(command -v nanalogue 2>/dev/null) || true
    if [ -z "$existing_path" ]; then
        return 0
    fi

    existing_dir=$(dirname "$existing_path")
    if [ "$existing_dir" = "$INSTALL_DIR" ]; then
        return 0
    fi

    error "nanalogue is already installed at: $existing_path
  To update it, run:
    NANALOGUE_INSTALL_DIR=$existing_dir curl -fsSL https://raw.githubusercontent.com/DNAReplicationLab/nanalogue/main/install.sh | sh"
}

install_binary() {
    $need_sudo cp "$1" "$INSTALL_DIR/$2"
    $need_sudo chmod +x "$INSTALL_DIR/$2"
    success "Installed: $2"
}

main() {
    printf "\n"
    printf "${BLUE}╔═══════════════════════════════════════╗${NC}\n"
    printf "${BLUE}║${NC}      Nanalogue Installer              ${BLUE}║${NC}\n"
    printf "${BLUE}╚═══════════════════════════════════════╝${NC}\n"
    printf "\n"

    info "Detecting platform..."
    os=$(detect_os)
    if [ "$os" = "windows" ]; then
        show_windows_message
    fi

    check_dependencies
    prompt_install_dir
    check_existing_install

    arch=$(detect_arch)
    info "Detected: $os / $arch"

    tmp_dir=$(mktemp -d)
    trap 'rm -rf "$tmp_dir"' EXIT
    release_json_file="${tmp_dir}/release.json"

    info "Fetching latest release info..."
    fetch_release_json "$release_json_file"
    validate_json "$release_json_file"
    check_github_error "$release_json_file"
    version=$(parse_json_field "$release_json_file" ".tag_name")
    test -n "$version" || error "Failed to parse version from GitHub API response."
    info "Latest version: $version"

    asset_name=$(get_asset_name "$os" "$arch")
    info "Selected binary: $asset_name"

    ensure_asset_digest "$release_json_file" "$asset_name"
    expected_digest=$(get_asset_digest "$release_json_file" "$asset_name")

    download_url="https://github.com/${REPO}/releases/download/${version}/${asset_name}"
    zip_file="${tmp_dir}/${asset_name}"

    info "Downloading from: $download_url"
    download "$download_url" "$zip_file" || error "Failed to download binary."
    success "Download complete."

    info "Verifying checksum..."
    verify_checksum "$zip_file" "$expected_digest"

    info "Extracting archive..."
    unzip -q "$zip_file" -d "$tmp_dir"

    nanalogue_bin=$(find "$tmp_dir" -name "nanalogue-*" -type f | head -1)
    nanalogue_sim_bin=$(find "$tmp_dir" -name "nanalogue_sim_bam-*" -type f | head -1)
    test -n "$nanalogue_bin" || error "Could not find nanalogue binary in archive."

    need_sudo=""
    check_path="${INSTALL_DIR}"
    if [ ! -d "$INSTALL_DIR" ]; then
        check_path=$(dirname "$INSTALL_DIR")
    fi
    if [ ! -w "$check_path" ]; then
        has_cmd sudo || error "Cannot write to $INSTALL_DIR and sudo is not available."
        need_sudo="sudo"
        info "Installation requires sudo privileges."
    fi

    if [ ! -d "$INSTALL_DIR" ]; then
        info "Creating directory: $INSTALL_DIR"
        $need_sudo mkdir -p "$INSTALL_DIR"
    fi

    info "Installing to $INSTALL_DIR..."
    install_binary "$nanalogue_bin" "nanalogue"
    if [ -n "$nanalogue_sim_bin" ]; then
        install_binary "$nanalogue_sim_bin" "nanalogue_sim_bam"
    fi

    printf "\n"
    if has_cmd nanalogue; then
        installed_version=$(nanalogue --version 2>/dev/null || echo "unknown")
        success "Installation complete!"
        printf "\n  ${GREEN}nanalogue${NC} is now available at: ${BLUE}$INSTALL_DIR/nanalogue${NC}\n"
        printf "  Version: ${BLUE}$installed_version${NC}\n\n"
        printf "  Run ${YELLOW}nanalogue --help${NC} to get started.\n\n"
    else
        warn "Installation complete, but nanalogue is not in your PATH."
        printf "\n  Add ${BLUE}$INSTALL_DIR${NC} to your PATH, or run directly:\n"
        printf "  ${YELLOW}$INSTALL_DIR/nanalogue --help${NC}\n\n"
    fi
}

main "$@"
