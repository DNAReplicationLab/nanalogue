//! System certificate location probing adapted from `openssl-probe` version 0.2.1.

use std::env;
use std::path::{Path, PathBuf};

/// Find certificate locations configured in the environment or known to the platform.
pub(crate) fn probe() -> (Option<PathBuf>, Vec<PathBuf>) {
    let existing_path = |name| {
        env::var_os(name)
            .map(PathBuf::from)
            .filter(|path| path.exists())
    };
    let cert_file = existing_path("SSL_CERT_FILE").or_else(|| {
        CERTIFICATE_FILE_NAMES
            .iter()
            .find_map(|path| Path::new(path).exists().then(|| PathBuf::from(path)))
    });

    let mut cert_dir: Vec<PathBuf> = existing_path("SSL_CERT_DIR").into_iter().collect();
    cert_dir.extend(
        CERTIFICATE_DIRS
            .iter()
            .map(PathBuf::from)
            .filter(|path| path.exists()),
    );
    (cert_file, cert_dir)
}

/// Known certificate directories on Linux.
#[cfg(target_os = "linux")]
const CERTIFICATE_DIRS: &[&str] = &[
    "/etc/ssl/certs",
    "/etc/pki/tls/certs",
    "/etc/security/certificates",
];

/// Known certificate directories on FreeBSD.
#[cfg(target_os = "freebsd")]
const CERTIFICATE_DIRS: &[&str] = &["/etc/ssl/certs", "/usr/local/share/certs"];

/// Known certificate directories on illumos and Solaris.
#[cfg(any(target_os = "illumos", target_os = "solaris"))]
const CERTIFICATE_DIRS: &[&str] = &["/etc/certs/CA"];

/// Known certificate directories on NetBSD.
#[cfg(target_os = "netbsd")]
const CERTIFICATE_DIRS: &[&str] = &["/etc/openssl/certs"];

/// Known certificate directories on AIX.
#[cfg(target_os = "aix")]
const CERTIFICATE_DIRS: &[&str] = &["/var/ssl/certs"];

/// Default certificate directory on other platforms, including macOS.
#[cfg(not(any(
    target_os = "linux",
    target_os = "freebsd",
    target_os = "illumos",
    target_os = "solaris",
    target_os = "netbsd",
    target_os = "aix"
)))]
const CERTIFICATE_DIRS: &[&str] = &["/etc/ssl/certs"];

/// Known certificate bundle files on Linux.
#[cfg(target_os = "linux")]
const CERTIFICATE_FILE_NAMES: &[&str] = &[
    "/etc/ssl/certs/ca-certificates.crt",
    "/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem",
    "/etc/pki/tls/certs/ca-bundle.crt",
    "/etc/ssl/ca-bundle.pem",
    "/etc/pki/tls/cacert.pem",
    "/etc/ssl/cert.pem",
    "/opt/etc/ssl/certs/ca-certificates.crt",
    "/etc/ssl/certs/cacert.pem",
];

/// Known certificate bundle files on FreeBSD.
#[cfg(target_os = "freebsd")]
const CERTIFICATE_FILE_NAMES: &[&str] = &["/usr/local/etc/ssl/cert.pem"];

/// Known certificate bundle files on DragonFly BSD.
#[cfg(target_os = "dragonfly")]
const CERTIFICATE_FILE_NAMES: &[&str] = &["/usr/local/share/certs/ca-root-nss.crt"];

/// Known certificate bundle files on NetBSD.
#[cfg(target_os = "netbsd")]
const CERTIFICATE_FILE_NAMES: &[&str] = &["/etc/openssl/certs/ca-certificates.crt"];

/// Known certificate bundle files on OpenBSD.
#[cfg(target_os = "openbsd")]
const CERTIFICATE_FILE_NAMES: &[&str] = &["/etc/ssl/cert.pem"];

/// Known certificate bundle files on Solaris.
#[cfg(target_os = "solaris")]
const CERTIFICATE_FILE_NAMES: &[&str] = &["/etc/certs/ca-certificates.crt"];

/// Known certificate bundle files on illumos.
#[cfg(target_os = "illumos")]
const CERTIFICATE_FILE_NAMES: &[&str] = &["/etc/ssl/cacert.pem", "/etc/certs/ca-certificates.crt"];

/// Known certificate bundle files on Android under Termux.
#[cfg(target_os = "android")]
const CERTIFICATE_FILE_NAMES: &[&str] = &["/data/data/com.termux/files/usr/etc/tls/cert.pem"];

/// Known certificate bundle files on Haiku.
#[cfg(target_os = "haiku")]
const CERTIFICATE_FILE_NAMES: &[&str] = &["/boot/system/data/ssl/CARootCertificates.pem"];

/// Default certificate bundle file on other platforms, including macOS.
#[cfg(not(any(
    target_os = "linux",
    target_os = "freebsd",
    target_os = "dragonfly",
    target_os = "netbsd",
    target_os = "openbsd",
    target_os = "solaris",
    target_os = "illumos",
    target_os = "android",
    target_os = "haiku",
)))]
const CERTIFICATE_FILE_NAMES: &[&str] = &["/etc/ssl/certs/ca-certificates.crt"];
