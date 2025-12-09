//! `PathOrURLOrStdin` enum for handling input sources
//! Represents stdin, file paths, or URLs as input sources

use crate::Error;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::path::PathBuf;
use std::str::FromStr;
use url::Url;

/// Represents different input sources: stdin, file path, or URL.
///
/// This enum allows flexible input handling where data can come from:
/// - Standard input (represented by "-")
/// - A local file path (existence is not validated during parsing)
/// - A remote URL (only http, https, and ftp schemes are recognized)
///
/// # Examples
///
/// ```
/// use nanalogue_core::PathOrURLOrStdin;
/// use std::str::FromStr;
///
/// // Parse from stdin marker
/// let stdin = PathOrURLOrStdin::from_str("-")?;
/// assert!(matches!(stdin, PathOrURLOrStdin::Stdin));
///
/// // Parse from URL
/// let url = PathOrURLOrStdin::from_str("https://example.com/data.bam")?;
/// assert!(matches!(url, PathOrURLOrStdin::URL(_)));
///
/// // Parse from file path (existence not checked during parsing)
/// let path = PathOrURLOrStdin::from_str("examples/example_1.bam")?;
/// assert!(matches!(path, PathOrURLOrStdin::Path(_)));
///
/// // Non-existent paths are also accepted (I/O errors occur at use-time)
/// let path = PathOrURLOrStdin::from_str("/nonexistent/file.txt")?;
/// assert!(matches!(path, PathOrURLOrStdin::Path(_)));
///
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
#[non_exhaustive]
pub enum PathOrURLOrStdin {
    /// Standard input
    #[default]
    Stdin,
    /// A local file path
    Path(PathBuf),
    /// A URL
    URL(Url),
}

impl FromStr for PathOrURLOrStdin {
    type Err = Error;

    /// Parses a string into a `PathOrURLOrStdin` variant.
    ///
    /// The parsing logic follows this order:
    /// 1. If the string is "-", returns `Stdin`
    /// 2. If the string is a valid URL with an allowed network scheme (http, https, ftp), returns `URL(parsed_url)`
    /// 3. Otherwise, returns `Path(parsed_path)` for any syntactically valid path
    ///
    /// **Note**: This method performs parsing only and does not validate file existence.
    /// I/O errors will be surfaced when the path is actually used. This avoids TOCTOU
    /// (Time-of-check to time-of-use) race conditions.
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::PathOrURLOrStdin;
    /// use std::str::FromStr;
    ///
    /// // Stdin
    /// let input = PathOrURLOrStdin::from_str("-")?;
    /// assert!(matches!(input, PathOrURLOrStdin::Stdin));
    ///
    /// // URL
    /// let input = PathOrURLOrStdin::from_str("https://example.com/file.txt")?;
    /// assert!(matches!(input, PathOrURLOrStdin::URL(_)));
    ///
    /// // Path (even if it doesn't exist yet)
    /// let input = PathOrURLOrStdin::from_str("/path/to/file.txt")?;
    /// assert!(matches!(input, PathOrURLOrStdin::Path(_)));
    ///
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// # Errors
    ///
    /// This method should not fail for typical input strings. Parsing is lenient
    /// and treats most inputs as valid paths.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Check for stdin marker
        if s == "-" {
            return Ok(PathOrURLOrStdin::Stdin);
        }

        // Try to parse as URL with allowed network schemes
        if let Ok(parsed_url) = Url::parse(s) {
            // Only accept known network schemes to avoid misclassifying local paths
            const ALLOWED_SCHEMES: &[&str] = &["http", "https", "ftp"];
            if ALLOWED_SCHEMES.contains(&parsed_url.scheme()) {
                return Ok(PathOrURLOrStdin::URL(parsed_url));
            }
            // If it's a valid URL but with an unsupported scheme, fall through to treat as path
        }

        // Otherwise, treat as path (don't check existence to avoid TOCTOU)
        let path = PathBuf::from(s);
        Ok(PathOrURLOrStdin::Path(path))
    }
}

impl fmt::Display for PathOrURLOrStdin {
    /// display "-" or the underlying path or URL.
    /// If the path contains non-UTF-8 characters, you may not get a valid display
    #[expect(
        clippy::pattern_type_mismatch,
        reason = "&self/self/etc. does not make a difference to readability here"
    )]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            PathOrURLOrStdin::Stdin => String::from("-"),
            PathOrURLOrStdin::Path(v) => v.to_string_lossy().to_string(),
            PathOrURLOrStdin::URL(v) => v.to_string(),
        }
        .fmt(f)
    }
}

#[cfg(test)]
#[expect(
    clippy::panic,
    reason = "panic is acceptable in tests for assertion failures"
)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write as _;
    use uuid::Uuid;

    #[test]
    fn from_str_parses_stdin() {
        let result = PathOrURLOrStdin::from_str("-").expect("should parse stdin");
        assert!(matches!(result, PathOrURLOrStdin::Stdin));
    }

    #[test]
    fn from_str_parses_url() {
        let result =
            PathOrURLOrStdin::from_str("https://example.com/file.bam").expect("should parse URL");
        match result {
            PathOrURLOrStdin::URL(u) => {
                assert_eq!(u.scheme(), "https");
                assert_eq!(u.host_str(), Some("example.com"));
                assert_eq!(u.path(), "/file.bam");
            }
            PathOrURLOrStdin::Stdin | PathOrURLOrStdin::Path(_) => {
                panic!("Expected URL variant")
            }
        }
    }

    #[test]
    fn from_str_parses_http_url() {
        let result =
            PathOrURLOrStdin::from_str("http://example.com/data").expect("should parse URL");
        match result {
            PathOrURLOrStdin::URL(u) => {
                assert_eq!(u.scheme(), "http");
            }
            PathOrURLOrStdin::Stdin | PathOrURLOrStdin::Path(_) => {
                panic!("Expected URL variant")
            }
        }
    }

    #[test]
    fn from_str_parses_existing_path() {
        // Create a temporary file with a random UUID name in platform temp directory
        let temp_dir = std::env::temp_dir();
        let temp_filename = temp_dir.join(format!("nanalogue_test_{}.txt", Uuid::new_v4()));
        {
            let mut file = File::create(&temp_filename).expect("should create temp file");
            file.write_all(b"test content")
                .expect("should write to file");
        }

        let result = PathOrURLOrStdin::from_str(
            temp_filename
                .to_str()
                .expect("temp path should be valid UTF-8"),
        )
        .expect("should parse path");
        match result {
            PathOrURLOrStdin::Path(p) => {
                assert_eq!(p, temp_filename);
            }
            PathOrURLOrStdin::Stdin | PathOrURLOrStdin::URL(_) => {
                panic!("Expected Path variant")
            }
        }

        // Clean up
        std::fs::remove_file(&temp_filename).expect("should remove temp file");
    }

    #[test]
    fn from_str_accepts_nonexistent_path() {
        // Use a random UUID to ensure the path doesn't exist
        let nonexistent_path = format!("/nonexistent/path/to/{}.txt", Uuid::new_v4());
        let result =
            PathOrURLOrStdin::from_str(&nonexistent_path).expect("should accept nonexistent path");
        assert!(
            matches!(result, PathOrURLOrStdin::Path(_)),
            "Expected Path variant for non-existent path"
        );
    }

    #[test]
    fn from_str_accepts_any_string_as_path() {
        // Any string that's not "-" and not a valid network URL should be treated as a path
        let result =
            PathOrURLOrStdin::from_str("not a url or valid path").expect("should accept as path");
        assert!(
            matches!(result, PathOrURLOrStdin::Path(_)),
            "Expected Path variant for arbitrary string"
        );
    }

    #[test]
    fn default_is_stdin() {
        let default_val = PathOrURLOrStdin::default();
        assert!(matches!(default_val, PathOrURLOrStdin::Stdin));
    }

    #[test]
    fn url_scheme_variants() {
        // Test allowed network schemes
        let ftp_result = PathOrURLOrStdin::from_str("ftp://example.com/file.txt");
        assert!(
            matches!(ftp_result, Ok(PathOrURLOrStdin::URL(_))),
            "ftp:// should be recognized as URL"
        );

        // Test that non-network schemes like file:// are treated as paths
        let file_result = PathOrURLOrStdin::from_str("file:///path/to/file");
        assert!(
            matches!(file_result, Ok(PathOrURLOrStdin::Path(_))),
            "file:// scheme should be treated as Path, not URL"
        );

        // Test other unsupported schemes are also treated as paths
        let data_result = PathOrURLOrStdin::from_str("data:text/plain,hello");
        assert!(
            matches!(data_result, Ok(PathOrURLOrStdin::Path(_))),
            "data: scheme should be treated as Path"
        );

        // Windows-like paths should also be treated as paths
        let windows_path = PathOrURLOrStdin::from_str("C:/path/to/file.txt");
        assert!(
            matches!(windows_path, Ok(PathOrURLOrStdin::Path(_))),
            "Windows-like paths should be treated as Path"
        );
    }

    #[test]
    fn display_stdin() {
        let stdin = PathOrURLOrStdin::Stdin;
        assert_eq!(stdin.to_string(), "-");
    }

    #[test]
    fn display_path() {
        let path = PathOrURLOrStdin::Path("/some/path/to/file.bam".into());
        assert_eq!(path.to_string(), "/some/path/to/file.bam");
    }

    #[test]
    fn display_url() {
        let url = PathOrURLOrStdin::URL(Url::parse("https://example.com/data.bam").unwrap());
        assert_eq!(url.to_string(), "https://example.com/data.bam");
    }
}
