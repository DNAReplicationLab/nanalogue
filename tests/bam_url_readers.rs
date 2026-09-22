#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! URL-reader scheme-policy coverage.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{
        Error, nanalogue_bam_reader_from_url, nanalogue_indexed_bam_reader_from_url,
    };
    use rust_htslib::bam::FetchDefinition;
    use url::Url;

    /// Both public URL readers reject unsupported schemes before `HTSlib` can
    /// open a local path or invoke an unconfigured transport handler.
    #[test]
    fn url_readers_reject_disallowed_schemes() {
        for literal in [
            "file:///etc/passwd",
            "data:text/plain,not-a-bam",
            "ssh://example.com/input.bam",
            "s3://bucket/input.bam",
            "gs://bucket/input.bam",
        ] {
            let url = Url::parse(literal).expect("test URL must parse");
            let expected = format!(
                "URL scheme `{}` is not in the allow-list (http, https, ftp)",
                url.scheme()
            );
            for error in [
                nanalogue_bam_reader_from_url(&url)
                    .expect_err("ordinary URL reader must reject a disallowed scheme"),
                nanalogue_indexed_bam_reader_from_url(&url, FetchDefinition::All)
                    .expect_err("indexed URL reader must reject a disallowed scheme"),
            ] {
                assert!(
                    matches!(&error, Error::InvalidState(message) if message == &expected),
                    "URL reader must report its allow-list error for `{literal}`, got {error:?}"
                );
            }
        }
    }
}
