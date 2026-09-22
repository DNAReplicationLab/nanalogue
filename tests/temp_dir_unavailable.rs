#![cfg_attr(coverage_nightly, feature(coverage_attribute))]
#![cfg(unix)]

//! An inaccessible `TMPDIR` is rejected before simulation output starts.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::Error;
    use nanalogue_core::simulate_mod_bam::{AlignmentFormat, SimulationConfig, TempBamSimulation};
    use std::ffi::OsString;
    use std::path::Path;

    /// Unavailable `TMPDIR`, asserted below not to be a directory.
    const UNAVAILABLE_TMPDIR: &str = "/definitely/not/a/dir";

    /// Minimal valid config, mirroring other simulation fixtures.
    const MINIMAL_CONFIG: &str = r#"{
        "contigs": { "number": 1, "len_range": [40, 40] },
        "reads": [{ "number": 2, "mapq_range": [10, 20], "base_qual_range": [10, 20],
            "len_range": [0.5, 0.5] }]
    }"#;

    /// Restores the saved `TMPDIR` value on drop, even if an assertion panics.
    #[derive(Debug)]
    struct TmpdirGuard {
        /// Value to restore, or `None` if `TMPDIR` was originally unset.
        saved: Option<OsString>,
    }

    impl Drop for TmpdirGuard {
        fn drop(&mut self) {
            match self.saved.take() {
                // SAFETY: this single-test binary has no concurrent env reader.
                Some(value) => unsafe { std::env::set_var("TMPDIR", value) },
                // SAFETY: as above; removes only this binary's override.
                None => unsafe { std::env::remove_var("TMPDIR") },
            }
        }
    }

    /// An inaccessible `TMPDIR` reports `InvalidState` rather than being
    /// silently ignored or replaced by another directory.
    #[test]
    fn inaccessible_tmpdir_is_rejected() {
        assert!(
            !Path::new(UNAVAILABLE_TMPDIR).is_dir(),
            "{UNAVAILABLE_TMPDIR} must not be a directory"
        );
        let _guard = TmpdirGuard {
            saved: std::env::var_os("TMPDIR"),
        };
        // SAFETY: the only test in this binary sets this override, so no other
        // thread reads the environment; `_guard` restores it even on panic.
        unsafe {
            std::env::set_var("TMPDIR", UNAVAILABLE_TMPDIR);
        }
        let config: SimulationConfig =
            serde_json::from_str(MINIMAL_CONFIG).expect("minimal config should deserialize");
        let error = TempBamSimulation::new(config, AlignmentFormat::Bam)
            .expect_err("an inaccessible TMPDIR must be rejected");
        assert!(
            matches!(&error, Error::InvalidState(message)
                if message == "the temp directory environmental variable is inaccessible"),
            "expected the inaccessible-temp-directory InvalidState, got {error:?}"
        );
    }
}
