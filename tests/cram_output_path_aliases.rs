#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Output aliases must be rejected before CRAM or reference files are created.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{Error, file_utils::write_cram_denovo, uuid};
    use rust_htslib::bam::Record;
    use std::{fs, num::NonZeroU32, path::Path};

    /// No records or valid reference are needed: path validation must run first.
    fn write_empty_cram(output: &Path, reference: &Path) -> Result<(), Error> {
        write_cram_denovo(
            Vec::<Record>::new(),
            [("chr1".to_owned(), 12)],
            Vec::<String>::new(),
            Vec::<String>::new(),
            output,
            reference,
            NonZeroU32::MIN,
        )
    }

    #[test]
    fn nonexistent_reference_aliases_future_crai_through_parent_components() {
        let root = std::env::temp_dir().join(uuid::v4_random());
        let child = root.join("child");
        fs::create_dir_all(&child).expect("create temporary directory");
        let output = root.join("output.cram");
        let reference = child.join("..").join("output.cram.crai");

        // Neither inode exists yet, and the strings are distinct. Only resolving
        // the parent directory reveals that the reference aliases the CRAI.
        assert!(!output.exists());
        assert!(!reference.exists());
        let error = write_empty_cram(&output, &reference).expect_err("aliased outputs");
        assert!(matches!(
            error,
            Error::InvalidState(message)
                if message == "CRAM, CRAI, and FASTA outputs must use different paths"
        ));
        assert!(!output.exists(), "validation must precede opening the CRAM");
        assert!(!reference.exists(), "no CRAI or reference may be created");
        assert!(!root.join("output.cram.crai.fai").exists());
        assert_eq!(fs::read_dir(&root).expect("read directory").count(), 1);
        fs::remove_dir_all(root).expect("remove temporary directory");
    }

    #[cfg(unix)]
    #[test]
    fn nonexistent_output_aliases_future_fai_through_symlinked_parent() {
        let root = std::env::temp_dir().join(uuid::v4_random());
        let real = root.join("real");
        let alias = root.join("alias");
        fs::create_dir_all(&real).expect("create temporary directory");
        std::os::unix::fs::symlink(&real, &alias).expect("create directory alias");
        let reference = real.join("reference.fa");
        let output = alias.join("reference.fa.fai");
        let sentinel = b">chr1\nACGTACGTACGT\n";
        fs::write(&reference, sentinel).expect("write reference sentinel");

        // The reference exists, but the colliding output and FASTA index do not.
        // This must not truncate the reference or create either sidecar.
        assert!(!output.exists());
        let error = write_empty_cram(&output, &reference).expect_err("aliased outputs");
        assert!(matches!(
            error,
            Error::InvalidState(message)
                if message == "CRAM, CRAI, and FASTA outputs must use different paths"
        ));
        assert_eq!(fs::read(&reference).expect("read reference"), sentinel);
        assert!(!output.exists());
        assert!(!real.join("reference.fa.fai.crai").exists());
        assert_eq!(fs::read_dir(&real).expect("read directory").count(), 1);
        fs::remove_dir_all(root).expect("remove temporary directory");
    }

    #[cfg(unix)]
    #[test]
    #[expect(clippy::panic, reason = "unexpected error variants fail the test")]
    fn metadata_failure_is_not_mistaken_for_a_missing_output() {
        let root = std::env::temp_dir().join(uuid::v4_random());
        fs::create_dir_all(&root).expect("create temporary directory");
        let output = root.join("loop.cram");
        let reference = root.join("reference.fa");
        std::os::unix::fs::symlink(&output, &output).expect("create self-referential link");
        let sentinel = b"reference must remain unchanged";
        fs::write(&reference, sentinel).expect("write reference sentinel");

        // A symlink loop is a metadata error, not NotFound. Compare against the
        // OS error directly rather than depending on a platform's ELOOP number.
        let expected = fs::metadata(&output).expect_err("link cannot resolve");
        let error = write_empty_cram(&output, &reference).expect_err("metadata failure");
        let Error::InputOutputError(actual) = error else {
            panic!("expected filesystem error, got {error:?}");
        };
        assert_eq!(actual.raw_os_error(), expected.raw_os_error());
        assert_eq!(fs::read_link(&output).expect("link remains"), output);
        assert_eq!(fs::read(&reference).expect("read reference"), sentinel);
        assert!(!root.join("loop.cram.crai").exists());
        assert!(!root.join("reference.fa.fai").exists());
        fs::remove_dir_all(root).expect("remove temporary directory");
    }
}
