#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Output aliases must be rejected before CRAM or reference files are created.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{
        Error,
        file_utils::{write_bam_denovo, write_cram_denovo},
        uuid,
    };
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

    /// A regular file cannot be traversed as a directory. Canonicalizing an
    /// output below such a component must report the I/O condition instead of
    /// treating the parent as missing and resolving the output beside the file.
    #[test]
    #[expect(clippy::panic, reason = "unexpected error variants fail the test")]
    fn file_parent_component_is_an_io_error_not_a_missing_output() {
        let root = std::env::temp_dir().join(uuid::v4_random());
        fs::create_dir_all(&root).expect("create temporary directory");
        let blocker = root.join("example_1.bam");
        fs::write(&blocker, b"regular file, not a directory").expect("write blocking file");
        // `examples/example_1.bam/out.cram` has this shape: the parent component
        // exists, but it is a regular file rather than a directory.
        let output = blocker.join("out.cram");
        let reference = root.join("reference.fa");

        // Derive the expected condition from the OS for the same path so the
        // assertion does not depend on a platform's raw ENOTDIR number.
        let expected = fs::metadata(&output).expect_err("a file cannot be traversed");
        let error = write_empty_cram(&output, &reference).expect_err("output cannot be created");
        let Error::InputOutputError(actual) = error else {
            panic!("expected a filesystem error, got {error:?}");
        };
        assert_eq!(
            actual.raw_os_error(),
            expected.raw_os_error(),
            "raw OS condition must match a direct metadata call"
        );
        assert_eq!(actual.kind(), expected.kind(), "error kind must match");
        assert_ne!(
            actual.kind(),
            std::io::ErrorKind::NotFound,
            "a regular file is not a missing output"
        );
        #[cfg(unix)]
        assert_eq!(
            actual.kind(),
            std::io::ErrorKind::NotADirectory,
            "Unix reports ENOTDIR when a path component is a regular file"
        );
        assert_eq!(
            fs::read(&blocker).expect("read blocking file"),
            b"regular file, not a directory",
            "the blocking file must remain unchanged"
        );
        assert!(!reference.exists(), "the reference must not be created");
        assert!(!root.join("reference.fa.fai").exists());
        fs::remove_dir_all(root).expect("remove temporary directory");
    }

    /// An interior NUL byte cannot name a file on Unix. The output-identity
    /// check rejects it first, so the reachable failure is `InvalidInput` from
    /// the metadata call; `path_to_c_string`'s own NUL rejection is never
    /// reached through this entry point.
    #[cfg(unix)]
    #[test]
    #[expect(clippy::panic, reason = "unexpected error variants fail the test")]
    fn interior_nul_byte_in_output_is_invalid_input() {
        use std::ffi::OsStr;
        use std::os::unix::ffi::OsStrExt as _;

        let root = std::env::temp_dir().join(uuid::v4_random());
        fs::create_dir_all(&root).expect("create temporary directory");
        let output = root.join(OsStr::from_bytes(b"out\0.cram"));
        let reference = root.join("reference.fa");

        let error = write_empty_cram(&output, &reference).expect_err("NUL cannot name a file");
        let Error::InputOutputError(actual) = error else {
            panic!("expected a filesystem error, got {error:?}");
        };
        assert_eq!(
            actual.kind(),
            std::io::ErrorKind::InvalidInput,
            "an interior NUL byte must be reported as invalid input"
        );
        assert!(!reference.exists(), "the reference must not be created");
        assert!(!root.join("reference.fa.fai").exists());
        fs::remove_dir_all(root).expect("remove temporary directory");
    }

    /// The dangling-symlink diagnostic also covers the output path itself: the
    /// link must be left alone and no CRAM or sidecar may be created.
    #[cfg(unix)]
    #[test]
    fn dangling_output_symlink_is_rejected_before_any_output_is_created() {
        let root = std::env::temp_dir().join(uuid::v4_random());
        fs::create_dir_all(&root).expect("create temporary directory");
        let output = root.join("output.cram");
        let missing_target = root.join("missing.cram");
        std::os::unix::fs::symlink(&missing_target, &output).expect("create dangling output link");
        let reference = root.join("reference.fa");
        let sentinel = b">chr1\nACGTACGTACGT\n";
        fs::write(&reference, sentinel).expect("write reference sentinel");

        let error = write_empty_cram(&output, &reference).expect_err("dangling output link");
        assert!(matches!(
            error,
            Error::InvalidState(message)
                if message == "output paths must not be symbolic links to files that do not exist"
        ));
        assert_eq!(
            fs::read_link(&output).expect("link remains"),
            missing_target,
            "the dangling link must not be replaced"
        );
        assert_eq!(fs::read(&reference).expect("read reference"), sentinel);
        assert!(!root.join("output.cram.crai").exists());
        assert!(!root.join("reference.fa.fai").exists());
        fs::remove_dir_all(root).expect("remove temporary directory");
    }

    /// The same diagnostic covers a dangling reference-FASTA symlink, which is
    /// also the base of the future FASTA index.
    #[cfg(unix)]
    #[test]
    fn dangling_reference_symlink_is_rejected_before_any_output_is_created() {
        let root = std::env::temp_dir().join(uuid::v4_random());
        fs::create_dir_all(&root).expect("create temporary directory");
        let output = root.join("output.cram");
        let reference = root.join("reference.fa");
        let missing_target = root.join("missing.fa");
        std::os::unix::fs::symlink(&missing_target, &reference)
            .expect("create dangling reference link");

        let error = write_empty_cram(&output, &reference).expect_err("dangling reference link");
        assert!(matches!(
            error,
            Error::InvalidState(message)
                if message == "output paths must not be symbolic links to files that do not exist"
        ));
        assert!(!output.exists(), "no CRAM output may be created");
        assert!(!root.join("output.cram.crai").exists());
        assert_eq!(
            fs::read_link(&reference).expect("link remains"),
            missing_target,
            "the dangling link must not be replaced"
        );
        assert!(!root.join("reference.fa.fai").exists());
        fs::remove_dir_all(root).expect("remove temporary directory");
    }

    /// `write_bam_denovo` opens its output directly, so a missing parent
    /// directory is reported by the BAM writer. `write_cram_denovo` checks
    /// output identities first and reports the same missing directory as an
    /// I/O `NotFound` instead.
    #[test]
    #[expect(clippy::panic, reason = "unexpected error variants fail the test")]
    fn write_bam_denovo_missing_directory_reports_the_writer_error() {
        let root = std::env::temp_dir().join(uuid::v4_random());
        fs::create_dir_all(&root).expect("create temporary directory");
        let output = root.join("missing").join("output.bam");

        let result = write_bam_denovo(
            Vec::<Record>::new(),
            vec![("chr1".to_owned(), 12)],
            Vec::<String>::new(),
            Vec::<String>::new(),
            &output,
        );
        assert!(
            matches!(result, Err(Error::RustHtslibError(_))),
            "a missing directory must surface as the writer error, got {result:?}"
        );

        // The identity check runs before any CRAM is opened, so the same
        // missing directory surfaces as an I/O `NotFound` there.
        let reference = root.join("reference.fa");
        let cram_output = root.join("missing").join("output.cram");
        let error = write_empty_cram(&cram_output, &reference).expect_err("missing directory");
        let Error::InputOutputError(source) = error else {
            panic!("expected a filesystem error, got {error:?}");
        };
        assert_eq!(
            source.kind(),
            std::io::ErrorKind::NotFound,
            "the identity check reports the missing parent as not found"
        );
        assert!(!output.exists(), "no BAM output may be created");
        assert!(
            !root.join("missing").exists(),
            "no directory may be created"
        );
        fs::remove_dir_all(root).expect("remove temporary directory");
    }
}
