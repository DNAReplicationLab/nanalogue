# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 0.x     | :x: Alpha          |

**Note**: This project is currently in alpha development.
No versions are officially supported at this time.
Use at your own risk.

## Reporting a Vulnerability

If you discover a security vulnerability in `nanalogue`, please report it responsibly:

1. **Do NOT open a public issue** for security vulnerabilities
2. Email the maintainer directly at: mail AT unintegrable dot com
3. Include:
   - Description of the vulnerability
   - Steps to reproduce
   - Potential impact
   - Any suggested fixes (optional)

We aim to respond to security reports as soon as possible.

## Security Considerations

### Data Sensitivity

This tool processes BAM files which may contain:
- Sensitive genetic/genomic data
- Personal health information
- Research data subject to confidentiality agreements

**Users are responsible for**:
- Ensuring proper access controls on input/output files
- Complying with data protection regulations (GDPR, HIPAA, etc.)
- Maintaining appropriate file permissions in multi-user environments

### File Access

When processing BAM files:
- Ensure proper file permissions are set
- Be aware of where output files are written
- Consider the privacy implications of processing sensitive genomic data
- Verify integrity of input files from untrusted sources

### Dependencies

This package relies on:
- `rust-htslib` - for BAM file I/O
- `fibertools-rs` - for mod BAM parsing
- `bio` and `bio-types` - for bioinformatics operations
- Various other Rust crates for data processing

For a complete list of dependencies, see `Cargo.toml`.
Security of these dependencies is maintained upstream.
Users should keep dependencies updated.

## Disclosure Policy

Once a security issue is fixed:
- We will release a patch as soon as possible
- A security advisory will be published
- Credit will be given to the reporter (unless anonymity is requested)

## Contact

For security concerns, contact: mail AT unintegrable dot com

For general bugs and features, use the GitHub issue tracker.
