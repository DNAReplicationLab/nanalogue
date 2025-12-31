# Contributing to nanalogue

Thank you for your interest in contributing to `nanalogue`! This project is largely maintained by one developer, so your contributions are especially valuable.

## Getting Started

1. Fork the repository
2. Clone your fork: `git clone https://github.com/YOUR_USERNAME/nanalogue.git`
3. Create a branch for your changes: `git checkout -b feature/your-feature-name`

## Development Setup

### Prerequisites

- Rust toolchain (rustc, cargo) - Install from [rustup.rs](https://rustup.rs/)
- Git for version control

### Building the Project

```bash
cargo build
```

### Running Tests

```bash
cargo test
```

## Code Quality Standards

This project maintains high code quality standards. Before submitting a PR, ensure:

1. **All tests pass**: Run `cargo test`
2. **Clippy lints pass**: Run `cargo clippy --all-features --all-targets -- -D warnings`
3. **Code is formatted**: Run `cargo fmt`
4. **Code coverage is maintained**: We aim for >92% test coverage

The project uses extensive clippy lints (see `Cargo.toml` for the full list). Your code must pass all linting checks.

## Contribution Guidelines

### Code Style

- Follow Rust naming conventions and idioms
- Add documentation comments (`///`) for public APIs
- Include inline comments for complex logic
- All code files should start with a brief 2-line comment explaining what the file does
- Match the style and formatting of surrounding code for consistency

### Testing

- Write tests for new functionality
- Update existing tests if behavior changes
- Include both unit tests and integration tests where appropriate
- Test edge cases and error conditions

### Documentation

- Update README.md if adding new commands or features
- Add docstrings to new public functions, structs, and modules
- Update CHANGELOG.md with your changes

### Commit Messages

- Write clear, descriptive commit messages
- Use present tense ("Add feature" not "Added feature")
- Reference issue numbers when applicable

## Pull Request Process

1. **Before submitting**:
   - Ensure all tests pass
   - Run clippy and fix all warnings
   - Run `cargo fmt` to format your code
   - Update documentation as needed

2. **Submit your PR**:
   - Provide a clear description of the changes
   - Reference any related issues
   - Explain why the change is needed

3. **Review process**:
   - The maintainer will review your PR when available
   - Be patient - this is a solo-maintained project
   - Be open to feedback and suggested changes
   - Once approved, your PR will be merged

## Types of Contributions

We welcome various types of contributions:

- **Bug fixes**: Found a bug? Fix it and submit a PR!
- **New features**: Have an idea? Open an issue first to discuss
- **Documentation**: Improvements to docs, examples, or comments
- **Tests**: Additional test coverage is always appreciated
- **Performance improvements**: Optimizations with benchmarks
- **Code quality**: Refactoring to improve maintainability

## Reporting Bugs

When reporting bugs, please include:

- A clear description of the issue
- Steps to reproduce the problem
- Expected vs actual behavior
- Your environment (OS, Rust version, etc.)
- Sample BAM files or data (if applicable and not sensitive)

**Note**: For security vulnerabilities, see SECURITY.md

## Feature Requests

For feature requests:

- Open an issue describing the feature
- Explain the use case and why it's valuable
- Discuss implementation approach if you have ideas
- Wait for maintainer feedback before starting work

## Questions?

- Open an issue for questions about the codebase
- Check existing issues and documentation first
- Be respectful and patient

## Code of Conduct

- Be respectful and inclusive
- Focus on constructive feedback
- Help create a welcoming environment for all contributors

## License

By contributing to nanalogue, you agree that your contributions will be licensed under the MIT License.

## Contact

Maintainer: Sathish Thiyagarajan (mail AT unintegrable dot com)

## Acknowledgments

All contributors will be acknowledged in the project. Thank you for helping make `nanalogue` better!
