# Third-Party Notices

This project includes code adapted from third-party open source software.

## Vendored HTSlib-related crates

Two crates are vendored locally to keep the build reproducible and to avoid
regenerating incompatible HTSlib bindings with the current toolchain. The
exact deltas from upstream are recorded in `third_party/patches/`:

- `third_party/patches/hts-sys-2.2.0.patch`
- `third_party/patches/rust-htslib-1.0.0.patch`

The vendored crates are:

- `hts-sys` v2.2.0
  - Local path: `vendor/hts-sys-2.2.0/`
- `rust-htslib` v1.0.0
  - Local path: `vendor/rust-htslib-1.0.0/`

These crates are vendored because the current Rust/clang/libclang environment
produces bindings that are not compatible with upstream `rust-htslib` as
published on crates.io. The vendored copies pin the known-good prebuilt
bindings path and the small compatibility fixes required for this repository.

The respective license files are stored under each vendored crate directory:

- `vendor/hts-sys-2.2.0/LICENSE.md`
- `vendor/rust-htslib-1.0.0/LICENSE.md`

## fibertools-rs

Portions of this project are adapted from the published crate `fibertools-rs` v0.8.2.

- Crate page: <https://crates.io/crates/fibertools-rs>
- Repository referenced by the crate: <https://github.com/fiberseq/fibertools-rs>
- License declared by the published crate metadata: MIT

This notice reflects the license declared in the published crate metadata for `fibertools-rs` v0.8.2.

```text
MIT License

Copyright (c) the fibertools-rs authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## rust-bio

This project includes code adapted from the published crate `bio` v3.0.0,
specifically the DNA complement and reverse-complement helpers now vendored in
`src/utils/complement.rs`.

- Crate page: <https://crates.io/crates/bio>
- Repository: <https://github.com/rust-bio/rust-bio>
- License declared by the published crate metadata: MIT
- Source copied on: 2026-04-08
- Source commit referenced by user: `aec47df`

```text
The MIT License (MIT)

Copyright (c) 2016 Johannes Köster, the Rust-Bio team, Google Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
