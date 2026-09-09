# Third-Party Notices

This project includes code adapted from third-party open source software.

## bedrs

This project includes a minimal adaptation of the published `bedrs` v0.2.26
crate in `src/bedrs.rs`. Unused functionality was removed.

- Crate page: <https://crates.io/crates/bedrs>
- Repository: <https://github.com/noamteyssier/bedrs>
- Reference crate checksum: `e80a9ee52ad2ad5233b261be535926bea4a3fdfd068fd7d4cb48c8edd7518173`
- Reference source commit: `f71d798d5a7dc21a7b731a1de6f0858fc37ee667`
- License: MIT; see `licenses/bedrs-LICENSE`

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

## openssl-probe

This project includes code adapted from the published crate `openssl-probe`
v0.2.1, specifically the system certificate location probing now vendored in
`src/utils/openssl_probe.rs`.

- Crate page: <https://crates.io/crates/openssl-probe>
- Repository: <https://github.com/rustls/openssl-probe>
- License declared by the published crate metadata: MIT OR Apache-2.0
- Copied license: `licenses/openssl-probe-LICENSE`
- Source copied on: 2026-07-29

```text
The MIT License (MIT)

Copyright (c) 2014 Alex Crichton

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

## static_assertions

This project includes a minimal adaptation of the published crate
`static_assertions` v1.1.0, specifically the `assert_not_impl_any!` macro
inlined into the `temp_bam_simulation_ownership` test module in
`src/simulate_mod_bam.rs`. The macro is reproduced verbatim except for removing
the upstream `#[allow(dead_code)]` attribute (this repository denies
`allow_attributes`) and relocating the upstream explanatory inline comments
into the macro's doc comment.

- Crate page: <https://crates.io/crates/static_assertions>
- Repository: <https://github.com/nvzqz/static-assertions-rs>
- License declared by the published crate metadata: MIT OR Apache-2.0
- Copied license: `licenses/static_assertions-LICENSE` (MIT)
- Source copied on: 2026-09-08
- Source file: <https://crates.io/crates/static_assertions/1.1.0/code/src/assert_impl.rs>
- Reference source tag: `v1.1.0`

```text
MIT License

Copyright (c) 2017 Nikolai Vazquez

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
