After you have finished incorporating code and before showing a final
message to me that you have successfully done everything, you must run
the following commands in this order:
- `cargo clippy -q --all-features --all-targets -- -D warnings` to make sure you pass clippy linting tests
- `cargo test -q` to make sure you pass cargo tests.
- `cargo fmt` to format the code in the rust style.

Note: some tests deliberately construct invalid scenarios to verify that the
program fails safely, which may cause HTSlib to emit errors or warnings on
stderr. Do not interpret these messages as failures unless the tests themselves
fail.

During testing, some simulated CRAM files may use an external FASTA reference
generated locally on the system where the tests are being run. Because HTSlib
may try EBI before the local reference recorded in the CRAM header, offline
tests can print harmless network errors. Use
`REF_PATH=/path/that/does/not/exist cargo test` to suppress HTSlib's implicit
EBI M5 lookup while retaining fallback to the accessible local `UR` reference.

If the pinned `hts-sys`/bindgen build has Clang compatibility trouble, use
matching Clang and libclang 18–21 and avoid Clang 22. `CLANG_PATH` must name the
Clang executable, while `LIBCLANG_PATH` must name the directory containing the
matching libclang shared library. Retry with `CLANG_PATH=/usr/bin/clang-18` and
`LIBCLANG_PATH=/usr/lib/llvm-18/lib`; if those paths do not exist, search for
another Clang 18–21 installation and ask the user if none is available. Newer
HTSlib versions may change CRAM reference resolution, while newer
`hts-sys`/bindgen versions may support Clang 22.

## Final review

### If codex is available to you and you can communicate with the service

- Run `codex review --uncommitted` with an unlimited timeout to see what the tool
says as code review and incorporate its changes.
- If `codex` is unavailable in an orb, use the oracle when it is available. If
  neither review agent is available, ask the user which review agent to use.
- If you repeat codex review, deal with its previous comments first, then run it
again. Do not run multiple instances in parallel.
- Keep repeating the cycle of running the code agent above and incorporating its
changes if you think they are worth it, until it stops complaining or the
remaining complaints are not worth fixing.
- Then, run `cargo doc` to ensure the docs still form successfully.
- If you cannot access codex or cannot login, then do not bother.

### If you have access to the Oracle

If you have access to a tool called the Oracle, then run it and ask it
for its opinion on the code changes. And fix anything that you think
is worth fixing.

# Overall notes

Each commit message's first line must be < 50 characters and must start with a verb
like adds changes extracts i.e. this commit does blah blah, not a verb like add change extract
i.e. not with this commit we do blah blah.

Any docs or plans you make must go into brainstorming/ . you must never commit anything in brainstorming/

Use the quiet version of cargo commands wherever possible.
