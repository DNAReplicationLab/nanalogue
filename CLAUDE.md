After you have finished incorporating code and before showing a final
message to me that you have successfully done everything, you must run
the following commands in this order:
- `cargo clippy -q --all-features --all-targets -- -D warnings` to make sure you pass clippy linting tests
- `cargo test -q` to make sure you pass cargo tests.
- `cargo fmt` to format the code in the rust style.

## Final review
- Run `codex review --uncommitted` with an unlimited timeout to see what the tool
says as code review and incorporate its changes.
- If you repeat codex review, deal with its previous comments first, then run it
again. Do not run multiple instances in parallel.
- Keep repeating the cycle of running the code agent above and incorporating its
changes if you think they are worth it, until it stops complaining or the
remaining complaints are not worth fixing.
- Then, run `cargo doc` to ensure the docs still form successfully.

# Overall notes

Each commit message's first line must be < 50 characters and must start with a verb
like adds changes extracts i.e. this commit does blah blah, not a verb like add change extract
i.e. not with this commit we do blah blah.

Any docs or plans you make must go into brainstorming/ . you must never commit anything in brainstorming/

Use the quiet version of cargo commands wherever possible.
