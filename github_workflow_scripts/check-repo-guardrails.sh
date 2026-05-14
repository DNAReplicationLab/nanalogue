#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [ -e CLAUDE.md ] && { [ ! -L AGENTS.md ] || [ ! AGENTS.md -ef CLAUDE.md ]; }; then
  echo "pre-commit: AGENTS.md must resolve to CLAUDE.md" >&2
  echo "  Fix with: rm AGENTS.md && ln -s CLAUDE.md AGENTS.md" >&2
  exit 1
fi

if ! git check-ignore -q brainstorming/sentinel 2>/dev/null; then
  echo "pre-commit: brainstorming/ is not covered by .gitignore" >&2
  echo "  Add 'brainstorming*' (or equivalent) to .gitignore." >&2
  exit 1
fi

if git diff --cached --name-only --diff-filter=ACMR | grep -q '^brainstorming/'; then
  echo "pre-commit: staged files found under brainstorming/ — this directory must stay local" >&2
  git diff --cached --name-only --diff-filter=ACMR | grep '^brainstorming/' | sed 's/^/  /' >&2
  exit 1
fi
