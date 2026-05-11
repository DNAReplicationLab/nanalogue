#!/usr/bin/env node

import fs from 'node:fs';
import path from 'node:path';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
process.chdir(repoRoot);

const cachedOnly = process.argv.includes('--cached-only');
const changelogPath = 'CHANGELOG.md';

function isStaged(file) {
  try {
    const out = execFileSync('git', ['diff', '--cached', '--name-only', '--diff-filter=ACMR', '--', file], { encoding: 'utf8' }).trim();
    return out.length > 0;
  } catch {
    return false;
  }
}

function readChangelog() {
  if (cachedOnly) {
    if (!isStaged(changelogPath)) {
      return null;
    }
    return execFileSync('git', ['show', `:${changelogPath}`], { encoding: 'utf8' });
  }
  if (!fs.existsSync(changelogPath)) {
    return null;
  }
  return fs.readFileSync(changelogPath, 'utf8');
}

const text = readChangelog();
if (text === null) {
  process.exit(0);
}
const lines = text.split(/\r?\n/);
const errors = [];

const firstNonEmpty = lines.find((line) => line.trim().length > 0);
if (firstNonEmpty !== '# Changelog') {
  errors.push('CHANGELOG.md must start with "# Changelog"');
}

let sawUnreleased = false;
const allowedSections = new Set(['Added', 'Changed', 'Deprecated', 'Removed', 'Fixed', 'Security']);

for (const [idx, line] of lines.entries()) {
  if (line.startsWith('## ')) {
    if (line === '## [Unreleased]') {
      sawUnreleased = true;
      continue;
    }

    if (!/^## \[[^\]]+\] - \d{4}-\d{2}-\d{2}$/.test(line)) {
      errors.push(`${changelogPath}:${idx + 1}: invalid version heading: ${line}`);
    }
    continue;
  }

  if (line.startsWith('### ')) {
    const section = line.slice(4).trim();
    if (!allowedSections.has(section)) {
      errors.push(`${changelogPath}:${idx + 1}: invalid section heading: ${line}`);
    }
  }
}

if (!sawUnreleased) {
  errors.push('CHANGELOG.md must contain a "## [Unreleased]" section');
}

if (errors.length > 0) {
  for (const error of errors) console.error(error);
  process.exit(1);
}
