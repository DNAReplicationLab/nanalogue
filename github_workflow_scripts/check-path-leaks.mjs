#!/usr/bin/env node

import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
process.chdir(repoRoot);

const cachedOnly = process.argv.includes('--cached-only');
const explicitFiles = process.argv.slice(2).filter((arg) => !arg.startsWith('--'));

let username = process.env.USER;
if (!username) {
  try {
    username = os.userInfo().username;
  } catch {
    username = '';
  }
}

if (!username) {
  console.warn('check-path-leaks: unable to determine a username; skipping scan');
  process.exit(0);
}

const needle = `/${username}`;

function stagedFiles() {
  try {
    const out = execFileSync('git', ['diff', '--cached', '--name-only', '--diff-filter=ACMR'], { encoding: 'utf8' }).trim();
    return out ? out.split('\n').filter(Boolean) : [];
  } catch {
    return [];
  }
}

function readCandidate(file) {
  try {
    if (cachedOnly) {
      return execFileSync('git', ['show', `:${file}`], { encoding: 'buffer' });
    }
    return fs.readFileSync(file);
  } catch {
    return null;
  }
}

const files = explicitFiles.length
  ? explicitFiles
  : cachedOnly
    ? stagedFiles()
    : execFileSync('git', ['ls-files', '--cached', '--others', '--exclude-standard'], { encoding: 'utf8' }).trim().split('\n').filter(Boolean);

const uniqueFiles = [...new Set(files)];
const errors = [];
const scopeLabel = cachedOnly ? 'staged files' : 'tracked files';

for (const file of uniqueFiles) {
  const blob = readCandidate(file);
  if (!blob) continue;
  if (blob.includes(0)) continue;

  const text = blob.toString('utf8');
  if (!text.includes(needle)) continue;

  const matches = [];
  for (const [idx, line] of text.split(/\r?\n/).entries()) {
    let searchFrom = 0;
    while (true) {
      const pos = line.indexOf(needle, searchFrom);
      if (pos === -1) break;
      matches.push(`${file}:${idx + 1}: ${line}`);
      searchFrom = pos + needle.length;
    }
  }

  if (matches.length > 0) {
    errors.push(`pre-commit: blocked — ${scopeLabel} contain a hardcoded path with '${username}':`);
    errors.push(...matches.map((match) => `  ${match}`));
  }
}

if (errors.length > 0) {
  for (const error of errors) console.error(error);
  console.error('');
  console.error('  Remove or replace any hardcoded paths containing your local username before committing.');
  process.exit(1);
}
