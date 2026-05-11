#!/usr/bin/env node

import fs from 'node:fs';
import path from 'node:path';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
process.chdir(repoRoot);

const args = new Set(process.argv.slice(2));
const cachedOnly = args.has('--cached-only');
const explicitFiles = [...process.argv.slice(2).filter((arg) => !arg.startsWith('--'))];
function gitList(commandArgs) {
  try {
    const out = execFileSync('git', commandArgs, { encoding: 'utf8' }).trim();
    return out ? out.split('\n').filter(Boolean) : [];
  } catch {
    return [];
  }
}

function readMarkdown(file) {
  if (cachedOnly) {
    try {
      return execFileSync('git', ['show', `:${file}`], { encoding: 'utf8' });
    } catch {
      // fall through to the working tree for untracked files
    }
  }
  return fs.readFileSync(file, 'utf8');
}

function isExternal(target) {
  return /^[a-z][a-z0-9+.-]*:/i.test(target);
}

function slugifyHeading(text) {
  return text
    .replace(/<[^>]+>/g, '')
    .toLowerCase()
    .replace(/[^\p{L}\p{N}\s-]/gu, '')
    .trim()
    .replace(/\s+/g, '-');
}

function headingsFor(fileText) {
  const headings = new Set();
  const counts = new Map();

  for (const line of fileText.split(/\r?\n/)) {
    const match = /^(#{1,6})\s+(.+?)\s*$/.exec(line);
    if (!match) continue;

    const slug = slugifyHeading(match[2]);
    const count = counts.get(slug) ?? 0;
    counts.set(slug, count + 1);
    headings.add(count === 0 ? slug : `${slug}-${count}`);
  }

  return headings;
}

function parseLinkDestination(rawTarget) {
  const trimmed = rawTarget.trim();
  if (trimmed.startsWith('<')) {
    const closing = trimmed.indexOf('>');
    return closing === -1 ? trimmed.slice(1) : trimmed.slice(1, closing);
  }
  return trimmed.split(/\s+/)[0] ?? '';
}

function resolveTarget(currentFile, rawTarget) {
  const [targetPath, fragment = ''] = rawTarget.split('#');
  let normalizedPath = '';
  if (targetPath) {
    try {
      normalizedPath = decodeURIComponent(targetPath);
    } catch {
      normalizedPath = targetPath;
    }
  }
  const baseDir = path.dirname(currentFile);
  const filePath = !normalizedPath
    ? currentFile
    : normalizedPath.startsWith('/')
      ? path.join(repoRoot, normalizedPath.slice(1))
      : path.resolve(baseDir, normalizedPath);

  return {
    filePath,
    fragment,
    rawPath: normalizedPath,
  };
}

function scanFile(file) {
  const text = readMarkdown(file);
  let inFence = false;
  const errors = [];

  for (const [idx, line] of text.split(/\r?\n/).entries()) {
    if (/^```|^~~~/.test(line.trim())) {
      inFence = !inFence;
      continue;
    }
    if (inFence) continue;

    const linkRegex = /!?\[[^\]]*\]\(([^)]+)\)/g;
    let match;
    while ((match = linkRegex.exec(line)) !== null) {
      const rawTarget = parseLinkDestination(match[1]);
      if (!rawTarget) continue;
      if (isExternal(rawTarget)) continue;

      const { filePath, fragment, rawPath } = resolveTarget(file, rawTarget);
      let targetText;
      try {
        targetText = readMarkdown(filePath);
      } catch {
        errors.push(`${file}:${idx + 1}: missing linked file: ${rawTarget}`);
        continue;
      }

      if (fragment) {
        const targetHeadings = headingsFor(targetText);
        if (!targetHeadings.has(slugifyHeading(fragment))) {
          errors.push(`${file}:${idx + 1}: missing heading anchor: ${rawTarget}`);
        }
      } else if (!rawPath.startsWith('#') && rawPath.endsWith('.md')) {
        // no-op; the file exists and has already been read successfully.
      }
    }
  }

  return errors;
}

const files = explicitFiles.length
  ? explicitFiles
  : cachedOnly
    ? gitList(['diff', '--cached', '--name-only', '--diff-filter=ACMR', '--', '*.md'])
    : gitList(['ls-files', '--', '*.md']);

const uniqueFiles = [...new Set(files.filter((file) => file.endsWith('.md') && fs.existsSync(file)))];
const errors = [];
for (const file of uniqueFiles) {
  errors.push(...scanFile(file));
}

if (errors.length > 0) {
  for (const error of errors) console.error(error);
  process.exit(1);
}
