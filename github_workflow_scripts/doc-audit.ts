#!/usr/bin/env node

const { spawnSync } = require('node:child_process');
const path = require('node:path');

const repoRoot = path.resolve(__dirname, '..');
process.chdir(repoRoot);

const apiKey = process.env.OPENAI_API_KEY;
if (!apiKey) {
  console.error('OPENAI_API_KEY is required');
  process.exit(1);
}

const model = process.env.OPENAI_MODEL || 'gpt-5.4-mini';
const prompt = [
  'You are auditing this repository for documentation problems.',
  'Review README files, Markdown docs, Rust doc comments, and workflow comments for stale claims, broken instructions, mismatched examples, missing warnings, and inconsistencies with the code.',
  'Use the repository files as the source of truth.',
  'Make the necessary documentation edits directly in the working tree.',
  'Keep the changes minimal and do not refactor code unless required to correct documentation.',
  'When you are done editing, stop.',
].join('\n');

const result = spawnSync(
  'npx',
  [
    '-y',
    '-p',
    '@mariozechner/pi-coding-agent@0.72.1',
    'pi',
    '--mode',
    'json',
    '--provider',
    'openai',
    '--model',
    model,
    '--no-session',
    '--no-config',
    '--no-extensions',
    '--no-skills',
    '--no-prompt-templates',
    '--no-themes',
    '--tools',
    'read,edit,grep,find,ls',
    prompt,
  ],
  {
    cwd: repoRoot,
    env: {
      HOME: process.env.HOME || '',
      PATH: process.env.PATH || '',
      OPENAI_API_KEY: apiKey,
      OPENAI_MODEL: model,
    },
    stdio: 'inherit',
    shell: false,
  },
);

if (result.error) {
  console.error(result.error);
  process.exit(1);
}

if (result.status !== 0) {
  process.exit(result.status ?? 1);
}
