#!/usr/bin/env node

const { spawnSync } = require('node:child_process');
const fs = require('node:fs');
const path = require('node:path');

const repoRoot = path.resolve(__dirname, '..');
process.chdir(repoRoot);

// Verify this npm build recognizes the min-release-age config key before
// relying on it for the npx package install cooldown, and reject any
// inherited min-release-age setting so this script controls the value.
function ensureSupportedNpmConfig() {
  const minReleaseAgeResult = spawnSync('npm', ['config', 'get', 'min-release-age'], {
    encoding: 'utf8',
    stdio: ['ignore', 'pipe', 'inherit'],
    shell: false,
  });

  if (minReleaseAgeResult.error) {
    console.error(minReleaseAgeResult.error);
    process.exit(1);
  }

  if (minReleaseAgeResult.status !== 0) {
    process.exit(minReleaseAgeResult.status ?? 1);
  }

  if (minReleaseAgeResult.stdout.trim() !== 'null') {
    console.error(
      `npm must support min-release-age config; expected ${JSON.stringify('null')}, got ${JSON.stringify(minReleaseAgeResult.stdout.trim())}`,
    );
    process.exit(1);
  }
}

ensureSupportedNpmConfig();

const apiKey = process.env.OPENAI_API_KEY;
if (!apiKey) {
  console.error('OPENAI_API_KEY is required');
  process.exit(1);
}

const model = process.env.OPENAI_MODEL || 'gpt-5.4-mini';
const minimumReleaseAge = '10080';
const outputPath = '/tmp/doc-audit.jsonl';
if (fs.existsSync(outputPath)) {
  console.error(`Refusing to overwrite existing JSONL log: ${outputPath}`);
  process.exit(1);
}
const outputFd = fs.openSync(outputPath, 'w');
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
    '@mariozechner/pi-coding-agent@0.70.2',
    'pi',
    '--mode',
    'json',
    '--provider',
    'openai',
    '--model',
    model,
    '--no-session',
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
      npm_config_min_release_age: minimumReleaseAge,
    },
    stdio: ['ignore', outputFd, 'inherit'],
    shell: false,
  },
);

fs.closeSync(outputFd);

if (result.error) {
  console.error(result.error);
  process.exit(1);
}

if (result.status !== 0) {
  process.exit(result.status ?? 1);
}
