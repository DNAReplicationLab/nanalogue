#!/usr/bin/env node

const fs = require('node:fs');
const path = require('node:path');

const cargoTomlPath = path.resolve(__dirname, '..', 'Cargo.toml');
const content = fs.readFileSync(cargoTomlPath, 'utf8');
const lines = content.split(/\r?\n/u);

const sectionHeaderPattern = /^\s*\[([^\]]+)\]\s*$/u;
const lintSectionPattern = /^lints\.[^\]]+$/u;
const keyValuePattern = /^\s*([A-Za-z0-9_-]+)\s*=/u;

/** @type {{sectionName: string, entries: Array<{key: string, lineNumber: number}>} | null} */
let currentSection = null;
let failed = false;

function compareKeys(left, right) {
  return left < right ? -1 : left > right ? 1 : 0;
}

function validateCurrentSection() {
  if (currentSection === null) {
    return;
  }

  const keys = currentSection.entries.map((entry) => entry.key);
  const sortedKeys = [...keys].sort(compareKeys);

  for (let index = 0; index < keys.length; index += 1) {
    if (keys[index] !== sortedKeys[index]) {
      failed = true;
      console.error(
        `${cargoTomlPath}: [${currentSection.sectionName}] lint keys must be alphabetic.`,
      );
      console.error(`  First out-of-order key: ${JSON.stringify(keys[index])} on line ${currentSection.entries[index].lineNumber}`);
      console.error(`  Expected order: ${sortedKeys.join(', ')}`);
      return;
    }
  }
}

for (let index = 0; index < lines.length; index += 1) {
  const line = lines[index];
  const sectionMatch = line.match(sectionHeaderPattern);
  if (sectionMatch) {
    validateCurrentSection();

    const sectionName = sectionMatch[1];
    currentSection = lintSectionPattern.test(sectionName)
      ? { sectionName, entries: [] }
      : null;
    continue;
  }

  if (currentSection === null) {
    continue;
  }

  if (/^\s*(#.*)?$/u.test(line)) {
    continue;
  }

  const keyMatch = line.match(keyValuePattern);
  if (keyMatch) {
    currentSection.entries.push({ key: keyMatch[1], lineNumber: index + 1 });
  }
}

validateCurrentSection();

if (failed) {
  process.exit(1);
}
