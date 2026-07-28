#!/usr/bin/env node

// Scans tracked and untracked non-ignored UTF-8 text files for characters
// outside printable ASCII. Horizontal tabs and line feeds are allowed;
// carriage returns and all other control or non-ASCII characters are reported.
// Pass --cached-only to scan staged Git blobs in a pre-commit hook. Exits 1
// when findings are present and 0 otherwise.

const { execFileSync } = require("node:child_process");
const { readFileSync } = require("node:fs");
const { resolve } = require("node:path");
const { argv } = require("node:process");

const ROOT = resolve(__dirname, "..");
const UTF8_DECODER = new TextDecoder("utf-8", { fatal: true });
const CACHED_ONLY = argv.includes("--cached-only");
const MAX_GIT_BLOB_BYTES = 64 * 1024 * 1024;

/** Visible Unicode symbols intentionally permitted in repository text. */
const ALLOWED_UNICODE_CODE_POINTS = new Set([
    0x00a9, // copyright sign
    0x00ae, // registered sign
    0x00b0, // degree sign
    0x00b2, // superscript two
    0x00b7, // middle dot
    0x00d7, // multiplication sign
    0x00e9, // latin small letter e with acute
    0x00f6, // latin small letter o with diaeresis
    0x2013, // en dash
    0x2014, // em dash
    0x2018, // left single quotation mark
    0x2019, // right single quotation mark
    0x201c, // left double quotation mark
    0x201d, // right double quotation mark
    0x2022, // bullet
    0x2026, // ellipsis
    0x20ac, // euro sign
    0x2122, // trademark sign
    0x2248, // almost equal to
    0x2264, // less-than or equal to
    0x2265, // greater-than or equal to
    0x25b6, // black right-pointing triangle
    0x25b8, // black right-pointing small triangle
    0x25ba, // black right-pointing pointer
    0x25bc, // black down-pointing triangle
    0x25c0, // black left-pointing triangle
    0x26a0, // warning sign
    0x2713, // check mark
    0x2705, // white heavy check mark
    0x2717, // ballot x
    0x754c, // cjk ideograph meaning "world/boundary"
    0x1c20, // lepcha letter la
    0x1f389, // party popper
    0x1f49a, // green heart
    0x1f4f8, // camera with flash
    0x1f50d, // left-pointing magnifying glass
    0x1f9ec, // DNA double helix
]);

/** Inclusive ranges of visible Unicode symbols allowed in repository text. */
const ALLOWED_UNICODE_RANGES = [
    [0x2190, 0x2193], // common directional arrows
    [0x2500, 0x257f], // box drawing
];

/** Returns repository files that could be included by a future git add. */
function listRepositoryFiles() {
    const args = CACHED_ONLY
        ? ["ls-files", "--cached", "-z"]
        : ["ls-files", "--cached", "--others", "--exclude-standard", "-z"];
    const output = execFileSync("git", args, { cwd: ROOT });
    return output
        .toString("utf8")
        .split("\0")
        .filter((file) => file.length > 0);
}

/** Decodes a UTF-8 text file, returning null for likely binary content. */
function readTextFile(file) {
    const content = CACHED_ONLY
        ? execFileSync("git", ["show", `:${file}`], {
              cwd: ROOT,
              maxBuffer: MAX_GIT_BLOB_BYTES,
          })
        : readFileSync(resolve(ROOT, file));
    if (content.includes(0)) return null;
    try {
        return UTF8_DECODER.decode(content);
    } catch {
        return null;
    }
}

/** Returns whether a Unicode code point is allowed in repository text. */
function isAllowedCodePoint(codePoint) {
    return (
        codePoint === 9 ||
        codePoint === 10 ||
        (codePoint >= 32 && codePoint <= 126) ||
        ALLOWED_UNICODE_CODE_POINTS.has(codePoint) ||
        ALLOWED_UNICODE_RANGES.some(
            ([start, end]) => codePoint >= start && codePoint <= end,
        )
    );
}

/** Formats a Unicode code point for searchable output. */
function formatCodePoint(codePoint) {
    return `U+${codePoint.toString(16).toUpperCase().padStart(4, "0")}`;
}

const findings = [];
let scannedFiles = 0;
let skippedBinaryFiles = 0;

for (const file of listRepositoryFiles()) {
    const content = readTextFile(file);
    if (content === null) {
        skippedBinaryFiles++;
        continue;
    }
    scannedFiles++;

    let line = 1;
    let column = 1;
    for (const character of content) {
        const codePoint = character.codePointAt(0);
        if (codePoint === undefined) continue;
        if (!isAllowedCodePoint(codePoint)) {
            findings.push({
                file,
                line,
                column,
                codePoint: formatCodePoint(codePoint),
            });
        }
        if (codePoint === 10) {
            line++;
            column = 1;
        } else {
            column++;
        }
    }
}

if (findings.length > 0) {
    console.error("Strange characters found:\n");
    for (const finding of findings) {
        console.error(
            `  ${finding.file}:${finding.line}:${finding.column} (${finding.codePoint})`,
        );
    }
    console.error(
        `\nFound ${findings.length} disallowed character(s) across ` +
            `${scannedFiles} text file(s); ${skippedBinaryFiles} binary file(s) skipped.`,
    );
    process.exitCode = 1;
} else {
    console.log(
        `No strange characters found across ${scannedFiles} text file(s) ` +
            `(${skippedBinaryFiles} binary file(s) skipped).`,
    );
}
