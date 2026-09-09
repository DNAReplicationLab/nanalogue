#!/usr/bin/env node

// Scans tracked and untracked non-ignored text files for credentials that may
// have been committed accidentally. Reports locations without printing secret
// values. Pass --cached-only to scan staged Git blobs in a pre-commit hook.
// Exits 1 when findings are present and 0 otherwise.

const { execFileSync } = require("node:child_process");
const { readFileSync } = require("node:fs");
const { resolve } = require("node:path");
const { argv } = require("node:process");

const ROOT = resolve(__dirname, "..");
const UTF8_DECODER = new TextDecoder("utf-8", { fatal: true });
const CACHED_ONLY = argv.includes("--cached-only");
const MAX_GIT_BLOB_BYTES = 64 * 1024 * 1024;

const SECRET_PATTERNS = [
    { name: "AWS access key ID", expression: /\b(?:AKIA|ASIA)[A-Z0-9]{16}\b/g },
    {
        name: "GitHub token",
        expression: /\bgh[pousr]_[A-Za-z0-9]{36,255}\b/g,
    },
    {
        name: "GitHub fine-grained token",
        expression: /\bgithub_pat_[A-Za-z0-9_]{40,255}\b/g,
    },
    { name: "GitLab token", expression: /\bglpat-[A-Za-z0-9_-]{20,}\b/g },
    { name: "npm token", expression: /\bnpm_[A-Za-z0-9]{36}\b/g },
    {
        name: "PyPI token",
        expression: /\bpypi-AgEIcHlwaS5vcmc[A-Za-z0-9_-]{20,}\b/g,
    },
    { name: "Google API key", expression: /\bAIza[A-Za-z0-9_-]{35}\b/g },
    {
        name: "OpenAI API key",
        expression: /\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{20,}\b/g,
    },
    {
        name: "Anthropic API key",
        expression: /\bsk-ant-[A-Za-z0-9_-]{20,}\b/g,
    },
    {
        name: "Slack token",
        expression: /\bxox[baprs]-[A-Za-z0-9-]{20,}\b/g,
    },
    {
        name: "Stripe live secret key",
        expression: /\bsk_live_[A-Za-z0-9]{16,}\b/g,
    },
];

const PRIVATE_KEY_MARKER = new RegExp(
    [
        "-----BEGIN ",
        "(?:(?:RSA|DSA|EC|OPENSSH|ENCRYPTED) )?PRIVATE KEY-----",
        "|-----BEGIN PGP PRIVATE",
        " KEY BLOCK-----",
    ].join(""),
    "g",
);
const EMBEDDED_URL_CREDENTIALS = /\bhttps?:\/\/[^/\s:@]+:[^@\s/]+@[^\s"'`<>]+/g;
const ALLOWED_LOCAL_TEST_USERNAME = "example-user-test";
const ALLOWED_LOCAL_TEST_PASSWORD = "example-password-test";
const CONTEXTUAL_SECRET =
    /(["']?)\b(?:api[_-]?key|access[_-]?key|client[_-]?secret|private[_-]?key|password|passwd|secret|token)\b\1\s*[:=]\s*(?:"([^"\n]+)"|'([^'\n]+)'|([A-Za-z0-9_+/@:=-]+))/gi;
const PLACEHOLDER_VALUES = new Set([
    "changeme",
    "dummy",
    "example",
    "fake",
    "placeholder",
    "redacted",
    "sample",
    "test",
    ALLOWED_LOCAL_TEST_PASSWORD,
]);
const DESCRIPTIVE_PLACEHOLDER =
    /^(?:your|example|test)[-_](?:api[-_]?key|access[-_]?key|client[-_]?secret|private[-_]?key|password|passwd|secret|token)(?:[-_]here)?$/;

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

/** Returns the one-based line containing a character offset. */
function lineAt(content, offset) {
    let line = 1;
    for (let index = 0; index < offset; index++) {
        if (content.charCodeAt(index) === 10) line++;
    }
    return line;
}

/** Estimates Shannon entropy in bits per character. */
function calculateEntropy(value) {
    const counts = new Map();
    for (const character of value) {
        counts.set(character, (counts.get(character) ?? 0) + 1);
    }
    let entropy = 0;
    for (const count of counts.values()) {
        const probability = count / value.length;
        entropy -= probability * Math.log2(probability);
    }
    return entropy;
}

/** Returns whether a context-matched value resembles a real credential. */
function resemblesSecret(value) {
    const normalized = value.trim().toLowerCase();
    if (value.trim().length < 12) return false;
    if (
        normalized.startsWith("${") ||
        normalized.startsWith("process.env") ||
        normalized.startsWith("<") ||
        PLACEHOLDER_VALUES.has(normalized) ||
        DESCRIPTIVE_PLACEHOLDER.test(normalized)
    ) {
        return false;
    }
    if (/^(.)\1+$/.test(value.trim())) return false;
    return calculateEntropy(value.trim()) >= 3;
}

/** Returns whether a URL uses the one approved localhost test credential pair. */
function isAllowedLocalTestCredentialUrl(value) {
    try {
        const parsed = new URL(value);
        return (
            parsed.hostname === "localhost" &&
            parsed.username === ALLOWED_LOCAL_TEST_USERNAME &&
            parsed.password === ALLOWED_LOCAL_TEST_PASSWORD
        );
    } catch {
        return false;
    }
}

/** Adds a finding once per detector and source line. */
function addFinding(findings, seen, finding) {
    const key = `${finding.file}:${finding.line}:${finding.detector}`;
    if (seen.has(key)) return;
    seen.add(key);
    findings.push(finding);
}

const findings = [];
const seen = new Set();
let scannedFiles = 0;
let skippedBinaryFiles = 0;

for (const file of listRepositoryFiles()) {
    const content = readTextFile(file);
    if (content === null) {
        skippedBinaryFiles++;
        continue;
    }
    scannedFiles++;

    for (const pattern of SECRET_PATTERNS) {
        pattern.expression.lastIndex = 0;
        for (const match of content.matchAll(pattern.expression)) {
            addFinding(findings, seen, {
                file,
                line: lineAt(content, match.index),
                detector: pattern.name,
            });
        }
    }

    EMBEDDED_URL_CREDENTIALS.lastIndex = 0;
    for (const match of content.matchAll(EMBEDDED_URL_CREDENTIALS)) {
        if (isAllowedLocalTestCredentialUrl(match[0])) continue;
        addFinding(findings, seen, {
            file,
            line: lineAt(content, match.index),
            detector: "credential in URL",
        });
    }

    PRIVATE_KEY_MARKER.lastIndex = 0;
    for (const match of content.matchAll(PRIVATE_KEY_MARKER)) {
        addFinding(findings, seen, {
            file,
            line: lineAt(content, match.index),
            detector: "private key block",
        });
    }

    CONTEXTUAL_SECRET.lastIndex = 0;
    for (const match of content.matchAll(CONTEXTUAL_SECRET)) {
        const value = match[2] ?? match[3] ?? match[4];
        if (!resemblesSecret(value)) continue;
        addFinding(findings, seen, {
            file,
            line: lineAt(content, match.index),
            detector: "credential-like assignment",
        });
    }
}

if (findings.length > 0) {
    console.error("Possible committed secrets found:\n");
    for (const finding of findings) {
        console.error(
            `  ${finding.file}:${finding.line} (${finding.detector})`,
        );
    }
    console.error(
        `\nFound ${findings.length} possible secret(s). Values are intentionally redacted.`,
    );
    process.exitCode = 1;
} else {
    console.log(
        `No possible secrets found across ${scannedFiles} text file(s) ` +
            `(${skippedBinaryFiles} binary file(s) skipped).`,
    );
}
