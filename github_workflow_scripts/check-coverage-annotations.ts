#!/usr/bin/env node

// Ensures test modules are excluded from nightly coverage instrumentation.
// Pass --cached-only to validate staged Git blobs from the pre-commit hook.

const { execFileSync } = require("node:child_process");
const { readFileSync } = require("node:fs");
const { resolve } = require("node:path");
const { argv } = require("node:process");

const ROOT = resolve(__dirname, "..");
const CACHED_ONLY = argv.includes("--cached-only");
const TEST_CFG = "#[cfg(test)]";
const COVERAGE_OFF = "#[cfg_attr(coverage_nightly, coverage(off))]";
const COVERAGE_FEATURE =
    "#![cfg_attr(coverage_nightly, feature(coverage_attribute))]";
const MAX_DISTANCE = 5;
const CRATE_ROOT_PREFIX_LINES = 5;

/** Returns Rust files in the selected repository view. */
function listRustFiles() {
    const args = CACHED_ONLY
        ? ["ls-files", "--cached", "-z", "--", "*.rs"]
        : ["ls-files", "--cached", "--others", "--exclude-standard", "-z", "--", "*.rs"];
    return execFileSync("git", args, { cwd: ROOT })
        .toString("utf8")
        .split("\0")
        .filter((file) => file.length > 0);
}

/** Reads a file from the selected repository view. */
function readFile(file) {
    if (CACHED_ONLY) {
        return execFileSync("git", ["show", `:${file}`], { cwd: ROOT }).toString("utf8");
    }
    return readFileSync(resolve(ROOT, file), "utf8");
}

/** Replaces comments and strings with spaces while preserving line numbers. */
function codeOnly(content) {
    let output = "";
    let blockCommentDepth = 0;
    let inString = false;
    let rawStringTerminator = null;

    for (let index = 0; index < content.length; index += 1) {
        const character = content[index];
        const next = content[index + 1];

        if (blockCommentDepth > 0) {
            if (character === "/" && next === "*") {
                blockCommentDepth += 1;
                output += "  ";
                index += 1;
            } else if (character === "*" && next === "/") {
                blockCommentDepth -= 1;
                output += "  ";
                index += 1;
            } else {
                output += character === "\n" ? "\n" : " ";
            }
            continue;
        }

        if (rawStringTerminator !== null) {
            if (content.startsWith(rawStringTerminator, index)) {
                output += " ".repeat(rawStringTerminator.length);
                index += rawStringTerminator.length - 1;
                rawStringTerminator = null;
            } else {
                output += character === "\n" ? "\n" : " ";
            }
            continue;
        }

        if (inString) {
            if (character === "\\") {
                output += " ";
                if (next !== undefined) {
                    output += next === "\n" ? "\n" : " ";
                    index += 1;
                }
            } else if (character === '"') {
                output += " ";
                inString = false;
            } else {
                output += character === "\n" ? "\n" : " ";
            }
            continue;
        }

        if (character === "/" && next === "/") {
            const newline = content.indexOf("\n", index);
            if (newline === -1) {
                output += " ".repeat(content.length - index);
                break;
            }
            output += " ".repeat(newline - index);
            index = newline - 1;
            continue;
        }
        if (character === "/" && next === "*") {
            blockCommentDepth = 1;
            output += "  ";
            index += 1;
            continue;
        }

        const characterLiteralMatch = content
            .slice(index)
            .match(/^(?:b)?'(?:\\x[0-9A-Fa-f]{2}|\\u\{[0-9A-Fa-f_]+\}|\\.|[^'\\\r\n])'/u);
        if (characterLiteralMatch !== null) {
            output += " ".repeat(characterLiteralMatch[0].length);
            index += characterLiteralMatch[0].length - 1;
            continue;
        }

        const rawStringMatch = content
            .slice(index)
            .match(/^(?:br|r)(#*)"/u);
        if (rawStringMatch !== null) {
            output += " ".repeat(rawStringMatch[0].length);
            index += rawStringMatch[0].length - 1;
            rawStringTerminator = `"${rawStringMatch[1]}`;
            continue;
        }
        if (character === '"') {
            output += " ";
            inString = true;
            continue;
        }

        output += character;
    }

    return output;
}

/** Returns whether an active-code line starts with an exact attribute. */
function hasAttribute(line, attribute) {
    const trimmed = line.trimStart();
    if (!trimmed.startsWith(attribute)) return false;
    const following = trimmed.slice(attribute.length);
    return following.length === 0 || /^\s/u.test(following);
}

/** Returns whether a path is a conventional Cargo target root. */
function isCrateRoot(file) {
    return (
        file === "src/lib.rs" ||
        file === "src/main.rs" ||
        /^src\/bin\/(?:[^/]+\.rs|[^/]+\/main\.rs)$/u.test(file) ||
        /^(?:tests|benches|examples)\/(?:[^/]+\.rs|[^/]+\/main\.rs)$/u.test(file)
    );
}

const rustFiles = new Set(listRustFiles());
const contents = new Map();
const errors = [];

for (const file of rustFiles) {
    const lines = codeOnly(readFile(file)).split(/\r?\n/u);
    contents.set(file, lines);

    for (let index = 0; index < lines.length; index += 1) {
        if (!hasAttribute(lines[index], TEST_CFG)) continue;

        const followingLines = lines.slice(index + 1, index + 1 + MAX_DISTANCE);
        if (!followingLines.some((line) => hasAttribute(line, COVERAGE_OFF))) {
            errors.push(
                `${file}:${index + 1}: ${COVERAGE_OFF} must appear within ` +
                    `${MAX_DISTANCE} lines after ${TEST_CFG}`,
            );
        }
    }
}

const crateRoots = [...rustFiles].filter(isCrateRoot);
for (const crateRoot of crateRoots) {
    const lines = contents.get(crateRoot);
    if (
        lines === undefined ||
        !lines
            .slice(0, CRATE_ROOT_PREFIX_LINES)
            .some((line) => hasAttribute(line, COVERAGE_FEATURE))
    ) {
        errors.push(
            `${crateRoot}: ${COVERAGE_FEATURE} must appear within the first ` +
                `${CRATE_ROOT_PREFIX_LINES} lines`,
        );
    }
}

if (errors.length > 0) {
    console.error("Nightly coverage annotation checks failed:\n");
    for (const error of errors) console.error(`  ${error}`);
    process.exit(1);
}

console.log(
    `Nightly coverage annotations valid for ${crateRoots.length} crate root(s).`,
);
