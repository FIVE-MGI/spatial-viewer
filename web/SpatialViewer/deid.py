r"""De-identify sample IDs before anything is published (HIPAA).

Sample IDs such as INT9901, INT9902b, T9903, I994, TFI95 or IO96 must never  (fake examples; deid:ok)
reach the viewer, the blob, config.js or GitHub. Every ID is replaced by a
plain number (1, 2, 3 ...). The same ID always gets the same number, across
every dataset and every run, because the mapping is kept in ONE private key
file that is never published:

    sample_id_key.csv     (set SPATIALVIZ_KEY to keep it outside the working tree)

Key columns: number, original_id, first_seen, context.

Forms produced
    value : "1"          obs/metadata cell values ("sample" column of an RNA-seq table)
    text  : "Sample 1"   free text such as dataset labels
    slug  : "s1"         dataset ids / folder names

Shorthand is expanded before numbering: "INT9911_9912" = INT9911 + INT9912,   (deid:ok)
"INT9921-2" = INT9921 + INT9922. A letter suffix (INT9931b) keeps the patient   (deid:ok)
number and the letter: "1b".

Command line
    python deid.py text  "T9941 / INT9942 / 9943 - Akoya"   -> Samples 1 / 2 / 3 - Akoya   (deid:ok)
    python deid.py check web\SpatialViewer D:\SpatialVizStatic\x\manifest.json
                        (exit 1 and list every leftover ID; scans names and text files)
    python deid.py key  (print the key - keep it off screens you share)
"""
from __future__ import annotations

import csv
import datetime as _dt
import json
import os
import re
import sys
import threading
from pathlib import Path

KEY_PATH = Path(os.environ.get("SPATIALVIZ_KEY", "sample_id_key.csv"))

# prefix + 2-4 digits, not glued to a preceding letter/digit; optional single-letter suffix
#
# UTD and NM are two further external id prefixes. NM is safe next to RefSeq
# transcript accessions only because those carry 6 or 9 digits (NM_000123) and this
# pattern caps at 4 with (?![0-9]) after it, so an accession never matches. Do not
# widen the digit count without checking that again - it would rewrite gene names.
_ID = re.compile(r"(?<![A-Za-z0-9])(INT|TFI|UTD|IO|NM|T|I)[-_ ]?(\d{2,4})([A-Za-z](?![A-Za-z]))?(?![0-9])", re.I)
# continuation right after an ID: "_9912", " / 9912", "-5"
_CONT = re.compile(r"^([\s/_,&+]*|-)(\d{1,4})(?![0-9A-Za-z])")

_lock = threading.Lock()
_key: dict[str, int] | None = None


def _load() -> dict[str, int]:
    global _key
    if _key is None:
        _key = {}
        if KEY_PATH.exists():
            with open(KEY_PATH, newline="", encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    _key[row["original_id"]] = int(row["number"])
    return _key


def _canon(prefix: str, digits: str) -> str:
    """INT95, int0095 and INT_0095 are the same sample (deid:ok)."""
    return f"{prefix.upper()}{int(digits):04d}"


def number_for(original: str, context: str = "") -> int:
    """Stable number for a canonical ID, appending to the private key if new."""
    with _lock:
        key = _load()
        if original in key:
            return key[original]
        n = max(key.values(), default=0) + 1
        key[original] = n
        KEY_PATH.parent.mkdir(parents=True, exist_ok=True)
        new = not KEY_PATH.exists()
        with open(KEY_PATH, "a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["number", "original_id", "first_seen", "context"])
            w.writerow([n, original, _dt.datetime.now().isoformat(timespec="seconds"), context])
        return n


def _render(n: int, suffix: str, form: str) -> str:
    tag = f"{n}{suffix.lower()}"
    return {"value": tag, "slug": f"s{tag}"}.get(form, f"Sample {tag}")


def anonymize(text: str, form: str = "text", context: str = "") -> str:
    """Replace every sample ID in `text`. `form` = text | value | slug."""
    if not isinstance(text, str) or not text:
        return text
    out, pos, found = [], 0, 0
    for m in _ID.finditer(text):
        if m.start() < pos or m.group(0) in _DTYPES:
            continue
        out.append(text[pos:m.start()])
        prefix, digits, suffix = m.group(1), m.group(2), m.group(3) or ""
        canon = _canon(prefix, digits)
        out.append(_render(number_for(canon, context), suffix, form))
        found += 1
        pos = m.end()
        while True:                                   # shorthand continuations
            c = _CONT.match(text[pos:])
            if not c:
                break
            sep, more = c.group(1), c.group(2)
            if len(more) == 4 and sep != "-":
                digits = more
            elif sep == "-" and len(more) <= 2:
                digits = digits[: len(digits) - len(more)] + more
            else:
                break
            canon = _canon(prefix, digits)
            out.append(sep + _render(number_for(canon, context), "", form))
            found += 1
            pos += c.end()
    out.append(text[pos:])
    s = "".join(out)
    if form == "text" and found > 1:
        s = re.sub(r"\bSample (\d+[a-z]?)((?:\s*[/_,&+-]\s*Sample \d+[a-z]?)+)",
                   lambda mm: "Samples " + re.sub(r"Sample ", "", mm.group(1) + mm.group(2)), s)
    return s


_DTYPES = {"int8", "int16", "int32", "int64"}      # numpy/C type names, never sample IDs

# Columns that ARE the sample/donor identity. Every value in them is replaced by a number,
# whatever its format - donor labels come in many shapes the ID pattern cannot know about.
SAMPLE_COLUMN = re.compile(r"^(sample|samples|sample_?id|sample_?name|donor|donor_?id|patient|patient_?id|subject|"
                           r"individual|orig\.?ident|library|library_?id|batch|specimen)$", re.I)


def is_sample_column(name: str, categories=()) -> bool:
    """A column is treated as sample identity if its name says so, or if most of its
    values contain a sample ID."""
    if SAMPLE_COLUMN.match(str(name).strip()) or "sample" in str(name).lower():
        return True
    cats = [str(c) for c in categories]
    return bool(cats) and sum(bool(find_ids(c)) for c in cats) >= max(1, len(cats) // 2)


def sample_value(value: str, context: str = "") -> str:
    """Whole-value replacement for a sample-identity column: "<anything>" -> "3".
    Keyed on the full original value, so the same sample gets the same number in every dataset."""
    v = str(value)
    if v in ("", "NA", "nan", "None"):
        return "NA"
    return str(number_for(f"SAMPLE:{v}", context))


def find_ids(text: str) -> list[str]:
    return [m.group(0) for m in _ID.finditer(text or "") if m.group(0) not in _DTYPES]


def find_in_obj(obj, path="") -> list[tuple[str, str]]:
    """Every (json-path, id) left inside a nested JSON-like object (keys and values)."""
    hits = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            hits += [(f"{path}/{k}", i) for i in find_ids(str(k))]
            hits += find_in_obj(v, f"{path}/{k}")
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            hits += find_in_obj(v, f"{path}[{i}]")
    elif isinstance(obj, str):
        hits += [(path, i) for i in find_ids(obj)]
    return hits


def assert_clean(obj, what: str):
    hits = find_in_obj(obj)
    if hits:
        raise SystemExit(f"[deid] REFUSING to publish {what}: sample IDs still present: {hits[:10]}")


_TEXT_EXT = {".js", ".json", ".html", ".txt", ".md", ".css", ".geojson", ".csv", ".py", ".ps1", ".yml", ".yaml"}
_ALLOW = "deid:ok"          # a line containing this marker is a documented FAKE example
_SKIP_DIRS = {".venv", "__pycache__", ".git", "tiles", "gene"}


def check_paths(paths) -> list[str]:
    """Scan file/folder NAMES and text-file CONTENTS for sample IDs."""
    problems = []
    for p in map(Path, paths):
        items = [p] if p.is_file() else [q for q in p.rglob("*") if not (set(q.parts) & _SKIP_DIRS)]
        for q in items:
            for i in find_ids(q.name):
                problems.append(f"name    {q}  ->  {i}")
            if q.is_file() and q.suffix.lower() in _TEXT_EXT and q.stat().st_size < 50_000_000:
                if q.name.endswith(".link.txt"):
                    continue                      # encrypted payloads are random base64
                try:
                    txt = q.read_text(encoding="utf-8", errors="ignore")
                except OSError:
                    continue
                for ln, line in enumerate(txt.splitlines(), 1):
                    if _ALLOW in line:
                        continue
                    for i in find_ids(line):
                        problems.append(f"content {q}:{ln}  ->  {i}")
    return problems


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd in ("text", "value", "slug"):
        print(anonymize(" ".join(rest), form=cmd, context="cli"))
    elif cmd == "check":
        probs = check_paths(rest or ["."])
        for p in probs:
            print(p)
        print(f"[deid] {len(probs)} sample ID(s) found" if probs else "[deid] clean")
        return 1 if probs else 0
    elif cmd == "key":
        print(KEY_PATH.read_text(encoding="utf-8") if KEY_PATH.exists() else "(no key yet)")
    else:
        print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
