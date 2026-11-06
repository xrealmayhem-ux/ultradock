#!/usr/bin/env python3
"""Verify files in the transferred 5-HT2A recovery snapshot."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results" / "docking_recovery" / "transfer_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    failures = []
    for item in data["files"]:
        path = ROOT / item["path"]
        if not path.is_file():
            failures.append(f"missing: {item['path']}")
            continue
        if path.stat().st_size != item["size_bytes"]:
            failures.append(f"size mismatch: {item['path']}")
        elif sha256(path) != item["sha256"]:
            failures.append(f"sha256 mismatch: {item['path']}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"Verified {len(data['files'])} files in the recovery snapshot.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
