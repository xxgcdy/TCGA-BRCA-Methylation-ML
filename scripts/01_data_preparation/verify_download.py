#!/usr/bin/env python3
import csv
import hashlib
import json
from datetime import datetime, timezone
import os
from pathlib import Path

PROJECT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
EXPECTED = PROJECT / "manifests" / "expected_md5_tcga_brca_hm450.tsv"
ROOT = PROJECT / "raw_beta"
REPORT = PROJECT / "results" / "md5_verification.tsv"
SUMMARY = PROJECT / "results" / "md5_verification_summary.json"


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


rows = []
with EXPECTED.open(newline="", encoding="utf-8-sig") as handle:
    for item in csv.DictReader(handle, delimiter="\t"):
        relative = item["relative_path"]
        path = ROOT / relative
        expected_md5 = item["expected_md5"].lower()
        expected_bytes = int(item["expected_bytes"])
        if not path.is_file():
            rows.append((relative, expected_md5, "", expected_bytes, "", "MISSING"))
            continue
        observed_bytes = path.stat().st_size
        observed_md5 = md5_file(path)
        status = "PASS" if observed_bytes == expected_bytes and observed_md5 == expected_md5 else "FAIL"
        rows.append((relative, expected_md5, observed_md5, expected_bytes, observed_bytes, status))

REPORT.parent.mkdir(parents=True, exist_ok=True)
with REPORT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["relative_path", "expected_md5", "observed_md5", "expected_bytes", "observed_bytes", "status"])
    writer.writerows(rows)

counts = {status: sum(row[-1] == status for row in rows) for status in ("PASS", "FAIL", "MISSING")}
summary = {
    "checked_at_utc": datetime.now(timezone.utc).isoformat(),
    "expected_files": len(rows),
    **{key.lower(): value for key, value in counts.items()},
}
SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
raise SystemExit(0 if counts["FAIL"] == 0 and counts["MISSING"] == 0 else 1)

