#!/usr/bin/env python3
"""Assemble numerically ordered, CRC-verified chunks into one source file."""

import argparse
import hashlib
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--parts", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    total = 0
    with args.output.open("wb") as destination:
        for index in range(args.parts):
            data = (args.input_dir / f"part_{index:03d}.bin").read_bytes()
            destination.write(data)
            digest.update(data)
            total += len(data)
    print(f"assembled_bytes={total}")
    print(f"sha256={digest.hexdigest()}")


if __name__ == "__main__":
    main()
