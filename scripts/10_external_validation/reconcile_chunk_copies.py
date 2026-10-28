#!/usr/bin/env python3
"""Reconcile two FTP chunk copies by a third independent transfer on disagreement."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from download_ftp_ranges import download_part


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="ftp.ncbi.nlm.nih.gov")
    parser.add_argument("--remote", required=True)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--parts", type=int, required=True)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--copy", type=Path, required=True)
    parser.add_argument("--tie-break", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    args.tie_break.mkdir(parents=True, exist_ok=True)
    part_size = (args.size + args.parts - 1) // args.parts
    disagreements = []
    for index in range(args.parts):
        name = f"part_{index:03d}.bin"
        primary = (args.primary / name).read_bytes()
        copy = (args.copy / name).read_bytes()
        if primary == copy:
            continue
        disagreements.append(index)

    def reconcile(index):
        name = f"part_{index:03d}.bin"
        primary_path = args.primary / name
        copy_path = args.copy / name
        primary = primary_path.read_bytes()
        copy = copy_path.read_bytes()
        start = index * part_size
        end = min(args.size - 1, (index + 1) * part_size - 1)
        target = args.tie_break / name
        for attempt in range(1, 11):
            target.unlink(missing_ok=True)
            download_part(args.host, args.remote, target, start, end)
            third = target.read_bytes()
            if third == primary:
                break
            if third == copy:
                primary_path.write_bytes(copy)
                return index
            print(f"part {index}: tie-break attempt {attempt} matched neither copy", flush=True)
        else:
            raise RuntimeError(f"No majority for part {index}")
        return index

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for index in pool.map(reconcile, disagreements):
            print(f"part {index}: majority resolved", flush=True)
    print(f"reconciled {len(disagreements)} disagreements: {disagreements}", flush=True)


if __name__ == "__main__":
    main()
