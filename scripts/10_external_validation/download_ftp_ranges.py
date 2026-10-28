#!/usr/bin/env python3
"""Resumable parallel FTP downloader for a single public file."""

import argparse
import ftplib
import os
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


class RangeComplete(Exception):
    pass


def download_part(host, remote, output, start, end, retries=30):
    output.parent.mkdir(parents=True, exist_ok=True)
    target_size = end - start + 1
    for attempt in range(1, retries + 1):
        current = output.stat().st_size if output.exists() else 0
        if current == target_size:
            return str(output), current
        # NCBI's FTP REST behavior was not byte-stable after interrupted gzip
        # transfers in testing. A failed small part is therefore restarted.
        if output.exists(): output.unlink()
        ftp = None
        try:
            ftp = ftplib.FTP(host, timeout=60)
            ftp.login(); ftp.set_pasv(True)
            remaining = target_size
            with output.open("wb") as handle:
                def write(block):
                    nonlocal remaining
                    if len(block) >= remaining:
                        handle.write(block[:remaining]); remaining = 0
                        raise RangeComplete
                    handle.write(block); remaining -= len(block)
                try:
                    ftp.retrbinary(f"RETR {remote}", write, blocksize=1024 * 1024, rest=start)
                except RangeComplete:
                    pass
            if output.stat().st_size == target_size:
                return str(output), target_size
        except Exception as exc:
            print(f"{output.name}: attempt {attempt} interrupted at {output.stat().st_size if output.exists() else 0}: {exc}", flush=True)
            time.sleep(min(attempt, 10))
        finally:
            if ftp is not None:
                try: ftp.close()
                except Exception: pass
    raise RuntimeError(f"Failed after {retries} attempts: {output}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--host", default="ftp.ncbi.nlm.nih.gov")
    p.add_argument("--remote", required=True)
    p.add_argument("--size", type=int, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--parts", type=int, default=8)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--first-part", type=int, default=0)
    p.add_argument("--last-part", type=int)
    p.add_argument("--skip-validation", action="store_true")
    args = p.parse_args()
    part_size = (args.size + args.parts - 1) // args.parts
    jobs = []
    last_part = args.parts - 1 if args.last_part is None else args.last_part
    if not (0 <= args.first_part <= last_part < args.parts):
        p.error("part range must satisfy 0 <= first-part <= last-part < parts")
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i in range(args.first_part, last_part + 1):
            start = i * part_size; end = min(args.size - 1, (i + 1) * part_size - 1)
            jobs.append(pool.submit(download_part, args.host, args.remote,
                                    args.output_dir / f"part_{i:03d}.bin", start, end))
        for job in jobs:
            print("complete", job.result(), flush=True)

    if args.skip_validation:
        return

    # A full streaming gzip check detects rare same-size corrupt FTP replies.
    # Redownload the first failing small part and restart validation until the
    # complete gzip member, including its CRC trailer, is valid.
    for validation_round in range(1, 31):
        decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
        failed = None
        try:
            for i in range(args.parts):
                decoder.decompress((args.output_dir / f"part_{i:03d}.bin").read_bytes())
            decoder.flush()
            if not decoder.eof:
                failed = args.parts - 1
        except zlib.error:
            failed = i
        if failed is None:
            print(f"gzip stream validation PASS after {validation_round} round(s)", flush=True)
            return
        path = args.output_dir / f"part_{failed:03d}.bin"
        print(f"gzip stream validation failed at {path.name}; redownloading", flush=True)
        path.unlink(missing_ok=True)
        start = failed * part_size; end = min(args.size - 1, (failed + 1) * part_size - 1)
        download_part(args.host, args.remote, path, start, end)
    raise RuntimeError("gzip stream validation did not converge")


if __name__ == "__main__":
    main()
