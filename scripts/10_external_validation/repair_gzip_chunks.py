#!/usr/bin/env python3
"""Repair rare corrupt FTP range chunks until the complete gzip CRC passes."""

import argparse
import zlib
from pathlib import Path

from download_ftp_ranges import download_part


def first_failure(output_dir: Path, parts: int):
    decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
    for index in range(parts):
        try:
            decoder.decompress((output_dir / f"part_{index:03d}.bin").read_bytes())
        except zlib.error as exc:
            return index, str(exc)
    decoder.flush()
    return (None, None) if decoder.eof else (parts - 1, "gzip stream did not reach EOF")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="ftp.ncbi.nlm.nih.gov")
    parser.add_argument("--remote", required=True)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--parts", type=int, required=True)
    parser.add_argument("--repair-radius", type=int, default=4)
    parser.add_argument("--max-rounds", type=int, default=100)
    args = parser.parse_args()
    part_size = (args.size + args.parts - 1) // args.parts

    for repair_round in range(1, args.max_rounds + 1):
        failed, error = first_failure(args.output_dir, args.parts)
        if failed is None:
            print(f"gzip stream and CRC PASS after {repair_round - 1} repair rounds", flush=True)
            return
        first = max(0, failed - args.repair_radius)
        print(f"round {repair_round}: failure at part {failed}: {error}; redownloading {first}-{failed}", flush=True)
        for index in range(first, failed + 1):
            path = args.output_dir / f"part_{index:03d}.bin"
            path.unlink(missing_ok=True)
            start = index * part_size
            end = min(args.size - 1, (index + 1) * part_size - 1)
            download_part(args.host, args.remote, path, start, end)
    raise RuntimeError("gzip repair did not converge")


if __name__ == "__main__":
    main()
