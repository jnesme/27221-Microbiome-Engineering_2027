#!/usr/bin/env python3
"""Download the Module 2 sample subset (Module2_MetaG/M2_sample_selection.tsv)
from ENA, verifying each file's MD5 against the manifest. Run via
submit_download.sh (bsub) rather than directly.
"""
import csv
import hashlib
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "M2_sample_selection.tsv"
OUT_DIR = HERE / "raw_reads"
N_WORKERS = 4


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, expect_md5, dest):
    if dest.exists() and dest.exists() and md5(dest) == expect_md5:
        return dest.name, "already-ok"
    rc = subprocess.run(["curl", "-sS", "-o", str(dest), url]).returncode
    if rc != 0:
        return dest.name, f"curl-error({rc})"
    got = md5(dest)
    if got != expect_md5:
        return dest.name, f"md5-mismatch({got}!={expect_md5})"
    return dest.name, "ok"


def main():
    OUT_DIR.mkdir(exist_ok=True)
    jobs = []
    with open(MANIFEST) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for mate in ("1", "2"):
                url = row[f"fastq_r{mate}_ftp"]
                exp_md5 = row[f"fastq_r{mate}_md5"]
                fname = f"{row['run_accession']}_{mate}.fastq.gz"
                jobs.append((url, exp_md5, OUT_DIR / fname))

    print(f"{len(jobs)} files to fetch ({N_WORKERS} concurrent)", flush=True)
    ok, bad = 0, []
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = [ex.submit(fetch, u, m, d) for u, m, d in jobs]
        for fut in as_completed(futures):
            name, status = fut.result()
            print(f"  {status:16s} {name}", flush=True)
            if status in ("ok", "already-ok"):
                ok += 1
            else:
                bad.append((name, status))

    print(f"\nDone: {ok}/{len(jobs)} ok", flush=True)
    if bad:
        print(f"{len(bad)} FAILED:", flush=True)
        for name, status in bad:
            print(f"  {name}: {status}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
