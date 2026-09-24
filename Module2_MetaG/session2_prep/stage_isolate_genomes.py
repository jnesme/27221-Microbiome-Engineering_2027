#!/usr/bin/env python3
"""Download the ~963 known SynCom isolate genome assemblies (AtSC/HvSC/LjSC,
from isolate_genomes_combined.tsv) and build an anvi'o contigs database for
each, for use as ANI ground-truth references in Session 2 (external-genomes.txt).

Resumable: skips any genome whose contigs-db already exists, so a killed/
requeued job picks up where it left off rather than starting over.
"""
import csv
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "isolate_genomes_combined.tsv"
FASTA_DIR = HERE / "isolate_genomes" / "fasta"
CONTIGS_DIR = HERE / "isolate_genomes" / "contigs_dbs"
EXTERNAL_GENOMES = HERE / "external-genomes.txt"
N_WORKERS = 8
ANVIO_BIN = "/work3/josne/miniconda3/envs/anvio-9/bin"


def download_fasta(accession, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return True
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}"
    rc = subprocess.run(["curl", "-sS", "-o", str(dest), url]).returncode
    return rc == 0 and dest.exists() and dest.stat().st_size > 0


def build_contigs_db(accession, fasta_path, db_path, project_name):
    if db_path.exists():
        return True
    reformatted = fasta_path.with_suffix(".reformatted.fa")
    r1 = subprocess.run(
        [f"{ANVIO_BIN}/anvi-script-reformat-fasta", str(fasta_path),
         "-o", str(reformatted), "--simplify-names", "--min-len", "0", "--seq-type", "NT"],
        capture_output=True, text=True,
    )
    if r1.returncode != 0:
        print(f"  reformat FAILED {accession}: {r1.stderr[-300:]}", flush=True)
        return False
    r2 = subprocess.run(
        [f"{ANVIO_BIN}/anvi-gen-contigs-database", "-f", str(reformatted),
         "-o", str(db_path), "-n", project_name, "-T", "1"],
        capture_output=True, text=True,
    )
    reformatted.unlink(missing_ok=True)
    if r2.returncode != 0:
        print(f"  gen-contigs-db FAILED {accession}: {r2.stderr[-300:]}", flush=True)
        db_path.unlink(missing_ok=True)
        return False
    return True


def process_one(row):
    acc = row["accession"]
    fasta_path = FASTA_DIR / f"{acc}.fa"
    db_path = CONTIGS_DIR / f"{acc}-contigs.db"
    if not download_fasta(acc, fasta_path):
        return acc, "download-failed"
    project_name = f"{row['collection']}_{acc}".replace(".", "_").replace("-", "_")
    if not build_contigs_db(acc, fasta_path, db_path, project_name):
        return acc, "contigs-db-failed"
    return acc, "ok"


def main():
    FASTA_DIR.mkdir(parents=True, exist_ok=True)
    CONTIGS_DIR.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(MANIFEST), delimiter="\t"))
    print(f"{len(rows)} genomes to stage ({N_WORKERS} concurrent)", flush=True)

    ok, failed = [], []
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(process_one, row): row["accession"] for row in rows}
        for i, fut in enumerate(as_completed(futures), 1):
            acc, status = fut.result()
            if i % 25 == 0 or status != "ok":
                print(f"  [{i}/{len(rows)}] {acc}: {status}", flush=True)
            if status == "ok":
                ok.append(acc)
            else:
                failed.append((acc, status))

    print(f"\nDone: {len(ok)}/{len(rows)} ok, {len(failed)} failed", flush=True)
    for acc, status in failed:
        print(f"  FAILED {acc}: {status}", flush=True)

    # external-genomes.txt: name<TAB>contigs_db_path (anvi'o expects this exact format)
    with open(EXTERNAL_GENOMES, "w") as f:
        f.write("name\tcontigs_db_path\n")
        for row in rows:
            acc = row["accession"]
            db_path = CONTIGS_DIR / f"{acc}-contigs.db"
            if db_path.exists():
                name = f"{row['collection']}_{acc}".replace(".", "_").replace("-", "_")
                f.write(f"{name}\t{db_path}\n")
    print(f"external-genomes.txt written with {sum(1 for _ in open(EXTERNAL_GENOMES)) - 1} entries", flush=True)

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
