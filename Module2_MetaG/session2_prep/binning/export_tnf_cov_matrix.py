#!/usr/bin/env python3
"""Export a TNF + coverage + GC feature matrix for splits >= a length
threshold, replicating anvi'o's own default 'tnf-cov' merged-profile
clustering recipe (anvio/data/clusterconfigs/merged/tnf-cov), for a subset
small enough to cluster independently of any binning tool (see plan:
cluster the merged profile independently of the binning tools).

Columns: TNF (kmer_contigs, CONTIGS.db) + raw mean_coverage_contigs
(PROFILE.db, per sample) + log mean_coverage_Q2Q3_contigs (PROFILE.db, per
sample, log-transformed not normalized) + gc_content_parent (CONTIGS.db
splits_basic_info).
"""
import sys
import math
import anvio.tables as t
import anvio.dbops as dbops

contigs_db_path, profile_db_path, min_len, out_path = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]

contigs_db = dbops.ContigsDatabase(contigs_db_path)
profile_db = dbops.ProfileDatabase(profile_db_path)
sample_names = sorted(list(profile_db.samples))

splits_basic_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name)
keep = {s for s, info in splits_basic_info.items() if info['length'] >= min_len}
print(f"{len(keep)} of {len(splits_basic_info)} splits >= {min_len} bp", file=sys.stderr)

kmer = contigs_db.db.get_table_as_dict(t.kmer_contigs_table_name if hasattr(t, 'kmer_contigs_table_name') else 'kmer_contigs')
kmer_cols = None

cov, _ = profile_db.db.get_view_data('mean_coverage_contigs', splits_basic_info=splits_basic_info)
cov_q2q3, _ = profile_db.db.get_view_data('mean_coverage_Q2Q3_contigs', splits_basic_info=splits_basic_info, log_norm_numeric_values=True)

with open(out_path, "w") as f:
    header = None
    n_written = 0
    for s in keep:
        if s not in kmer or s not in cov or s not in cov_q2q3:
            continue
        kmer_entry = {k: v for k, v in kmer[s].items() if k != 'contig'}
        if kmer_cols is None:
            kmer_cols = sorted(kmer_entry.keys())
        row = [str(kmer_entry[k]) for k in kmer_cols]
        row += [str(cov[s][sm]) for sm in sample_names]
        row += [str(cov_q2q3[s][sm]) for sm in sample_names]
        gc = splits_basic_info[s]['gc_content']
        row.append(str(gc))
        if header is None:
            header = ["item"] + [f"tnf_{k}" for k in kmer_cols] + [f"cov_{sm}" for sm in sample_names] + [f"covq23_{sm}" for sm in sample_names] + ["gc"]
            f.write("\t".join(header) + "\n")
        f.write(s + "\t" + "\t".join(row) + "\n")
        n_written += 1

print(f"wrote {n_written} rows to {out_path}", file=sys.stderr)
