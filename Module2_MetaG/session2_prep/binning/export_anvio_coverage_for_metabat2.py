#!/usr/bin/env python3
"""Export per-contig, per-sample mean coverage from the merged anvi'o profile
DB (the same computation CONCOCT's anvi-cluster-contigs driver already used
successfully: anvio/cli/cluster_contigs.py's prepare_input_files - take the
first split's 'mean_coverage_contigs' view entry per parent contig), and
write it as a MetaBAT2 --cvExt depth file (contigName, contigLen, <sample>...,
no totalAvgDepth, no variance columns - avoids jgi_summarize_bam_contig_depths
entirely, which was confirmed to produce wrong depth values on several
contigs for this dataset).
"""
import sys
import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.fastalib as fastalib

contigs_db_path = sys.argv[1]
profile_db_path = sys.argv[2]
fasta_path = sys.argv[3]
out_path = sys.argv[4]

# MetaBAT2 requires the abundance file's contig order to exactly match the
# FASTA's order ("[Error!] the order of contigs in abundance file is not the
# same as the assembly file") - get that order directly from the FASTA
# actually being fed to metabat2, not from whatever order the DB query returns.
fasta_order = []
fa = fastalib.SequenceSource(fasta_path)
while next(fa):
    fasta_order.append(fa.id)
fa.close()

contigs_db = dbops.ContigsDatabase(contigs_db_path)
profile_db = dbops.ProfileDatabase(profile_db_path)
sample_names = sorted(list(profile_db.samples))

splits_basic_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name)
split_coverages, _ = profile_db.db.get_view_data('mean_coverage_contigs', splits_basic_info=splits_basic_info)

contig_coverages = {}
for split_name, entry in split_coverages.items():
    c = entry['__parent__']
    if c not in contig_coverages:
        contig_coverages[c] = entry

contig_lengths = contigs_db.db.get_table_as_dict(t.contigs_info_table_name)

with open(out_path, "w") as f:
    f.write("contigName\tcontigLen\t" + "\t".join(f"{s}.bam" for s in sample_names) + "\n")
    missing_len = missing_cov = 0
    for c in fasta_order:
        entry = contig_coverages.get(c)
        length = contig_lengths.get(c, {}).get('length')
        if entry is None:
            missing_cov += 1
            continue
        if length is None:
            missing_len += 1
            continue
        vals = "\t".join(str(entry[s]) for s in sample_names)
        f.write(f"{c}\t{length}\t{vals}\n")

print(f"wrote {out_path}: {len(fasta_order)} contigs in FASTA order, "
      f"{len(sample_names)} samples, {missing_len} missing length, {missing_cov} missing coverage")
