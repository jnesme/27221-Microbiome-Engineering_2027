#!/usr/bin/env python3
"""Read a SAM stream on stdin. For every aligned read (mate-aware key
QNAME_1 / QNAME_2) write the key to <prefix>.<class>.keys, where class is
'plastid' if the reference name matches --plastid-regex, else 'other'.
Prints "<class>\t<count>" lines to stdout when done."""
import argparse
import re
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--prefix", required=True)
ap.add_argument("--plastid-regex", default=None)
args = ap.parse_args()

rx = re.compile(args.plastid_regex) if args.plastid_regex else None
files = {"plastid": open(f"{args.prefix}.plastid.keys", "w"),
         "other": open(f"{args.prefix}.other.keys", "w")}
counts = {"plastid": 0, "other": 0}

for line in sys.stdin:
    if line.startswith("@"):
        continue
    f = line.split("\t", 4)
    flag = int(f[1])
    if flag & 4 or flag & 256 or flag & 2048:
        continue
    mate = 1 if flag & 64 else 2
    cls = "plastid" if (rx and rx.search(f[2])) else "other"
    files[cls].write(f"{f[0]}_{mate}\n")
    counts[cls] += 1

for fh in files.values():
    fh.close()
for k, v in counts.items():
    print(f"{k}\t{v}")
