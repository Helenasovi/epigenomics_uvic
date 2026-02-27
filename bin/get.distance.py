#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Find closest gene to a coordinate")
parser.add_argument("--input", required=True, help="gene.starts.tsv")
parser.add_argument("--start", required=True, type=int, help="Regulatory element start")

args = parser.parse_args()

query_start = args.start

closest_gene = None
closest_start = None
min_distance = None

with open(args.input) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        gene, start = line.split("\t")
        start = int(start)

        distance = abs(start - query_start)

        if min_distance is None or distance < min_distance:
            min_distance = distance
            closest_gene = gene
            closest_start = start

print(f"{closest_gene}\t{closest_start} {min_distance}")
