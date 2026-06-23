#!/usr/bin/env python3
"""
blast_to_count_table.py

Create an OTU count/taxonomy table from:
1. headerless BLAST tabular output
2. taxonomy TSV
3. OTU count file

Supports:
    -d unite
        taxonomy columns: id1 id2 kingdom phylum class order family genus species
        BLAST ref_id expected as: id1|id2

    -d silva
        taxonomy columns: sequence_id domain phylum class order family genus species
        BLAST ref_id expected as: sequence_id

    -d amf
        taxonomy columns: id1 id2 family genus description vtx
        BLAST ref_id expected as: id1|id2

Expected BLAST columns:
    query_id ref_id pident length mismatch gapopen qstart qend sstart send evalue bitscore

Expected OTU count file:
    #OTU ID    SAMPLE_NAME
    OTU_1      170545

Usage:
    python blast_to_count_table.py \
        -d unite \
        -b blast_results.tsv \
        -t unite_taxonomy.tsv \
        -c otu_counts.txt \
        -o otu_count_taxonomy.tsv

    python blast_to_count_table.py \
        -d silva \
        -b blast_results.tsv \
        -t silva_taxonomy.tsv \
        -c otu_counts.txt \
        -o otu_count_taxonomy.tsv
        
    python blast_to_count_table.py \
        -d amf \
        -b blast_results.tsv \
        -t amf_taxonomy.tsv \
        -c otu_counts.txt \
        -o otu_count_taxonomy.tsv
"""

import argparse
import csv
import sys
import re
from collections import defaultdict


BLAST_COLUMNS = [
    "query_id",
    "ref_id",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

UNITE_TAX_COLUMNS = [
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

SILVA_TAX_COLUMNS = [
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

AMF_TAX_COLUMNS = [
    "family",
    "genus",
    "description",
    "vtx",
]

NGSPECIESID_PATTERN = re.compile(
    r"consensus_cl_id_(\d+)_total_supporting_reads_(\d+)"
)

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create OTU count/taxonomy table from BLAST output, "
            "taxonomy TSV, and OTU count file."
        )
    )
    parser.add_argument(
        "-d", "--database",
        required=True,
        choices=["unite", "silva", "amf"],
        help="Database/taxonomy format: unite, silva or amf."
    )
    parser.add_argument(
        "-b", "--blast",
        required=True,
        help="Input BLAST tabular TSV without header."
    )
    parser.add_argument(
        "-t", "--taxonomy",
        required=True,
        help="Input taxonomy TSV."
    )
    parser.add_argument(
        "-c", "--counts",
        required=False,
        help="Input OTU count table. Required unless --ngspeciesid is used."
    )
    parser.add_argument(
        "--ngspeciesid",
        action="store_true",
        help=(
            "Extract cluster ID and read count from NGSpeciesID-style query names, "
            "e.g. consensus_cl_id_0_total_supporting_reads_138085. "
            "If used, -c/--counts is not required."
        )
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file."
    )
    return parser.parse_args()


def extract_ngspeciesid_info(query_id):
    match = NGSPECIESID_PATTERN.search(query_id)

    if not match:
        raise ValueError(
            f"Could not parse NGSpeciesID query ID: {query_id}"
        )

    cluster_id = match.group(1)
    read_count = int(match.group(2))

    return cluster_id, read_count

def read_counts(count_file):
    counts = {}

    with open(count_file, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")

        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue

            if row[0].startswith("#"):
                continue

            if len(row) < 2:
                raise ValueError(
                    f"Count table row {line_number} has fewer than 2 columns: {row}"
                )

            otu_id = row[0].strip()

            try:
                read_count = int(row[1].strip())
            except ValueError as e:
                raise ValueError(
                    f"Could not parse read count in count table row {line_number}: {row}"
                ) from e

            counts[otu_id] = read_count

    return counts


def read_taxonomy(taxonomy_file, database):
    taxonomy = {}

    if database == "unite":
        tax_columns = UNITE_TAX_COLUMNS
        required = {"id1", "id2"} | set(tax_columns)
    elif database == "silva":
        tax_columns = SILVA_TAX_COLUMNS
        required = {"sequence_id"} | set(tax_columns)
    elif database == "amf":
        tax_columns = AMF_TAX_COLUMNS
        required = {"id1", "id2"} | set(tax_columns)
    else:
        raise ValueError(f"Unsupported database: {database}")

    with open(taxonomy_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("Taxonomy file appears to be empty.")

        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"{database.upper()} taxonomy file is missing required columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            if database in ["unite", "amf"]:
                ref_id = f"{row['id1']}|{row['id2']}"
            elif database == "silva":
                ref_id = row["sequence_id"]
            else:
                raise ValueError(
                    f"Unsupported database type: {database}"
                )

            taxonomy[ref_id] = {
                col: row.get(col, "") for col in tax_columns
            }

    return taxonomy, tax_columns


def aggregate_blast(blast_file):
    aggregated = defaultdict(lambda: {
        "ref_id": None,
        "sum_bitscore": 0.0,
        "sum_length": 0,
        "min_evalue": float("inf"),
        "max_pident": 0.0,
    })

    with open(blast_file, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")

        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue

            if row[0].startswith("#"):
                continue

            if len(row) < len(BLAST_COLUMNS):
                raise ValueError(
                    f"BLAST row {line_number} has {len(row)} columns, "
                    f"expected at least {len(BLAST_COLUMNS)}:\n{row}"
                )

            rec = dict(zip(BLAST_COLUMNS, row[:len(BLAST_COLUMNS)]))

            query_id = rec["query_id"].strip()
            ref_id = rec["ref_id"].strip()

            try:
                bitscore = float(rec["bitscore"])
                align_len = int(float(rec["length"]))
                evalue = float(rec["evalue"])
                pident = float(rec["pident"])
            except ValueError as e:
                raise ValueError(
                    f"Could not parse numeric BLAST values in row {line_number}: {row}"
                ) from e

            key = (query_id, ref_id)

            aggregated[key]["ref_id"] = ref_id
            aggregated[key]["sum_bitscore"] += bitscore
            aggregated[key]["sum_length"] += align_len
            aggregated[key]["min_evalue"] = min(aggregated[key]["min_evalue"], evalue)
            aggregated[key]["max_pident"] = max(aggregated[key]["max_pident"], pident)

    by_query = defaultdict(list)

    for (query_id, ref_id), values in aggregated.items():
        by_query[query_id].append({
            "query_id": query_id,
            "ref_id": ref_id,
            "sum_bitscore": values["sum_bitscore"],
            "sum_length": values["sum_length"],
            "min_evalue": values["min_evalue"],
            "max_pident": values["max_pident"],
        })

    return by_query


def choose_best_hit(hit_list):
    return sorted(
        hit_list,
        key=lambda x: (
            -x["sum_bitscore"],
            -x["sum_length"],
            x["min_evalue"],
            -x["max_pident"],
            x["ref_id"],
        )
    )[0]


def otu_sort_key(otu_id):
    if otu_id.startswith("OTU_"):
        try:
            return int(otu_id.replace("OTU_", ""))
        except ValueError:
            pass
    return otu_id


def main():
    args = parse_args()

    if args.ngspeciesid:
        counts = None
    else:
        if args.counts is None:
            raise ValueError(
                "-c/--counts is required unless --ngspeciesid is used."
            )
        counts = read_counts(args.counts)
    taxonomy, tax_columns = read_taxonomy(args.taxonomy, args.database)
    by_query = aggregate_blast(args.blast)

    out_columns = [
        "otu_id",
        "read_count",
        "ref_id",
        "sum_bitscore",
        "sum_length",
        "min_evalue",
        "max_pident",
    ] + tax_columns

    n_written = 0
    n_missing_counts = 0
    n_missing_tax = 0

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=out_columns, delimiter="\t")
        writer.writeheader()

        for query_id in sorted(by_query.keys(), key=otu_sort_key):
            best_hit = choose_best_hit(by_query[query_id])
            ref_id = best_hit["ref_id"]

            if args.ngspeciesid:
                output_otu_id, read_count = extract_ngspeciesid_info(query_id)
            else:
                output_otu_id = query_id
                read_count = counts.get(query_id)
                if read_count is None:
                    n_missing_counts += 1
                    read_count = 0

            tax = taxonomy.get(ref_id)
            if tax is None:
                n_missing_tax += 1
                tax = {col: "" for col in tax_columns}

            writer.writerow({
                "otu_id": output_otu_id,
                "read_count": read_count,
                "ref_id": ref_id,
                "sum_bitscore": best_hit["sum_bitscore"],
                "sum_length": best_hit["sum_length"],
                "min_evalue": best_hit["min_evalue"],
                "max_pident": best_hit["max_pident"],
                **tax,
            })

            n_written += 1

    sys.stderr.write(f"Done. Wrote {n_written} OTUs to {args.output}\n")

    if n_missing_counts:
        sys.stderr.write(f"Warning: {n_missing_counts} OTUs had no count entry.\n")

    if n_missing_tax:
        sys.stderr.write(f"Warning: {n_missing_tax} best-hit Ref_ID values had no taxonomy match.\n")


if __name__ == "__main__":
    main()
