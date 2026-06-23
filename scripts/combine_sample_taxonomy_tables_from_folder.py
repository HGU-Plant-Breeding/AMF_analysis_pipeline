#!/usr/bin/env python3
"""
combine_sample_taxonomy_tables_from_folder.py

Combine multiple per-sample OTU taxonomy tables from one folder into a
single taxonomy-by-sample count matrix.

Supports two compatible taxonomy formats:

SILVA:
    domain phylum class order family genus species

UNITE:
    kingdom phylum class order family genus species

The script checks all input files and stops if taxonomy headers are inconsistent.

Expected input files:
    <sample>.cluster_count_taxonomy.tsv

Usage:
    python combine_sample_taxonomy_tables_from_folder.py \
        -i per_sample_tables \
        -o combined_species_matrix.tsv
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


SILVA_TAXONOMY_COLUMNS = [
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

UNITE_TAXONOMY_COLUMNS = [
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

AMF_TAXONOMY_COLUMNS = [
    "family",
    "genus",
    "description",
    "vtx",
]

BASE_REQUIRED_COLUMNS = [
    "otu_id",
    "read_count",
    "ref_id",
]

FILENAME_PATTERN = re.compile(r"^(.+)\.cluster_count_taxonomy\.tsv$")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Combine per-sample taxonomy tables from one folder into "
            "a taxonomy-by-sample count matrix. Automatically detects "
            "SILVA-style domain or UNITE-style kingdom taxonomy headers."
        )
    )
    parser.add_argument(
        "-i", "--input-folder",
        required=True,
        help="Folder containing *.cluster_count_taxonomy.tsv files."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV matrix file."
    )
    return parser.parse_args()


def natural_sample_sort_key(sample_name):
    if sample_name.isdigit():
        return (0, int(sample_name))
    return (1, sample_name)


def find_input_files(input_folder):
    sample_files = []

    for entry in os.listdir(input_folder):
        match = FILENAME_PATTERN.match(entry)
        if match:
            sample_name = match.group(1)
            filepath = os.path.join(input_folder, entry)
            if os.path.isfile(filepath):
                sample_files.append((sample_name, filepath))

    if not sample_files:
        raise ValueError(
            f"No files matching '*.cluster_count_taxonomy.tsv' found in {input_folder}"
        )

    sample_files.sort(key=lambda x: natural_sample_sort_key(x[0]))
    return sample_files


def detect_taxonomy_columns(fieldnames, filepath):
    fieldset = set(fieldnames)

    has_silva = all(col in fieldset for col in SILVA_TAXONOMY_COLUMNS)
    has_unite = all(col in fieldset for col in UNITE_TAXONOMY_COLUMNS)
    has_amf = all(col in fieldset for col in AMF_TAXONOMY_COLUMNS) and "description" in fieldset and "vtx" in fieldset

    if has_silva and has_unite:
        raise ValueError(
            f"File {filepath} contains both 'domain' and 'kingdom'. "
            "Cannot decide whether this is SILVA or UNITE format."
        )

    if has_silva:
        return "silva", SILVA_TAXONOMY_COLUMNS

    if has_unite:
        return "unite", UNITE_TAXONOMY_COLUMNS

    if has_amf:
        return "amf", AMF_TAXONOMY_COLUMNS

    raise ValueError(
        f"File {filepath} does not match SILVA, UNITE, or AMF taxonomy format.\n"
        f"Expected SILVA columns: {', '.join(SILVA_TAXONOMY_COLUMNS)}\n"
        f"Expected UNITE columns: {', '.join(UNITE_TAXONOMY_COLUMNS)}\n"
        f"Expected AMF columns: {', '.join(AMF_TAXONOMY_COLUMNS)}"
    )

def read_sample_table(filepath, expected_database=None, expected_taxonomy_columns=None):
    counts = defaultdict(int)

    with open(filepath, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"No header found in file: {filepath}")

        database, taxonomy_columns = detect_taxonomy_columns(reader.fieldnames, filepath)

        if expected_database is not None and database != expected_database:
            raise ValueError(
                f"Inconsistent taxonomy format detected.\n"
                f"Expected {expected_database.upper()} based on previous files, "
                f"but file {filepath} looks like {database.upper()}."
            )

        required_columns = BASE_REQUIRED_COLUMNS + taxonomy_columns
        missing = set(required_columns) - set(reader.fieldnames)

        if missing:
            raise ValueError(
                f"File {filepath} is missing required columns: "
                f"{', '.join(sorted(missing))}"
            )

        for row_num, row in enumerate(reader, start=2):
            try:
                reads = int(row["read_count"])
            except ValueError as e:
                raise ValueError(
                    f"Invalid read_count value in {filepath} line {row_num}: "
                    f"{row['read_count']}"
                ) from e

            tax_key = tuple(row[col] for col in taxonomy_columns)
            counts[tax_key] += reads

    return database, taxonomy_columns, counts


def main():
    args = parse_args()

    sample_files = find_input_files(args.input_folder)

    expected_database = None
    expected_taxonomy_columns = None

    per_sample_counts = {}
    all_taxa = set()

    for sample_name, filepath in sample_files:
        database, taxonomy_columns, sample_counts = read_sample_table(
            filepath,
            expected_database=expected_database,
            expected_taxonomy_columns=expected_taxonomy_columns,
        )

        if expected_database is None:
            expected_database = database
            expected_taxonomy_columns = taxonomy_columns
            sys.stderr.write(
                f"Detected {expected_database.upper()} taxonomy format.\n"
            )

        per_sample_counts[sample_name] = sample_counts
        all_taxa.update(sample_counts.keys())

        sys.stderr.write(f"Read sample {sample_name}: {filepath}\n")

    sample_names = [sample_name for sample_name, _ in sample_files]
    sorted_taxa = sorted(all_taxa)

    output_columns = expected_taxonomy_columns + sample_names

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(output_columns)

        for tax_key in sorted_taxa:
            row = list(tax_key)

            for sample_name in sample_names:
                row.append(per_sample_counts[sample_name].get(tax_key, 0))

            writer.writerow(row)

    sys.stderr.write(
        f"Done. Wrote {len(sorted_taxa)} taxa across {len(sample_names)} samples "
        f"to {args.output}\n"
    )


if __name__ == "__main__":
    main()
