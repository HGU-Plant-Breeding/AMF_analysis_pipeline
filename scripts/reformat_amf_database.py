#!/usr/bin/env python3
"""
reformat_amf_fasta.py

Reformat AMF FASTA database headers.

Input header:
>gi|ESA00914|gb|AJ854100| Paraglomeraceae Paraglomus sp. 18S rDNA TYPE VTX00001

Output FASTA header:
>ESA00914|AJ854100

Output taxonomy TSV columns:
ref_id    family    genus    description    vtx

Usage:
    python reformat_amf_fasta.py \
        -i amf_database.fasta \
        -o amf_database_reformatted.fasta \
        -t amf_taxonomy.tsv
"""

import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Reformat AMF FASTA headers and write taxonomy TSV."
    )
    parser.add_argument("-i", "--input", required=True, help="Input AMF FASTA")
    parser.add_argument("-o", "--output-fasta", required=True, help="Output FASTA with clean ref_id headers")
    parser.add_argument("-t", "--taxonomy-tsv", required=True, help="Output taxonomy TSV")
    return parser.parse_args()


def fasta_records(handle):
    header = None
    seq_chunks = []

    for line in handle:
        line = line.rstrip("\n")
        if not line:
            continue

        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_chunks)
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())

    if header is not None:
        yield header, "".join(seq_chunks)


def wrap_sequence(seq, width=80):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def parse_header(header):
    """
    Expected:
    gi|ESA00914|gb|AJ854100| Paraglomeraceae Paraglomus sp. 18S rDNA TYPE VTX00001
    """
    match = re.match(r"^gi\|([^|]+)\|gb\|([^|]+)\|\s*(.+)$", header)

    if not match:
        raise ValueError(f"Header does not match expected AMF format: {header}")

    id1 = match.group(1).strip()
    id2 = match.group(2).strip()
    rest = match.group(3).strip()

    ref_id = f"{id1}|{id2}"

    fields = rest.split()
    family = fields[0] if len(fields) >= 1 else ""
    genus = fields[1] if len(fields) >= 2 else ""

    vtx_match = re.search(r"(VTX\d+)", rest)
    vtx = vtx_match.group(1) if vtx_match else ""

    return id1, id2, ref_id, family, genus, rest, vtx


def main():
    args = parse_args()

    with open(args.input) as infile, \
         open(args.output_fasta, "w") as fasta_out, \
         open(args.taxonomy_tsv, "w") as tax_out:

        tax_out.write("id1\tid2\tfamily\tgenus\tdescription\tvtx\n")

        for header, sequence in fasta_records(infile):
            id1, id2, ref_id, family, genus, description, vtx = parse_header(header)

            fasta_out.write(f">{ref_id}\n")
            for chunk in wrap_sequence(sequence):
                fasta_out.write(chunk + "\n")

            tax_out.write(
                f"{id1}\t{id2}\t{family}\t{genus}\t{description}\t{vtx}\n"
            )


if __name__ == "__main__":
    main()
