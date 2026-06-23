#!/usr/bin/env python3

import argparse
import gzip
import bz2
import sys
from collections import defaultdict

try:
    import edlib
except ImportError:
    sys.stderr.write(
        "ERROR: Python module 'edlib' is required.\n"
        "Install with: mamba install -c conda-forge edlib\n"
        "or: pip install edlib\n"
    )
    sys.exit(1)


def rc(seq):
    """Return reverse complement of DNA sequence."""
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def open_in_any(path):
    """Open plain, gzipped, bzipped FASTQ, or stdin."""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".bz2"):
        return bz2.open(path, "rt")
    return open(path, "rt")


def open_out_any(path):
    """Open plain or gzipped output by filename suffix."""
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def iter_fastq(path):
    """
    Stream FASTQ records.

    Yields:
        read_id, sequence, quality
    """
    fh = open_in_any(path)
    close = path != "-"

    try:
        while True:
            header = fh.readline().rstrip()
            if not header:
                break

            seq = fh.readline().rstrip()
            plus = fh.readline().rstrip()
            qual = fh.readline().rstrip()

            if not seq or not plus or not qual:
                raise ValueError(f"Incomplete FASTQ record after header: {header}")

            if not header.startswith("@"):
                raise ValueError(f"FASTQ header does not start with '@': {header}")

            read_id = header[1:].split()[0]
            yield read_id, seq, qual

    finally:
        if close:
            fh.close()


def load_primers(tsv):
    """
    Load primer panel TSV.

    Required header:
        marker_id    fwd_primer    rev_primer    min_bp    max_bp
    """
    panel = []

    with open(tsv, "r") as f:
        header = f.readline().strip().split("\t")
        if len(header) == 1:
            # fallback for whitespace-separated header
            header = header[0].split()

        required = ["marker_id", "fwd_primer", "rev_primer", "min_bp", "max_bp"]
        missing = [x for x in required if x not in header]

        if missing:
            raise ValueError(
                "Primer panel is missing required header columns: "
                + ", ".join(missing)
                + "\nRequired header:\n"
                + "\t".join(required)
            )

        cols = {name: i for i, name in enumerate(header)}

        for line_no, line in enumerate(f, start=2):
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) == 1:
                # fallback for whitespace-separated rows
                parts = line.strip().split()

            try:
                marker_id = parts[cols["marker_id"]]
                fwd = parts[cols["fwd_primer"]].upper()
                rev = parts[cols["rev_primer"]].upper()
                min_bp = int(parts[cols["min_bp"]])
                max_bp = int(parts[cols["max_bp"]])
            except Exception as e:
                raise ValueError(f"Could not parse primer panel line {line_no}: {line.rstrip()}") from e

            panel.append((marker_id, fwd, rev, min_bp, max_bp))

    if not panel:
        raise ValueError("Primer panel contains zero primer pairs.")

    return panel


def find_primer_hit(seq, primer, max_errors, clamp3):
    """
    Find approximate primer match in read sequence.

    Uses edlib semi-global mode:
      query = primer
      target = read sequence

    Returns best hit as (start, end), 0-based inclusive coordinates,
    or None if no acceptable hit is found.

    The last clamp3 bases of the primer must match exactly at the
    3' end of the primer.
    """
    res = edlib.align(primer, seq, mode="HW", task="locations", k=max_errors)

    if res["editDistance"] < 0:
        return None

    best = None

    for start, end in res["locations"]:
        read_seg = seq[start:end + 1].upper()

        if clamp3 > 0:
            if len(primer) < clamp3 or len(read_seg) < clamp3:
                continue
            if read_seg[-clamp3:] != primer[-clamp3:].upper():
                continue

        # Prefer the shortest/tightest alignment among accepted hits
        if best is None or (end - start) < (best[1] - best[0]):
            best = (start, end)

    return best


def write_fastq_record(out, read_id, seq, qual):
    out.write(f"@{read_id}\n")
    out.write(f"{seq}\n")
    out.write("+\n")
    out.write(f"{qual}\n")


def products_string(sizes):
    """Return NA or semicolon-separated sorted unique product sizes."""
    if not sizes:
        return "NA"
    return ";".join(str(x) for x in sorted(set(sizes)))


def add_candidate_if_valid(
    candidates,
    orientation,
    left_hit,
    right_hit,
    min_bp,
    max_bp,
    remove_primers,
):
    """
    Add a candidate amplicon if two primer hits are in correct order.

    Coordinates are 0-based, end-inclusive for primer hits.
    Slices use Python half-open coordinates.
    """
    left_start, left_end = left_hit
    right_start, right_end = right_hit

    if left_end >= right_start:
        return

    full_size = right_end - left_start + 1

    if not (min_bp <= full_size <= max_bp):
        return

    if remove_primers == "Y":
        amp_start = left_end + 1
        amp_end = right_start
    else:
        amp_start = left_start
        amp_end = right_end + 1

    candidates.append({
        "orientation": orientation,
        "full_size": full_size,
        "amp_start": amp_start,
        "amp_end": amp_end,
    })


def main():
    parser = argparse.ArgumentParser(
        description="""
Detect primer-defined amplicons in Nanopore FASTQ reads and write matching
amplicon sequences to FASTQ/FASTQ.GZ.

The output size used for filtering and summaries is the full amplicon size
including both primer sequences.

Both read orientations are detected:
  forward: FWD ... REVRC
  reverse: REV ... FWDRC

Amplicon FASTQ output is kept in the original read orientation.
Reverse hits are NOT reverse-complemented.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i", "--reads",
        required=True,
        help="Input FASTQ/FASTQ.GZ/FASTQ.BZ2, or '-' for stdin. "
             "If using '-' with gzipped input, decompress before piping, e.g. gzip -dc *.fastq.gz | script -i - ..."
    )

    parser.add_argument(
        "-p", "--primers",
        required=True,
        help="Primer panel TSV with header: marker_id fwd_primer rev_primer min_bp max_bp"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output amplicon FASTQ. Use .gz suffix for FASTQ.GZ output."
    )

    parser.add_argument(
        "--remove-primers",
        choices=["Y", "N"],
        default="N",
        help="Remove primer sequences from exported amplicon FASTQ records: Y or N. Default: N"
    )

    parser.add_argument(
        "--max-errors",
        type=int,
        default=2,
        help="Maximum edit distance allowed per primer. Default: 2"
    )

    parser.add_argument(
        "--clamp3",
        type=int,
        default=8,
        help="Require the last N bases of each primer to match exactly. Default: 8"
    )

    parser.add_argument(
        "--calls",
        default=None,
        help="Optional read-level calls TSV output: read_id, marker_id, size_bp, orientation."
    )

    parser.add_argument(
        "--summary",
        default=None,
        help="Optional per-marker summary TSV output: marker_id/products."
    )

    parser.add_argument(
        "--combined",
        default=None,
        help="Optional combined matrix TSV. Appends one headerless line per run, one column per marker."
    )

    args = parser.parse_args()

    panel = load_primers(args.primers)

    # Precompute sequences for both possible read orientations.
    panel2 = []
    for marker_id, fwd, rev, min_bp, max_bp in panel:
        panel2.append({
            "marker_id": marker_id,
            "fwd": fwd,
            "rev": rev,
            "fwdrc": rc(fwd),
            "revrc": rc(rev),
            "min_bp": min_bp,
            "max_bp": max_bp,
        })

    calls = []
    per_marker_sizes = defaultdict(list)

    n_reads = 0
    n_amplicons = 0
    n_forward = 0
    n_reverse = 0

    with open_out_any(args.output) as amp_out:
        for read_id, seq, qual in iter_fastq(args.reads):
            n_reads += 1

            for marker in panel2:
                marker_id = marker["marker_id"]
                min_bp = marker["min_bp"]
                max_bp = marker["max_bp"]

                candidates = []

                # Forward-orientation amplicon:
                # FWD ... REVRC
                fwd_hit = find_primer_hit(seq, marker["fwd"], args.max_errors, args.clamp3)
                revrc_hit = find_primer_hit(seq, marker["revrc"], args.max_errors, args.clamp3)

                if fwd_hit and revrc_hit:
                    add_candidate_if_valid(
                        candidates=candidates,
                        orientation="forward",
                        left_hit=fwd_hit,
                        right_hit=revrc_hit,
                        min_bp=min_bp,
                        max_bp=max_bp,
                        remove_primers=args.remove_primers,
                    )

                # Reverse-orientation amplicon:
                # REV ... FWDRC
                # Output is kept exactly as it appears in the read.
                rev_hit = find_primer_hit(seq, marker["rev"], args.max_errors, args.clamp3)
                fwdrc_hit = find_primer_hit(seq, marker["fwdrc"], args.max_errors, args.clamp3)

                if rev_hit and fwdrc_hit:
                    add_candidate_if_valid(
                        candidates=candidates,
                        orientation="reverse",
                        left_hit=rev_hit,
                        right_hit=fwdrc_hit,
                        min_bp=min_bp,
                        max_bp=max_bp,
                        remove_primers=args.remove_primers,
                    )

                for hit in candidates:
                    amp_seq = seq[hit["amp_start"]:hit["amp_end"]]
                    amp_qual = qual[hit["amp_start"]:hit["amp_end"]]

                    if not amp_seq:
                        continue

                    if hit["orientation"] == "reverse":
                        n_reverse += 1
                    else:
                        n_forward += 1

                    amp_id = (
                        f"{read_id}|marker={marker_id}|size={hit['full_size']}"
                        f"|orientation={hit['orientation']}"
                        f"|amplicon_start={hit['amp_start'] + 1}"
                        f"|amplicon_end={hit['amp_end']}"
                        f"|remove_primers={args.remove_primers}"
                    )

                    write_fastq_record(amp_out, amp_id, amp_seq, amp_qual)

                    calls.append((
                        read_id,
                        marker_id,
                        hit["full_size"],
                        hit["orientation"],
                    ))

                    per_marker_sizes[marker_id].append(hit["full_size"])
                    n_amplicons += 1

    if args.calls is not None:
        with open(args.calls, "w") as out:
            out.write("read_id\tmarker_id\tsize_bp\torientation\n")
            for row in calls:
                out.write("\t".join(map(str, row)) + "\n")

    if args.summary is not None:
        with open(args.summary, "w") as out:
            out.write("marker_id\tproducts\n")
            for marker_id, *_ in panel:
                out.write(f"{marker_id}\t{products_string(per_marker_sizes.get(marker_id, []))}\n")

    if args.combined is not None:
        row_values = [
            products_string(per_marker_sizes.get(marker_id, []))
            for marker_id, *_ in panel
        ]
        with open(args.combined, "a") as out:
            out.write("\t".join(row_values) + "\n")

    sys.stderr.write(f"Processed reads: {n_reads}\n")
    sys.stderr.write(f"Amplicons written: {n_amplicons}\n")
    sys.stderr.write(f"Forward-orientation amplicons: {n_forward}\n")
    sys.stderr.write(f"Reverse-orientation amplicons: {n_reverse}\n")


if __name__ == "__main__":
    main()
