#!/usr/bin/env python3
"""Command-line interface for the bioinformatics toolkit."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from biotoolkit.fasta import load_fasta_records
from biotoolkit.fastq import summarize_fastq, write_summary_tsv
from biotoolkit.orf import (
    find_main_start_codon,
    find_start_codon_positions,
    get_reading_frame,
    translate_to_protein,
)
from biotoolkit.protein import (
    amino_acid_composition,
    find_cysteines,
    find_signal_peptide_candidate,
    hydrophobicity_profile,
    protein_aromaticity,
    protein_isoelectric_point,
    protein_molecular_weight,
)
from biotoolkit.stats import gc_content, nucleotide_composition


def cmd_fasta_stats(args: argparse.Namespace) -> None:
    for record in load_fasta_records(args.fasta):
        composition = nucleotide_composition(record.seq)
        gc = gc_content(record.seq)

        print(record.id)
        print(record.description)
        if record.annotations:
            print(record.annotations)
        print(f"Length: {len(record.seq)}")
        print(f"GC content: {gc:.2f}%")
        print("Nucleotide counts:", composition["counts"])
        print("Nucleotide percentages:")
        for base, pct in composition["percentages"].items():
            print(f"  {base}: {pct:.2f}%")
        print()


def cmd_fasta_orf(args: argparse.Namespace) -> None:
    for record in load_fasta_records(args.fasta):
        main_start = find_main_start_codon(record.seq)
        reading_frame = get_reading_frame(main_start) if main_start else None
        start_codons = find_start_codon_positions(
            record.seq, reading_frame=reading_frame
        )

        print(record.id)
        if start_codons:
            print(f"Reading frame: {reading_frame}")
            print(f"Start codon (ATG) in reading frame at position(s): {start_codons}")
            print(f"Translation start (longest ORF) at position: {main_start}")
        else:
            print("No start codon (ATG) found")

        protein = translate_to_protein(record.seq)
        if protein:
            print(f"Protein length: {len(protein)} amino acids")
            print(f"Protein sequence: {protein}")
        else:
            print("No protein translated (no start codon found)")
        print()


def cmd_fasta_protein(args: argparse.Namespace) -> None:
    for record in load_fasta_records(args.fasta):
        protein = translate_to_protein(record.seq)
        if not protein:
            print(f"{record.id}: no protein translated")
            print()
            continue

        aa_composition = amino_acid_composition(protein)
        mw = protein_molecular_weight(protein)
        pi = protein_isoelectric_point(protein)
        aromaticity = protein_aromaticity(protein)

        print(record.id)
        print(f"Protein length: {len(protein)} amino acids")
        print(f"Molecular weight: {mw:.2f} Da")
        print(f"Isoelectric point (pI): {pi:.2f}")
        print(f"Aromaticity: {aromaticity:.3f} ({aromaticity * 100:.1f}%)")
        print("Amino acid counts:")
        for aa, count in sorted(aa_composition["counts"].items()):
            if count > 0:
                print(f"  {aa}: {count}")

        profile = hydrophobicity_profile(protein, window_size=10)
        signal_peptide = find_signal_peptide_candidate(profile)
        if signal_peptide:
            print("Potential signal peptide (highest hydrophobicity window):")
            print(
                f"  Location: positions {signal_peptide['window_start']}-"
                f"{signal_peptide['window_end']}"
            )
            print(f"  Sequence: {signal_peptide['window_sequence']}")
            print(f"  Average hydrophobicity: {signal_peptide['hydrophobicity']:.3f}")

        cysteines = find_cysteines(protein)
        print(f"Cysteine residues found: {len(cysteines)}")
        for cysteine in cysteines:
            print(
                f"  Position {cysteine['position']}: {cysteine['environment']} "
                f"(environment positions {cysteine['environment_start']}-"
                f"{cysteine['environment_end']})"
            )
        print()


def cmd_fastq_stats(args: argparse.Namespace) -> None:
    stats = summarize_fastq(args.fastq)
    if args.output:
        write_summary_tsv(args.fastq, args.output)
        print(f"Wrote summary to {args.output}")

    for row in stats:
        print(
            f"{row['id']}\tlength={row['length']}\t"
            f"mean_Q={row['mean_quality']:.2f}\t"
            f"median_Q={row['median_quality']:.2f}\t"
            f"pct_Q30={row['pct_q30']:.2f}%"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Python bioinformatics toolkit for FASTA/FASTQ analysis"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    fasta_parser = subparsers.add_parser("fasta", help="Analyze FASTA files")
    fasta_sub = fasta_parser.add_subparsers(dest="fasta_command", required=True)

    stats_parser = fasta_sub.add_parser("stats", help="Nucleotide statistics")
    stats_parser.add_argument("fasta", type=Path, help="Path to FASTA file")
    stats_parser.set_defaults(func=cmd_fasta_stats)

    orf_parser = fasta_sub.add_parser("orf", help="ORF detection and translation")
    orf_parser.add_argument("fasta", type=Path, help="Path to FASTA file")
    orf_parser.set_defaults(func=cmd_fasta_orf)

    protein_parser = fasta_sub.add_parser(
        "protein", help="Protein property analysis"
    )
    protein_parser.add_argument("fasta", type=Path, help="Path to FASTA file")
    protein_parser.set_defaults(func=cmd_fasta_protein)

    fastq_parser = subparsers.add_parser("fastq", help="Analyze FASTQ files")
    fastq_sub = fastq_parser.add_subparsers(dest="fastq_command", required=True)

    fastq_stats = fastq_sub.add_parser("stats", help="Per-read quality statistics")
    fastq_stats.add_argument("fastq", type=Path, help="Path to FASTQ file")
    fastq_stats.add_argument(
        "-o", "--output", type=Path, help="Optional TSV output path"
    )
    fastq_stats.set_defaults(func=cmd_fastq_stats)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
