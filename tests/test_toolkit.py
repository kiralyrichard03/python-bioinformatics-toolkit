import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biotoolkit.fasta import load_fasta_records, parse_header_attributes
from biotoolkit.fastq import read_quality_stats, summarize_fastq
from biotoolkit.orf import find_main_start_codon, translate_to_protein
from biotoolkit.stats import gc_content, nucleotide_composition

DATA = ROOT / "data"


def test_parse_header_attributes():
    description = "INS [organism=Homo sapiens] [GeneID=3630] [transcript=1]"
    attrs = parse_header_attributes(description)
    assert attrs["organism"] == "Homo sapiens"
    assert attrs["GeneID"] == "3630"


def test_gc_content():
    assert gc_content("ATGC") == 50.0
    assert gc_content("AAAA") == 0.0


def test_nucleotide_composition():
    result = nucleotide_composition("ATGC")
    assert result["counts"] == {"A": 1, "T": 1, "G": 1, "C": 1}
    assert result["percentages"]["A"] == 25.0


def test_load_fasta_records():
    records = load_fasta_records(DATA / "sample.fasta")
    assert len(records) == 1
    assert records[0].id == "NM_000207.3"
    assert records[0].annotations["GeneID"] == "3630"


def test_orf_translation_on_sample():
    records = load_fasta_records(DATA / "sample.fasta")
    protein = translate_to_protein(records[0].seq)
    assert len(protein) == 110
    assert protein.startswith("MALWMR")


def test_find_main_start_codon():
    records = load_fasta_records(DATA / "sample.fasta")
    start = find_main_start_codon(records[0].seq)
    assert start is not None
    assert start > 0


def test_fastq_summary():
    stats = summarize_fastq(DATA / "sample.fastq")
    assert len(stats) == 5
    assert all(row["length"] == 20 for row in stats)


def test_read_quality_stats_q30():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(
        Seq("ATGC"),
        id="test",
        letter_annotations={"phred_quality": [10, 20, 30, 40]},
    )
    stats = read_quality_stats(record)
    assert stats["mean_quality"] == 25.0
    assert stats["pct_q30"] == 50.0
