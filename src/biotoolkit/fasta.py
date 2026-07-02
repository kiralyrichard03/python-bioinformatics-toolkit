"""FASTA parsing utilities."""

import re
from pathlib import Path

from Bio import SeqIO


def parse_header_attributes(description: str) -> dict:
    """Extract [key=value] pairs from NCBI FASTA headers into a dict."""
    return dict(re.findall(r"\[(\w+)=([^\]]+)\]", description))


def load_fasta_records(path: str | Path):
    """Load FASTA records and attach parsed NCBI header attributes."""
    records = []
    for record in SeqIO.parse(path, "fasta"):
        record.annotations.update(parse_header_attributes(record.description))
        records.append(record)
    return records
