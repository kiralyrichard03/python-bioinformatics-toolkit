# Python Bioinformatics Toolkit

A lightweight Python toolkit for parsing FASTA/FASTQ files, computing sequence statistics, finding open reading frames, and analyzing translated proteins.

## Biological question

Given nucleotide sequences from public databases, what are their compositional properties, coding potential, and derived protein characteristics?

## Features

- **FASTA parser** with NCBI `[key=value]` header attribute extraction
- **FASTQ parser** with per-read quality statistics (mean/median Phred, % Q30)
- **GC content** and nucleotide composition
- **ORF finder** with reading-frame filtering and longest-ORF translation
- **Protein analysis**: molecular weight, pI, aromaticity, hydrophobicity profile, cysteine environments

## Prerequisites

- Python 3.10+
- pip

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/python-bioinformatics-toolkit.git
cd python-bioinformatics-toolkit
pip install -r requirements.txt
```

## Quick start

```bash
# Nucleotide statistics
python cli.py fasta stats data/sample.fasta

# ORF detection and translation
python cli.py fasta orf data/sample.fasta

# Protein properties (MW, pI, hydrophobicity, cysteines)
python cli.py fasta protein data/sample.fasta

# FASTQ quality summary
python cli.py fastq stats data/sample.fastq

# Export FASTQ stats to TSV
python cli.py fastq stats data/sample.fastq -o results/fastq_summary.tsv
```

## Sample data

| File | Source |
|------|--------|
| `data/sample.fasta` | Human insulin transcript `NM_000207.3` from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets) |
| `data/sample.fastq` | Synthetic demo reads for quality-statistics testing |

## Example output

```
NM_000207.3
Length: 465
GC content: 58.28%
...
Protein length: 110 amino acids
Molecular weight: 11980.91 Da
Isoelectric point (pI): 5.42
```

## Input / output

| Command | Input | Output |
|---------|-------|--------|
| `fasta stats` | `.fasta` | Length, GC%, nucleotide counts |
| `fasta orf` | `.fasta` | Start codons, translated protein |
| `fasta protein` | `.fasta` | Protein physicochemical properties |
| `fastq stats` | `.fastq` | Per-read quality metrics (stdout or TSV) |

## Run tests

```bash
pip install pytest
pytest tests/ -v
```

## Project structure

```
src/biotoolkit/   # Core library modules
cli.py            # Command-line interface
data/             # Demo sequences
tests/            # pytest suite
```

## Limitations

This is a portfolio/educational toolkit, not a production clinical pipeline. ORF detection uses a longest-ORF heuristic; signal peptide prediction is a simple hydrophobicity window scan.

## License

MIT
