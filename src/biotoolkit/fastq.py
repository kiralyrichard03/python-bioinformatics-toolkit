"""FASTQ parsing and quality statistics."""

from pathlib import Path

from Bio import SeqIO


def load_fastq_records(path: str | Path):
    """Load FASTQ records from a file."""
    return list(SeqIO.parse(path, "fastq"))


def read_quality_stats(record) -> dict:
    """Compute per-read quality statistics from a Biopython SeqRecord."""
    qualities = record.letter_annotations.get("phred_quality", [])
    length = len(record.seq)

    if not qualities:
        return {
            "id": record.id,
            "length": length,
            "mean_quality": 0.0,
            "median_quality": 0.0,
            "pct_q30": 0.0,
        }

    sorted_qualities = sorted(qualities)
    mid = len(sorted_qualities) // 2
    if len(sorted_qualities) % 2:
        median = float(sorted_qualities[mid])
    else:
        median = (sorted_qualities[mid - 1] + sorted_qualities[mid]) / 2

    return {
        "id": record.id,
        "length": length,
        "mean_quality": sum(qualities) / len(qualities),
        "median_quality": median,
        "pct_q30": sum(q >= 30 for q in qualities) / len(qualities) * 100,
    }


def summarize_fastq(path: str | Path) -> list[dict]:
    """Return quality statistics for every read in a FASTQ file."""
    return [read_quality_stats(record) for record in load_fastq_records(path)]


def write_summary_tsv(path: str | Path, output_path: str | Path) -> None:
    """Write per-read FASTQ statistics to a TSV file."""
    stats = summarize_fastq(path)
    lines = [
        "read_id\tlength\tmean_quality\tmedian_quality\tpct_q30",
    ]
    for row in stats:
        lines.append(
            f"{row['id']}\t{row['length']}\t{row['mean_quality']:.2f}"
            f"\t{row['median_quality']:.2f}\t{row['pct_q30']:.2f}"
        )
    Path(output_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
