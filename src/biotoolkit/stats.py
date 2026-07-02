"""Nucleotide sequence statistics."""


def gc_content(seq) -> float:
    """Return GC content as a percentage (0-100)."""
    sequence = str(seq).upper()
    gc_count = sequence.count("G") + sequence.count("C")
    total = len(sequence)
    return (gc_count / total * 100) if total else 0.0


def nucleotide_composition(seq) -> dict:
    """Return nucleotide counts and percentages for DNA sequences."""
    sequence = str(seq).upper()
    total = len(sequence)
    bases = ("A", "T", "G", "C")
    counts = {base: sequence.count(base) for base in bases}
    other = total - sum(counts.values())
    if other:
        counts["other"] = other
    percentages = {
        base: (count / total * 100) if total else 0.0
        for base, count in counts.items()
    }
    return {"counts": counts, "percentages": percentages}
