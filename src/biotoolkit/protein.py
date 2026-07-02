"""Protein sequence analysis."""

from Bio.SeqUtils.ProtParam import ProteinAnalysis

AMINO_ACID_WEIGHTS = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1311,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
}
WATER_MASS = 18.01528

HYDROPHOBICITY_SCALE = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}


def amino_acid_composition(protein: str) -> dict:
    """Return amino acid counts and percentages for a protein sequence."""
    sequence = protein.upper()
    total = len(sequence)
    counts = {aa: sequence.count(aa) for aa in AMINO_ACID_WEIGHTS}
    other = total - sum(counts.values())
    if other:
        counts["other"] = other
    percentages = {
        aa: (count / total * 100) if total else 0.0
        for aa, count in counts.items()
        if count > 0
    }
    return {"counts": counts, "percentages": percentages}


def protein_molecular_weight(protein: str) -> float:
    """Calculate protein molecular weight as the sum of amino acid residue masses."""
    sequence = protein.upper()
    if not sequence:
        return 0.0

    residue_mass = sum(AMINO_ACID_WEIGHTS.get(aa, 0.0) for aa in sequence)
    return residue_mass + WATER_MASS


def protein_isoelectric_point(protein: str) -> float:
    """Return the isoelectric point (pI) using Biopython ProtParam."""
    return ProteinAnalysis(protein.upper()).isoelectric_point()


def protein_aromaticity(protein: str) -> float:
    """Return aromaticity (fraction of Phe, Trp, Tyr) using Biopython ProtParam."""
    return ProteinAnalysis(protein.upper()).aromaticity()


def hydrophobicity_profile(protein: str, window_size: int = 10) -> list[dict]:
    """Calculate average hydrophobicity over a sliding window."""
    sequence = protein.upper()
    if len(sequence) < window_size:
        return []

    profile = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size]
        scores = [HYDROPHOBICITY_SCALE.get(aa, 0.0) for aa in window]
        profile.append(
            {
                "window_start": i + 1,
                "window_end": i + window_size,
                "window_sequence": window,
                "hydrophobicity": sum(scores) / window_size,
            }
        )
    return profile


def find_signal_peptide_candidate(profile: list[dict]) -> dict | None:
    """Return the window with the highest average hydrophobicity."""
    if not profile:
        return None
    return max(profile, key=lambda point: point["hydrophobicity"])


def find_cysteines(protein: str, context_size: int = 3) -> list[dict]:
    """Find cysteine residues and their local amino acid environment."""
    sequence = protein.upper()
    cysteines = []

    for i, amino_acid in enumerate(sequence):
        if amino_acid != "C":
            continue
        env_start = max(0, i - context_size)
        env_end = min(len(sequence), i + context_size + 1)
        cysteines.append(
            {
                "position": i + 1,
                "environment": sequence[env_start:env_end],
                "environment_start": env_start + 1,
                "environment_end": env_end,
            }
        )

    return cysteines
