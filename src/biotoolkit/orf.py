"""Open reading frame detection and translation."""

from Bio.Seq import Seq


def get_reading_frame(position: int) -> int:
    """Return reading frame (0, 1, or 2) for a 1-based nucleotide position."""
    return (position - 1) % 3


def find_main_start_codon(seq) -> int | None:
    """Return the 1-based ATG position that yields the longest protein."""
    sequence = str(seq).upper()
    best_position = None
    best_length = 0

    for i in range(len(sequence) - 2):
        if sequence[i : i + 3] != "ATG":
            continue
        protein = str(Seq(sequence[i:]).translate(to_stop=True, table=1))
        if len(protein) > best_length:
            best_length = len(protein)
            best_position = i + 1

    return best_position


def find_start_codon_positions(seq, reading_frame: int | None = None) -> list[int]:
    """Find ATG start codon positions in the correct reading frame.

    Returns 1-based positions. ATGs in other reading frames are omitted.
    If reading_frame is not given, it is inferred from the main ORF start.
    """
    sequence = str(seq).upper()
    all_positions = [
        i + 1
        for i in range(len(sequence) - 2)
        if sequence[i : i + 3] == "ATG"
    ]

    if not all_positions:
        return []

    if reading_frame is None:
        main_start = find_main_start_codon(seq)
        if main_start is None:
            return []
        reading_frame = get_reading_frame(main_start)

    return [
        position
        for position in all_positions
        if get_reading_frame(position) == reading_frame
    ]


def translate_to_protein(seq) -> str:
    """Translate the longest open reading frame starting from an ATG codon."""
    sequence = str(seq).upper()
    longest_protein = ""

    for i in range(len(sequence) - 2):
        if sequence[i : i + 3] != "ATG":
            continue
        protein = str(Seq(sequence[i:]).translate(to_stop=True, table=1))
        if len(protein) > len(longest_protein):
            longest_protein = protein

    return longest_protein
