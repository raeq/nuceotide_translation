from classes.codons import Codons
from classes.amino_acids import UnspecifiedOrUnknown


def translate_fasta(fasta: str) -> str:
    """
    Translate a FASTA file into a protein sequence.
    :param fasta: a FASTA file
    :return: a protein sequence
    """

    c = Codons()
    for i in range(0, len(fasta), 3):
        codon = fasta[i:i + 3]
        yield c.get(codon, UnspecifiedOrUnknown())