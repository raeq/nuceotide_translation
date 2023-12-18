from collections import UserDict
from typing import Dict
from src.exceptions import invalidcodonexception
from src.classes.amino_acids import *


class Codons(UserDict):
    """A codon is a sequence of three DNA or RNA nucleotides that corresponds with a specific amino acid or stop signal
    during protein synthesis. DNA and RNA molecules are written in a language of four nucleotides; meanwhile, the
    language of proteins includes 20 amino acids. The correspondence between a DNA or RNA sequence and an amino acid
    sequence is known as the genetic code. The genetic code consists of three-letter 'words' called codons formed from
    a sequence of three nucleotides (e.g. ACT, CAG, TTT).

    acids: a dictionary of known amino acids
    """

    _acids: dict = Dict[str, AbstractAminoAcid]
    data: dict = Dict[str, AbstractAminoAcid]

    def __init__(self):
        super().__init__()
        self._acids = {Alanine().symbol: Alanine(), Arginine().symbol: Arginine(), Asparagine().symbol: Asparagine(),
                       AsparticAcid().symbol: AsparticAcid(), Cysteine().symbol: Cysteine(),
                       GlutamicAcid().symbol: GlutamicAcid(), Glutamine().symbol: Glutamine(),
                       Glycine().symbol: Glycine(), Histidine().symbol: Histidine(), Isoleucine().symbol: Isoleucine(),
                       Leucine().symbol: Leucine(), Lysine().symbol: Lysine(), Methionine().symbol: Methionine(),
                       Phenylalanine().symbol: Phenylalanine(), Proline().symbol: Proline(), Serine().symbol: Serine(),
                       Threonine().symbol: Threonine(), Tryptophan().symbol: Tryptophan(),
                       Tyrosine().symbol: Tyrosine(), Valine().symbol: Valine(),
                       Selenocysteine().symbol: Selenocysteine(), Pyrrolysine().symbol: Pyrrolysine(),
                       UnspecifiedOrUnknown().symbol: UnspecifiedOrUnknown()}

        self.data = {}

        for k, v in self.acids.items():
            for codon in v.codons:
                self.data[codon.lower()] = v

    def __getitem__(self, key):
        """Returns the amino acid corresponding to the codon.
        :param key: the codon
        :return: the amino acid corresponding to the codon
        """

        if any(c not in {'a', 'c', 'g', 't'} for c in key.lower()):
            raise invalidcodonexception.InvalidCodonException(f"Key must be acgt {key=}")

        return self.data[key.lower()]


    @property
    def acids(self):
        """Returns the amino acids dictionary.
        :return: the amino acids dictionary
        """
        return self._acids
