from src.classes.amino_acids import AbstractAminoAcid, Methionine
from src.classes.codons import Codons
from src.api.translate import translate_fasta

def test_fasta_01():
    c: Codons = Codons()
    assert c['atg'] == Methionine()


def test_api():
    a: AbstractAminoAcid
    input_fasta = 'ggtaagtcctctagtacaaacacccccaatattgtgatataattaaaattatattcatattctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc'

    assert ''.join([a.symbol for a in translate_fasta(input_fasta)]) == 'GKSSSTNTPNIVI*LKLYSYSVARKNTFRLYOSHLLUSVV'

