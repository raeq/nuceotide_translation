from src.classes.amino_acids import Codons, AbstractAminoAcid, Methionine
from src.api.translate import translate_fasta

def test_fasta_01():
    c: Codons = Codons()
    assert c['atg'] == Methionine()


def test_api():
    a: AbstractAminoAcid
    input_fasta = 'ggtaagtcctctagtacaaacacccccaatattgtgatataattaaaattatattcatattctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc'

    assert ''.join([a.symbol for a in translate_fasta(input_fasta)]) == 'GKSSSTNTPNIVI*LKLYSYSVARKNTFRLYOSHLLUSVV'

