import pytest

from src.classes.amino_acids import AbstractAminoAcid, Methionine
from src.classes.codons import Codons
from src.api.translate import translate_fasta
from src.exceptions.invalidcodonexception import InvalidCodonException


def test_fasta_01():
    c: Codons = Codons()
    assert c['atg'] == Methionine()


def test_api():
    a: AbstractAminoAcid
    input_fasta = 'ggtaagtcctctagtacaaacacccccaatattgtgatataattaaaattatattcatattctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc'

    assert ''.join([a.symbol for a in translate_fasta(input_fasta)]) == 'GKSSSTNTPNIVI*LKLYSYSVARKNTFRLYOSHLLUSVV'

def test_mixedcase():
    a: AbstractAminoAcid
    input_fasta = 'gXtaagtcctctagtacaaacacccccaatattgtgatataattaaaattatattcatattctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc'

    with pytest.raises(InvalidCodonException):
        assert ''.join([a.symbol for a in translate_fasta(input_fasta)]) == 'GKSSSTNTPNIVI*LKLYSYSVARKNTFRLYOSHLLUSVV'
