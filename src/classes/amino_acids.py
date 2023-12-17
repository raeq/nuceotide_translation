from abc import ABC
from collections import UserDict
from dataclasses import dataclass
from enum import Enum
from typing import List, Dict


class Charge(Enum):
    POSITIVE = "+"
    NEGATIVE = "-"
    NEUTRAL = "0"
    UNKNOWN = "Unknown"


class Polarity(Enum):
    NON_POLAR = "Non-polar"
    POLAR = "Polar"
    BASIC_POLAR = "Basic Polar"
    BRONSTED_A = "Bronsted Acid"
    BRONSTED_B = "Bronsted Base"
    UNKNOWN = "Unknown"


@dataclass(init=False, repr=True, eq=True, order=True, unsafe_hash=False, frozen=False)
class AbstractAminoAcid(ABC):
    name: str
    symbol: str
    symbol_3: str
    codons: List[str]
    polarity: Polarity
    charge: Charge
    smiles: str
    inchi_key: str
    IUPAC_name: str
    Monoisotopic_mass: float


class Alanine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Alanine"
        self.symbol = "A"
        self.symbol_3 = "Ala"
        self.codons = ["GCT", "GCC", "GCA", "GCG"]
        self.polarity = Polarity.NON_POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "CC(C)(C(=O)O)N"
        self.inchi_key = "QNAYBMKLOCPYGJ-UHFFFAOYSA-N"
        self.IUPAC_name = "2-aminopropanoic acid"
        self.Monoisotopic_mass = 89.0477


class Arginine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Arginine"
        self.symbol = "R"
        self.symbol_3 = "Arg"
        self.codons = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
        self.polarity = Polarity.BASIC_POLAR
        self.charge = Charge.POSITIVE
        self.smiles = "NCCCCNC(=N)N"
        self.inchi_key = "NKRVEGMKGKLLBR-UHFFFAOYSA-N"
        self.IUPAC_name = "2-amino-5-(diaminomethylideneamino)pentanoic acid"
        self.Monoisotopic_mass = 174.1117


class Asparagine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Asparagine"
        self.symbol = "N"
        self.symbol_3 = "Asn"
        self.codons = ["AAT", "AAC"]
        self.polarity = Polarity.POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "NC(C(=O)O)C"
        self.inchi_key = "DFPAKSUCGFBDDF-UHFFFAOYSA-N"
        self.IUPAC_name = "2-aminobutanedioic acid"
        self.Monoisotopic_mass = 132.0535


class AsparticAcid(AbstractAminoAcid):

    def __init__(self):
        self.name = "Aspartic Acid"
        self.symbol = "D"
        self.symbol_3 = "Asp"
        self.codons = ["GAT", "GAC"]
        self.polarity = Polarity.BRONSTED_A
        self.charge = Charge.NEGATIVE
        self.smiles = "NC(C(=O)O)C(=O)O"
        self.inchi_key = "CSCPPACGZOOCGX-UHFFFAOYSA-N"
        self.IUPAC_name = "2-aminobutanedioic acid"
        self.Monoisotopic_mass = 133.0375


class Cysteine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Cysteine"
        self.symbol = "C"
        self.symbol_3 = "Cys"
        self.codons = ["TGT", "TGC"]
        self.polarity = Polarity.BRONSTED_B
        self.charge = Charge.NEUTRAL
        self.smiles = "SCC(C(=O)O)N"
        self.inchi_key = "RWSXRVCMGQZWBV-UHFFFAOYSA-N"
        self.IUPAC_name = "2-amino-3-sulfanylpropanoic acid"
        self.Monoisotopic_mass = 121.0197


class GlutamicAcid(AbstractAminoAcid):

    def __init__(self):
        self.name = "Glutamic Acid"
        self.symbol = "E"
        self.symbol_3 = "Glu"
        self.codons = ["GAA", "GAG"]
        self.polarity = Polarity.BRONSTED_A
        self.charge = Charge.NEGATIVE
        self.smiles = "N[C@@H](CCC(=O)O)C(=O)O"
        self.inchi_key = "FEWJPZIEWOKRBE-GASJEMHNSA-N"
        self.IUPAC_name = "2-aminopentanedioic acid"
        self.Monoisotopic_mass = 147.0532


class Glutamine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Glutamine"
        self.symbol = "Q"
        self.symbol_3 = "Gln"
        self.codons = ["CAA", "CAG"]
        self.polarity = Polarity.POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "N[C@@H](CCC(N)=O)C(=O)O"
        self.inchi_key = "ZDXPYRJPNDTMRX-VKHMYHEASA-N"
        self.IUPAC_name = "2-aminopentanedioic acid"
        self.Monoisotopic_mass = 146.0691


class Glycine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Glycine"
        self.symbol = "G"
        self.symbol_3 = "Gly"
        self.codons = ["GGT", "GGC", "GGA", "GGG"]
        self.polarity = Polarity.NON_POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "NCC(=O)O"
        self.inchi_key = "DHMQDGOQFOQNFH-UHFFFAOYSA-N"
        self.IUPAC_name = "2-aminoacetic acid"
        self.Monoisotopic_mass = 75.0320


class Histidine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Histidine"
        self.symbol = "H"
        self.symbol_3 = "His"
        self.codons = ["CAT", "CAC"]
        self.polarity = Polarity.BRONSTED_A
        self.charge = Charge.NEUTRAL
        self.smiles = "N[C@@H](Cc1cnc[nH]1)C(=O)O"
        self.inchi_key = "NMBLJWWWAFQDNQ-UHFFFAOYSA-N"
        self.IUPAC_name = "2-amino-3-(1H-imidazol-5-yl)propanoic acid"
        self.Monoisotopic_mass = 155.0694


class Isoleucine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Isoleucine"
        self.symbol = "I"
        self.symbol_3 = "Ile"
        self.codons = ["ATT", "ATC", "ATA"]
        self.polarity = Polarity.NON_POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "CC[C@H](C)[C@@H](C(=O)O)N"
        self.inchi_key = "QBLZDYSMRJIZSP-VKHMYHEASA-N"
        self.IUPAC_name = "2-amino-3-methylpentanoic acid"
        self.Monoisotopic_mass = 131.0946


class Leucine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Leucine"
        self.symbol = "L"
        self.symbol_3 = "Leu"
        self.codons = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
        self.polarity = Polarity.NON_POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "CC[C@@H](C)[C@@H](C(=O)O)N"
        self.inchi_key = "CZMRCDWAGMRECN-UGDNZRGBSA-N"
        self.IUPAC_name = "2-amino-4-methylpentanoic acid"
        self.Monoisotopic_mass = 131.0946


class Lysine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Lysine"
        self.symbol = "K"
        self.symbol_3 = "Lys"
        self.codons = ["AAA", "AAG"]
        self.polarity = Polarity.BRONSTED_A
        self.charge = Charge.POSITIVE
        self.smiles = "NCCCC[C@H](N)C(=O)O"
        self.inchi_key = "FKERNNLJNOPAEC-QMMMGPOBSA-N"
        self.IUPAC_name = "2,6-diaminohexanoic acid"
        self.Monoisotopic_mass = 146.1055


class Methionine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Methionine"
        self.symbol = "M"
        self.symbol_3 = "Met"
        self.codons = ["ATG"]
        self.polarity = Polarity.NON_POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "SC[C@H](C(=O)O)N"
        self.inchi_key = "VNWKTOKETHGBQD-UHFFFAOYSA-N"
        self.IUPAC_name = "2-amino-4-(methylsulfanyl)butanoic acid"
        self.Monoisotopic_mass = 149.0477


class Phenylalanine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Phenylalanine"
        self.symbol = "F"
        self.symbol_3 = "Phe"
        self.codons = ["TTT", "TTC"]
        self.polarity = Polarity.NON_POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "N[C@@H](Cc1ccccc1)C(=O)O"
        self.inchi_key = "ZDXPYRJPNDTMRX-VKHMYHEASA-N"
        self.IUPAC_name = "2-amino-3-phenylpropanoic acid"
        self.Monoisotopic_mass = 165.0789


class Proline(AbstractAminoAcid):

    def __init__(self):
        self.name = "Proline"
        self.symbol = "P"
        self.symbol_3 = "Pro"
        self.codons = ["CCT", "CCC", "CCA", "CCG"]
        self.polarity = Polarity.NON_POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "N1[C@@H](CCC1)C(=O)O"
        self.inchi_key = "YQWVGJAVOOVJII-UHFFFAOYSA-N"
        self.IUPAC_name = "2-aminopropanoic acid"
        self.Monoisotopic_mass = 115.0633


class Serine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Serine"
        self.symbol = "S"
        self.symbol_3 = "Ser"
        self.codons = ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]
        self.polarity = Polarity.POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "NC(CO)C(=O)O"
        self.inchi_key = "QNAYBMKLOCPYGJ-UHFFFAOYSA-N"
        self.IUPAC_name = "2-amino-3-hydroxypropanoic acid"
        self.Monoisotopic_mass = 105.0426


class Threonine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Threonine"
        self.symbol = "T"
        self.symbol_3 = "Thr"
        self.codons = ["ACT", "ACC", "ACA", "ACG"]
        self.polarity = Polarity.POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "NC([C@@H](C)O)C(=O)O"
        self.inchi_key = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        self.IUPAC_name = "2-amino-3-hydroxybutanoic acid"
        self.Monoisotopic_mass = 119.0582


class Tryptophan(AbstractAminoAcid):

    def __init__(self):
        self.name = "Tryptophan"
        self.symbol = "W"
        self.symbol_3 = "Trp"
        self.codons = ["TGG"]
        self.polarity = Polarity.NON_POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"
        self.inchi_key = "QIVBCDIJIAJPQS-VKHMYHEASA-N"
        self.IUPAC_name = "2-amino-3-(1H-indol-3-yl)propanoic acid"
        self.Monoisotopic_mass = 204.0899


class Tyrosine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Tyrosine"
        self.symbol = "Y"
        self.symbol_3 = "Tyr"
        self.codons = ["TAT", "TAC"]
        self.polarity = Polarity.BRONSTED_A
        self.charge = Charge.NEUTRAL
        self.smiles = "N[C@@H](Cc1ccc(O)cc1)C(=O)O"
        self.inchi_key = "QTBSBXVTEAMEQO-VKHMYHEASA-N"
        self.IUPAC_name = "2-amino-3-(4-hydroxyphenyl)propanoic acid"
        self.Monoisotopic_mass = 181.0739


class Valine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Valine"
        self.symbol = "V"
        self.symbol_3 = "Val"
        self.codons = ["GTT", "GTC", "GTA", "GTG"]
        self.polarity = Polarity.NON_POLAR
        self.charge = Charge.NEUTRAL
        self.smiles = "CC(C)[C@@H](C(=O)O)N"
        self.inchi_key = "CZMRCDWAGMRECN-UGDNZRGBSA-N"
        self.IUPAC_name = "2-amino-3-methylbutanoic acid"
        self.Monoisotopic_mass = 117.0789


class Selenocysteine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Selenocysteine"
        self.symbol = "U"
        self.symbol_3 = "Sec"
        self.codons = ["TGA"]
        self.polarity = Polarity.UNKNOWN
        self.charge = Charge.UNKNOWN
        self.smiles = "SC[C@H](C(=O)O)N"
        self.inchi_key = "VNWKTOKETHGBQD-UHFFFAOYSA-N"
        self.IUPAC_name = "2-amino-3-selanylpropanoic acid"
        self.Monoisotopic_mass = 168.0535


class Pyrrolysine(AbstractAminoAcid):

    def __init__(self):
        self.name = "Pyrrolysine"
        self.symbol = "O"
        self.symbol_3 = "Pyl"
        self.codons = ["TAG"]
        self.polarity = Polarity.UNKNOWN
        self.charge = Charge.UNKNOWN
        self.smiles = "NCCCC[C@H](N)C(=O)O"
        self.inchi_key = "FKERNNLJNOPAEC-QMMMGPOBSA-N"
        self.IUPAC_name = "2,6-diaminohexanoic acid"
        self.Monoisotopic_mass = 255.1583


class UnspecifiedOrUnknown(AbstractAminoAcid):

    def __init__(self):
        self.name = "Unspecified or Unknown"
        self.symbol = "*"
        self.symbol_3 = "Xaa"
        self.codons = ["NNN"]
        self.polarity = Polarity.UNKNOWN
        self.charge = Charge.UNKNOWN
        self.smiles = "Unknown"
        self.inchi_key = "Unknown"
        self.IUPAC_name = "Unknown"
        self.Monoisotopic_mass = 0.0000


class Codons(UserDict):
    """AminoAcids class contains all the amino acids as attributes.
    The attributes are instances of the AbstractAminoAcid class.
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

    @property
    def acids(self):
        """Returns the amino acids dictionary.
        :return: the amino acids dictionary
        """
        return self._acids
