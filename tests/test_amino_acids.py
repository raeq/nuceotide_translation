from src.classes import amino_acids


def test_arginine():
    assert amino_acids.Arginine().symbol == 'R'
    assert amino_acids.Arginine().symbol_3 == 'Arg'