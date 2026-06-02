from pathlib import Path

TEST_DATA = Path(__file__).parent / "data"

def test_sbml_import():
    from fluxbound import load_model

    model = load_model(TEST_DATA / "e_coli_core.xml")
    assert model is not None
    assert len(model.compartments) == 2
    assert len(model.metabolites) == 72
    assert len(model.genes) == 137
    assert len(model.reactions) == 95
    assert len(model.objective) == 1