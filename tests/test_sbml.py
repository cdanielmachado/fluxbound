from pathlib import Path

from fluxbound import load_model, save_model

TEST_DATA = str(Path(__file__).parent) + "/data/"


def test_sbml_small():

    model = load_model(TEST_DATA + "e_coli_core.xml")
    assert model is not None
    assert len(model.compartments) == 2
    assert len(model.metabolites) == 72
    assert len(model.genes) == 137
    assert len(model.reactions) == 95
    assert len(model.objective) == 1

    save_model(model, TEST_DATA + "e_coli_core_copy.xml")


def test_sbml_large():
    model = load_model(TEST_DATA + "Recon3D.xml.gz")
    assert model is not None
    assert len(model.compartments) == 9
    assert len(model.metabolites) == 5835
    assert len(model.genes) == 2248
    assert len(model.reactions) == 10600
    assert len(model.objective) == 1

    save_model(model, TEST_DATA + "Recon3D_copy.xml.gz")
