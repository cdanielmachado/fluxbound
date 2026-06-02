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

    model2 = load_model(TEST_DATA + "e_coli_core_copy.xml")
    assert model2 is not None
    assert len(model2.compartments) == 2
    assert len(model2.metabolites) == 72
    assert len(model2.genes) == 137
    assert len(model2.reactions) == 95
    assert len(model2.objective) == 1

    save_model(model, TEST_DATA + "e_coli_core_copy.xml")


def test_sbml_large():
    model = load_model(TEST_DATA + "iML1515.xml.gz")
    assert model is not None
    assert len(model.compartments) == 3
    assert len(model.metabolites) == 1877
    assert len(model.genes) == 1516
    assert len(model.reactions) == 2712
    assert len(model.objective) == 1

    save_model(model, TEST_DATA + "iML1515_copy.xml.gz")

    model2 = load_model(TEST_DATA + "iML1515_copy.xml.gz")
    assert model2 is not None
    assert len(model2.compartments) == 3
    assert len(model2.metabolites) == 1877
    assert len(model2.genes) == 1516
    assert len(model2.reactions) == 2712
    assert len(model2.objective) == 1
