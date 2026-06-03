from pathlib import Path
from pytest import approx

from fluxbound import load_model, FBA

TEST_DATA = str(Path(__file__).parent) + "/data/"


def test_fba_core():

    model = load_model(TEST_DATA + "e_coli_core.xml")
    sol = FBA(model)
    assert sol.status.value == "Optimal"
    assert approx(sol.fobj, 0.001) == 0.874


def test_fba_gem():
    model = load_model(TEST_DATA + "iML1515.xml.gz")
    sol = FBA(model)
    assert sol.status.value == "Optimal"
    assert approx(sol.fobj, 0.001) == 0.877
