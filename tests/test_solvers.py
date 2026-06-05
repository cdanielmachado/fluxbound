from pathlib import Path

from pytest import approx
from warnings import warn

from fluxbound import FBA, load_model, set_default_solver

TEST_DATA = str(Path(__file__).parent) + "/data/"


def test_scip():
    try:
        set_default_solver('scip')
    except Exception:
        warn('SCIP not found, skipping test.')
    else:
        model = load_model(TEST_DATA + "e_coli_core.xml")
        sol = FBA(model)
        assert sol.status.value == "Optimal"
        assert approx(sol.fobj, 0.001) == 0.874


def test_gurobi():
    try:
        set_default_solver('gurobi')
    except Exception:
        warn('Gurobi not found, skipping test.')
    else:
        model = load_model(TEST_DATA + "e_coli_core.xml")
        sol = FBA(model)
        assert sol.status.value == "Optimal"
        assert approx(sol.fobj, 0.001) == 0.874


def test_cplex():
    try:
        set_default_solver('cplex')
    except Exception:
        warn('CPLEX not found, skipping test.')
    else:
        model = load_model(TEST_DATA + "e_coli_core.xml")
        sol = FBA(model)
        assert sol.status.value == "Optimal"
        assert approx(sol.fobj, 0.001) == 0.874