from ..core.model import Model
from .solver import Solver

available_solvers = dict()

try:
    from .scip_wrapper import SCIPSolver

    available_solvers["scip"] = SCIPSolver
except ImportError:
    pass

try:
    from .gurobi_wrapper import GurobiSolver

    available_solvers["gurobi"] = GurobiSolver
except ImportError:
    pass


try:
    from .cplex_wrapper import CplexSolver

    available_solvers["cplex"] = CplexSolver
except ImportError:
    pass


default_solver: str | None = None


def get_default_solver() -> str:

    global default_solver

    if default_solver:
        return default_solver

    solver_order = ["gurobi", "cplex", "scip"]

    for solver in solver_order:
        if solver in list(available_solvers.keys()):
            default_solver = solver
            break

    if not default_solver:
        raise RuntimeError("No solver available.")

    return default_solver


def set_default_solver(solvername: str) -> None:

    global default_solver

    if solvername.lower() in list(available_solvers.keys()):
        default_solver = solvername.lower()
    else:
        available = ", ".join(list(available_solvers.keys()))
        raise RuntimeError(f"{solvername} not in available solvers: {available}")


def solver_instance(model: Model | None = None) -> Solver:
    return available_solvers[get_default_solver()](model)
