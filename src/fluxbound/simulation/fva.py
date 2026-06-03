from math import inf
from warnings import warn

from ..core.model import Model
from ..solvers import solver_instance
from ..solvers.solution import Status
from .fba import FBA


def FVA(
    model: Model,
    obj_frac: float = 0,
    reactions: list | None = None,
    constraints: dict | None = None,
):
    """Run Flux Variability Analysis (FVA).

    Arguments:
        model a constraint-based model
        obj_frac: minimum fraction of the default objective
        reactions: list of reactions to analyze (default: all)
        constraints: additional constraints

    Returns:
        flux variation
    """

    solver = solver_instance(model)

    if obj_frac > 0:
        solution = FBA(model, constraints=constraints, solver=solver)
        if solution.status != Status.OPTIMAL:
            raise RuntimeError("Unable to calculate objective: " + str(solution.status))
        obj_min = obj_frac * solution.fobj  # pyright: ignore[reportOperatorIssue]
        solver.add_constraint("obj_min", model.objective, ">", obj_min)

    if not reactions:
        reactions = list(model.reactions.keys())

    lower = {}
    for r_id in reactions:
        solution = FBA(
            model, r_id, True, constraints=constraints, solver=solver, get_values=False
        )

        if solution.status == Status.OPTIMAL:
            lower[r_id] = solution.fobj
        elif (
            solution.status == Status.UNBOUNDED or solution.status == Status.INF_OR_UNB
        ):
            lower[r_id] = -inf
        elif solution.status == Status.INFEASIBLE:
            lower[r_id] = None
            warn("Infeasible solution")
        else:
            lower[r_id] = None
            warn("Unknown solution status")

    upper = {}
    for r_id in reactions:
        solution = FBA(
            model, r_id, False, constraints=constraints, solver=solver, get_values=False
        )

        if solution.status == Status.OPTIMAL:
            upper[r_id] = solution.fobj
        elif (
            solution.status == Status.UNBOUNDED or solution.status == Status.INF_OR_UNB
        ):
            upper[r_id] = inf
        elif solution.status == Status.INFEASIBLE:
            upper[r_id] = None
            warn("Infeasible solution")
        else:
            upper[r_id] = None
            warn("Unknown solution status")

    return {(lower[r_id], upper[r_id]) for r_id in reactions}
