from math import inf
from warnings import warn

from ..solvers import solver_instance
from ..solvers.solver import Solution, Solver, Status


def FBA(
    model,
    objective: str | dict | None = None,
    minimize: bool = False,
    parsimonious: bool = False,
    obj_frac: float = 1.0,
    constraints: dict | None = None,
    solver: Solver | None = None,
    get_values: bool | list = True,
    shadow_prices: bool = False,
) -> Solution:
    """Run a Flux Balance Analysis (FBA) simulation:

    Arguments:
        model: a constraint-based model
        objective: objective coefficients or single reaction to optimize
        minimize: whether to minimize the objective function
        parsimonious: run parsimonious FBA (pFBA)
        obj_frac: minimum fraction of FBA objective (for pFBA)
        constraints: environmental or additional constraints
        solver: re-use a solver pre-initialized with the model (for speed)
        get_values: set to false if you only care about the objective value (for speed)
        shadow_prices: calculate shadow prices

    Returns:
        Solution
    """

    if not objective:
        objective = model.objective

    if not solver:
        solver = solver_instance(model)

    if not parsimonious:
        solution = solver.solve(
            objective,
            minimize=minimize,
            constraints=constraints,
            get_values=get_values,
            shadow_prices=shadow_prices,
        )
        return solution
    else:
        pre_solution = solver.solve(
            objective,
            minimize=minimize,
            constraints=constraints,
            get_values=False,
        )

        if pre_solution.status != Status.OPTIMAL:
            warn("Failed to find an optimal solution for initial problem.")
            return pre_solution

        if isinstance(objective, str):
            objective = {objective: 1}

        solver.add_constraint("obj", objective, ">", obj_frac * pre_solution.fobj)  # pyright: ignore

        if solver.reuse_for is None:
            solver.reuse_for = "pFBA"
            for r_id, rxn in model.reactions.items():
                if rxn.lb < 0:
                    solver.add_variable(f"{r_id}_abs", 0, inf)
            solver.update()
            for r_id, rxn in model.reactions.items():
                if rxn.lb < 0:
                    solver.add_constraint(
                        f"c_{r_id}_lower", {r_id: -1, f"{r_id}_abs": 1}, ">", 0
                    )
                    solver.add_constraint(
                        f"c_{r_id}_upper", {r_id: 1, f"{r_id}_abs": 1}, ">", 0
                    )
            solver.update()
        elif solver.reuse_for != "pFBA":
            raise RuntimeError(
                "Attempting to reuse a solver for pFBA previously used for: "
                + solver.reuse_for
            )

        objective = dict()
        for r_id, rxn in model.reactions.items():
            if rxn.lb < 0:
                objective[f"{r_id}_abs"] = 1
            else:
                objective[r_id] = 1

        solution = solver.solve(
            objective,
            minimize=True,
            constraints=constraints,
            get_values=list(model.reactions),
            shadow_prices=shadow_prices,
        )
        solver.remove_constraint("obj")

        return solution
