from ..solvers import solver_instance
from ..solvers.solver import Solution, Solver


def FBA(
    model,
    objective: str | dict | None = None,
    minimize: bool = False,
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

    solution = solver.solve(
        objective,
        minimize=minimize,
        constraints=constraints,
        get_values=get_values,
        shadow_prices=shadow_prices,
    )
    return solution
