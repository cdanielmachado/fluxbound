from enum import Enum
from math import inf

from ..core.model import Model
from .solution import Solution, Status


class VarType(Enum):
    """Enumeration of possible variable types."""

    BINARY = "binary"
    INTEGER = "integer"
    CONTINUOUS = "continuous"


class Parameter(Enum):
    """Enumeration of parameters common to all solvers."""

    TIME_LIMIT = 0
    FEASIBILITY_TOL = 1
    INT_FEASIBILITY_TOL = 2
    OPTIMALITY_TOL = 3
    MIP_REL_GAP = 4
    MIP_ABS_GAP = 5
    POOL_SIZE = 6
    POOL_GAP = 7


class Solver:
    """Abstract class representing a generic solver.

    All solver interfaces should implement the methods defined in this class.
    """

    def __init__(self, model: Model | None = None):
        self.problem = None
        self.model: Model | None = model
        self.variables: list = []
        self.constraints: list = []
        self.objective: dict = {}
        self.minimize: bool = True
        self.reuse_for: str | None = None
        self._cached_vars: dict = {}
        self._cached_constrs: dict = {}

    def add_variable(
        self,
        var_id: str,
        lb: float = -inf,
        ub: float = inf,
        vartype: VarType = VarType.CONTINUOUS,
    ) -> None:
        self._cached_vars[var_id] = (lb, ub, vartype)

    def add_variables(self, var_dict: dict):
        pass

    def add_constraint(
        self, constr_id: str, lhs: dict, sense: str = "=", rhs: float = 0
    ) -> None:
        """Add a constraint to the current problem.

        Arguments:
            constr_id: constraint identifier
            lhs: variables and respective coefficients
            sense: constraint sense ('<', '=', '>')
            rhs: right-hand side of equation
        """

        self._cached_constrs[constr_id] = (lhs, sense, rhs)

    def add_constraints(self, constr_dict: dict) -> None:
        """Solver specific implementation"""
        pass

    def remove_constraint(self, constr_id: str) -> None:
        """Solver specific implementation"""
        pass

    def update(self) -> None:
        """Update internal structure. Used for efficient lazy updating."""

        if len(self._cached_vars) > 0:
            self.add_variables(self._cached_vars)
            self._cached_vars = {}

        if len(self._cached_constrs) > 0:
            self.add_constraints(self._cached_constrs)
            self._cached_constrs = {}

    def build_problem(self, model: Model) -> None:

        for r_id, reaction in model.reactions.items():
            self.add_variable(r_id, reaction.lb, reaction.ub)

        table = model.metabolite_reaction_lookup()

        for m_id in model.metabolites:
            self.add_constraint(m_id, table[m_id])

        self.update()

    def solve(
        self,
        objective: str | dict | None = None,
        minimize: bool = True,
        model: Model | None = None,
        constraints: dict | None = None,
        get_values: bool | list = True,
        shadow_prices: bool = False,
    ) -> Solution:
        """Solve the optimization problem.

        Arguments:
            objective: linear objective
            minimize: solve a minimization problem (default: True)
            model : model (optional, leave blank to reuse previous model structure)
            constraints : additional constraints (optional)
            get_values: yes/no or list of variables to return values for
            shadow_prices: return shadow prices if available (default: False)

        Returns:
            solution
        """

        if model:
            self.build_problem(model)
        else:
            self.update()

        if objective is not None:
            self.set_objective(objective, minimize)

        if constraints:
            old_bounds = self.set_temporary_bounds(constraints)

        status = self.internal_solve()

        if status == Status.OPTIMAL or status == Status.SUBOPTIMAL:
            solution = self.get_solution(status, get_values, shadow_prices)
        else:
            solution = Solution(status)

        if constraints:
            self.reset_bounds(old_bounds)  # pyright: ignore

        return solution

    def set_objective(self, objective: str | dict, minimize: bool = True):
        pass

    def set_temporary_bounds(self, bounds: dict) -> dict:
        return {}

    def set_bounds(self, bounds: dict) -> None:
        self.set_temporary_bounds(bounds)

    def reset_bounds(self, bounds: dict) -> None:
        pass

    def internal_solve(self) -> Status:
        return Status.UNKNOWN

    def get_solution(
        self,
        status: Status,
        get_values: bool | list = True,
        shadow_prices: bool = False,
    ) -> Solution:
        return Solution(status)

    def set_parameter(self, parameter: Parameter, value: float) -> None:
        raise Exception("Not implemented for this solver.")

    def set_logging(self, enabled: bool = False) -> None:
        raise Exception("Not implemented for this solver.")

    def write_to_file(self, filename: str) -> None:
        raise Exception("Not implemented for this solver.")
