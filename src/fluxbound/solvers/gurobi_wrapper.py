from math import inf
from warnings import warn

from gurobipy import GRB, quicksum
from gurobipy import Model as GRBModel

from ..core.model import Model
from .solution import Solution, Status
from .solver import Parameter, Solver, VarType


def infinity_fix(val):
    if val == inf:
        return GRB.INFINITY
    elif val == -inf:
        return -GRB.INFINITY
    else:
        return val


status_mapping = {
    GRB.OPTIMAL: Status.OPTIMAL,
    GRB.UNBOUNDED: Status.UNBOUNDED,
    GRB.INFEASIBLE: Status.INFEASIBLE,
    GRB.INF_OR_UNBD: Status.INF_OR_UNB,
}

vartype_mapping = {
    VarType.BINARY: GRB.BINARY,
    VarType.INTEGER: GRB.INTEGER,
    VarType.CONTINUOUS: GRB.CONTINUOUS,
}

parameter_mapping = {
    Parameter.TIME_LIMIT: GRB.Param.TimeLimit,
    Parameter.FEASIBILITY_TOL: GRB.Param.FeasibilityTol,
    Parameter.INT_FEASIBILITY_TOL: GRB.Param.IntFeasTol,
    Parameter.OPTIMALITY_TOL: GRB.Param.OptimalityTol,
    Parameter.MIP_ABS_GAP: GRB.Param.MIPGapAbs,
    Parameter.MIP_REL_GAP: GRB.Param.MIPGap,
    Parameter.POOL_SIZE: GRB.Param.PoolSolutions,
    Parameter.POOL_GAP: GRB.Param.PoolGap,
}

sense_mapping = {"=": GRB.EQUAL, "<": GRB.LESS_EQUAL, ">": GRB.GREATER_EQUAL}


class GurobiSolver(Solver):
    def __init__(self, model: Model | None = None) -> None:
        Solver.__init__(self)
        self.problem: GRBModel = GRBModel()  # pyright: ignore[reportGeneralTypeIssues]
        self.set_logging(False)
        if model:
            self.build_problem(model)

    def add_variables(self, var_dict: dict) -> None:

        for var_id, (lb, ub, vartype) in var_dict.items():
            self.problem.addVar(
                name=var_id,
                lb=infinity_fix(lb),
                ub=infinity_fix(ub),
                vtype=vartype_mapping[vartype],
            )  # pyright: ignore[reportAttributeAccessIssue]

        self.variables.extend(var_dict.keys())
        self.problem.update()  # pyright: ignore[reportAttributeAccessIssue]

    def add_constraints(self, constr_dict: dict) -> None:

        for constr_id, (lhs, sense, rhs) in constr_dict.items():
            expr = quicksum(
                coeff * self.problem.getVarByName(r_id)
                for r_id, coeff in lhs.items()
                if coeff
            )  # pyright: ignore[reportAttributeAccessIssue]
            self.problem.addLConstr(expr, sense_mapping[sense], rhs, constr_id)  # pyright: ignore[reportAttributeAccessIssue]

        self.constraints.extend(constr_dict.keys())
        self.problem.update()  # pyright: ignore[reportAttributeAccessIssue]

    def remove_constraint(self, constr_id: str) -> None:
        if constr_id in self.constraints:
            self.problem.remove(self.problem.getConstrByName(constr_id))  # pyright: ignore[reportArgumentType, reportCallIssue, reportAttributeAccessIssue]
        else:
            warn("Constraint not in problem: " + constr_id)

    def set_objective(self, objective: str | dict, minimize: bool = True) -> None:

        if isinstance(objective, str):
            objective = {objective: 1.0}

        obj_expr = quicksum(
            [
                coeff * self.problem.getVarByName(r_id)
                for r_id, coeff in objective.items()
                if coeff != 0
            ]
        )  # pyright: ignore[reportAttributeAccessIssue]
        sense = GRB.MINIMIZE if minimize else GRB.MAXIMIZE

        self.problem.setObjective(obj_expr, sense)  # pyright: ignore[reportAttributeAccessIssue]

        self.objective = objective
        self.minimize = minimize

    def internal_solve(self) -> Status:
        self.problem.optimize()  # pyright: ignore[reportAttributeAccessIssue]
        status = status_mapping.get(self.problem.status, Status.UNKNOWN)  # pyright: ignore[reportAttributeAccessIssue]
        return status

    def get_solution(
        self,
        status: Status,
        get_values: bool | list = True,
        shadow_prices: bool = False,
    ) -> Solution:

        fobj = self.problem.ObjVal  # pyright: ignore[reportAttributeAccessIssue]

        if get_values:
            if isinstance(get_values, list):
                var_list = get_values
            else:
                var_list = self.variables

            values = {}
            for r_id in var_list:
                var = self.problem.getVarByName(r_id)
                value = var.X  # pyright: ignore[reportOptionalMemberAccess]
                values[r_id] = value
        else:
            values = None

        if shadow_prices:
            s_prices = {}
            for m_id in self.constraints:
                constr = self.problem.getConstrByName(m_id)
                value = constr.Pi  # pyright: ignore[reportOptionalMemberAccess]
                s_prices[m_id] = value
        else:
            s_prices = None

        return Solution(status, fobj=fobj, values=values, shadow_prices=s_prices)

    def set_temporary_bounds(self, bounds: dict) -> dict:

        old_constraints = {}
        for r_id, x in bounds.items():
            lb, ub = x if isinstance(x, tuple) else (x, x)
            if r_id in self.variables:
                lpvar = self.problem.getVarByName(r_id)  # pyright: ignore[reportAttributeAccessIssue, reportOptionalMemberAccess]
                old_constraints[r_id] = (lpvar.lb, lpvar.ub)  # pyright: ignore[reportAttributeAccessIssue, reportOptionalMemberAccess]
                lpvar.lb = infinity_fix(lb)  # pyright: ignore[reportAttributeAccessIssue, reportOptionalMemberAccess]
                lpvar.ub = infinity_fix(ub)  # pyright: ignore[reportAttributeAccessIssue, reportOptionalMemberAccess]
            else:
                warn(f"Constrained variable '{r_id}' not previously declared")

        self.problem.update()  # pyright: ignore[reportAttributeAccessIssue]

        return old_constraints

    def reset_bounds(self, bounds: dict) -> None:

        for r_id, (lb, ub) in bounds.items():
            lpvar = self.problem.getVarByName(r_id)  # pyright: ignore[reportAttributeAccessIssue]
            lpvar.lb, lpvar.ub = lb, ub  # pyright: ignore[reportAttributeAccessIssue, reportOptionalMemberAccess]
        self.problem.update()  # pyright: ignore[reportAttributeAccessIssue]

    def set_parameter(self, parameter: Parameter, value: float) -> None:

        if parameter in parameter_mapping:
            grb_param = parameter_mapping[parameter]
            self.problem.setParam(grb_param, value)  # pyright: ignore[reportAttributeAccessIssue]
        else:
            raise Exception(f"Unknown parameter: {str(parameter)}")

    def set_logging(self, enabled: bool = False) -> None:
        self.problem.setParam("OutputFlag", 1 if enabled else 0)  # pyright: ignore[reportAttributeAccessIssue]

    def write_to_file(self, filename: str) -> None:
        self.problem.write(filename)  # pyright: ignore[reportAttributeAccessIssue]
