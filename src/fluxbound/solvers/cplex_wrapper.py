import sys
from math import inf
from warnings import warn

from cplex import Cplex, SparsePair, infinity  # pyright: ignore

from ..core.model import Model
from .solution import Solution, Status
from .solver import Parameter, Solver, VarType


def infinity_fix(val: float) -> float:
    if val == inf:
        return infinity
    elif val == -inf:
        return -infinity
    else:
        return float(val)


class CplexSolver(Solver):
    def __init__(self, model: Model | None = None) -> None:
        Solver.__init__(self)
        self.problem: Cplex = Cplex()
        self._cached_lower_bounds: dict = {}
        self._cached_upper_bounds: dict = {}

        status = self.problem.solution.status
        self.status_mapping: dict = {
            status.optimal: Status.OPTIMAL,
            status.optimal_tolerance: Status.OPTIMAL,
            status.unbounded: Status.UNBOUNDED,
            status.infeasible: Status.INFEASIBLE,
            status.infeasible_or_unbounded: Status.INF_OR_UNB,
            status.MIP_optimal: Status.OPTIMAL,
            status.MIP_unbounded: Status.UNBOUNDED,
            status.MIP_infeasible: Status.INFEASIBLE,
            status.MIP_infeasible_or_unbounded: Status.INF_OR_UNB,
        }

        self.vartype_mapping: dict = {
            VarType.BINARY: self.problem.variables.type.binary,
            VarType.INTEGER: self.problem.variables.type.integer,
            VarType.CONTINUOUS: self.problem.variables.type.continuous,
        }

        parameters = self.problem.parameters
        self.parameter_mapping: dict = {
            Parameter.TIME_LIMIT: parameters.timelimit,  # pyright: ignore
            Parameter.FEASIBILITY_TOL: parameters.simplex.tolerances.feasibility,  # pyright: ignore
            Parameter.OPTIMALITY_TOL: parameters.simplex.tolerances.optimality,  # pyright: ignore
            Parameter.INT_FEASIBILITY_TOL: parameters.mip.tolerances.integrality,  # pyright: ignore
            Parameter.MIP_ABS_GAP: parameters.mip.tolerances.mipgap,  # pyright: ignore
            Parameter.MIP_REL_GAP: parameters.mip.tolerances.absmipgap,  # pyright: ignore
            Parameter.POOL_SIZE: parameters.mip.limits.populate,  # pyright: ignore
            Parameter.POOL_GAP: parameters.mip.pool.relgap,  # pyright: ignore
        }

        self.set_logging(False)

        if model:
            self.build_problem(model)

    def add_variables(self, var_dict: dict) -> None:

        var_ids = list(var_dict.keys())

        lbs = [infinity_fix(lb) for (lb, _, _) in var_dict.values()]
        ubs = [infinity_fix(ub) for (_, ub, _) in var_dict.values()]
        vartypes = [
            self.vartype_mapping[vartype] for (_, _, vartype) in var_dict.values()
        ]

        self.problem.variables.add(names=var_ids, lb=lbs, ub=ubs, types=vartypes)  # pyright: ignore

        self.variables.extend(var_ids)
        self._cached_lower_bounds.update(dict(zip(var_ids, lbs)))
        self._cached_upper_bounds.update(dict(zip(var_ids, ubs)))

    def add_constraints(self, constr_dict: dict) -> None:

        constr_ids = list(constr_dict.keys())
        lhs_all = [
            SparsePair(ind=list(lhs.keys()), val=list(lhs.values()))
            for (lhs, _, _) in constr_dict.values()
        ]
        map_sense = {"=": "E", "<": "L", ">": "G"}
        sense_all = [map_sense[sense] for (_, sense, _) in constr_dict.values()]
        rhs_all = [rhs for (_, _, rhs) in constr_dict.values()]

        self.problem.linear_constraints.add(
            lin_expr=lhs_all,
            senses=sense_all,  # pyright: ignore
            rhs=rhs_all,
            names=constr_ids,
        )
        self.constraints.extend(constr_ids)

    def remove_constraint(self, constr_id: str) -> None:
        if constr_id in self.constraints:
            self.problem.linear_constraints.delete(constr_id)
        else:
            warn("Constraint not in problem: " + constr_id)

    def set_objective(self, objective: str | dict, minimize: bool = True) -> None:

        if isinstance(objective, str):
            objective = {objective: 1.0}

        updated_coeffs = {}

        for var_id in self.variables:
            if objective.get(var_id, 0) != self.objective.get(var_id, 0):
                updated_coeffs[var_id] = objective.get(var_id, 0)

        if updated_coeffs:
            self.problem.objective.set_linear(list(updated_coeffs.items()))

        self.objective = objective

        if minimize != self.minimize:
            sense = (
                self.problem.objective.sense.minimize
                if minimize
                else self.problem.objective.sense.maximize
            )
            self.problem.objective.set_sense(sense)
            self.minimize = minimize

    def internal_solve(self) -> Status:
        self.problem.solve()
        status = self.status_mapping.get(
            self.problem.solution.get_status(), Status.UNKNOWN
        )
        return status

    def get_solution(
        self,
        status: Status,
        get_values: bool | list = True,
        shadow_prices: bool = False,
    ) -> Solution:

        fobj = self.problem.solution.get_objective_value()

        if get_values:
            if isinstance(get_values, list):
                values = dict(
                    zip(get_values, self.problem.solution.get_values(get_values))
                )
            else:
                values = dict(zip(self.variables, self.problem.solution.get_values()))
        else:
            values = None

        if shadow_prices:
            s_prices = dict(
                zip(
                    self.constraints,
                    self.problem.solution.get_dual_values(self.constraints),
                )
            )
        else:
            s_prices = None

        return Solution(
            status,
            fobj=fobj,
            values=values,
            shadow_prices=s_prices,
        )

    def set_temporary_bounds(self, bounds: dict) -> dict:

        lb_new, ub_new = {}, {}
        reset = {"lb": [], "ub": []}

        for r_id, x in bounds.items():
            if r_id in self.variables:
                if isinstance(x, tuple):
                    lb, ub = infinity_fix(x[0]), infinity_fix(x[1])
                else:
                    lb, ub = infinity_fix(x), infinity_fix(x)

                if lb != self._cached_lower_bounds[r_id]:
                    lb_new[r_id] = lb
                    reset["lb"].append(r_id)
                if ub != self._cached_upper_bounds[r_id]:
                    ub_new[r_id] = ub
                    reset["ub"].append(r_id)
            else:
                warn(f"Constrained variable not previously declared: {r_id}")

        if len(lb_new) > 0:
            self.problem.variables.set_lower_bounds(lb_new.items())

        if len(ub_new) > 0:
            self.problem.variables.set_upper_bounds(ub_new.items())

        return reset

    def reset_bounds(self, bounds):

        lb_old = [(r_id, self._cached_lower_bounds[r_id]) for r_id in bounds["lb"]]
        if len(lb_old) > 0:
            self.problem.variables.set_lower_bounds(lb_old)

        ub_old = [(r_id, self._cached_upper_bounds[r_id]) for r_id in bounds["ub"]]
        if len(ub_old) > 0:
            self.problem.variables.set_upper_bounds(ub_old)

    def set_parameter(self, parameter: Parameter, value: float) -> None:

        if parameter in self.parameter_mapping:
            self.parameter_mapping[parameter].set(value)
        else:
            raise Exception(f"Unknown parameter: {str(parameter)}")

    def set_logging(self, enabled: bool = False) -> None:

        if enabled:
            self.problem.set_log_stream(sys.stdout)
            self.problem.set_error_stream(sys.stderr)
            self.problem.set_warning_stream(sys.stderr)
            self.problem.set_results_stream(sys.stdout)
        else:
            self.problem.set_log_stream(None)
            self.problem.set_error_stream(None)
            self.problem.set_warning_stream(None)
            self.problem.set_results_stream(None)

    def write_to_file(self, filename: str) -> None:
        self.problem.write(filename)
