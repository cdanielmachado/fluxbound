from math import inf
from types import FunctionType
from warnings import warn

from .model import Model
from .utils import AttrDict


class Environment(AttrDict):
    def __init__(self) -> None:
        AttrDict.__init__(self)

    def __str__(self) -> str:
        entries = (f"{r_id}\t{lb}\t{ub}" for r_id, (lb, ub) in self.items())
        return "\n".join(entries)

    def __repr__(self) -> str:
        return str(self)

    @staticmethod
    def from_model(model: Model) -> "Environment":

        env = Environment()

        for r_id in model.get_exchange_reactions():
            rxn = model.reactions[r_id]
            env[r_id] = rxn.lb, rxn.ub

        return env

    @staticmethod
    def complete(
        model: Model, max_uptake: float = 10.0, inplace: bool = False
    ) -> "Environment | None":
        """
        Generate a complete growth medium for a given model

        Arguments:
            model (to extract a list of all exchange reactions)
            max_uptake: maximum uptake rate
            inplace: apply to model or return environment object

        Returns:
            Environment (if not inplace)

        """

        env = Environment()

        for r_id in model.get_exchange_reactions():
            env[r_id] = (-max_uptake, inf)

        if inplace:
            env.apply(model, exclusive=False, inplace=True)
        else:
            return env

    @staticmethod
    def empty(model: Model, inplace: bool = False) -> "Environment | None":
        """
        Generate an empty growth medium for a given model

        Arguments:
            model (to extract a list of all exchange reactions)
            inplace: apply to model or return environment object

        Returns:
            Environment (if not inplace)

        """

        env = Environment()

        for r_id in model.get_exchange_reactions():
            env[r_id] = (0, inf)

        if inplace:
            env.apply(model, exclusive=False, inplace=True)
        else:
            return env

    def apply(
        self,
        model: Model,
        exclusive: bool = True,
        inplace: bool = True,
        warning: bool = True,
    ) -> dict | None:
        """
        Apply environmental conditions to a given model

        Args:
            model (to extract a list of all exchange reactions)
            exclusive: block uptake compounds not specified in this environment
            warning: print warning for exchange reactions not found in the model
            inplace: apply to model, otherwise return a constraints dict
        """

        if exclusive:
            env = Environment.empty(model, inplace=False)
            env.update(self)  # pyright: ignore[reportOptionalMemberAccess]
        else:
            env = self

        constraints = {}

        for r_id, (lb, ub) in env.items():  # pyright: ignore[reportOptionalMemberAccess]
            if r_id in model.reactions:
                if inplace:
                    model.set_flux_bounds(r_id, lb, ub)
                else:
                    constraints[r_id] = (lb, ub)
            elif warning:
                warn(f"Exchange reaction not in model: {r_id}")

        if not inplace:
            return constraints

    def simplify(self, inplace: bool = False) -> "Environment | None":
        """Remove reactions with blocked uptake."""

        if inplace:
            env = self
        else:
            env = self.copy()

        to_remove = []

        for r_id, (lb, _) in env.items():
            if lb >= 0:
                to_remove.append(r_id)

        for r_id in to_remove:
            del env[r_id]

        if not inplace:
            return env  # pyright: ignore[reportReturnType]

    @staticmethod
    def from_compounds(
        compounds: list, fmt_func: FunctionType, max_uptake: float = 10.0
    ) -> "Environment":
        """
        Initialize environment from list of medium compounds

        Arguments:
            compounds: list of compounds in the medium
            fmt_func: function to convert metabolite ids to exchange reaction ids
            max_uptake: maximum uptake rate for given compounds

        """

        env = Environment()

        for cpd in compounds:
            r_id = fmt_func(cpd)
            env[r_id] = (-max_uptake, inf)

        return env
