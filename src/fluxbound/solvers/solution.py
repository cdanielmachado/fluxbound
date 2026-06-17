import re
from enum import Enum
from typing import Self
from warnings import warn

import pandas as pd


class Status(Enum):
    """Enumeration of possible solution status."""

    OPTIMAL = "Optimal"
    UNKNOWN = "Unknown"
    SUBOPTIMAL = "Suboptimal"
    UNBOUNDED = "Unbounded"
    INFEASIBLE = "Infeasible"
    INF_OR_UNB = "Infeasible or Unbounded"


def pretty_print_values(
    values: dict, pattern: str | None = None, sort: bool = False, abstol: float = 1e-9
):

    values = {key: value for key, value in values.items() if abs(value) > abstol}

    if pattern:
        re_expr = re.compile(pattern)
        values = {
            key: value
            for key, value in values.items()
            if re_expr.search(key) is not None
        }

    if sort:
        values = dict(sorted(values.items(), key=lambda x: x[1]))

    entries = (f"{key:<12} {value: .4g}" for key, value in values.items())

    print("\n".join(entries))


class Solution:
    def __init__(
        self,
        status: Status = Status.UNKNOWN,
        message: str | None = None,
        fobj: float | None = None,
        values: dict | None = None,
        shadow_prices: dict | None = None,
    ) -> None:
        self.status: Status = status
        self.message: str | None = message
        self.fobj: float | None = fobj
        self.values: dict | None = values
        self.shadow_prices: dict | None = shadow_prices

    def __str__(self) -> str:
        if self.fobj is None:
            return f"Status: {self.status.value}\n"
        else:
            return f"Objective: {self.fobj:.4g}\nStatus: {self.status.value}\n"

    def __repr__(self) -> str:
        return str(self)

    def show_values(
        self, pattern: str | None = None, sort: bool = False, abstol: float = 1e-9
    ) -> None:

        if self.values is None:
            warn("No solution to show")
        else:
            pretty_print_values(self.values, pattern=pattern, sort=sort, abstol=abstol)

    def compare(
        self,
        other: Self,
        intersect=True,
        pattern: str | None = None,
        sort: bool = False,
        abstol: float = 1e-9,
    ) -> None:
        if self.values is None or other.values is None:
            warn("One or both solutions are empty")
        else:
            values = {}
            for key, val in self.values.items():
                if key in other.values:
                    values[key] = val - other.values[key]
                elif not intersect:
                    values[key] = val

            pretty_print_values(values, pattern=pattern, sort=sort, abstol=abstol)

    def to_dataframe(self) -> pd.DataFrame | None:
        if self.values is None:
            return None
        else:
            return pd.DataFrame(
                list(self.values.values()),
                columns=["value"],
                index=list(self.values.keys()),
            )
