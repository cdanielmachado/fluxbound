import re
from enum import Enum

import pandas as pd


class Status(Enum):
    """Enumeration of possible solution status."""

    OPTIMAL = "Optimal"
    UNKNOWN = "Unknown"
    SUBOPTIMAL = "Suboptimal"
    UNBOUNDED = "Unbounded"
    INFEASIBLE = "Infeasible"
    INF_OR_UNB = "Infeasible or Unbounded"


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

    def __str__(self):
        return f"Objective: {self.fobj}\nStatus: {self.status.value}\n"

    def __repr__(self):
        return str(self)

    def show_values(
        self, pattern: str | None = None, sort: bool = False, abstol: float = 1e-9
    ):

        if self.values is None:
            return

        values = [
            (key, value) for key, value in self.values.items() if abs(value) > abstol
        ]

        if pattern:
            re_expr = re.compile(pattern)
            values = [x for x in values if re_expr.search(x[0]) is not None]

        if sort:
            values.sort(key=lambda x: x[1])

        entries = (f"{r_id:<12} {val: .6g}" for (r_id, val) in values)

        print("\n".join(entries))

    def to_dataframe(self):
        return pd.DataFrame(
            self.values.values(), columns=["value"], index=self.values.keys()
        )
