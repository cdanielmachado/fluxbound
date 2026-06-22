

# fluxbound

A minimalistic python package for flux balance analysis


## Installation

`pip install fluxbound`

## Solvers

SCIP solver is installed by default. You can optionally install GUROBI or CPLEX. 

If multiple solvers are installed, the selection order is: GUROBI > CPLEX > SCIP.

You can change this with `set_default_solver(solvername)`. 

## Basic usage


### FBA

Load a model and run flux balance analysis:

```python
from fluxbound import load_model, FBA

model = load_model("e_coli_core.xml")
sol = FBA(model)
print(sol)
```

You can switch to pFBA:

```python
sol = FBA(model, parsimonious=True)
```

You can add additional constraints:

```python
constr = {
    'R_EX_glc__D_e': (-5, 0),
    'R_EX_o2_e': 0,
}

sol = FBA(model, constraints=constr)
```

And easily change the objective:

```python
sol = FBA(model, objective='R_ATPM', minimize=True)
```

### FVA

To run flux variability analysis:

```python
from fluxbound import FVA

ranges = FVA(model)
```

FVA at 90% maximum growth in anaerobic conditions:

```python
ranges = FVA(
    model, 
    obj_frac = 0.9,
    constraints = {'R_EX_o2_e': 0},
)
```