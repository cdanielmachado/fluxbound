# fluxbound

A minimalistic python package for flux balance analysis

![logo](logo.png)



## Installation

`pip install fluxbound`

### Solvers

SCIP solver is installed by default. You can optionally install GUROBI or CPLEX. 

If multiple solvers are installed, the selection order is: GUROBI > CPLEX > SCIP.

But you can change your prefered solver:

```python
from fluxbound import set_default_solver

set_default_solver('scip')
```

## Basic usage

### Loading

Loading (and inspecting) a model: 

```python
from fluxbound import load_model

model = load_model("e_coli_core.xml")
model

Model: e_coli_core
  - Genes: 137
  - Metabolites: 72
      c: 52
      e: 20
  - Reactions: 95
      internal: 50
      transport: 25
      exchange: 20
  - Obj: R_BIOMASS_Ecoli_core_w_GAM
```



### FBA

Running a basic FBA simulation:

```python
from fluxbound import FBA

model = load_model("e_coli_core.xml")
sol = FBA(model)
print(sol)

Objective: 0.8739
Status: Optimal
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

### Solution

The solution object has such utilities such a filtering using regular expressions: 

```python
sol.show_values(pattern='R_EX', sort=True, abstol=0.1)

R_EX_o2_e    -21.8
R_EX_glc__D_e -10
R_EX_nh4_e   -4.765
R_EX_pi_e    -3.215
R_EX_h_e      17.53
R_EX_co2_e    22.81
R_EX_h2o_e    29.18
```

And converting to pandas DataFrames

```python
df = sol.to_dataframe()
df.query("value > 30")

              value
R_ATPS4r  45.514010
R_CYTBD   43.598985
R_NADH16  38.534610
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

### Deletion

You can simulate the deletion of multiple genes:

```python
from fluxbound import gene_deletion

sol = gene_deletion(model, ['G_b0008', 'G_b1241', 'G_b4301'])
```

Or reactions:

```python
from fluxbound import reaction_deletion

sol = reaction_deletion(model, ['R_PGI', 'R_PFK'])
```

### Essentiality 

Calculate all essential genes or reactions: 

```python
from fluxbound import essential_genes

genes = essential_genes(model)
```

Remember that you can always use the constraints argument to change conditions:

```python
from fluxbound import essential_reactions

essential_anaerobic = essential_reactions(model, constraints={'R_EX_o2_e': 0})
```

## Advanced

### Environments

The Environment class is an utility for easily creating and changing environmental conditions.


You can create a temporary environment to use in a simulation:

```python
from fluxbound import Environment

env = Environment.complete(model, inplace=False)
sol = FBA(model, constraints=env)
```

Or permanently change a model:

```python
Environment.complete(model, inplace=True)
sol = FBA(model)
```

Environments are dictionary objects that you can easily manipulate:

```python
env = Environment.from_model(model)
env['R_EX_o2_e'] = (-20, 0)
env.simplify()

{'R_EX_co2_e': (-1000.0, 1000.0),
 'R_EX_glc__D_e': (-10.0, 1000.0),
 'R_EX_h_e': (-1000.0, 1000.0),
 'R_EX_h2o_e': (-1000.0, 1000.0),
 'R_EX_nh4_e': (-1000.0, 1000.0),
 'R_EX_o2_e': (-20, 0),
 'R_EX_pi_e': (-1000.0, 1000.0)}
```

### Model manipulation

If instead of just simulating some environmental or genetic perturbations you want to
actually modify a model, you can play around with the model object.

Make changes:

```python
model.reactions.R_PGI.lb = 0
model.set_flux_bounds('R_ATPM', lb=1, ub=1)
```

Add elements:

```python
from fluxbound import Metabolite, Reaction

x = Metabolite("M_x_e", name="compound x", compartment='e')
r = Reaction("R_EX_x_e", stoichiometry={'M_x_e': -1}, lb=0, ub=1000)
model.add_metabolite(x)
model.add_reaction(r)
print(model.reactions.R_EX_y_e)

R_EX_y_e: M_x_e -->  [0, 1000]
```

Remove elements:

```python
model.remove_reactions(
    ['R_EX_fru_e', 'R_FRUpts2'],
    clean_orphans=True
) # removes orphaned metabolites and genes 
```

Finally you can save your new model as SBML:

```python
from fluxbound import save_model

save_model(model, 'new_model.xml')
```
