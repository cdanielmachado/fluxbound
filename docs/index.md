

# fluxbound

A minimalistic python package for flux balance analysis


## Installation

`pip install fluxbound`

## Basic usage

Load a model and run flux balance analysis:

```python
from fluxbound import load_model, FBA

model = load_model("e_coli_core.xml")
sol = FBA(model)
print(sol)
```
