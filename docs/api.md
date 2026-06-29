# API Reference


## SBML

::: fluxbound.io.sbml
    options:
      show_if_no_docstring: false
      members:
        - load_model
        - save_model


## Simulation

::: fluxbound.simulation.fba.FBA

::: fluxbound.simulation.fva.FVA


## Deletion

::: fluxbound.simulation.deletion
    options:
      members:
        - gene_deletion
        - reaction_deletion
        - essential_genes
        - essential_reactions

## Base classes

::: fluxbound.core.model
    options:
      merge_init_into_class: true
      show_bases: false
      show_if_no_docstring: false
      members:
        - Gene
        - Metabolite
        - Compartment
        - Reaction
        - Model

::: fluxbound.core.environment.Environment
    options:
      merge_init_into_class: true
      show_bases: false
      show_if_no_docstring: false
      members:
        - empty
        - complete
        - from_compounds
        - from_model
        - apply
        - simplify