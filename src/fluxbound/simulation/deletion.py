from ..core.model import Model
from ..solvers import solver_instance
from ..solvers.solution import Solution
from ..solvers.solver import Solver, Status
from .fba import FBA


def gene_deletion(
    model: Model,
    genes: str | list,
    parsimonious: bool = False,
    constraints: dict | None = None,
    skip_silent: bool = False,
    solver: Solver | None = None,
) -> Solution | None:

    if isinstance(genes, str):
        genes = [genes]

    del_rxns = deleted_genes_to_reactions(model, genes)

    if len(del_rxns) == 0 and skip_silent:
        return

    return reaction_deletion(model, del_rxns, parsimonious, constraints, solver)


def deleted_genes_to_reactions(model: Model, genes: list) -> list:

    gene_bool = {gene: gene not in genes for gene in model.genes}

    del_rxns = []
    for r_id, rxn in model.reactions.items():
        if rxn.gpr is not None:
            if not rxn.gpr.eval_rule(gene_bool):
                del_rxns.append(r_id)

    return del_rxns


def reaction_deletion(
    model: Model,
    reactions: str | list,
    parsimonious: bool = False,
    constraints: dict | None = None,
    solver: Solver | None = None,
) -> Solution:

    if isinstance(reactions, str):
        reactions = [reactions]

    del_and_env = {}

    if constraints is not None:
        del_and_env.update(constraints)

    for r_id in reactions:
        del_and_env[r_id] = 0

    return FBA(model, parsimonious=parsimonious, constraints=del_and_env, solver=solver)


def essential_genes(
    model: Model, growth_frac: float = 0.01, constraints: dict | None = None
) -> list:

    solver = solver_instance(model)
    sol0 = FBA(model, constraints=constraints, solver=solver)

    if sol0.status != Status.OPTIMAL:
        raise RuntimeError("Reference condition is infeasible.")

    essential = []

    for gene in model.genes.keys():
        sol = gene_deletion(
            model, gene, constraints=constraints, solver=solver, skip_silent=True
        )
        if sol is None:
            continue
        elif sol.status == Status.INFEASIBLE:
            essential.append(gene)
        elif sol.status == Status.OPTIMAL:
            if sol.fobj < growth_frac * sol0.fobj:  # pyright: ignore
                essential.append(gene)

    return essential


def essential_reactions(
    model: Model, growth_frac: float = 0.01, constraints: dict | None = None
) -> list:

    solver = solver_instance(model)
    sol0 = FBA(model, constraints=constraints, solver=solver)

    if sol0.status != Status.OPTIMAL:
        raise RuntimeError("Reference condition is infeasible.")

    essential = []

    for r_id in model.reactions.keys():
        sol = reaction_deletion(model, r_id, constraints=constraints, solver=solver)
        if sol.status == Status.INFEASIBLE:
            essential.append(r_id)
        elif sol.status == Status.OPTIMAL:
            if sol.fobj < growth_frac * sol0.fobj:  # pyright: ignore
                essential.append(r_id)

    return essential
