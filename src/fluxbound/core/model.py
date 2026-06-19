from collections import Counter
from enum import Enum
from math import inf
from warnings import warn

from .utils import AttrDict, print_linear_expr, valid_sbml_id


class Base:
    def __init__(self, elem_id: str, name: str | None = None) -> None:

        if not valid_sbml_id(elem_id):
            raise ValueError(f"'{elem_id}' is not a valid SBML identifier")

        self.id: str = elem_id
        self.name: str = name if name is not None else elem_id
        self.metadata: AttrDict = AttrDict()

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return str(self)


class Compartment(Base):
    def __init__(
        self,
        comp_id: str,
        name: str | None = None,
        external: bool = False,
        size: float = 1.0,
    ) -> None:

        super().__init__(elem_id=comp_id, name=name)
        self.size: float = size
        self.external = external


class Metabolite(Base):
    def __init__(self, met_id: str, compartment: str, name: str | None = None) -> None:
        super().__init__(elem_id=met_id, name=name)
        self.compartment: str = compartment


class Gene(Base):
    def __init__(self, gene_id: str, name: str | None = None) -> None:
        super().__init__(elem_id=gene_id, name=name)


class Protein:
    def __init__(self):
        self.genes = []

    def __str__(self):
        protein_str = " and ".join(self.genes)

        if len(self.genes) > 1:
            protein_str = "(" + protein_str + ")"

        return protein_str

    def __repr__(self):
        return str(self)


class GPR:
    def __init__(self):
        self.proteins = []

    def __str__(self):

        gpr_str = " or ".join(map(str, self.proteins))

        if len(self.proteins) > 1:
            gpr_str = "(" + gpr_str + ")"

        return gpr_str

    def __repr__(self):
        return str(self)

    def get_genes(self):
        return {gene for protein in self.proteins for gene in protein.genes}


class ReactionType(Enum):
    INTERNAL = "internal"
    TRANSPORT = "transport"
    EXCHANGE = "exchange"
    SINK_OR_DEMAND = "sink/demand"
    OTHER = "other"


class Reaction(Base):
    def __init__(
        self,
        rxn_id: str,
        name: str | None = None,
        stoichiometry: AttrDict | None = None,
        lb: float = -inf,
        ub: float = inf,
        gpr: GPR | None = None,
        rtype: ReactionType = ReactionType.OTHER,
    ) -> None:

        super().__init__(rxn_id, name)
        self.stoichiometry: AttrDict = (
            stoichiometry if stoichiometry is not None else AttrDict()
        )

        if ub < lb:
            raise ValueError(
                f"Upper bound ({ub}) cannot be less than lower bound ({lb})"
            )

        self.lb: float = lb
        self.ub: float = ub
        self.gpr: GPR | None = gpr
        self.rtype: ReactionType = rtype

    def get_substrates(self):
        return [m_id for m_id, coeff in self.stoichiometry.items() if coeff < 0]

    def get_products(self):
        return [m_id for m_id, coeff in self.stoichiometry.items() if coeff > 0]

    def get_genes(self):
        if self.gpr is not None:
            return self.gpr.get_genes()
        else:
            return []

    def to_string(self) -> str:

        lhs = {m_id: -coeff for m_id, coeff in self.stoichiometry.items() if coeff < 0}
        lhs = print_linear_expr(lhs)

        rhs = {m_id: coeff for m_id, coeff in self.stoichiometry.items() if coeff > 0}
        rhs = print_linear_expr(rhs)

        if self.lb < 0:
            if self.ub <= 0:
                arrow = "<--"
            else:
                arrow = "<=>"
        else:
            arrow = "-->"

        if self.lb not in {-inf, 0} or self.ub not in {inf, 0}:
            bounds = f" [{self.lb}, {self.ub}]"
        else:
            bounds = ""

        return f"{self.id}: {lhs} {arrow} {rhs}{bounds}"

    def __str__(self) -> str:
        return self.to_string()


class Model(Base):
    def __init__(self, model_id: str, name: str | None = None) -> None:

        super().__init__(elem_id=model_id, name=name)
        self.compartments: AttrDict = AttrDict()
        self.metabolites: AttrDict = AttrDict()
        self.genes: AttrDict = AttrDict()
        self.reactions: AttrDict = AttrDict()
        self.objective: AttrDict = AttrDict()
        self._updated: bool = True
        self._met_rxn_lookup: dict | None = None
        self._gene_rxn_lookup: dict | None = None

    def add_compartment(self, compartment: Compartment, replace: bool = False) -> None:
        self._updated = False

        if compartment.id in self.compartments and not replace:
            raise RuntimeError(f"Compartment {compartment.id} already exists.")

        self.compartments[compartment.id] = compartment

    def add_metabolite(self, met: Metabolite, replace: bool = False) -> None:
        self._updated = False

        if met.id in self.metabolites and not replace:
            raise RuntimeError(f"Metabolite {met.id} already exists.")

        if met.compartment not in self.compartments:
            raise RuntimeError(
                f"Metabolite {met.id} has invalid compartment {met.compartment}."
            )

        self.metabolites[met.id] = met

    def add_gene(self, gene: Gene, replace: bool = False) -> None:
        self._updated = False

        if gene.id in self.genes and not replace:
            raise RuntimeError(f"Gene {gene.id} already exists.")

        self.genes[gene.id] = gene

    def add_reaction(self, rxn: Reaction, replace: bool = False) -> None:
        self._updated = False

        if rxn.id in self.reactions and not replace:
            raise RuntimeError(f"Reaction {rxn.id} already exists.")

        for m_id in rxn.stoichiometry.keys():
            if m_id not in self.metabolites:
                raise RuntimeError(f"Reaction {rxn.id} has invalid metabolite {m_id}.")

        self.reactions[rxn.id] = rxn

    def remove_metabolites(self, m_ids: list, include_reactions: bool = False) -> None:

        m_r_lookup = self.metabolite_reaction_lookup()
        self._updated = False

        deleted = []
        for m_id in m_ids:
            if m_id in self.metabolites:
                del self.metabolites[m_id]
                deleted.append(m_id)
            else:
                warn(f"Metabolite {m_id} not in model")

        if include_reactions:
            r_ids = {r_id for m_id in deleted for r_id in m_r_lookup[m_id].keys()}
            if len(r_ids) > 0:
                self.remove_reactions(list(r_ids))

    def remove_reactions(self, r_ids: list, clean_orphans: bool = False) -> None:
        self._updated = False

        deleted = []
        for r_id in list(r_ids):
            if r_id in self.reactions:
                del self.reactions[r_id]
                deleted.append(r_id)
            else:
                warn(f"Reaction {r_id} not in model")

        if clean_orphans:
            m_r_lookup = self.metabolite_reaction_lookup()
            orphan_mets = {
                m_id for m_id, r_ids in m_r_lookup.items() if len(r_ids) == 0
            }
            if len(orphan_mets) > 0:
                self.remove_metabolites(list(orphan_mets), include_reactions=False)

            g_r_lookup = self.gene_reaction_lookup()
            for g_id, rxns in g_r_lookup.items():
                if len(rxns) == 0:
                    del self.genes[g_id]

    def print(self) -> None:
        for rxn in self.reactions.values():
            print(rxn.to_string())

    def summary(self) -> str:
        repr = f"Model: {self.id}\n"
        repr += f"  - Genes: {len(self.genes)}\n"
        repr += f"  - Metabolites: {len(self.metabolites)}\n"

        count_mets = Counter()
        for met in self.metabolites.values():
            count_mets[met.compartment] += 1
        for c_id, total in count_mets.most_common():
            repr += f"      {c_id}: {total}\n"

        repr += f"  - Reactions: {len(self.reactions)}\n"

        count_rxns = Counter()
        for rxn in self.reactions.values():
            count_rxns[rxn.rtype.value] += 1
        for rtype, total in count_rxns.most_common():
            repr += f"      {rtype}: {total}\n"

        repr += f"  - Obj: {print_linear_expr(self.objective)}\n"
        return repr

    def __str__(self) -> str:
        return self.summary()

    def update(self) -> None:
        self._updated = True
        self._met_rxn_lookup = None
        self._gene_rxn_lookup = None

    def metabolite_reaction_lookup(self) -> dict:

        if not self._updated:
            self.update()

        if self._met_rxn_lookup is None:
            self._met_rxn_lookup = {m_id: {} for m_id in self.metabolites}

            for r_id, reaction in self.reactions.items():
                for m_id, coeff in reaction.stoichiometry.items():
                    self._met_rxn_lookup[m_id][r_id] = coeff

        return self._met_rxn_lookup

    def gene_reaction_lookup(self) -> dict:

        if not self._updated:
            self.update()

        if self._gene_rxn_lookup is None:
            self._gene_rxn_lookup = {g_id: [] for g_id in self.genes}

            for r_id, reaction in self.reactions.items():
                if reaction.gpr is not None:
                    for g_id in reaction.gpr.get_genes():
                        self._gene_rxn_lookup[g_id].append(r_id)

        return self._gene_rxn_lookup

    def get_reactions_by_type(self, reaction_type: ReactionType) -> list:
        return [rxn.id for rxn in self.reactions.values() if rxn.rtype == reaction_type]

    def get_exchange_reactions(self) -> list:
        return self.get_reactions_by_type(ReactionType.EXCHANGE)

    def get_biomass(self) -> str:

        # 1st heuristic: objective function
        if len(self.objective) == 1:
            return list(self.objective.keys())[0]

        # 2nd heuristic: most substrates
        n_subs = [
            (r_id, len(rxn.get_substrates())) for r_id, rxn in self.reactions.items()
        ]
        n_subs.sort(key=lambda x: -x[1])
        return n_subs[0][0]

    def set_flux_bounds(
        self, r_id, lb: float | None = None, ub: float | None = None
    ) -> None:

        if r_id not in self.reactions:
            warn(f"Reaction {r_id} not found")
            return

        if lb is not None:
            self.reactions[r_id].lb = lb

        if ub is not None:
            self.reactions[r_id].ub = ub
