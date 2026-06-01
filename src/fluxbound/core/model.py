from enum import Enum

from .utils import valid_sbml_id, AttrDict
from math import inf

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
    
    def __init__(self, comp_id: str, name: str | None = None, external: bool = False, 
                 size: float = 1.0) -> None:

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


class Protein():

    def __init__(self):
        self.genes = []

    def __str__(self):
        protein_str = ' and '.join(self.genes)

        if len(self.genes) > 1:
            protein_str = '(' + protein_str + ')'

        return protein_str

    def __repr__(self):
        return str(self)
    

class GPR:
    def __init__(self):
        self.proteins = []

    def __str__(self):

        gpr_str = ' or '.join(map(str, self.proteins))

        if len(self.proteins) > 1:
            gpr_str = '(' + gpr_str + ')'

        return gpr_str

    def __repr__(self):
        return str(self)

    def get_genes(self):
        return [gene for protein in self.proteins for gene in protein.genes]


class ReactionType(Enum):
    """ Enumeration of possible reaction types. """
    ENZYMATIC = 'enzymatic'
    TRANSPORT = 'transport'
    EXCHANGE = 'exchange'
    UNBALANCED = 'unbalanced'
    OTHER = 'other'


class Reaction(Base):
    
    def __init__(self, rxn_id: str, name: str | None = None, stoichiometry: AttrDict | None = None, 
                 lb: float = -inf, ub: float = inf, gpr: GPR | None = None, rtype: ReactionType = ReactionType.OTHER) -> None:
        
        super().__init__(rxn_id, name)
        self.stoichiometry: AttrDict = stoichiometry if stoichiometry is not None else AttrDict()
        
        if ub < lb:
            raise ValueError(f"Upper bound ({ub}) cannot be less than lower bound ({lb})")
        
        self.lb: float = lb
        self.ub: float = ub
        self.gpr: GPR | None  = gpr
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

        def format_coeff(coeff, m_id):
            return m_id if coeff == 1.0 else f"{coeff} {m_id}"
            
        lhs = ' + '.join(format_coeff(-coeff, m_id) 
                         for m_id, coeff in self.stoichiometry.items() if coeff < 0)
        rhs = ' + '.join(format_coeff(coeff, m_id) 
                         for m_id, coeff in self.stoichiometry.items() if coeff > 0)
        
        if self.lb < 0:
            if self.ub <= 0:
                arrow = '<--'
            else:
                arrow = '<=>'
        else:
            arrow = '-->'
        
        if self.lb not in {-inf, 0} or self.ub not in {inf, 0}:
            bounds = f" [{self.lb}, {self.ub}]"
        else:
            bounds = ""

        return f"{lhs} {arrow} {rhs}{bounds}"

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

