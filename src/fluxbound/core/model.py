from .utils import valid_sbml_id


class Base:

    def __init__(self, elem_id: str, name: str | None = None) -> None:

        if not valid_sbml_id(elem_id):
            raise ValueError(f"'{elem_id}' is not a valid SBML identifier")
        
        self.id: str = elem_id
        self.name: str = name if name is not None else elem_id
        self.metadata: dict = {}

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return str(self)
    

class Compartment(Base):
    
    def __init__(self, comp_id: str, name: str | None = None, external: bool = False, 
                 size: float = 1.0) -> None:
        """
        Arguments:
            comp_id (str): a valid unique identifier
            name (str): compartment name (optional)
            external (bool): is external (default: false)
            size (float): compartment size (default: 1.0)
        """
        super().__init__(elem_id=comp_id, name=name)
        self.size: float = size
        self.external = external


class Metabolite(Base):

    def __init__(self, met_id: str, compartment: str, name: str | None = None,) -> None:
        """
        Arguments:
            met_id (str): a valid unique identifier
            name (str): metabolite name (optional)
            compartment (str): the compartment the metabolite is in (optional)
        """
        super().__init__(elem_id=met_id, name=name)
        self.compartment: str = compartment




