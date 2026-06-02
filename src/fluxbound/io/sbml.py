import os
import re
from math import inf, isfinite
from warnings import warn

import libsbml as sb
from sympy import Symbol
from sympy.logic.boolalg import And, Boolean, Or, is_dnf, to_dnf

from ..core.model import (
    GPR,
    Base,
    Compartment,
    Gene,
    Metabolite,
    Model,
    Protein,
    Reaction,
    ReactionType,
)
from ..core.utils import AttrDict

IDENTIFIERS_PATTERN = re.compile(r"/([^/]+)/([^/]+)$")


def load_sbml(filename: str) -> sb.Model:

    if not os.path.exists(filename):
        raise IOError("Model file not found")

    reader = sb.SBMLReader()
    document = reader.readSBML(filename)
    sbml_model = document.getModel()

    if sbml_model is None:
        document.printErrors()
        raise IOError(f"Failed to load {filename}.")

    return sbml_model


def load_model(filename: str) -> Model:

    sbml_model = load_sbml(filename)
    fbc_model = sbml_model.getPlugin("fbc")
    is_fbc = fbc_model is not None

    params = {p.getId(): p.getValue() for p in sbml_model.getListOfParameters()}

    model = Model(sbml_model.getId())
    load_compartments(sbml_model, model)
    load_metabolites(sbml_model, model, is_fbc)
    if is_fbc:
        load_genes(sbml_model, model)
    load_reactions(sbml_model, model, is_fbc, params)
    if is_fbc:
        load_objective(sbml_model, model)

    return model


def load_compartments(sbml_model: sb.Model, model: Model) -> None:
    for compartment in sbml_model.getListOfCompartments():
        size = compartment.getSize()
        if not isfinite(size):
            size = 1.0
        external = False  # TODO: determine if compartment is external
        comp = Compartment(compartment.getId(), name=compartment.getName(), external=external, size=size)
        extract_metadata(compartment, comp)
        model.add_compartment(comp)


def load_metabolites(sbml_model: sb.Model, model: Model, is_fbc: bool) -> None:
    for species in sbml_model.getListOfSpecies():
        met = Metabolite(species.getId(), name=species.getName(), compartment=species.getCompartment())
        extract_metadata(species, met)
        if is_fbc:
            fbc_species = species.getPlugin("fbc")
            if fbc_species.isSetChemicalFormula():
                formula = fbc_species.getChemicalFormula()
                met.metadata["FORMULA"] = formula
            if fbc_species.isSetCharge():
                charge = fbc_species.getCharge()
                met.metadata["CHARGE"] = str(charge)
        model.add_metabolite(met)


def load_genes(sbml_model: sb.Model, model: Model) -> None:
    fbc_model = sbml_model.getPlugin("fbc")
    for fbc_gene in fbc_model.getListOfGeneProducts():
        gene = Gene(fbc_gene.getId(), fbc_gene.getName())
        extract_metadata(fbc_gene, gene)
        model.add_gene(gene)


def load_reactions(
    sbml_model: sb.Model, model: Model, is_fbc: bool, params: dict) -> None:
    for reaction in sbml_model.getListOfReactions():
        stoichiometry = load_stoichiometry(reaction)
        lb, ub = load_bounds(reaction, is_fbc, params)
        gpr = load_gpr(reaction) if is_fbc else None
        rtype = ReactionType.OTHER  # TODO: determine reaction type
        rxn = Reaction(
            reaction.getId(), name=reaction.getName(), stoichiometry=stoichiometry, lb=lb, ub=ub, gpr=gpr, rtype=rtype
        )
        extract_metadata(reaction, rxn)
        model.add_reaction(rxn)


def load_stoichiometry(reaction: sb.Reaction) -> AttrDict:
    stoichiometry = AttrDict()
    for reactant in reaction.getListOfReactants():
        m_id = reactant.getSpecies()
        coeff = -reactant.getStoichiometry()

        if m_id not in stoichiometry:
            stoichiometry[m_id] = coeff
        else:
            stoichiometry[m_id] += coeff

    for product in reaction.getListOfProducts():
        m_id = product.getSpecies()
        coeff = product.getStoichiometry()

        if m_id not in stoichiometry:
            stoichiometry[m_id] = coeff
        else:
            stoichiometry[m_id] += coeff
        if stoichiometry[m_id] == 0.0:
            del stoichiometry[m_id]

    return stoichiometry


def load_bounds(reaction: sb.Reaction, is_fbc: bool, params: dict) -> tuple[float, float]:
    if is_fbc:
        fbc_reaction = reaction.getPlugin("fbc")
        param_lb = fbc_reaction.getLowerFluxBound()
        param_ub = fbc_reaction.getUpperFluxBound()
        lb = params[param_lb]
        ub = params[param_ub]
    else:
        lb = -inf if reaction.getReversible() else 0.0
        ub = inf

    return lb, ub


def load_gpr(reaction: sb.Reaction) -> str:
    fbc_reaction = reaction.getPlugin("fbc")
    fbc_gpr = fbc_reaction.getGeneProductAssociation()

    if fbc_gpr is None:
        return None
    else:
        fbc_gpr = fbc_gpr.getAssociation() 

    try:
        gpr = easy_gpr_parse(fbc_gpr)
    except RuntimeError:
        gpr = hard_gpr_parse(fbc_gpr)

    return gpr


def easy_gpr_parse(fbc_gpr: sb.Association) -> GPR:
    gpr = GPR()

    if fbc_gpr.isFbcOr():
        for item in fbc_gpr.getListOfAssociations():
            protein = Protein()
            if item.isFbcAnd():
                for subitem in item.getListOfAssociations():
                    if subitem.isGeneProductRef():
                        protein.genes.append(subitem.getGeneProduct())
                    else:
                        raise RuntimeError("Unsupported GPR structure")
            elif item.isGeneProductRef():
                protein.genes.append(item.getGeneProduct())
            else:
                raise RuntimeError("Unsupported GPR structure")
            gpr.proteins.append(protein)
    elif fbc_gpr.isFbcAnd():
        protein = Protein()
        for item in fbc_gpr.getListOfAssociations():
            if item.isGeneProductRef():
                protein.genes.append(item.getGeneProduct())
            else:
                raise RuntimeError("Unsupported GPR structure")
        gpr.proteins = [protein]
    elif fbc_gpr.isGeneProductRef():
        protein = Protein()
        protein.genes = [fbc_gpr.getGeneProduct()]
        gpr.proteins = [protein]
    else:
        raise RuntimeError("Unsupported GPR structure")

    return gpr


def hard_gpr_parse(fbc_gpr: sb.Association) -> GPR:
    sympy_expr = fbc_association_to_sympy(fbc_gpr)

    if not is_dnf(sympy_expr):
        sympy_expr = to_dnf(sympy_expr)

    gpr = GPR()

    if type(sympy_expr) is Or:
        for sub_expr in sympy_expr.args:
            protein = Protein()
            if type(sub_expr) is And:
                protein.genes = [str(gene) for gene in sub_expr.args]
            else:
                protein.genes = [str(sub_expr)]
            gpr.proteins.append(protein)
    elif type(sympy_expr) is And:
        protein = Protein()
        protein.genes = [str(gene) for gene in sympy_expr.args]
        gpr.proteins = [protein]
    else:
        protein = Protein()
        protein.genes = [str(sympy_expr)]
        gpr.proteins = [protein]

    return gpr


def fbc_association_to_sympy(node: sb.Association) -> Boolean:

    # GeneProductRef
    if node.isGeneProductRef():
        gene_id = node.getGeneProduct()
        return Symbol(gene_id)

    # And
    if node.isFbcAnd():
        children = [
            fbc_association_to_sympy(node.getAssociation(i))
            for i in range(node.getNumAssociations())
        ]
        return And(*children)

    # Or
    if node.isFbcOr():
        children = [
            fbc_association_to_sympy(node.getAssociation(i))
            for i in range(node.getNumAssociations())
        ]
        return Or(*children)

    raise TypeError(f"Unsupported FBC association type: {node.getElementName()}")


def load_objective(sbml_model: sb.Model, model: Model) -> None:
    fbc_model = sbml_model.getPlugin("fbc")
    objective = fbc_model.getActiveObjective()

    for rxn_obj in objective.getListOfFluxObjectives():
        r_id = rxn_obj.getReaction()
        coeff = rxn_obj.getCoefficient()
        if coeff:
            model.objective[r_id] = coeff


def extract_metadata(sbml_elem: sb.SBase, elem: Base) -> None:
    # TODO: based on old code, check if refactoring is needed

    sboterm = sbml_elem.getSBOTermID()
    if sboterm:
        elem.metadata["SBOTerm"] = sboterm

    notes = sbml_elem.getNotes()
    if notes:
        recursive_node_parser(notes, elem.metadata)

    parse_annotations(sbml_elem, elem)


def parse_annotations(sbml_elem: sb.SBase, elem: Base) -> None:
    # TODO: implemented by Snorre for reframed, check if refactoring is needed

    """
    Parse the annotations found in XML Annotations like
    http://identifiers.org/reactome/R-ALL-419151

    There is a question whether all annotations should be added as lists, or
    only in those cases where you have multiple annotations from the same db,
    and then keep the rest as strings. Now, all annotations gets added as lists.

    """
    for term in sbml_elem.getCVTerms():
        n_resources = term.getNumResources()
        for i in range(n_resources):
            match = IDENTIFIERS_PATTERN.search(term.getResourceURI(i))
            if match:
                db = match.group(1)
                db_id = match.group(2)

                try:
                    elem.metadata[db]
                except KeyError:
                    elem.metadata[db] = [db_id]
                else:
                    elem.metadata[db].append(db_id)

            else:
                warn(f"Could not extract annotation from {term.getResourceURI(i)}")


def recursive_node_parser(node: sb.SBase, cache: dict) -> None:
    node_data = node.getCharacters()
    if ":" in node_data:
        key, value = node_data.split(":", 1)
        cache[key.strip()] = value.strip()

    for i in range(node.getNumChildren()):
        recursive_node_parser(node.getChild(i), cache)
