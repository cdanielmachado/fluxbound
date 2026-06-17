import os
import re
from collections import Counter
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
    document = reader.readSBML(filename)  # pyright: ignore[reportArgumentType]
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
    apply_heuristics(model)

    return model


def load_compartments(sbml_model: sb.Model, model: Model) -> None:
    for compartment in sbml_model.getListOfCompartments():
        size = compartment.getSize()
        if not isfinite(size):
            size = 1.0
        external = False
        comp = Compartment(
            compartment.getId(),
            name=compartment.getName(),
            external=external,
            size=size,
        )
        extract_metadata(compartment, comp)
        model.add_compartment(comp)


def load_metabolites(sbml_model: sb.Model, model: Model, is_fbc: bool) -> None:
    for species in sbml_model.getListOfSpecies():
        met = Metabolite(
            species.getId(),
            name=species.getName(),
            compartment=species.getCompartment(),
        )
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
    sbml_model: sb.Model, model: Model, is_fbc: bool, params: dict
) -> None:
    for reaction in sbml_model.getListOfReactions():
        stoichiometry = load_stoichiometry(reaction)
        lb, ub = load_bounds(reaction, is_fbc, params)
        gpr = load_gpr(reaction) if is_fbc else None
        rtype = ReactionType.OTHER
        rxn = Reaction(
            reaction.getId(),
            name=reaction.getName(),
            stoichiometry=stoichiometry,
            lb=lb,
            ub=ub,
            gpr=gpr,
            rtype=rtype,
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


def load_bounds(
    reaction: sb.Reaction, is_fbc: bool, params: dict
) -> tuple[float, float]:
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


def load_gpr(reaction: sb.Reaction) -> GPR | None:
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
        for item in fbc_gpr.getListOfAssociations():  # pyright: ignore[reportAttributeAccessIssue]
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
        for item in fbc_gpr.getListOfAssociations():  # pyright: ignore[reportAttributeAccessIssue]
            if item.isGeneProductRef():
                protein.genes.append(item.getGeneProduct())
            else:
                raise RuntimeError("Unsupported GPR structure")
        gpr.proteins = [protein]
    elif fbc_gpr.isGeneProductRef():
        protein = Protein()
        protein.genes = [fbc_gpr.getGeneProduct()]  # pyright: ignore[reportAttributeAccessIssue]
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
        gene_id = node.getGeneProduct()  # pyright: ignore[reportAttributeAccessIssue]
        return Symbol(gene_id)

    # And
    if node.isFbcAnd():
        children = [
            fbc_association_to_sympy(node.getAssociation(i))  # pyright: ignore[reportAttributeAccessIssue]
            for i in range(node.getNumAssociations())
        ]
        return And(*children)

    # Or
    if node.isFbcOr():
        children = [
            fbc_association_to_sympy(node.getAssociation(i))  # pyright: ignore[reportAttributeAccessIssue]
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


def recursive_node_parser(node: sb.XMLNode, cache: dict) -> None:
    node_data = node.getCharacters()
    if ":" in node_data:
        key, value = node_data.split(":", 1)
        cache[key.strip()] = value.strip()

    for i in range(node.getNumChildren()):
        recursive_node_parser(node.getChild(i), cache)


def apply_heuristics(model: Model) -> None:
    counter = Counter()

    for rxn in model.reactions.values():
        if len(rxn.stoichiometry) == 1:
            m_id = list(rxn.stoichiometry.keys())[0]
            comp = model.metabolites[m_id].compartment
            counter[comp] += 1

    ext_comp = counter.most_common(1)[0][0]

    for c_id, comp in model.compartments.items():
        comp.external = bool(c_id == ext_comp)

    for rxn in model.reactions.values():
        subs = [m_id for m_id, coeff in rxn.stoichiometry.items() if coeff < 0]
        prods = [m_id for m_id, coeff in rxn.stoichiometry.items() if coeff > 0]
        if len(subs) == 0 or len(prods) == 0:
            if len(rxn.stoichiometry) == 1:
                m_id = list(rxn.stoichiometry.keys())[0]
                comp = model.metabolites[m_id].compartment
                if model.compartments[comp].external:
                    rxn.rtype = ReactionType.EXCHANGE
                else:
                    rxn.rtype = ReactionType.SINK_OR_DEMAND
        else:
            comps = {model.metabolites[m_id].compartment for m_id in rxn.stoichiometry}
            if len(comps) == 1:
                rxn.rtype = ReactionType.INTERNAL
            else:
                rxn.rtype = ReactionType.TRANSPORT


def save_model(model: Model, filename: str) -> None:

    document = sb.SBMLDocument()
    sbml_model = document.createModel(model.id)
    document.enablePackage(sb.FbcExtension.getXmlnsL3V1V2(), "fbc", True)
    fbc_model = sbml_model.getPlugin("fbc")
    fbc_model.setStrict(True)
    document.setPackageRequired("fbc", False)

    save_compartments(model, sbml_model)
    save_metabolites(model, sbml_model)
    save_genes(model, sbml_model)
    save_reactions(model, sbml_model)
    save_objective(model, sbml_model)
    save_fluxbounds(model, sbml_model)
    save_metadata(model, sbml_model)
    writer = sb.SBMLWriter()
    writer.writeSBML(document, filename)


def save_compartments(model: Model, sbml_model: sb.Model) -> None:
    for comp in model.compartments.values():
        sbml_comp = sbml_model.createCompartment()
        sbml_comp.setId(comp.id)
        sbml_comp.setName(comp.name)
        sbml_comp.setSize(comp.size)
        sbml_comp.setConstant(True)
        save_metadata(comp, sbml_comp)


def save_metabolites(model: Model, sbml_model: sb.Model) -> None:
    for met in model.metabolites.values():
        sbml_met = sbml_model.createSpecies()
        sbml_met.setId(met.id)
        sbml_met.setName(met.name)
        sbml_met.setCompartment(met.compartment)
        sbml_met.setConstant(False)
        sbml_met.setHasOnlySubstanceUnits(True)
        sbml_met.setBoundaryCondition(False)

        fbc_species = sbml_met.getPlugin("fbc")

        if "FORMULA" in met.metadata:
            formula = met.metadata["FORMULA"]
            try:
                fbc_species.setChemicalFormula(formula)
            except Exception as e:
                warn(f"Failed to set formula '{formula}' for {met.id}: {e}")
        if "CHARGE" in met.metadata:
            try:
                charge = int(met.metadata["CHARGE"])
                fbc_species.setCharge(charge)
            except Exception as e:
                warn(
                    f"Failed to set charge '{met.metadata['CHARGE']}' for {met.id}: {e}"
                )

        save_metadata(met, sbml_met)


def save_genes(model: Model, sbml_model: sb.Model) -> None:
    fbc_model = sbml_model.getPlugin("fbc")
    for gene in model.genes.values():
        fbc_gene = fbc_model.createGeneProduct()
        fbc_gene.setId(gene.id)
        fbc_gene.setName(gene.name)
        fbc_gene.setLabel(gene.name)
        save_metadata(gene, fbc_gene)


def save_reactions(model: Model, sbml_model: sb.Model) -> None:
    for rxn in model.reactions.values():
        sbml_rxn = sbml_model.createReaction()
        sbml_rxn.setId(rxn.id)
        sbml_rxn.setName(rxn.name)
        sbml_rxn.setReversible(rxn.lb < 0)
        sbml_rxn.setFast(False)

        for m_id, coeff in rxn.stoichiometry.items():
            if coeff < 0:
                reactant = sbml_rxn.createReactant()
                reactant.setSpecies(m_id)
                reactant.setStoichiometry(-coeff)
                reactant.setConstant(True)
            elif coeff > 0:
                product = sbml_rxn.createProduct()
                product.setSpecies(m_id)
                product.setStoichiometry(coeff)
                product.setConstant(True)

        if rxn.gpr is not None:
            fbcrxn = sbml_rxn.getPlugin("fbc")
            gpr_assoc = fbcrxn.createGeneProductAssociation()

            if len(rxn.gpr.proteins) > 1:
                gpr_assoc = gpr_assoc.createOr()

            for protein in rxn.gpr.proteins:
                if len(protein.genes) > 1:
                    protein_assoc = gpr_assoc.createAnd()
                else:
                    protein_assoc = gpr_assoc

                for gene in protein.genes:
                    gene_ref = protein_assoc.createGeneProductRef()
                    gene_ref.setGeneProduct(gene)

        save_metadata(rxn, sbml_rxn)


def save_objective(model: Model, sbml_model: sb.Model) -> None:
    fbcmodel = sbml_model.getPlugin("fbc")
    obj = fbcmodel.createObjective()
    obj.setId("objective")
    fbcmodel.setActiveObjectiveId("objective")
    obj.setType("maximize")
    for r_id, coeff in model.objective.items():
        if coeff != 0:
            r_obj = obj.createFluxObjective()
            r_obj.setReaction(r_id)
            r_obj.setCoefficient(coeff)


def save_fluxbounds(model: Model, sbml_model: sb.Model) -> None:
    lb_inf = sbml_model.createParameter()
    lb_inf.setId("LB_INF")
    lb_inf.setValue(-inf)
    lb_inf.setConstant(True)

    ub_inf = sbml_model.createParameter()
    ub_inf.setId("UB_INF")
    ub_inf.setValue(inf)
    ub_inf.setConstant(True)

    zero = sbml_model.createParameter()
    zero.setId("ZERO")
    zero.setValue(0.0)
    zero.setConstant(True)

    for r_id, rxn in model.reactions.items():
        fbcrxn = sbml_model.getReaction(r_id).getPlugin("fbc")

        if rxn.lb == -inf:
            fbcrxn.setLowerFluxBound("LB_INF")
        elif rxn.lb == 0:
            fbcrxn.setLowerFluxBound("ZERO")
        else:
            lb_id = f"{r_id}_lb"
            lb_param = sbml_model.createParameter()
            lb_param.setId(lb_id)
            lb_param.setValue(rxn.lb)
            lb_param.setConstant(True)
            fbcrxn.setLowerFluxBound(lb_id)

        if rxn.ub == inf:
            fbcrxn.setUpperFluxBound("UB_INF")
        elif rxn.ub == 0:
            fbcrxn.setUpperFluxBound("ZERO")
        else:
            ub_id = f"{r_id}_ub"
            ub_param = sbml_model.createParameter()
            ub_param.setId(ub_id)
            ub_param.setValue(rxn.ub)
            ub_param.setConstant(True)
            fbcrxn.setUpperFluxBound(ub_id)


def save_metadata(elem: Base, sbml_elem: sb.SBase) -> None:
    # TODO: based on old code, check if refactoring is needed

    meta_id = f"meta_{sbml_elem.getId()}"
    sbml_elem.setMetaId(meta_id)
    note_keys = ["CHARGE", "FORMULA"]
    notes_dict = {}

    if elem.metadata:
        for key, annotations in elem.metadata.items():
            if key == "SBOTerm":
                sbml_elem.setSBOTerm(annotations)
            elif key in note_keys:
                notes_dict[key] = annotations
            elif key == "XMLAnnotation":
                continue
            else:
                # Assume this is an annotation
                if not isinstance(annotations, list):
                    annotations = [annotations]
                for annotation in annotations:
                    if annotation:
                        cv = sb.CVTerm()
                        cv.setQualifierType(sb.BIOLOGICAL_QUALIFIER)
                        cv.setBiologicalQualifierType(sb.BQB_IS)
                        annotation_string = (
                            f"https://identifiers.org/{key}/{annotation}"
                        )
                        cv.addResource(annotation_string)
                        sbml_elem.addCVTerm(cv)

    if len(notes_dict):
        notes = [
            f"<p>{key}: {re.escape(value)}</p>"
            for key, value in notes_dict.items()
            if key != "XMLAnnotation"
        ]
        note_string = "<html>" + "".join(notes) + "</html>"
        note_xml = sb.XMLNode.convertStringToXMLNode(note_string)
        note_xml.getNamespaces().add("http://www.w3.org/1999/xhtml")
        sbml_elem.setNotes(note_xml)
