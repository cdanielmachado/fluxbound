from pathlib import Path

from fluxbound import (
    essential_genes,
    essential_reactions,
    load_model,
    set_default_solver,
)

TEST_DATA = str(Path(__file__).parent) + "/data/"

set_default_solver("scip")

ess_genes_ana = {
    "G_b0720",
    "G_b2779",
    "G_b2416",
    "G_b2415",
    "G_b1779",
    "G_b1136",
    "G_b4025",
    "G_b2926",
    "G_b3956",
    "G_b3919",
}

ess_rxns_ana = {
    "R_PFK",
    "R_PGI",
    "R_PGK",
    "R_PGM",
    "R_PIt2r",
    "R_PPC",
    "R_ACONTa",
    "R_ACONTb",
    "R_BIOMASS_Ecoli_core_w_GAM",
    "R_CS",
    "R_RPI",
    "R_ENO",
    "R_TPI",
    "R_EX_glc__D_e",
    "R_EX_h_e",
    "R_EX_nh4_e",
    "R_EX_pi_e",
    "R_FBA",
    "R_GAPD",
    "R_GLCpts",
    "R_GLNS",
    "R_ICDHyr",
    "R_NH4t",
}


def test_gene_essentiality():

    model = load_model(TEST_DATA + "e_coli_core.xml")
    anaerobic = {"R_EX_o2_e": 0}
    genes = essential_genes(model, constraints=anaerobic)
    assert set(genes) == ess_genes_ana


def test_reaction_essentiality():

    model = load_model(TEST_DATA + "e_coli_core.xml")
    anaerobic = {"R_EX_o2_e": 0}
    reactions = essential_reactions(model, constraints=anaerobic)
    assert set(reactions) == ess_rxns_ana
