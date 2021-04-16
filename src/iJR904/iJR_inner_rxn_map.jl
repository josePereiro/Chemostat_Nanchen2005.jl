# ## ------------------------------------------------------------------
# # Dev
# # LOAD RAW MODEL
# src_file = "/Users/Pereiro/University/Research/Metabolism/MaxEntEP2020/WIP/Chemostat_Nanchen2005/data/raw/iJR904/iJR904.mat"
# mat_model = MAT.matread(src_file)["model"]
# model = ChU.MetNet(mat_model; reshape=true)
# println(size(model))

# ## ------------------------------------------------------------------
# ChU.summary(model, "glx_c")

# ## ------------------------------------------------------------------
# ChU.search(model, "biomass")

# ## ------------------------------------------------------------------
# let
#     mets1 = [ "nadh_c"]
#     mets2 = ["nad_c"]
#     isempty(mets2) && (mets2 = mets1)

#     rxnis = []
#     for (met1, met2) in Iterators.product(mets1, mets2)
#         rxnis1 = ChU.met_rxns(model, met1)
#         rxnis2 = ChU.met_rxns(model, met2)
#         push!(rxnis, intersect(rxnis1, rxnis2)...)
#     end
#     rxnis = unique!(rxnis)

#     println("RXNS ------------------------------------------------- ")
#     isempty(rxnis) && @warn "rxnis is empty"
#     for rxni in rxnis
#         rxn = ChU.rxns(model, rxni)
#         ChU.summary(model, rxn)
#         println()
#     end
# end

## ------------------------------------------------------------------
TABLE2_INNER_RXN_MAP = Dict(
    # GLC => "glc_DASH_D_c"
    # ATP => "atp_c"
    # G6P => "g6p_c"
    # T3P => "g3p_c"
    # 6PG => "2dh3dgal6p_c", "2ddg6p_c", "6pgl_c", "6pgc_c" # Not sure some of them
    # P5P => "ru5p_DASH_D_c", "xu5p_DASH_D_c"
    # F6P => "f6p_c"
    # 2T3P => "fdp_c" # Not sure!!!
    # S7P => "s7p_c"
    # E4P => "e4p_c"
    # 2P5P => "xu5p_DASH_D_c", "r5p_c" # Not sure (maybe 2 is a stoi coef)
    # NADH => "nadh_c"
    # PGA => "13dpg_c" # no sure
    # PEP => "pep_c"
    # PYR => "pyr_c"
    # AcCoA => "accoa_c"
    # ICT => "cit_c", "icit_c"
    # OGA => "oaa_c", "akg_c"
    # CO2 => "co2_c"
    # OGA => "fum_c"
    # MAL => "mal_DASH_L_c"


    # DDPGALA: (-1.0) 2dh3dgal6p_c <==> (1.0) g3p_c + (1.0) pyr_c
    # EDA: (-1.0) 2ddg6p_c ==> (1.0) g3p_c + (1.0) pyr_c
    "6PG -> T3P + PYR" => ["DDPGALA", "EDA"], 

    # HEX1: (-1.0) atp_c + (-1.0) glc_DASH_D_c ==> (1.0) adp_c + (1.0) g6p_c + (1.0) h_c
    "GLC + ATP -> G6P" => ["HEX1"], # glucose-6-phosphate isomerase
    
    # G6PDH2r: "(-1.0) g6p_c + (-1.0) nadp_c <==> (1.0) 6pgl_c + (1.0) h_c + (1.0) nadph_c"
    "G6P -> 6PG + NADPH" => ["G6PDH2r"], 

    # (-1.0) 6pgc_c + (-1.0) nadp_c ==> (1.0) co2_c + (1.0) nadph_c + (1.0) ru5p_DASH_D_c
    "6PG -> P5P + CO2 + NADPH" => ["GND"], 

    # PGI: (-1.0) g6p_c <==> (1.0) f6p_c
    "G6P -> F6P" => ["PGI"],

    # Not sure!!!
    # PFK: (-1.0) atp_c + (-1.0) f6p_c ==> (1.0) adp_c + (1.0) fdp_c + (1.0) h_c
    "F6P + ATP -> 2T3P" => ["PFK"], 

    # Not sure!!!
    # TKT1: (-1.0) r5p_c + (-1.0) xu5p_DASH_D_c <==> (1.0) g3p_c + (1.0) s7p_c
    "2P5P -> S7P + T3P" => ["TKT1"], 

    # TKT2: (-1.0) e4p_c + (-1.0) xu5p_DASH_D_c <==> (1.0) f6p_c + (1.0) g3p_c
    "P5P + E4P -> F6P + T3P" => ["TKT2"],

    # (-1.0) g3p_c + (-1.0) s7p_c <==> (1.0) e4p_c + (1.0) f6p_c
    "S7P + T3P -> E4P + F6P" => ["TALA"], 

    # Not sure!!!
    # GAPD: (-1.0) g3p_c + (-1.0) nad_c + (-1.0) pi_c <==> (1.0) 13dpg_c + (1.0) h_c + (1.0) nadh_c
    "T3P -> PGA + ATP + NADH" => ["GAPD"],

    # Not Found
    "PGA -> PEP" => [""],

    # PYK: (-1.0) adp_c + (-1.0) h_c + (-1.0) pep_c ==> (1.0) atp_c + (1.0) pyr_c
    "PEP -> PYR + ATP" => ["PYK"],
    
    # PYR: (-1.0) coa_c + (-1.0) nad_c + (-1.0) pyr_c ==> (1.0) accoa_c + (1.0) co2_c + (1.0) nadh_c
    "PYR -> AcCoA + CO2 + NADH" => ["PYR"],
    
    # CS: (-1.0) accoa_c + (-1.0) h2o_c + (-1.0) oaa_c ==> (1.0) cit_c + (1.0) coa_c + (1.0) h_c
    "OAA + AcCoA -> ICT" => ["CS"],

    # ICDHyr: (-1.0) icit_c + (-1.0) nadp_c <==> (1.0) akg_c + (1.0) co2_c + (1.0) nadph_c
    "ICT -> OGA + CO2 + NADPH" => ["ICDHyr"],

    # Not found!!!
    "OGA -> FUM + CO2 + 1.5*ATP + 2NADH" => [""],

    # FUM: (-1.0) fum_c + (-1.0) h2o_c <==> (1.0) mal_DASH_L_c
    "FUM -> MAL" => ["FUM"],

    # MDH: (-1.0) mal_DASH_L_c + (-1.0) nad_c <==> (1.0) h_c + (1.0) nadh_c + (1.0) oaa_c
    "MAL -> OAA + NADH" => ["MDH"],

    # ME1: (-1.0) mal_DASH_L_c + (-1.0) nad_c ==> (1.0) co2_c + (1.0) nadh_c + (1.0) pyr_c
    "MAL -> PYR + CO2 + NADH" => ["ME1"],

    # PPCK: (-1.0) atp_c + (-1.0) oaa_c ==> (1.0) adp_c + (1.0) co2_c + (1.0) pep_c
    "OAA + ATP -> PEP + CO2" => ["PPCK"],

    # PPC: (-1.0) co2_c + (-1.0) h2o_c + (-1.0) pep_c ==> (1.0) h_c + (1.0) oaa_c + (1.0) pi_c
    "PEP + CO2 -> OAA" => ["PPC"],

    # ACS: (-1.0) ac_c + (-1.0) atp_c + (-1.0) coa_c ==> (1.0) accoa_c + (1.0) amp_c + (1.0) ppi_c
    "AcCoA -> Acetate + ATP" => ["ACS"],

    # MALS: (-1.0) accoa_c + (-1.0) glx_c + (-1.0) h2o_c ==> (1.0) coa_c + (1.0) h_c + (1.0) mal_DASH_L_c
    "ICT + AcCoA -> MAL + FUM + NADH" => ["MALS"],

    #  Not found
    "NADPH -> NADH" => [""],
    
    # Not found
    "O2 + 2NADH -> 2P/O x ATP" => [""],

    # BiomassEcoli
    "Biomass" => ["BiomassEcoli"],
)