# Here I include some maps between the experimental
# data ids and the model ids

## ------------------------------------------------------------------
# maps between Kayser2005 https://doi.org/10.1099/mic.0.27481-0.
# data and the model
function load_mets_map() 
    mets_map = Dict()
    mets_map["GLC"] = "glc_DASH_D_e"
    mets_map["SUCC"] = "succ_e"
    mets_map["FORM"] = "for_e"
    mets_map["THM"] = "thm_e"
    mets_map["NH4"] = "nh4_e"
    mets_map["CIT"] = "cit_e"
    mets_map["CO2"] = "co2_e"
    mets_map["O2"] = "o2_e"
    mets_map["AC"] = "ac_e"
    mets_map["PYR"] = "pyr_e"
    mets_map["LAC"] = "lac_DASH_D_e"
    mets_map["MAL"] = "mal_DASH_D_e";
    for (k, v) in mets_map
        mets_map[v] = k
    end
    return mets_map
end

# EX_glc_LPAREN_e_RPAREN_
function load_rxns_map()
    rxns_map = Dict()
    rxns_map["D"] = "BiomassEcoli"
    rxns_map["AC"] = "EX_ac_LPAREN_e_RPAREN_"
    rxns_map["NH4"] = "EX_nh4_LPAREN_e_RPAREN_"
    rxns_map["GLC"] = "EX_glc_LPAREN_e_RPAREN_"
    rxns_map["THM"] = "EX_thm_LPAREN_e_RPAREN_"
    rxns_map["CO2"] = "EX_co2_LPAREN_e_RPAREN_"
    rxns_map["FORM"] = "EX_for_LPAREN_e_RPAREN_"
    rxns_map["CIT"] = "EX_cit_LPAREN_e_RPAREN_"
    rxns_map["SUCC"] = "EX_succ_LPAREN_e_RPAREN_"
    rxns_map["O2"] = "EX_o2_LPAREN_e_RPAREN_"
    rxns_map["MAL"] = "EX_mal_L_LPAREN_e_RPAREN_"
    rxns_map["LAC"] = "EX_lac_D_LPAREN_e_RPAREN_"
    rxns_map["PYR"] = "EX_pyr_LPAREN_e_RPAREN_"
    for (k, v) in rxns_map
        rxns_map[v] = k
    end
    return rxns_map
end


## ------------------------------------------------------------------
# enzymatic costs
# from Beg, (2007) https://doi.org/10.1073/pnas.0609845104.

# A map between model ids and the reactions reported in Beg2007
function load_beg_rxns_map()
    beg_rxns_map = Dict(
        "carbonyl reductase (NADPH)" => ["P5CR"],
        "alcohol dehydrogenase (NADP+)" => ["ALCD19"],
        "quinate/shikimate dehydrogenase" => ["SHK3Dr"],
        "malate dehydrogenase (decarboxylating)" => ["MDH","MDH2", "MDH3"],
        "3alpha-hydroxysteroid dehydrogenase (B-specific)" => [""],
        "2-hydroxy-3-oxopropionate reductase" => [""],
        "glucose dehydrogenase (acceptor)" => ["G6PDH2r", "UDPGD"],
        "cellobiose dehydrogenase (acceptor)" => [""],
        "peroxidase" => [""],
        "catechol 2,3-dioxygenase" => [""],
        "arachidonate 8-lipoxygenase" => [""],
        "calcidiol 1-monooxygenase" => [""],
        "nitric-oxide synthase" => [""],
        "phenylalanine 4-monooxygenase" => [""],
        "tryptophan 5-monooxygenase" => [""],
        "Carboxylate reductase" => ["P5CR"],
        "arsenate reductase (donor)" => [""],
        "biliverdin reductase" => [""],
        "15-oxoprostaglandin 13-oxidase" => [""],
        "coproporphyrinogen oxidase" => ["CPPPGO","PPPGO"],
        "long-chain-acyl-CoA dehydrogenase" => [""],
        "butyryl-CoA dehydrogenase" => [""],
        "acyl-CoA dehydrogenase" => [""],
        "L-amino-acid oxidase" => ["ASPO3","ASPO4","ASPO5","ASPO6"],
        "amine oxidase (flavin-containing)" => ["PYAM5PO"],
        "methylenetetrahydrofolate reductase [NAD(P)H]" => ["MTHFR2"],
        "formyltetrahydrofolate dehydrogenase" => ["MTHFD"],
        "sarcosine oxidase" => [""],
        "nitrate reductase (NADH)" => [""],
        "nitrite reductase (NO-forming)" => [""],
        "nitrate reductase" => ["NO3R2"],
        "trypanothione-disulfide reductase" => [""],
        "glutathione-disulfide reductase" => [""],
        "thioredoxin-disulfide reductase" => [""],
        "thiol oxidase" => [""],
        "nitrate reductase (cytochrome)" => [""],
        "aspartate carbamoyltransferase" => ["ASPCT"],
        "serine O-acetyltransferase" => ["SERAT"],
        "protein-glutamine gamma-glutamyltransferase" => [""],
        "gamma-glutamyltransferase" => ["CRNBTCT"],
        "citrate (Si)-synthase" => ["CS"],
        "kaempferol 3-O-galactosyltransferase" => [""],
        "NAD+ ADP-ribosyltransferase" => ["NNDMBRT"],
        "di-trans,poly-cis-decaprenylcistransferase" => ["ACGAMT"],
        "cystathionine gamma-synthase" => [""],
        "adenosine kinase" => ["ADNK1"],
        "glycerate kinase" => ["GLYCK"],
        "galactokinase" => ["GALKr"],
        "[pyruvate dehydrogenase (acetyl-transferring)] kinase" => [""],
        "guanylate kinase" => ["GK1"],
        "FMN adenylyltransferase" => ["FMNAT"],
        "tRNA adenylyltransferase" => [""],
        "aryl sulfotransferase" => [""],
        "aminoacyl-tRNA hydrolase" => [""],
        "carboxymethylenebutenolidase" => [""],
        "ubiquitin thiolesterase" => [""],
        "fructose-bisphosphatase" => ["FBP","FBA"],
        "[phosphorylase] phosphatase" => [""],
        "phosphoglycolate phosphatase" => ["PGLYCP"],
        "protein-tyrosine-phosphatase" => [""],
        "inositol-polyphosphate 5-phosphatase" => ["MI1PP"],
        "3',5'-cyclic-GMP phosphodiesterase" => [""],
        "beta-glucosidase" => ["MLTG1","MLTG2","MLTG3","MLTG4","MLTG5"],
        "beta-glucuronidase" => [""],
        "glucosylceramidase" => [""],
        "cyclomaltodextrinase" => [""],
        "alpha-N-arabinofuranosidase" => [""],
        "purine nucleosidase" => ["AMPN","AHCYSNS","CMPN","MTAN","NMNN"],
        "rRNA N-glycosylase" => [""],
        "NAD+ nucleosidase" => [""],
        "Xaa-Pro aminopeptidase" => [""],
        "dipeptidyl-peptidase I" => [""],
        "peptidyl-dipeptidase A" => [""],
        "coagulation factor Xa" => [""],
        "t-Plasminogen activator" => [""],
        "cathepsin B" => [""],
        "envelysin" => [""],
        "amidase" => ["","GSPMDA","NMNDA","NNAM"],
        "formamidase" => [""],
        "arginase" => [""],
        "guanidinoacetase" => ["GUAD"],
        "apyrase" => [""],
        "phloretin hydrolase" => [""],
        "Orotidine-5'-phosphate decarboxylase" => ["OMPDC"],
        "4-Hydroxybenzoate decarboxylase" => ["OPHBDC"],
        "Threonine aldolase" => ["THRAr"],
        "enoyl-CoA hydratase" => [""],
        "Uroporphyrinogen-III synthase" => ["UPP3S"],
        "dihydroxy-acid dehydratase" => ["DHAD1","DHAD2"],
        "pectin lyase" => [""],
        "DNA-(apurinic or apyrimidinic site) lyase" => [""],
        "lactoylglutathione lyase" => ["LGTHL"],
        "guanylate cyclase" => [""],
        "dTDP-4-dehydrorhamnose 3,5-epimerase" => ["TDPDRE"],
        "UDP-glucose 4-epimerase" => ["UDPG4E"],
        "Triose-phosphate isomerase" => ["TPI"],
        "steroid DELTA-isomerase" => [""],
        "dodecenoyl-CoA isomerase" => [""],
        "Glutamate-1-semialdehyde 2,1-aminomutase" => [""],
        "Chalcone isomerase" => [""],
        "Chloromuconate cycloisomerase" => [""],
        "Tyrosine-tRNA ligase" => [""],
        "Threonine-tRNA ligase" => [""],
        "Isoleucine-tRNA ligase" => [""],
        "Lysine-tRNA ligase" => [""],
        "formate-tetrahydrofolate ligase" => [""],
        "Adenylosuccinate synthase" => ["ADSS"],
        "DNA ligase (NAD+)" => [""]
    )

    # inverse relation
    for (beg_id, model_ids) in beg_rxns_map
        for model_id in model_ids
            beg_ids = get!(beg_rxns_map, model_id, [])
            push!(beg_ids, beg_id)
        end
    end
    return beg_rxns_map
end



## ------------------------------------------------------------------
# base model exch met map
# A quick way to get exchages from mets and te other way around
load_exch_met_map() = UJL.load_data(procdir("exch_met_map.bson"); verbose=false)

## ------------------------------------------------------------------
# The intakes bounds of the network are determined by the 
# medium concentration in the Chemostat model (see Cossios paper)
# This is a base medium for modeling

function load_base_intake_info()
    return Dict(
        "EX_glc_LPAREN_e_RPAREN_" => Dict("c"=> maximum(Nd.val(:cGLC)), "lb"=> -ABS_MAX_BOUND),
        
        "EX_nh4_LPAREN_e_RPAREN_" => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_o2_LPAREN_e_RPAREN_"  => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_pi_LPAREN_e_RPAREN_"  => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_so4_LPAREN_e_RPAREN_" => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
    )
end

function intake_info(exp)
    return Dict(
        "EX_glc_LPAREN_e_RPAREN_" => Dict("c"=> Nd.val(:cGLC, exp), "lb"=> -ABS_MAX_BOUND),
        
        "EX_nh4_LPAREN_e_RPAREN_" => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_o2_LPAREN_e_RPAREN_"  => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_pi_LPAREN_e_RPAREN_"  => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_so4_LPAREN_e_RPAREN_" => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
    )
end