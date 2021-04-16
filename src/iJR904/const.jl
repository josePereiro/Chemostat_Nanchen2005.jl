const PROJ_IDER = "iJR904"
const BIOMASS_IDER = "BiomassEcoli"
const EX_GLC_IDER = "EX_glc_LPAREN_e_RPAREN_"
const ATPM_IDER = "ATPM"
const COST_IDER = "tot_cost"
const ABS_MAX_BOUND = 100.0
const MAX_CONC = 9999.0

krebs_iders = ["SUCD1", "SUCOAS", "AKGDH", "ICDHyr", 
    "ACONT", "CS", "MDH", "FUM", "MALS", "ICL"
]

kreps_idermap = Dict(
    "SUCD1" => ["SUCD1i"], 
    "SUCOAS" => ["SUCOAS_fwd", "SUCOAS_bkwd"], 
    "AKGDH" => ["AKGDH"],
    "ICDHyr" => ["ICDHyr_fwd", "ICDHyr_bkwd"],
    "ACONT" => ["ACONT_bkwd", "ACONT_fwd"],
    "CS" => ["CS"],
    "MDH" => ["MDH_fwd", "MDH_bkwd"],
    "FUM" => ["FUM_fwd", "FUM_bkwd"],
    "MALS" => ["MALS"],
    "ICL" => ["ICL"]
)