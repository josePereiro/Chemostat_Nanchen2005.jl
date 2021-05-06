import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Nanchen2006")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    #  ----------------------------------------------------------------------------
    import Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006

    const iJR = ChN.iJR904
    const Nd = ChN.NanchenData # experimental data
    const Bd = ChN.BegData    # cost data

    #  ----------------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChLP = Ch.LP

    using Serialization
    import UtilsJL
    const UJL = UtilsJL
end

## -----------------------------------------------------------------------------------------------
const FBA_Z_FIX_MIN_COST    = :FBA_Z_FIX_MIN_COST
const FBA_MAX_Z_MIN_COST    = :FBA_MAX_Z_MIN_COST
const FBA_Z_FIX_MAX_VG_MIN_COST = :FBA_Z_FIX_MAX_VG_MIN_COST
const FBA_Z_FIX_MIN_VG_COST = :FBA_Z_FIX_MIN_VG_COST
const FBA_Z_VG_FIX_MIN_COST = :FBA_Z_VG_FIX_MIN_COST
const EXPS = Nd.EXPS

## -----------------------------------------------------------------------------------------------
# Data container
LPDAT = UJL.DictTree()

## -------------------------------------------------------------------
# FBA
let
    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    max_sense = -1.0
    min_sense = 1.0

    iterator = Nd.val(:D) |> enumerate |> collect 
    for (exp, D) in iterator

        @info("Doing ", exp); println()
        model0 = iJR.load_model("fva_models", exp)

        # FBA_MAX_Z_MIN_COST
        let
            model = deepcopy(model0)
            fbaout = ChLP.fba(model, objider, costider)
            
            LPDAT[FBA_MAX_Z_MIN_COST, :model, exp] = model
            LPDAT[FBA_MAX_Z_MIN_COST, :fbaout, exp] = fbaout
        end
        
        # FBA_Z_FIX_MIN_COST
        let
            model = deepcopy(model0)
            exp_growth = Nd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout = ChLP.fba(model, objider, costider)

            LPDAT[FBA_Z_FIX_MIN_COST, :model, exp] = model
            LPDAT[FBA_Z_FIX_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MIN_VG_COST
        let
            model = deepcopy(model0)
            exp_growth = Nd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout1 = ChLP.fba(model, exglcider; sense = min_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            ChU.bounds!(model, exglcider, exglc, exglc)
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[FBA_Z_FIX_MIN_VG_COST, :model, exp] = model
            LPDAT[FBA_Z_FIX_MIN_VG_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MAX_VG_MIN_COST
        let
            model = deepcopy(model0)
            exp_growth = Nd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout1 = ChLP.fba(model, exglcider; sense = max_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            ChU.bounds!(model, exglcider, exglc, exglc)
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[FBA_Z_FIX_MAX_VG_MIN_COST, :model, exp] = model
            LPDAT[FBA_Z_FIX_MAX_VG_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_VG_FIX_MIN_COST
        let
            model = deepcopy(model0)
        
            exp_growth = Nd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            
            exp_exglc = Nd.uval("GLC", exp)
            # exp_exglc must be negative
            exp_exglc = exp_exglc > 0 ? -exp_exglc : exp_exglc 
            ChU.bounds!(model, exglcider, exp_exglc, exp_exglc)
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[FBA_Z_VG_FIX_MIN_COST, :model, exp] = model
            LPDAT[FBA_Z_VG_FIX_MIN_COST, :fbaout, exp] = fbaout
        end
    end
end

## -------------------------------------------------------------------
LP_DAT_FILE = iJR.procdir("lp_dat_file.bson")
ChU.save_data(LP_DAT_FILE, LPDAT)
