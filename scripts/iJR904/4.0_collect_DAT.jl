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

    using Base.Threads
end

# ----------------------------------------------------------------------------------
DAT = ChU.DictTree();

# ----------------------------------------------------------------------------------
FLX_IDERS = ["GLC", "SUCC", "AC", "FORM"]
DAT[:FLX_IDERS] = FLX_IDERS;

# -------------------------------------------------------------------
# ME methods
const ME_Z_OPEN_G_OPEN        = :ME_Z_OPEN_G_OPEN
const ME_MAX_POL              = :ME_MAX_POL
const ME_MAX_POL_B0           = :ME_MAX_POL_B0
const ME_Z_EXPECTED_G_MOVING  = :ME_Z_EXPECTED_G_MOVING
const ME_Z_EXPECTED_G_BOUNDED = :ME_Z_EXPECTED_G_BOUNDED
const ME_Z_FIXXED_G_BOUNDED   = :ME_Z_FIXXED_G_BOUNDED

# LP methods
const FBA_Z_FIX_MIN_COST      = :FBA_Z_FIX_MIN_COST
const FBA_MAX_BIOM_MIN_COST   = :FBA_MAX_BIOM_MIN_COST
const FBA_Z_FIX_MIN_VG_COST   = :FBA_Z_FIX_MIN_VG_COST
const FBA_Z_VG_FIX_MIN_COST   = :FBA_Z_VG_FIX_MIN_COST

LP_METHODS = [
    FBA_Z_FIX_MIN_COST, FBA_MAX_BIOM_MIN_COST, 
    FBA_Z_FIX_MIN_VG_COST, FBA_Z_VG_FIX_MIN_COST
]
DAT[:LP_METHODS] = LP_METHODS

ME_METHODS = [
    # ME_Z_OPEN_G_OPEN, 
    ME_MAX_POL,
    ME_MAX_POL_B0,
    # ME_Z_FIXXED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_MOVING
]
DAT[:ME_METHODS] = ME_METHODS

ALL_METHODS = [LP_METHODS; ME_METHODS]
DAT[:ALL_METHODS] = ALL_METHODS

EXPS = 1:4
DAT[:EXPS] = EXPS;

# -------------------------------------------------------------------
ME_INDEX_FILE = iJR.procdir("maxent_ep_index.bson")
ME_INDEX = ChU.load_data(ME_INDEX_FILE; verbose = false);

# -------------------------------------------------------------------
LP_DAT_FILE = iJR.procdir("lp_dat_file.bson")
LP_DAT = ChU.load_data(LP_DAT_FILE; verbose = false);

# ----------------------------------------------------------------------------------
Fd_mets_map = iJR.load_mets_map()
Fd_rxns_map = iJR.load_rxns_map();

# ----------------------------------------------------------------------------------
const DAT_FILE_PREFFIX =  "maxent_ep_dat"
function dat_file(;kwargs...)
    fname = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; kwargs...)
    iJR.procdir(fname)
end

## ----------------------------------------------------------------------------------
# COMMON DAT
let
    max_model = iJR.load_model("max_model"; uncompress = false)
    objider = iJR.BIOMASS_IDER

    for exp in EXPS
        # exp dat
        Fd_biom = Nd.val("D", exp)
        DAT[:exp, :flx, "D", exp] = Fd_biom
        DAT[:exp, :err, "D", exp] = 0.0

        # bounds
        fva_model = iJR.load_model("fva_models", exp; uncompress = false)
        max_lb, max_ub = ChU.bounds(max_model, objider)
        fva_lb, fva_ub = ChU.bounds(fva_model, objider)
        lb = max(max_lb, fva_lb)
        ub = min(max_ub, fva_ub)
        DAT[:bounds, "D", exp] = (lb, ub)

        for Fd_ider in FLX_IDERS
            model_exch = Fd_rxns_map[Fd_ider]

            # exp dat
            Fd_flx = Nd.uval(Fd_ider, exp)
            DAT[:exp, :flx, Fd_ider, exp] = Fd_flx
            DAT[:exp, :err, Fd_ider, exp] = 0.0

            # bounds
            max_lb, max_ub = ChU.bounds(max_model, model_exch)
            fva_lb, fva_ub = ChU.bounds(fva_model, model_exch)
            lb = max(max_lb, fva_lb)
            ub = min(max_ub, fva_ub)
            DAT[:bounds, Fd_ider, exp] = (lb, ub)

        end
    end
end

# ## ----------------------------------------------------------------------------------
# # Dev
# let
    
#     # In addition to a better quantitative match overall, the
#     # maximum entropy model correctly predicted non-zero flux
#     # through the glyoxylate shunt, i.e., for the isocitrate lyase (ICL)
#     # as well as ME1 reactions, which FBA misses qualitatively by
#     # setting them to zero. As a consequence, this also leads to a better
#     # match of our model with data for reactions isocitrate dehydrogenase 
#     # (ICDH) and alpha ketoglutarate dehydrogenase
#     # (AKGDH) that channel pyruvate through the Krebs cycle.
    
#     # ICL
#     # ME1

#     exp = 4
#     method = ME_MAX_POL
#     # method = ME_MAX_POL_B0
    
#     datfile = dat_file(;method, exp)
#     # datfile = dat_file(;method)
#     dat = deserialize(datfile)
            
#     model = dat[:model]
#     epouts = dat[:epouts]
#     exp_beta = maximum(keys(epouts)) # dat[:exp_beta]
#     epout = epouts[exp_beta]

#     exglc = ChU.av(model, epout, iJR.EX_GLC_IDER)
#     for (ider, model_iders) in iJR.kreps_idermap
        
#         flx = ChU.av(model, epout, model_iders[1])/ abs(exglc)
#         if length(model_iders) == 2 # reversible
#             flx -= ChU.av(model, epout, model_iders[2])/ abs(exglc)
#         end
#         @show ider flx
#         println()
#     end
# end

## ----------------------------------------------------------------------------------
# MAXENT DAT
let 
    WLOCK = ReentrantLock()
    objider = iJR.BIOMASS_IDER

    # util fun
    isexpdep(method) = (method != ME_MAX_POL_B0)
    depks(method, typ, Kd_met, exp) = 
        isexpdep(method) ? (method, typ, Kd_met, exp) : (method, typ, Kd_met)
    function dat_file(;method, exp)
        kwargs = isexpdep(method) ? (;method, exp) : (;method)
        fname = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; kwargs...)
        iJR.procdir(fname)
    end

    # Feed jobs
    nths = nthreads()
    Ch = Channel(nths) do ch
        for method in DAT[:ME_METHODS]
            for exp in DAT[:EXPS]
                put!(ch, (exp, method))
                !isexpdep(method) && break # Do only once
            end
        end
    end

    @threads for thid in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
        for (exp, method) in Ch
            
            # ME data
            datfile = dat_file(;method, exp)
            !isfile(datfile) && continue
            dat = deserialize(datfile)
            
            model = dat[:model]
            epouts = dat[:epouts]
            exp_beta = maximum(keys(epouts)) # dat[:exp_beta]
            epout = epouts[exp_beta]
            
            lock(WLOCK) do
                @info("Doing", 
                    exp, method, 
                    length(dat[:epouts]), 
                    epout.iter, thid
                ); println()
            end

            # Biomass
            ep_biom = ChU.av(model, epout, objider)
            ep_std = sqrt(ChU.va(model, epout, objider))
            
            # store
            lock(WLOCK) do
                DAT[depks(method, :flx, "D", exp)...] = ep_biom
                DAT[depks(method, :err, "D", exp)...] = ep_std
            end

            for Fd_met in FLX_IDERS
                model_exch = Fd_rxns_map[Fd_met]

                # flxs
                ep_av = ChU.av(model, epout, model_exch)
                ep_std = sqrt(ChU.va(model, epout, model_exch))

                # proj 2d
                proj = ChLP.projection2D(model, objider, model_exch; l = 50)
                        
                lock(WLOCK) do
                    DAT[depks(method, :proj, Fd_met, exp)...] = proj
                    DAT[depks(method, :flx, Fd_met, exp)...] = ep_av
                    DAT[depks(method, :err, Fd_met, exp)...] = ep_std
                end
            end

            # additional fluxs
            for (ider, model_iders) in iJR.kreps_idermap
                # flxs
                ep_av = ChU.av(model, epout, model_iders[1])
                ep_std = sqrt(ChU.va(model, epout, model_iders[1]))
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    ep_av -= ChU.av(model, epout, model_iders[2])
                    ep_std += sqrt(ChU.va(model, epout, model_iders[2]))
                end

                # proj 2d (fwd only)
                proj = ChLP.projection2D(model, objider, model_iders[1]; l = 50)
                        
                lock(WLOCK) do
                    DAT[depks(method, :proj, ider, exp)...] = proj
                    DAT[depks(method, :flx, ider, exp)...] = ep_av
                    DAT[depks(method, :err, ider, exp)...] = ep_std
                end
            end

        end # for (exp, method)
    end # for thid
end

## ----------------------------------------------------------------------------------
# LP DAT
let
    objider = iJR.BIOMASS_IDER

    for method in LP_METHODS
            
        for exp in EXPS

            model = LP_DAT[method, :model, exp]
            fbaout = LP_DAT[method, :fbaout, exp]

            # Biomass
            fba_flx = ChU.av(model, fbaout, objider)
            DAT[method, :flx, "D", exp] = fba_flx

            for Fd_ider in FLX_IDERS
                model_ider = Fd_rxns_map[Fd_ider]

                fba_flx = ChU.av(model, fbaout, model_ider)
                DAT[method, :flx, Fd_ider, exp] = fba_flx
            end

            # additional fluxs
            for (ider, model_iders) in iJR.kreps_idermap
                # flxs
                fba_flx = ChU.av(model, fbaout, model_iders[1])
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    fba_flx -= ChU.av(model, fbaout, model_iders[2])
                end
                        
                DAT[method, :flx, ider, exp] = fba_flx
            end
        end

    end # for method
end

## ----------------------------------------------------------------------------------
DAT_FILE = iJR.procdir("dat.bson")
UJL.save_data(DAT_FILE, DAT)