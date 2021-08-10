using ProjAssistant
@quickactivate

@time begin
    using MAT

    import Chemostat
    const ChLP = Chemostat.LP
    const ChSS = Chemostat.SteadyState
    const ChSU = Chemostat.SimulationUtils
    const ChU = Chemostat.Utils

    import Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006

    const iJR = ChN.iJR904
    const Nd = ChN.NanchenData # experimental data

    using ProgressMeter
    using Plots
    import SparseArrays
end

## ----------------------------------------------------------------------------
BASE_MODELS = ldat(iJR, "base_models.bson");

## ----------------------------------------------------------------------------
# Biomass medium sensibility
let
    model = BASE_MODELS["base_model"]
    obj_ider = iJR.BIOMASS_IDER
    xi = Nd.val("xi") |> minimum
    intake_info = iJR.load_base_intake_info()
    results = Dict()
    factors = 0.0:0.01:1.0
    prog = Progress(length(factors) * length(intake_info))
    for (exch, dat) in intake_info
        c0 = dat["c"]
        res = get!(results, exch, [])
        for f in factors
            dat["c"] = f * c0
            ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)
            growth = try 
                    fbaout = ChLP.fba(model, obj_ider, iJR.COST_IDER)
                    ChU.av(model, fbaout, obj_ider)
                catch err; 0.0 end
            push!(res, growth)

            next!(prog; showvalues = [
                    (:exch, exch),
                    (:f, f),
                    (:c0, c0),
                    (:c, dat["c"]),
                    (:growth, growth)
                ]
            )
        end 
        dat["c"] = c0 
    end
    finish!(prog)

    p = plot(
            title = "Biomass medium sensivility", 
            xlabel = "fraction of initial conc", 
            ylabel = "growth"
        )
    sresults = sort(collect(results); by = (x) -> sum(x[2]))
    lcount = 4
    for (exch, res) in sresults
        lb_ = lcount > 0 ? exch : ""
        plot!(p, factors, res; label = lb_, lw = 3)
        lcount -= 1
    end
    sfig(iJR, p, @fileid, "medium_sesitivity_study", ".png")
end

## ----------------------------------------------------------------------------
# Checking fba_obj_val < exp_obj_val
let 
    to_map = Nd.val("D") |> enumerate
    for (exp, D) in to_map

        model = BASE_MODELS["fva_models"][exp]
        # model = BASE_MODELS["max_model"]

        fbaout = ChLP.fba(model, iJR.BIOMASS_IDER, iJR.COST_IDER);
        fba_obj_val = ChU.av(model, fbaout, iJR.BIOMASS_IDER)
        fba_obj_val = ChU.av(model, fbaout, iJR.BIOMASS_IDER)
        fba_ex_glc_val = ChU.av(model, fbaout, iJR.EX_GLC_IDER)
        fba_ex_glc_b = ChU.bounds(model, iJR.EX_GLC_IDER)
        exp_obj_val = Nd.val("D", exp)
        fba_cost_val = ChU.av(model, fbaout, iJR.COST_IDER)

        @info("FBA SOLUTION", 
            iJR.BIOMASS_IDER, fba_ex_glc_val, fba_ex_glc_b, 
            fba_obj_val, exp_obj_val, 
            iJR.COST_IDER, fba_cost_val
        )

        (fba_obj_val < exp_obj_val) && @warn("fba objval < exp objval", fba_obj_val, exp_obj_val)
    end
end

## ----------------------------------------------------------------------------
# Find cGLC that fit experimental growth
let
    to_map = Nd.val("D") |> enumerate
    for (Di, D) in to_map
        Kd_cGLC = Nd.val(:cGLC, Di)
        Kd_growth = Nd.val(:D, Di)
        xi = Nd.val(:xi, Di)
        model = ChU.load_data(iJR.BASE_MODEL_FILE; verbose = false)
        intake_info = iJR.load_base_intake_info()
        for (exch, info) in intake_info
            info["c"] = iJR.MAX_CONC # open medium
        end
        
        function work_fun(cGLC)
            ## Open intakes except Glucose
            intake_info[iJR.EX_GLC_IDER]["c"] = first(cGLC)
            
            # impose constraint
            ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)

            ## fba
            fbaout = ChLP.fba(model, iJR.BIOMASS_IDER, iJR.COST_IDER)
            fba_growth = ChU.av(model, fbaout, iJR.BIOMASS_IDER)
            return [fba_growth]
        end

        cGLC = ChSU.grad_desc(work_fun; x0 = [Kd_cGLC], x1 = [Kd_cGLC * 0.9], th = 1e-5, 
            C = [Kd_cGLC * 0.1], target = [Kd_growth], maxiters = 500) |> first
        @info "Results" Di Kd_cGLC cGLC
    end

end

## ----------------------------------------------------------------------------
# Testing scaled model
let
    model = ChU.load_data(iJR.BASE_MODEL_FILE; verbose = false)
    fbaout = ChLP.fba(model, iJR.BIOMASS_IDER, iJR.COST_IDER);
    println("FBA SOLUTION", 
        "\nobj_ider:         ", iJR.BIOMASS_IDER,
        "\nsize:             ", size(model),
        "\nfba obj_val:      ", ChU.av(model, fbaout, iJR.BIOMASS_IDER),
        "\ncost_ider:        ", iJR.COST_IDER,
        "\nfba cost_val:     ", ChU.av(model, fbaout, iJR.COST_IDER),
        "\n\n"
    )
    model = ChU.well_scaled_model(model, 100.0; verbose = false)
    fbaout = ChLP.fba(model, iJR.BIOMASS_IDER, iJR.COST_IDER);
    println("FBA SOLUTION", 
        "\nobj_ider:         ", iJR.BIOMASS_IDER,
        "\nsize:             ", size(model),
        "\nfba obj_val:      ", ChU.av(model, fbaout, iJR.BIOMASS_IDER),
        "\ncost_ider:        ", iJR.COST_IDER,
        "\nfba cost_val:     ", ChU.av(model, fbaout, iJR.COST_IDER),
        "\n\n"
    )
end