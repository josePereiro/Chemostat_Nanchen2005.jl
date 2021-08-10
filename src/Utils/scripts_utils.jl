function test_fba(model, obj_ider, cost_ider = nothing; summary = false)
    hascost = cost_ider in model.rxns
    fbaout = !hascost ? ChLP.fba(model, obj_ider) : ChLP.fba(model, obj_ider, cost_ider)
    
    # info
    obj_val = ChU.av(model, fbaout, obj_ider)
    exp_obj_val = maximum(NanchenData.val("D"))
    cost_ider = !hascost ? "Not included" : cost_ider
    cost_val = !hascost ? "Not included" : ChU.av(model, fbaout, cost_ider)
    @info("FBA SOLUTION", obj_val, exp_obj_val, cost_ider, cost_val)
    
    summary && ChU.summary(model, fbaout)

    nothing
end

function test_fba(exp::Int, model, obj_ider, cost_ider = nothing; summary = false)
    hascost = cost_ider in model.rxns
    fbaout = !hascost ? ChLP.fba(model, obj_ider) : ChLP.fba(model, obj_ider, cost_ider)

    # info
    obj_val = ChU.av(model, fbaout, obj_ider)
    exp_obj_val = NanchenData.val("D", exp)
    cost_ider = !hascost ? "Not included" : cost_ider
    cost_val = !hascost ? "Not included" : ChU.av(model, fbaout, cost_ider)
    @info("FBA SOLUTION", obj_val, exp_obj_val, cost_ider, cost_val)

    summary && ChU.summary(model, fbaout)
    nothing
end