function test_fba(model, obj_ider, cost_ider = nothing; summary = false)
    hascost = cost_ider in model.rxns
    fbaout = !hascost ? ChLP.fba(model, obj_ider) : ChLP.fba(model, obj_ider, cost_ider)
    ChU.tagprintln_inmw("FBA SOLUTION", 
        "\nobj_ider:         ", obj_ider,
        "\nfba obj_val:      ", ChU.av(model, fbaout, obj_ider),
        "\nmax exp obj_val:  ", maximum(NanchenData.val("D")),
        "\ncost_ider:        ", !hascost ? "Not included" : cost_ider,
        "\nfba cost_val:     ", !hascost ? "Not included" : ChU.av(model, fbaout, cost_ider),
    )
    summary && ChU.summary(model, fbaout)
    nothing
end

function test_fba(exp::Int, model, obj_ider, cost_ider = nothing; summary = false)
    hascost = cost_ider in model.rxns
    fbaout = !hascost ? ChLP.fba(model, obj_ider) : ChLP.fba(model, obj_ider, cost_ider)
    ChU.tagprintln_inmw("FBA SOLUTION", 
        "\nobj_ider:         ", obj_ider,
        "\nfba obj_val:      ", ChU.av(model, fbaout, obj_ider),
        "\nexp obj_val:      ", NanchenData.val("D", exp),
        "\ncost_ider:        ", !hascost ? "Not included" : cost_ider,
        "\nfba cost_val:     ", !hascost ? "Not included" : ChU.av(model, fbaout, cost_ider),
    )
    summary && ChU.summary(model, fbaout)
    nothing
end