## -------------------------------------------------------------------
function load_model(modelkey::String; uncompress = true)
    if modelkey == "raw_model"
        model = _load_raw_model()
    else 
        models = lprocdat(iJR904, "base_models.bson"; verbose = false);
        model = models[modelkey]
    end
    return uncompress ? ChU.uncompressed_model(model) : model
end

function load_model(modelkey::String, exp::Int; uncompress = true)
    models = lprocdat(iJR904, "base_models.bson"; verbose = false);
    model = models[modelkey][exp]
    return uncompress ? ChU.uncompressed_model(model) : model
end

function _load_raw_model()
    src_file = rawdir(iJR904, "iJR904.mat")
    mat_model = MAT.matread(src_file)["model"]
    return ChU.MetNet(mat_model; reshape = true)
end