module NanchenData

    import ..Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006
    using Statistics

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_sub_proj(@__MODULE__)
    
    include("data.jl")

    function __init__()
        _populate_bundle()
        UJL.create_proj_dirs(@__MODULE__)
    end

end