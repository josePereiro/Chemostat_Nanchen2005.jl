module iJR904

    import BSON
    import ..Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006
    import ..BegData
    const Bd = BegData
    import ..NanchenData
    const Nd = NanchenData
    import Chemostat
    const ChU = Chemostat.Utils

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_sub_proj(@__MODULE__)

    include("const.jl")
    include("load_data.jl")
    include("beg_enz_cost.jl")
    include("load_model.jl")
    include("iJR_inner_rxn_map.jl")
    
    function __init__()
        UJL.create_proj_dirs(@__MODULE__)
    end

end