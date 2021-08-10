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
    import MAT

    using ProjAssistant
    @gen_sub_proj

    include("const.jl")
    include("load_data.jl")
    include("beg_enz_cost.jl")
    include("load_model.jl")
    include("ME_MODES.jl")
    
    function __init__()
        @create_proj_dirs
    end

end