module Chemostat_Nanchen2006

    import BSON
    import DrWatson
    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP

    using ProjAssistant
    @gen_top_proj

    include("NanchenData/NanchenData.jl")
    include("Utils/Utils.jl")
    include("BegData/BegData.jl")
    include("iJR904/iJR904.jl")

    function __init__()
        @create_proj_dirs
    end

end
