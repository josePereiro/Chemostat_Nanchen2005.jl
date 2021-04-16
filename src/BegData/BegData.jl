# Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
module BegData

    import ..Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006
    using DataFrames
    import CSV

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_sub_proj(@__MODULE__)

    ## ------------------------------------------------------------------------
    # The average crowding coefficient a was fitted to obtain the minimum square 
    # deviation between the measured and model predicted growth rates, 
    # resulting in a = 0.0040 ± 0.0005 h•g/mmol, in which g is grams dry weight. 
    # However, the maximum growth rates on glucose and glycerol are more consistent 
    # with a = 0.0031 ± 0.0001 h•g/mmol and a = 0.0053 ± 0.0001 h•g/mmol, respectively.
    const AVE_COST  = 0.0031 # I will use this to set the unknown costs

    ## ------------------------------------------------------------------------
    # enzymatic costs from Beg et al. (2007): 
    # https://doi.org/10.1073/pnas.0609845104.
    load_enz_data() = CSV.read(rawdir("beg2007___enzymatic_data.tsv"), DataFrame)
    
    function __init__()
        UJL.create_proj_dirs(@__MODULE__)
    end
end