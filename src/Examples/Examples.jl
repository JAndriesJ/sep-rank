module Examples

using LinearAlgebra
using Random
using DataFrames
using CSV

srcDir = dirname(@__FILE__)*"\\"

include(srcDir *"Utils_states.jl")
include(srcDir *"Examples_sep.jl")
include(srcDir *"Examples_ent.jl")

# using .Utils
using .Utils_states, .Examples_sep, .Examples_ent

export get_example,
       get_example_overview,
       gen_ρ_RAND

gen_ρ_RAND(d,r::Integer,isℂ= false) = Examples_sep.gen_ρ_RAND(d, r, isℂ )

"""Returns dict with three keys: ["rand", "sep", "ent"] each entry is also a dict."""
get_example() = merge(Examples_sep.get_sep_example(),Examples_ent.get_ent_example())
get_example(ex::String) = get_example()[ex]
get_example(ex::Int) = get_example()[[keys(get_example())...][ex]]

siz_dic = Dict(  4 => (2,2),
                 6 => (2,3),
                 8 => (2,4),
                 9 => (3,3),
                 10 => (2,5),
                 12 => (3,4),
                 16 => (4,4))

"""Returns a DataFrame object containing the"""
function get_example_meta(isWrite = false)
    df  = DataFrame( Example = String[], Size = Any[], isSeparable = String[], isReal = Bool[], Bi_rank = Any[], hasPPᵀ = Bool[]);
    for sep in ["sep","ent"]
        ρ_dict = (sep =="sep") ? Examples_sep.get_sep_example() : Examples_ent.get_ent_example()
        for k in keys(ρ_dict)
            ρ      =  ρ_dict[k]
            ex     = k
            siz    = siz_dic[size(ρ)[1]]
            isr    = isreal(ρ)
            rank1  = rank(ρ)
            rank2  = rank(take_Pᵀ(ρ, 2, siz))
            isPPt  = Utils_states.isPPᵀ(ρ,siz)

            push!(df, [ex,siz,sep ,isr,(rank1,rank2),isPPt] )
        end
    end
    sort!(df,[:Example,:isSeparable])
    isWrite ? (CSV.write(dirname(dirname(@__FILE__))*"\\ExSummary.csv", df , delim = "&")) : (return df)
end
#
get_example_meta(ex::Int) = get_example_meta()[ex,:]
function get_example_meta(ex::String)
    df = get_example_meta()
    return df[df.Example .== ex,:]
end



end
