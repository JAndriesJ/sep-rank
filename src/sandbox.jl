### Begin here
using Pkg
# Pkg.activate(".")
Pkg.status()

# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("JuMP")
# Pkg.add("MosekTools")
# Pkg.add("LinearAlgebra")
# Pkg.add("Test")

    srcDir  = dirname(@__FILE__)*"\\";
    include(srcDir*"Examples\\Examples.jl")
    include(srcDir*"Examples\\Utils_states.jl")
    include(srcDir*"Moments.jl")
    include(srcDir*"Constraints\\C_constraints.jl")
    include(srcDir*"Constraints\\R_constraints.jl")
    include(srcDir*"Constraints\\Utils_cons.jl")
    include(srcDir*"Model\\R_sep_Model.jl")
    include(srcDir*"Model\\C_sep_Model.jl")
    include(srcDir*"Model\\Utils_Model.jl")
    include(srcDir*"sep_Compute.jl")

    using .sep_Compute
    using .Examples
    using .Utils_states
    using .Moments ; const mom = Moments
    using .Utils_cons ; const uc = Utils_cons
    using .C_constraints ; const ccon = C_constraints
    using .R_constraints
    using .Utils_Model
    using .R_sep_Model
    using .C_sep_Model
    using .sep_Compute

    using LinearAlgebra ; const la = LinearAlgebra
    using JuMP

## Goal get the C block diag to get the same values as the non- block diag.

ρ  =  [1. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 1.]
d = (2,2) ; t = (2,2) ; cl = "S∞" #"S∞"

ℂb_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list=cl,noBlock=false)
ℂ_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list=cl,noBlock=true)
ℝ_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list=cl,isRe=true)

## Block diag
C,E = mom.get_ℂ_block_diag(d,t;noBlock = true)
Cb,Eb = Moments.get_ℂ_block_diag(d,t;noBlock = false)
## Moments
MM = Moments.get_xx̄yȳMM_blocks(d,t)


MM[0,-2]


## Utils_cons
