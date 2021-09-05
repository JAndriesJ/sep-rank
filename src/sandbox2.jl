    Pkg.activate(".")
    srcDir  = dirname(@__FILE__)*"\\";
    include(srcDir*"Examples\\Examples.jl")
    include(srcDir*"Examples\\Utils_states.jl")
    include(srcDir*"Moments.jl")
    include(srcDir*"Constraints\\Utils_cons.jl")
    include(srcDir*"Constraints\\R_constraints.jl")
    include(srcDir*"Constraints\\C_constraints.jl")
    include(srcDir*"Model\\Utils_Model.jl")
    include(srcDir*"Model\\R_sep_Model.jl")
    include(srcDir*"Model\\C_sep_Model.jl")
    include(srcDir*"sep_Compute.jl")

    using .Examples      ; const ex = Examples
    using .Utils_states  ; const us = Utils_states
    using .Moments       ; const mom = Moments
    using .Utils_cons    ; const uc = Utils_cons
    using .C_constraints ; const ccon = C_constraints
    using .R_constraints ; const rcon = C_constraints
    using .Utils_Model
    using .R_sep_Model
    using .C_sep_Model
    using .sep_Compute

    using LinearAlgebra ; const la = LinearAlgebra # Pkg.add("LinearAlgebra")
    using JuMP

## Try d = (3,3), t = 3 without block-diag
ρ = ex.get_example("RCai")
d = ex.get_example_meta("RCai").Size[1]
t = 3
ℂ1_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S1G",noBlock=true)

## Try d = (2,2), t = 4 with block-diag
d = (7,7) ; r = 5 ; t = 2
ρ = Examples.gen_ρ_RAND(d, r, true)
ℂ1_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S1G",noBlock=false)

## Get variable reduction
d = (3,3) ; t = 3
## Block
MMᴿcoef,MMᴿexᴿ = mom.get_ℂ_block_diag(d,t,noBlock=false)
for b in keys(MMᴿcoef)
    unique([MMᴿexᴿ[(-2, -1)]...])
end

B = [size(unique([MMᴿexᴿ[b]...]))[1] for b in keys(MMᴿcoef)]
sum(B)

## No block
# MMᴿcoef,MMᴿexᴿ = mom.get_ℂ_block_diag(d, 2*t,noBlock=true)  # Take too long
unique(MMᴿexᴿ["Default"])
