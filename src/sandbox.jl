### Begin here
# using Pkg

# Pkg.status()

# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("MosekTools")
# Pkg.add("Test")
# Pkg.add("JuMP")
# Pkg.add("LinearAlgebra")

    Pkg.activate(".")
    #Pkg.instantiate()
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

## Goal get the C block diag to get the same values as the non- block diag.
ρ = [1. 0. 0. 0.
     0. 0. 0. 0.
     0. 0. 0. 0.
     0. 0. 0. 1.] ;d = (2,2)
t = (3,3) ;


#df = Examples.get_example_meta()
ex1_m = ex.get_example_meta("RCai")
d = ex1_m.Size[1]
ρ = ex.get_example("RCai")


ℂ1_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S1G",noBlock=false)
ℂ2_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S2G",noBlock=false)
ℂ3_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S3G",noBlock=false)

objective_value(ℂ1_sep_mod)
objective_value(ℂ2_sep_mod)
objective_value(ℂ3_sep_mod)

ℝ1_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S1G",noBlock=false,isRe = true)
ℝ2_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S2G",noBlock=false,isRe = true)
ℝ3_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S3G",noBlock=false,isRe = true)

objective_value(ℝ1_sep_mod)
objective_value(ℝ2_sep_mod)
objective_value(ℝ3_sep_mod)


println("Primal: ", primal_status(ℝ3_sep_mod))
println("Dual: ", dual_status(ℝ3_sep_mod))
## Batch run
include(srcDir*"batch.jl")
using .batch

df     = Examples.get_example_meta()
ρ_dict = ex.get_example()
t      = (3,3)
batch.batch_model(t,ρ_dict,df)

batch.batch_Computeξₜˢᵉᵖ("C:\\Users\\andries\\all-my-codes\\sep-rank\\assets\\bounds\\")
nar = batch.unstack_constraints(df)


nat = DataFrames.select(nar,[:Example,:Size,:Bi_rank,:CSᵦₐₗₗsG, :CS₂sG, :CS∞sG, :RSᵦₐₗₗsG, :RS₂sG, :RS∞sG,:isSeparable])
nat.isSeparable[poes.isSeparable .== "ent"] .= "∞"

CSV.write("C:\\Users\\andries\\all-my-codes\\sep-rank\\assets\\bounds\\t=$(t[1]).csv", poes, delim="&")

##   Size check
t = (3,3)
exa = "RCd" # "RCe"
ex1_m = Examples.get_example_meta(exa)
d = ex1_m.Size[1]
ρ = ex.get_example(exa)

model = JuMP.Model()
@variable(model, Lx[ccon.make_mon_expo_keys(d,t[1])])## Create variables

# size(mom.make_mon_expo(d,t[1]-2))
# t = (2,2)
# binomial(sum(d .*2) + t[1] ,t[1])
##  blocks
Gᴿ_con = ccon.make_Gᴿ_con(ρ,d,t,Lx;noBlock=false)
for b in keys(Gᴿ_con)
    println(size(Gᴿ_con[b]))
    # prod(size(Gᴿ_con[b]))
end
ADSA1 = sum([prod(size(Gᴿ_con[b])) for b in keys(Gᴿ_con)])
#
loc_cons_S1 = ccon.make_loc_cons_S1(ρ,d,t,Lx;noBlock=false)
for b in keys(loc_cons_S1)
    contains(b[2],"x²") ? println(b[1],size(loc_cons_S1[b])) : nothing
end
asdfas = [size(loc_cons_S1[b]) for b in keys(loc_cons_S1) if contains(b[2],"y²")  ]

ADSA2 = sum([contains(b[2],"x²") ? prod(size(loc_cons_S1[b])) : 0 for b in keys(loc_cons_S1)])
#
MMᴿcoef,MMᴿexᴿ = mom.get_ℂ_block_diag(d, t,noBlock=false)
for b in keys(MMᴿcoef)
    println(size(MMᴿcoef[b]))
end

ADSA3 = sum([prod(size(MMᴿcoef[b])) for b in keys(MMᴿcoef)])
ADSA1 + ADSA2 + ADSA3
## No blocks
Gᴿ_con = ccon.make_Gᴿ_con(ρ,d,t,Lx;noBlock=true)
for b in keys(Gᴿ_con)
    println(size(Gᴿ_con[b]))
end
ADSA1 = sum([prod(size(Gᴿ_con[b])) for b in keys(Gᴿ_con)])
#
loc_cons_S1 = ccon.make_loc_cons_S1(ρ,d,t,Lx;noBlock=true)
for b in keys(loc_cons_S1)
    contains(b[2],"x²") ? println(b[1],size(loc_cons_S1[b])) : nothing
end
ADSA2 = sum([contains(b[2],"x²") ? prod(size(loc_cons_S1[b])) : 0 for b in keys(loc_cons_S1)])
#
MMᴿcoef,MMᴿexᴿ = mom.get_ℂ_block_diag(d, t,noBlock=true)
for b in keys(MMᴿcoef)
    println(size(MMᴿcoef[b]))
end
ADSA3 = sum([prod(size(MMᴿcoef[b])) for b in keys(MMᴿcoef)])
ADSA1 + ADSA2 + ADSA3


## Running many random matrices:

using Random

Random.seed!(345)
function gen__ρ_RAND(d, r::Integer)
    d₁,d₂ = d ; ρ = zeros(d₁*d₂,d₁*d₂)
    a = [randn(1,d₁) + im*randn(1,d₁) for i in 1:r ]
    b = [randn(1,d₂) + im*randn(1,d₂) for i in 1:r ]

    for i in 1:r
        A = transpose(a[i])*conj(a[i])
        B = transpose(b[i])*conj(b[i])
        ρ = kron(A,B) + ρ
    end
    return ρ
end

d = (3,3)
t = 3


file_loc = dirname(dirname(@__FILE__))*"\\assets\\bounds\\Randomt$(t)_d_$d.csv"
touch(file_loc)
open(file_loc,"a") do io
    write(io, "|S1|S2|S3|T1|T2|T3|P1|P2|P3|D1|D2|D3|\n")
end



for i in 1:30
    ρ = gen__ρ_RAND(d, 5)
    ρ =  ρ / sum([ρ[i,i] for i in 1:size(ρ)[1]])

    O1 = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S1G",noBlock=false)
    O2 = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S2G",noBlock=false)
    # O3 = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list="S3G",noBlock=false)

    P1 = string(JuMP.primal_status(O1))
    D1 = string(JuMP.dual_status(O1))
    P2 = string(JuMP.primal_status(O2))
    D2 = string(JuMP.dual_status(O2))
    P3 = "-"#string(JuMP.primal_status(O3))
    D3 = "-"#string(JuMP.dual_status(O3))

    S1   = round(JuMP.objective_value(O1),digits=3)
    S2   = round(JuMP.objective_value(O2),digits=3)
    S3   = 0#round(JuMP.objective_value(O3),digits=3)

    T1 = JuMP.solution_summary(O1).solve_time
    T2 = JuMP.solution_summary(O2).solve_time
    T3 = 0#JuMP.solution_summary(O3).solve_time

    open(file_loc,"a") do io
        write(io, "|$S1|$S2|$S3|$T1|$T2|$T3|$P1|$P2|$P3|$D1|$D2|$D3|\n")
    end
end
