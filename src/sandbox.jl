### Begin here
using Pkg
# Pkg.activate(".")
Pkg.status()

# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("MosekTools")
# Pkg.add("Test")

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

    using LinearAlgebra ; const la = LinearAlgebra # Pkg.add("LinearAlgebra")
    using JuMP # Pkg.add("JuMP")

## Goal get the C block diag to get the same values as the non- block diag.

ρ  =  [1. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 1.]
d = (2,2) ; t = (2,2) ; cl = "S∞" ; n = sum(2 .* d)  #"S∞"

ℂb_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list=cl,noBlock=false)
ℂ_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list=cl,noBlock=true)
# ℝ_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list=cl,isRe=true)


NBD_PSD_con = ccon.make_PSD_con(d,t,ℂb_sep_mod[:Lx];noBlock=true)
BD_PSD_con = C_constraints.make_PSD_con(d,t,ℂb_sep_mod[:Lx];noBlock=false)

NBD = NBD_PSD_con["Default"]
sep_Compute.get_sol_vals(NBD)
sep_Compute.get_sol_min_eigval(NBD)

BD = BD_PSD_con[0,0]
sep_Compute.get_sol_vals(BD)
sep_Compute.get_sol_min_eigval(BD)


uc.idx2varxx̄ᵀtyȳᵀ(Lx,mom.make_xx̄ᵀ⨂yȳᵀ(d))

L_xx̄ᵀ_tens_yȳᵀ = ccon.make_ord4_con(d,ℂb_sep_mod[:Lx])

sep_Compute.get_sol_vals(L_xx̄ᵀ_tens_yȳᵀ["imag"])

sep_Compute.get_sol_vals(L_xx̄ᵀ_tens_yȳᵀ["real"])



B = mom.make_xx̄ᵀ⨂yȳᵀ(d)
ord4_con = ccon.make_ord4_con(d,ℂb_sep_mod[:Lx])
ord4_con["real"]
ord4_con["real"][1,1]


split_expo(ααᶥββᶥ,d₁,d₂) = (ααᶥββᶥ[1:d₁],ααᶥββᶥ[d₁+1:2*d₁],
                            ααᶥββᶥ[1+2*d₁:2*d₁+d₂],ααᶥββᶥ[1+2*d₁+d₂:end])


## Approach 2

## Approach 2
d = (2,2)
n = sum(2 .* d); p_ex = make_mon_expo(n,1;isle = false)
x, x̄, y, ȳ = split_expo(p_ex,d...)
B = var_kron(x .+ reshape(x̄,1,:), y .+ reshape(ȳ,1,:))



M = Moments.make_mon_expo(d,2 .* t)
M_vec = M[:,1]


k = [0, 0, 0, 0, 0, 0, 3, 1]
XX = findall(M .== [k])
fd = hcat(map(x -> [x[1],x[2]],findall(M .== [k]))...)
γδ_dict[k] = [M_vec[fd[1,:]], M_vec[fd[2,:]]]

γδ_dict = Dict()
for k in Moments.make_mon_expo(d,2*t[1])
    fd = hcat(map(x -> [x[1],x[2]],findall(M .== [k]))...)
    γδ_dict[k] = [M_vec[fd[1,:]], M_vec[fd[2,:]]]
end


γδ_dict[B[1,1]]







γδ_dict = mom.get_γδ_dict(d,t)
(ℝcoef,ℝexpo),(ℂcoef,ℂexpo) = Moments.get_ℜBℑB(d,B,γδ_dict)







Cᴿ[1:4,1:4],Bᴿ[1:4,1:4]



## Block diag
# C,E = mom.get_ℂ_block_diag(d,t;noBlock = true)
# Cb,Eb = Moments.get_ℂ_block_diag(d,t;noBlock = false)
## Moments

MM_b = Moments.get_xx̄yȳMM_blocks(d,t)
b_ind =[j for k in keys(MM_b_ind) for j in MM_b_ind[k]]


BD_PSD_con = ccon.make_PSD_con(d,t,ℂb_sep_mod[:Lx];noBlock=true)
B = BD_PSD_con["Default"]
sep_Compute.get_sol_vals(B)

sep_Compute.get_sol_vals(ℂb_sep_mod[:Lx])


Moments.make_xxᵀ⨂yyᵀ(d)




# JuMP.value.(PSD_con[0,1])
# sep_Compute.get_sol_vals(PSD_con[0,1])
for b in keys(BD_PSD_con)
    sep_Compute.get_sol_min_eigval(BD_PSD_con[b]) < 0 ? println(b,sep_Compute.get_sol_min_eigval(BD_PSD_con[b])) : nothing
end

sep_Compute.get_sol_vals(BD_PSD_con[0, 1])

PSD_conNOBLOCK = ccon.make_PSD_con(d,t,ℂb_sep_mod[:Lx];noBlock=true)
SMM = sep_Compute.get_sol_vals(PSD_conNOBLOCK["Default"])
sep_Compute.get_sol_min_eigval(PSD_conNOBLOCK["Default"])
sep_Compute.get_sol_vals(PSD_conNOBLOCK["Default"][P,P])





Iᵗ = make_mon_expo(d,t[1])
Iᵗ_temp = map(x ->  sum.(split_expo(x,d...)),Iᵗ)
r_list  = map(v -> v[1]-v[2],Iᵗ_temp) ; s_list  = map(v -> v[3]-v[4],Iᵗ_temp)


r = 0 ; s = 1
T = findall((r_list .== r) .& (s_list .== s))

P = []
S = []
for r ∈ -t[1]:t[1], s ∈ -t[1]:t[1]
    T = findall((r_list .== r) .& (s_list .== s))
    isempty(T) ? nothing : push!(P,T...); push!(S,(r,s))
end

P













γδ_dict = get_γδ_dict(d,t)

γ,δ = γδ_dict[b]


model = JuMP.Model()
@variable(model, Lx[ccon.make_mon_expo_keys(d,t[1])] )## Create variables
PSD_con = ccon.make_PSD_con(d,t,Lx;noBlock=false)


function get_expo(p1,p2,d)
    γ,γᶥ,ζ,ζᶥ = split_expo(p1,d...) ; δ,δᶥ,η,ηᶥ = split_expo(p2,d...)
    return vcat(γ+γᶥ,δ+δᶥ,ζ+ζᶥ,η+ηᶥ)
end



get_expo(p1p2_arr,d) = map((x,y)->get_expo(x,y,d),p1p2_arr...)



γδ = γδ_dict[xx̄yȳMM_blocks[0,1][2,1]]

Eb[0,1][2,1]
Eb[0,1][4,1]

d = (2,2) ; t = (2,2)
γδ_dict = get_γδ_dict(d,t)
γδ = γδ_dict[[0, 0, 0, 0, 0, 1, 1, 0]]
@assert get_expo(γδ,d) == [[0, 0, 0, 0, 1, 1, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 1],
                         [0, 0, 0, 0, 0, 1, 1, 0],
                         [0, 0, 0, 0, 0, 0, 1, 1]]

Cb[0,1]
Eb[0,1]
PSD_con[0,1][1,2]
## Utils_cons
