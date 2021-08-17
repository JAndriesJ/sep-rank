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

    using .Examples ; const ex = Examples
    using .Utils_states ; const us = Utils_states
    using .Moments ; const mom = Moments
    using .Utils_cons ; const uc = Utils_cons
    using .C_constraints ; const ccon = C_constraints
    using .R_constraints ; const rcon = C_constraints
    using .Utils_Model
    using .R_sep_Model
    using .C_sep_Model
    using .sep_Compute

    using LinearAlgebra ; const la = LinearAlgebra # Pkg.add("LinearAlgebra")
    using JuMP

## Goal get the C block diag to get the same values as the non- block diag.

ρ =   [1. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 1.]
d = (2,2)
t = (3,3) ; cl = "S1sG"  #"S∞ S₂ Sᵦₐₗₗ sG"

# t = (3,3)

df = Examples.get_example_meta()
ex1_m = ex.get_example_meta("E7_i")
d = ex1_m.Size[1]
ρ = ex.get_example("E7_i")

ℂb_sep_mod = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list=cl,noBlock=false)
ℂ_sep_mod  = sep_Compute.Computeξₜˢᵉᵖ(ρ,d,t;con_list=cl,noBlock=true)
##  Size check
t = (3,3)
exa = "RCai"
ex1_m = ex.get_example_meta(exa)
d = ex1_m.Size[1]
ρ = ex.get_example(exa)

model = JuMP.Model()
@variable(model, Lx[ccon.make_mon_expo_keys(d,t[1])])## Create variables

##  blocks
Gᴿ_con = ccon.make_Gᴿ_con(ρ,d,t,Lx;noBlock=false)
for b in keys(Gᴿ_con)
    println(size(Gᴿ_con[b]))
end

loc_cons_S1 = ccon.make_loc_cons_S1(ρ,d,t,Lx;noBlock=false)
for b in keys(loc_cons_S1)
    contains(b[2],"x²") ? println(b[1],size(loc_cons_S1[b])) : nothing
end

[keys(loc_cons_S1)...]


MMᴿcoef,MMᴿexᴿ = mom.get_ℂ_block_diag(d, t,noBlock=false)

for b in keys(MMᴿcoef)
    println(size(MMᴿcoef[b]))
end



## No blocks

Gᴿ_con = ccon.make_Gᴿ_con(ρ,d,t,Lx;noBlock=true)
for b in keys(Gᴿ_con)
    println(size(Gᴿ_con[b]))
end

loc_cons_S1 = ccon.make_loc_cons_S1(ρ,d,t,Lx;noBlock=true)
for b in keys(loc_cons_S1)
    contains(b[2],"x²") ? println(b[1],size(loc_cons_S1[b])) : nothing
end

[keys(loc_cons_S1)...]


MMᴿcoef,MMᴿexᴿ = mom.get_ℂ_block_diag(d, t,noBlock=true)

for b in keys(MMᴿcoef)
    println(size(MMᴿcoef[b]))
end





## Batch run
include(srcDir*"batch.jl")
using .batch

df     = ex.get_example_meta()
ρ_dict = ex.get_example()
t      = (2,2)
batch.batch_model(t,ρ_dict,df)
batch.batch_Computeξₜˢᵉᵖ("C:\\Users\\andries\\all-my-codes\\sep-rank\\assets\\bounds\\")
nar = batch.unstack_constraints(df)

nar
nat = DataFrames.select(nar,[:Example,:Size,:Bi_rank,:CSᵦₐₗₗsG, :CS₂sG, :CS∞sG, :RSᵦₐₗₗsG, :RS₂sG, :RS∞sG,:isSeparable])
nat.isSeparable[poes.isSeparable .== "ent"] .= "∞"

CSV.write("C:\\Users\\andries\\all-my-codes\\sep-rank\\assets\\bounds\\t=$(t[1]).csv", poes, delim="&")
## The operating table:
using CSV
using DataFrames


temp_df  = CSV.read("C:\\Users\\andries\\all-my-codes\\sep-rank\\assets\\bounds\\"*"Summary.csv", DataFrame)
temp2_df = select(innerjoin(df, temp_df, on = :Example),[:Example,:isReal,:isSeparable,:Model,:Size,:Bi_rank,:Constraint,:obj_val])
temp4_df = unstack(temp2_df, :Model, :obj_val)

temp4R_df = select(temp4_df,[:Example,:isReal,:isSeparable,:Size,:Bi_rank,:Constraint,:R])
temp4Rus_df = unstack(temp4R_df, :Constraint, :R)
rename!(temp4Rus_df,Dict(:S₂sG => :RS₂sG, :SᵦₐₗₗsG => :RSᵦₐₗₗsG, :S∞sG=>:RS∞sG));

temp4C_df = select(temp4_df,[:Example,:isReal,:isSeparable,:Size,:Bi_rank,:Constraint,:C])
temp4Cus_df = unstack(temp4C_df, :Constraint, :C)
rename!(temp4Cus_df,Dict(:S₂sG => :CS₂sG, :SᵦₐₗₗsG => :CSᵦₐₗₗsG, :S∞sG=>:CS∞sG));

temp5_df = innerjoin(temp4Cus_df,select(temp4Rus_df,[:Example,:RS₂sG, :RSᵦₐₗₗsG, :RS∞sG]), on = :Example)

CSV.write(bDir*"SummaryUnstacked.csv", temp5_df, delim="&")
return temp5_df




## 3



## Old patients
MM = mom.make_mon_expo(d,t)
m = MM[44,43]


n = length(m)
nz_coord = findall(.!(m .== 0))
isempty(nz_coord) ? (return [m]) : nothing
c = pop!(nz_coord)
Nar_old = [ j*mom.eᵢ(n,c) for j in 0:m[c]]
global Nar_old
while  !isempty(nz_coord)
    c = pop!(nz_coord)
    Nar_new = [ j*mom.eᵢ(n,c) for j in 0:m[c]]
    Nar_old = [ [n + m] for n in Nar_old for m in Nar_new ]
end

[m] .- Moments.get_γδ(m)

MM_vec = make_mon_expo(d,t[1]*2)
Dict(zip(MM_vec,map(m ->([m] .- Moments.get_γδ(m),Moments.get_γδ(m)) ,MM_vec)))


γδ_dict = Moments.get_γδ_dict((2,2),(2,2))
γδ_dict_2 = Moments.get_γδ_dict_2((2,2),(2,2))

γδ = γδ_dict[[0, 0, 0, 1, 2, 0, 1, 0]]
γδ_2 = γδ_dict_2[[0, 0, 0, 1, 2, 0, 1, 0]]

γδ[1]
γδ_2[1]

using JLD2, FileIO ; # Pkg.add("JLD2") ; Pkg.add("FileIO")









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


## Approach 1
B1 = mom.make_xx̄ᵀ⨂yȳᵀ(d)
## Approach 2
d = (2,2)
p_ex = make_mon_expo(n,1;isle = false)
x, x̄, y, ȳ = split_expo(p_ex,d...)
B2temp = var_kron(x .+ reshape(x̄,1,:), y .+ reshape(ȳ,1,:))
γδ_dict = mom.get_γδ_dict(d,t)
B2c,B2e = Moments.get_ℜℑααᶥββᶥᴿ(d,B2temp,γδ_dict)





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
