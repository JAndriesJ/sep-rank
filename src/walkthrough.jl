## The real model.
println("""Specify: ρ, registers: d, constraint list: cl,level of hierarchy: t
ρ =  [1 0 0 0
      0 0 0 0
      0 0 0 0
      0 0 0 1]
d  = (2,2)
cl = "S∞ sG" or "S₂ sG" or "S₂₁ sG"
t  = (2,2)""")
ρ  =  [1. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 1.]
d  = (2,2)
cl = "S∞ sG"
t  = (3,3)
println("""Generate the model:
sep_mod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list = cl)
""")
sep_mod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list = cl)
display(sep_mod)

## Solve the OP.
println("We initialize the computation module")
sep_mod = sep_Compute.Computeξₜˢᵉᵖ(sep_mod)
Lx = sep_mod[:Lx]
## A closer look at the real constraints
println("""Let us take a deeper look into the constraints:""")
println("These are $(binomial(sum(d)+2*t[1],2*t[1])) the moments present in the hierachy at level t = (2,2)")
mon_expo_keys = R_constraints.make_mon_expo_keys(sum(d),t[1])
display(mon_expo_keys)

println("""These are the blocks of the 'L([x,y]ₜ[x,y]ₜᵀ)'.
 """)

zero_moms = R_constraints.zeroprop(ρ,d,(3,3),Lx)


PSD_con = R_constraints.make_PSD_con(d,t,Lx)
display(PSD_con)
println("""
"ee": moments of the form xᵅyᵝ ⋅ xᵞyᵟ where |α|,|β|,|γ|,|δ| are all even.
"eo": moments of the form xᵅyᵝ ⋅ xᵞyᵟ where |α|,|γ| even and |β|,|δ| odd.
"oe": moments of the form xᵅyᵝ ⋅ xᵞyᵟ where |α|,|γ| odd and |β|,|δ| even.
"oo": moments of the form xᵅyᵝ ⋅ xᵞyᵟ where |α|,|β|,|γ|,|δ| are all odd.
 """)
display(PSD_con["ee"])
display(PSD_con["eo"])
display(PSD_con["oe"])
display(PSD_con["oo"])

println("This is 'L(xxᵀ⊗yyᵀ)'.")
ord4_con = R_constraints.make_ord4_con(d,Lx)
display(ord4_con)

println("This is the 'L(( √ρₘₐₓ - xᵢ²)⋅η), L(( √ρₘₐₓ - yⱼ²)⋅η) ⪰ 0 ∀ i ∈ [d₁], j ∈ [d₂]' constraints.")
cons_S_inf = R_constraints.make_loc_cons_S_inf(ρ,d,t,Lx)
display(cons_S_inf)
display(cons_S_inf[("ee", "x²_1")])
display(cons_S_inf[("eo", "x²_1")])
display(cons_S_inf[("oe", "x²_1")])
# display(cons_S_inf[("oo", "x_1")])

println("This is the L(( (Tr(ρ)) - ∑xᵢ²)⋅η), L(( (Tr(ρ)) - ∑yⱼ²)⋅η) ⪰ 0' constraints.")
cons_S₂    = R_constraints.make_loc_cons_S₂(ρ,d,t,Lx)
display(cons_S₂)
display(cons_S₂[("ee", "x")])
display(cons_S₂[("eo", "x")])
display(cons_S₂[("oe", "x")])

println("This is the ' L(( (Tr(ρ)) - ∑xᵢ²)⋅η) ⪰ 0, L(( 1 - ∑yⱼ²)⋅η) = 0' constraints.")
cons_S₂₁x,cons_S₂₁y   = R_constraints.make_loc_cons_S₂₁(ρ,d,t,Lx)
display(cons_S₂₁x)
display(cons_S₂₁y)


println("This is the 'ρ⊗L(η) - L( (xxᵀ⊗yyᵀ) ⊗ η ) ⪰ 0 ∀ η ∈ even-degree-principle-submatrices of ([x,y]₌ₜ₋₂[x,y]₌ₜ₋₂ᵀ)' constraints.")
G_con = R_constraints.make_G_con(ρ,d,t,Lx)
display(G_con["ee"])
## Inspecting the solution
using JuMP
# println("We can recover the moments as with the following command: Dict(zip(Lx.data,value.(Lx.data)))")
# display(Dict(zip(Lx.data,value.(Lx.data))))
#
# display(value.(PSD_con["ee"]))
#
# println("The fourth order constraint should equal ρ")
# display(value.(ord4_con))
#
# println("One of the S_inf constraints is")
# display(value.(cons_S_inf[("oe", "x²_2")]))
#
# println("The G constriant blocks are")
# display(value.(G_con["ee"]))
## The Real model quick
ex = ["Eq22HHH96a" "RANDa3" "p12HK14b" "T1CD12i"][1]
ρ = Examples.get_example(7)
d = (3,3)#filter(:Example => ==(ex),df).Size[1]
t = (2,2)
cl = "S₂"#"S∞ sG" # "S₂ sG" # "S₂₁ sG"
sepr_mod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list = cl)
sepr_mod = sep_Compute.Computeξₜˢᵉᵖ(sepr_mod)


display(filter(:Example => ==(ex),df))

#Lrx = sepr_mod[:Lx]
## The complex model quick
include(srcDir*"Model\\C_sep_Model.jl")
using .C_sep_Model
ex = ["RANDℂb2" "Eq2-3BP00b" "RANDa3"][3]
ρ = Examples.get_example(ex)
d = filter(:Example => ==(ex),df).Size[1]
t = (2,2)
cl = "S∞sG" # "S₂ sG" # "S₂₁ sG"
sep_mod = C_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list = cl)
sep_mod = sep_Compute.Computeξₜˢᵉᵖ(sep_mod)
display(filter(:Example => ==(ex),df))

Lx = sep_mod[:Lx]
## Running things in batches
include(srcDir*"batch.jl")
using .batch
t = (2,2)
batch.batch_model(t,Examples.get_examples(),df)
batch.batch_Computeξₜˢᵉᵖ("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\assets\\bounds\\")
nar = batch.unstack_constraints(df)


## A closer look at the complex constraints
println("""Let us take a deeper look into the constraints:""")
include(srcDir*"Constraints\\C_constraints.jl")
using .C_constraints
#println("These are $(binomial(sum(d)+2*t[1],2*t[1])) the moments present in the hierachy at level t = (2,2)")
mon_expo_keys = C_constraints.make_mon_expo_keys(d,t[1])
display(mon_expo_keys)


println("""These are the blocks of the 'L([x,y]ₜ[x,y]ₜᵀ)'.
 """)

#zero_moms = C_constraints.zeroprop(ρ,d,(3,3),Lx)
PSD_con = C_constraints.make_PSD_con(d,t,Lx)
display(PSD_con)

println("This is 'L(xx*⊗yy*)'.")
ord4_con = C_constraints.make_ord4_con(d,Lx)
display(ord4_con)

println("This is the 'L(( √ρₘₐₓ - xᵢ²)⋅η), L(( √ρₘₐₓ - yⱼ²)⋅η) ⪰ 0 ∀ i ∈ [d₁], j ∈ [d₂]' constraints.")
cons_S_inf = C_constraints.make_loc_cons_S_inf(ρ,d,t,Lx)
display(cons_S_inf)


println("This is the L(( (Tr(ρ)) - ∑xᵢ²)⋅η), L(( (Tr(ρ)) - ∑yⱼ²)⋅η) ⪰ 0' constraints.")
cons_S₂    = C_constraints.make_loc_cons_S₂(ρ,d,t,Lx)
display(cons_S₂)


println("This is the ' L(( (Tr(ρ)) - ∑xᵢ²)⋅η) ⪰ 0, L(( 1 - ∑yⱼ²)⋅η) = 0' constraints.")
cons_S₂₁x,cons_S₂₁y   = C_constraints.make_loc_cons_S₂₁(ρ,d,t,Lx)
display(cons_S₂₁x)
display(cons_S₂₁y)


println("This is the 'ρ⊗L(η) - L( (xxᵀ⊗yyᵀ) ⊗ η ) ⪰ 0 ∀ η ∈ even-degree-principle-submatrices of ([x,y]₌ₜ₋₂[x,y]₌ₜ₋₂ᵀ)' constraints.")
G_con = C_constraints.make_G_con(ρ,d,t,Lx)
display(G_con)
