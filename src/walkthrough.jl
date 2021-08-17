# SEP-RANK
## Greeting and welcome to my code for computing lowerboounds on the separable rank.

Pkg.activate(".")
    srcDir  = dirname(@__FILE__)*"\\";
    include(srcDir*"Examples\\Examples.jl")
    include(srcDir*"Model\\R_sep_Model.jl")
    include(srcDir*"Model\\C_sep_Model.jl")
    include(srcDir*"sep_Compute.jl")

    using .Examples
    using .R_sep_Model
    using .C_sep_Model
    using .sep_Compute
### Quick start
#### Manually load data and run the model:
ρ  =  [1. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 1.]
d₁,d₂  = (2,2) ; cl = "S1sG" ; t  = (2,2)
sep_Cmod = C_sep_Model.Modelξₜˢᵉᵖ(ρ,(d₁,d₂),t;con_list = cl)
sep_Cmod_sol = sep_Compute.Computeξₜˢᵉᵖ(sep_Cmod)

Lx = sep_Cmod[:Lx]

#### Run the **ℝeal** model with the same data as above:
sep_Rmod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,(d₁,d₂),t;con_list = cl)
sep_Rmod_sol = sep_Compute.Computeξₜˢᵉᵖ(sep_Rmod)

#### Load example data and run complex model
ex = "E7_i"
d = ex.get_example_meta(ex).Size[1] ; ρ = ex.get_example(ex) ; cl = "S1sG" ; t  = (2,2)
sep_Cmod = C_sep_Model.Modelξₜˢᵉᵖ(ρ,(d₁,d₂),t;con_list = cl)
sep_Cmod_sol = sep_Compute.Computeξₜˢᵉᵖ(sep_Cmod)

### Details ***(optional)***
#### Preloaded Examples 
df = Examples.get_example_meta()
ρ_dict = Examples.get_example()

df["E7_i",:]
ρ_dict["E7_i"]
#Look up example in appendix.


#### Running a batch of examples 
include(srcDir*"batch.jl")
using .batch

batch.batch_model(t,ρ_dic,df;sdir = save_directory,con_list= ["S1sG","S2sG","S3sG"])


#### Loading the constraint matrices
include(srcDir*"Constraints\\R_constraints.jl")
include(srcDir*"Constraints\\C_constraints.jl")
using .C_constraints 
using .R_constraints 
using JuMP

model = JuMP.Model() ; @variable(model, Lx[ccon.make_mon_expo_keys(d,t[1])])

#### The **ℂomplex** moment matrix, 4ᵗʰ-order constraints, localizing maps and matrix localizing map.
C_PSD_con = C_constraints.make_PSD_con(d,t,Lx)
C_ord4_con = C_constraints.make_ord4_con(d,Lx)
C_loc_cons_S1 = C_constraints.make_loc_cons_S1(ρ,d,t,Lx)
C_loc_cons_S2 = C_constraints.make_loc_cons_S2(ρ,d,t,Lx)
C_loc_cons_S3 = C_constraints.make_loc_cons_S3(ρ,d,t,Lx)
C_Gᴿ_con = C_constraints.make_Gᴿ_con(ρ,d,t,Lx)

#### The **ℝeal** moment matrix, 4ᵗʰ-order constraints, localizing maps and matrix localizing map.
R_PSD_con = R_constraints.make_PSD_con(d,t,Lx)
R_ord4_con = R_constraints.make_ord4_con(d,Lx)
R_loc_cons_S1 = R_constraints.make_loc_cons_S1(ρ,d,t,Lx)
R_loc_cons_S2 = R_constraints.make_loc_cons_S2(ρ,d,t,Lx)
R_loc_cons_S3 = R_constraints.make_loc_cons_S3(ρ,d,t,Lx)
R_G_con = R_constraints.make_G_con(ρ,d,t,Lx)

