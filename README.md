# sep-rank

# SEP-RANK
Greeting and welcome to my code for computing lowerboounds on the separable rank.
see [my talk](https://www.youtube.com/watch?v=46ciiXD3zUk) for a 25 minute talk on the topic.

### Quick start
Start by loading the modules
```Julia
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
```

#### Option 1. Manually load data and run the model:
```Julia
ρ  =  [1. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 1.]
d₁,d₂  = (2,2) ; cl = "S1sG" ; t  = (2,2)
sep_Cmod = C_sep_Model.Modelξₜˢᵉᵖ(ρ,(d₁,d₂),t;con_list = cl)
sep_Cmod_sol = sep_Compute.Computeξₜˢᵉᵖ(sep_Cmod)
```


#### Run the **ℝeal** model with the same data as above:
```Julia
sep_Rmod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,(d₁,d₂),t;con_list = cl)
sep_Rmod_sol = sep_Compute.Computeξₜˢᵉᵖ(sep_Rmod)
```

#### Load example data and run complex model
```Julia
ex = "E7_i"
d = ex.get_example_meta(ex).Size[1] ; ρ = ex.get_example(ex) ; cl = "S1sG" ; t  = (2,2)
sep_Cmod = C_sep_Model.Modelξₜˢᵉᵖ(ρ,(d₁,d₂),t;con_list = cl)
sep_Cmod_sol = sep_Compute.Computeξₜˢᵉᵖ(sep_Cmod)
```

### Details ***(optional)***
#### Preloaded Examples (Look up example in appendix [][])
```Julia
df = Examples.get_example_meta()
ρ_dict = Examples.get_example()

df["E7_i",:]
ρ_dict["E7_i"]
```

#### Running a batch of examples 
```Julia
include(srcDir*"batch.jl")
using .batch

batch.batch_model(t,ρ_dic,df;sdir = save_directory,con_list= ["S1sG","S2sG","S3sG"])
```

#### Loading the constraint matrices
```Julia
include(srcDir*"Constraints\\R_constraints.jl")
include(srcDir*"Constraints\\C_constraints.jl")
using .C_constraints 
using .R_constraints 
using JuMP

model = JuMP.Model() ; @variable(model, Lx[ccon.make_mon_expo_keys(d,t[1])])
```

#### The **ℂomplex** moment matrix, 4ᵗʰ-order constraints, localizing maps and matrix localizing map.
```Julia
C_PSD_con = C_constraints.make_PSD_con(d,t,Lx)
C_ord4_con = C_constraints.make_ord4_con(d,Lx)
C_loc_cons_S1 = C_constraints.make_loc_cons_S1(ρ,d,t,Lx)
C_loc_cons_S2 = C_constraints.make_loc_cons_S2(ρ,d,t,Lx)
C_loc_cons_S3 = C_constraints.make_loc_cons_S3(ρ,d,t,Lx)
C_Gᴿ_con = C_constraints.make_Gᴿ_con(ρ,d,t,Lx)
```
#### The **ℝeal** moment matrix, 4ᵗʰ-order constraints, localizing maps and matrix localizing map.
```Julia
R_PSD_con = R_constraints.make_PSD_con(d,t,Lx)
R_ord4_con = R_constraints.make_ord4_con(d,Lx)
R_loc_cons_S1 = R_constraints.make_loc_cons_S1(ρ,d,t,Lx)
R_loc_cons_S2 = R_constraints.make_loc_cons_S2(ρ,d,t,Lx)
R_loc_cons_S3 = R_constraints.make_loc_cons_S3(ρ,d,t,Lx)
R_G_con = R_constraints.make_G_con(ρ,d,t,Lx)
```





## An example where normalizing ρ to have ones on the diagonal breaks the separability:

``` Julia
include(srcDir*"Examples\\Utils_states.jl")
using .Utils_states
using LinearAlgebra

ρ =    [2 3 4 3 4 5 4 5 6
        3 5 7 4 6 8 5 7 9
        4 7 10 5 8 11 6 9 12
        3 4 5 5 6 7 7 8 9
        4 6 8 6 8 10 8 10 12
        5 8 11 7 10 13 9 12 15
        4 5 6 7 8 9 10 11 12
        5 7 9 8 10 12 11 13 15
        6 9 12 9 12 15 12 15 18]
ρᵀᵇ =  Utils_states.take_Pᵀ(ρ,1,(3,3))

m = diagm(sqrt.(1 ./ diag(ρ)))
mρm = m*ρ*m
mρmᵀᵇ =  Utils_states.take_Pᵀ(mρm,1,(3,3))

eigvals(ρ)[1]
eigvals(ρᵀᵇ)[1]

eigvals(m*ρ*m)[1]
eigvals(m*ρᵀᵇ*m)[1]

eigvals(mρmᵀᵇ)[1]
```

