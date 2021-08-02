# sep-rank


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

