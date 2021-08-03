module sep_Compute
using MosekTools ; using JuMP
using LinearAlgebra ; const la = LinearAlgebra

srcDir = dirname(dirname(@__FILE__))*"\\src\\"
include(srcDir*"Model\\C_sep_Model.jl")
include(srcDir*"Model\\R_sep_Model.jl")
using .C_sep_Model
using .R_sep_Model

export Computeξₜˢᵉᵖ,
       get_sol_vals,
       get_sol_min_eigval

"""
Inpute: Jump model of separability problem.
Output: MOSEK solution response.
"""
function Computeξₜˢᵉᵖ(model)
    JuMP.set_optimizer(model, Mosek.Optimizer)
    # optimize
    optimize!(model)

    # output results
    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    return model
end
Computeξₜˢᵉᵖ(ρ,d,t;con_list ="S∞ sG",noBlock = false) = Computeξₜˢᵉᵖ(C_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list,noBlock))

get_sol_vals(arr) = JuMP.value.(arr)
get_sol_min_eigval(arr) = la.eigmin(JuMP.value.(arr))

end
