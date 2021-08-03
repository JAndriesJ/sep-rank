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
    using .Moments
    using .Utils_cons
    using .C_constraints
    using .R_constraints
    using .Utils_Model
    using .R_sep_Model
    using .C_sep_Model
    using .sep_Compute

    using LinearAlgebra
    using JuMP

##
