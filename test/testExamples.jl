module testExamples
using LinearAlgebra
using Test

srcDir = dirname(dirname(@__FILE__))*"\\src\\Examples\\"
include(srcDir *"Examples_sep.jl")
include(srcDir *"Examples_ent.jl")
include(srcDir *"Examples.jl")

using .Examples_sep
using .Examples_ent
using .Examples

@testset "Test the call examples" begin
    ρ_ex = Examples.get_example()
    @test ρ_ex["RANDℂb5"] == Examples.get_example("RANDℂb5")
    @test ρ_ex["RANDℂb5"] == Examples.get_example(5)
end

@testset "Test the call example meta data" begin
    ρ_ex = Examples.get_example()
    df = Examples.get_example_meta()
    @test df[df.Example .== "Eq14Hor97a",:] ==  Examples.get_example_meta("Eq14Hor97a")
    @test df[df.Example .== "Eq14Hor97a",:].Example[1] ==  Examples.get_example_meta(1).Example
    length(ρ_ex) == size(df)[1]
end

##
df = ex.get_example_meta()
entdf = df[df.isSeparable .== "ent",:]
ρ_dict = ex.get_example()

for n_ex in 1:23
    exa = entdf.Example[n_ex]
    ρ = ρ_dict[exa]
    d = entdf.Size[n_ex]
    println("$(entdf.Example[n_ex]) $(la.eigmin(ρ))")
    @assert us.isPPᵀ(ρ,d) == entdf.hasPPᵀ[n_ex]
end

##
df = ex.get_example_meta()
sepdf = df[df.isSeparable .== "sep",:]
ρ_dict = ex.get_example()

for n_ex in 1:36
    exa = sepdf.Example[n_ex]
    ρ = ρ_dict[exa]
    d = sepdf.Size[n_ex]
    println("$(sepdf.Example[n_ex]) $(la.eigvals(ρ)[1])")
    @assert us.isPPᵀ(ρ,d) == sepdf.hasPPᵀ[n_ex]
end



end
