module testExamples
using LinearAlgebra
using Test

srcDir = pwd()*"\\src\\Examples\\"
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
    df = Examples.get_example_meta()
    @test df[df.Example .== "Eq14Hor97a",:] ==  Examples.get_example_meta("Eq14Hor97a")
    @test df[df.Example .== "Eq14Hor97a",:] ==  Examples.get_example_meta(1)
    length(ρ_ex) == size(df)[1]
end


end
