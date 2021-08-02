module testUtils_states
using Test
using LinearAlgebra

include("src\\Examples\\Utils_states.jl")
using .Utils_states



@testset "take_Pᵀ" begin
    d = Dict()
    d[1] = (2,2)
    d[2] = (2,3)
    d[3] = (3,3)
    d[4] = (3,4)

    ρ_test = Dict()
    ρ_test[1] =  [1  2   9  10
                  5  6  13  14
                  3  4  11  12
                  7  8  15  16]

    ρ_test[2] =     [ 1   2  13  14  25  26
                      5   6  17  18  29  30
                      9  10  21  22  33  34
                      3   4  15  16  27  28
                      7   8  19  20  31  32
                     11  12  23  24  35  36]

    ρ_test[3] =    [ 1   2   3  28  29  30  55  56  57
                    10  11  12  37  38  39  64  65  66
                    19  20  21  46  47  48  73  74  75
                    4   5   6  31  32  33  58  59  60
                    13  14  15  40  41  42  67  68  69
                    22  23  24  49  50  51  76  77  78
                    7   8   9  34  35  36  61  62  63
                    16  17  18  43  44  45  70  71  72
                    25  26  27  52  53  54  79  80  81]

    ρ_test[4] =     [1   2   3  37  38  39   73   74   75  109  110  111
                     10  11  12  46  47  48   82   83   84  118  119  120
                     19  20  21  55  56  57   91   92   93  127  128  129
                     28  29  30  64  65  66  100  101  102  136  137  138
                      4   5   6  40  41  42   76   77   78  112  113  114
                     13  14  15  49  50  51   85   86   87  121  122  123
                     22  23  24  58  59  60   94   95   96  130  131  132
                     31  32  33  67  68  69  103  104  105  139  140  141
                      7   8   9  43  44  45   79   80   81  115  116  117
                     16  17  18  52  53  54   88   89   90  124  125  126
                     25  26  27  61  62  63   97   98   99  133  134  135
                     34  35  36  70  71  72  106  107  108  142  143  144]


    gen_ρ(d) =  reshape(1:(d[1]*d[2])^2,d[1]^2,d[2]^2)

    ρ = gen_ρ(d[1])
    Pᵀρ =  Utils_states.take_Pᵀ(ρ, 2, d[1])

    for i in 1:4
        ρ = gen_ρ(d[i])
        Pᵀρ =  Utils_states.take_Pᵀ(ρ, 1, d[i])
        @test Pᵀρ == ρ_test[i]
    end
    d = (3,2)

    ρ = reshape(1:36,6,6)
    @test ρ == take_Pᵀ(ρ,0,d)
    ρᵀ¹ = take_Pᵀ(ρ,1,d)
    ρᵀ² = take_Pᵀ(ρ,2,d)
    @test ρᵀ¹ == [1   2   3  19  20  21
                    7   8   9  25  26  27
                    13  14  15  31  32  33
                    4   5   6  22  23  24
                    10  11  12  28  29  30
                    16  17  18  34  35  36]
    @test ρᵀ² == [  1   7  13   4  10  16
                    2   8  14   5  11  17
                    3   9  15   6  12  18
                    19  25  31  22  28  34
                    20  26  32  23  29  35
                    21  27  33  24  30  36]

end

@testset "isPPᵀ" begin
    ρ_1=   [ 0.0833333  0.0        0.0        0.0        0.0833333  0.0        0.0       0.0        0.0833333
             0.0        0.0833333  0.0        0.0        0.0        0.0        0.0       0.0        0.0
             0.0        0.0        0.0833333  0.0        0.0        0.0        0.0       0.0        0.0
             0.0        0.0        0.0        0.0833333  0.0        0.0        0.0       0.0        0.0
             0.0833333  0.0        0.0        0.0        0.0833333  0.0        0.0       0.0        0.0833333
             0.0        0.0        0.0        0.0        0.0        0.0833333  0.0       0.0        0.0
             0.0        0.0        0.0        0.0        0.0        0.0        0.208333  0.0        0.161374
             0.0        0.0        0.0        0.0        0.0        0.0        0.0       0.0833333  0.0
             0.0833333  0.0        0.0        0.0        0.0833333  0.0        0.161374  0.0        0.208333]


    ρ_2=   [ 0.5  0.0  0.0  0.0
             0.0  0.0  0.0  0.0
             0.0  0.0  0.0  0.0
             0.0  0.0  0.0  0.5]

    @test !Utils_states.isPPᵀ(ρ_1,(3,3))
    @test Utils_states.isPPᵀ(ρ_2,(2,2))


    ρ_true = ones(4,4)
    ρ_false = reshape(1:16,4,4)
    d = (2,2)
    @assert isPPᵀ(ρ_true,d) == true
    @assert isPPᵀ(ρ_false,d) == false

end

@testset "maketraceone" begin
    @test tr(Utils_states.maketraceone(reshape(1:81,9,9))) == 1
    @test tr(Utils_states.maketraceone(reshape(1:36,6,6))) == 1
    @test tr(Utils_states.maketraceone(rand(4,4))) ≈ 1
    @test tr(Utils_states.maketraceone(rand(4,4))) ≈ 1
    @test tr(Utils_states.maketraceone(rand(5,5))) ≈ 1
end

end
