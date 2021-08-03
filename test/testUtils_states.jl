module testUtils_states
using Test
using LinearAlgebra

cd("..") 
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

    ρ_test[2] =     [   1   2   3  19  20  21
                        7   8   9  25  26  27
                        13  14  15  31  32  33
                        4   5   6  22  23  24
                        10  11  12  28  29  30
                        16  17  18  34  35  36]
    ρ_test[3] =    [ 1   2   3  28  29  30  55  56  57
                    10  11  12  37  38  39  64  65  66
                    19  20  21  46  47  48  73  74  75
                    4   5   6  31  32  33  58  59  60
                    13  14  15  40  41  42  67  68  69
                    22  23  24  49  50  51  76  77  78
                    7   8   9  34  35  36  61  62  63
                    16  17  18  43  44  45  70  71  72
                    25  26  27  52  53  54  79  80  81]

    ρ_test[4] =      [  1   2   3   4  49  50  51  52   97   98   99  100
                        13  14  15  16  61  62  63  64  109  110  111  112
                        25  26  27  28  73  74  75  76  121  122  123  124
                        37  38  39  40  85  86  87  88  133  134  135  136
                        5   6   7   8  53  54  55  56  101  102  103  104
                        17  18  19  20  65  66  67  68  113  114  115  116
                        29  30  31  32  77  78  79  80  125  126  127  128
                        41  42  43  44  89  90  91  92  137  138  139  140
                        9  10  11  12  57  58  59  60  105  106  107  108
                        21  22  23  24  69  70  71  72  117  118  119  120
                        33  34  35  36  81  82  83  84  129  130  131  132
                        45  46  47  48  93  94  95  96  141  142  143  144]


    gen_ρ(d) =  reshape(1:(d[1]*d[2])^2,d[1]*d[2],d[1]*d[2])

    ρ = gen_ρ(d[1])
    Pᵀρ =  Utils_states.take_Pᵀ(ρ, 2, d[1])

    for i in 1:4
        @test Utils_states.take_Pᵀ(gen_ρ(d[i]), 1, d[i]) == ρ_test[i]
    end

    d = (3,2)
    ρ = reshape(1:36,6,6)
    @test ρ == Utils_states.take_Pᵀ(ρ,0,d)
    ρᵀ¹ = Utils_states.take_Pᵀ(ρ,1,d)
    ρᵀ² = Utils_states.take_Pᵀ(ρ,2,d)
    @test ρᵀ¹ ==   [1   2  13  14  25  26
                    7   8  19  20  31  32
                    3   4  15  16  27  28
                    9  10  21  22  33  34
                    5   6  17  18  29  30
                    11  12  23  24  35  36]
    @test ρᵀ² == [  1   7   3   9   5  11
                    2   8   4  10   6  12
                    13  19  15  21  17  23
                    14  20  16  22  18  24
                    25  31  27  33  29  35
                    26  32  28  34  30  36]

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
    @test Utils_states.isPPᵀ(ρ_true,d) == true
    @test Utils_states.isPPᵀ(ρ_false,d) == false

end

@testset "maketraceone" begin
    @test tr(Utils_states.maketraceone(reshape(1:81,9,9))) == 1
    @test tr(Utils_states.maketraceone(reshape(1:36,6,6))) == 1
    @test tr(Utils_states.maketraceone(rand(4,4))) ≈ 1
    @test tr(Utils_states.maketraceone(rand(4,4))) ≈ 1
    @test tr(Utils_states.maketraceone(rand(5,5))) ≈ 1
end

end
