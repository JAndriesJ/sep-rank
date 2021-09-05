module Examples_sep
using LinearAlgebra
using Random

srcDir = dirname(@__FILE__)*"\\"
include(srcDir *"Utils_states.jl")
using .Utils_states ; const us = Utils_states

export get_sep_examples,
       gen_ρ_RAND

## Separable States
"""
Separable state examples:
"""
function get_sep_example()
    ρ = Dict()

##  T1CD12 i-ii
    ρ["S1a"] = get_ρ_T1CD12(1)
    ρ["S1b"] = get_ρ_T1CD12(2)
    # ρ["T1CD12iii"]  = get_ρ_T1CD12(3)

## T2CD12 i,ii, v
    ρ["S2a"]   = get_ρ_T2CD12(1)
    ρ["S2b"]   = get_ρ_T2CD12(2)
    # ρ["T2CD12iii"]  = get_ρ_T2CD12(3)
    # ρ["T2CD12iv"]  = get_ρ_T2CD12(4)
    # ρ["T2CD12v"]    = get_ρ_T2CD12(5)
## Ex26CD12
    ρ["S3"]   = get_ρ_Ex26CD12()

## Ex3-8DNY20
    ρ["S4"] = get_ρ_DNY20()

## DTT00a
    ρ["S5a"] = get_ρ_DTT00a(3)
    ρ["S5b"] = get_ρ_DTT00b()
## Random batch
## ℝeal
    ρ["RRai"]   = gen_ρ_RAND((3,3), 5)
    ρ["RRaii"]  = gen_ρ_RAND((3,3), 10)
    ρ["RRaiii"] = gen_ρ_RAND((3,3), 15)
    ρ["RRaiv"]  = gen_ρ_RAND((3,3), 20)
    ρ["RRav"]   = gen_ρ_RAND((3,3), 25)

    ρ["RRbi"]   = gen_ρ_RAND((2,4), 5)
    ρ["RRbii"]  = gen_ρ_RAND((2,4), 10)
    ρ["RRbiii"] = gen_ρ_RAND((2,4), 15)
    ρ["RRbiv"]  = gen_ρ_RAND((2,4), 20)
    ρ["RRbv"]   = gen_ρ_RAND((2,4), 25)

    ρ["RRci"]   = gen_ρ_RAND((2,5), 5)
    ρ["RRcii"]  = gen_ρ_RAND((2,5), 10)
    ρ["RRciii"] = gen_ρ_RAND((2,5), 15)
    ρ["RRciv"]  = gen_ρ_RAND((2,5), 20)
    ρ["RRcv"]   = gen_ρ_RAND((2,5), 25)

    ρ["RRdv"]   = gen_ρ_RAND((4,4), 25)

## ℂomplex
    ρ["RCai"]   = gen_ρ_RAND((3,3), 5, true)
    ρ["RCaii"]  = gen_ρ_RAND((3,3), 10, true)
    ρ["RCaiii"] = gen_ρ_RAND((3,3), 15, true)
    ρ["RCaiv"]  = gen_ρ_RAND((3,3), 20, true)
    ρ["RCav"]   = gen_ρ_RAND((3,3), 25, true)

    ρ["RCbi"]   = gen_ρ_RAND((2,4), 5, true)
    ρ["RCbii"]  = gen_ρ_RAND((2,4), 10, true)
    ρ["RCbiii"] = gen_ρ_RAND((2,4), 15, true)
    ρ["RCbiv"]  = gen_ρ_RAND((2,4), 20, true)
    ρ["RCbv"]   = gen_ρ_RAND((2,4), 25, true)

    ρ["RCci"]   = gen_ρ_RAND((2,5), 5, true)
    ρ["RCcii"]  = gen_ρ_RAND((2,5), 10, true)
    ρ["RCciii"] = gen_ρ_RAND((2,5), 15, true)
    ρ["RCciv"]  = gen_ρ_RAND((2,5), 20, true)
    ρ["RCcv"]   = gen_ρ_RAND((2,5), 25, true)


    ρ["RCd"]   = gen_ρ_RAND((3,4), 25, true)
    ρ["RCe"]   = gen_ρ_RAND((2,6), 25, true)
 ## Real
 #    ρ["RANDa1"]  = gen_ρ_RAND(3, 5)
 #    ρ["RANDa2"]  = gen_ρ_RAND(3, 10)
 #    ρ["RANDa3"]  = gen_ρ_RAND(3, 15)
 #    ρ["RANDa4"]  = gen_ρ_RAND(3, 20)
 #    ρ["RANDa5"]  = gen_ρ_RAND(3, 25)
 #
 #    ρ["RANDb1"]  = gen_ρ_RAND(4, 5)
 #    ρ["RANDb2"]  = gen_ρ_RAND(4, 10)
 #    ρ["RANDb3"]  = gen_ρ_RAND(4, 15)
 #    ρ["RANDb4"]  = gen_ρ_RAND(4, 20)
 #    ρ["RANDb5"]  = gen_ρ_RAND(4, 25)
 # ## Complex
 #    ρ["RANDℂa1"]  = gen_ρ_RAND(3,  5,true)
 #    ρ["RANDℂa2"]  = gen_ρ_RAND(3, 10,true)
 #    ρ["RANDℂa3"]  = gen_ρ_RAND(3, 15,true)
 #    ρ["RANDℂa4"]  = gen_ρ_RAND(3, 20,true)
 #    ρ["RANDℂa5"]  = gen_ρ_RAND(3, 25,true)
 #
 #    ρ["RANDℂb1"]  = gen_ρ_RAND(4,  5,true)
 #    ρ["RANDℂb2"]  = gen_ρ_RAND(4, 10,true)
 #    ρ["RANDℂb3"]  = gen_ρ_RAND(4, 15,true)
 #    ρ["RANDℂb4"]  = gen_ρ_RAND(4, 20,true)
 #    ρ["RANDℂb5"]  = gen_ρ_RAND(4, 25,true)


    for key in keys(ρ)
        ρ[key] = Utils_states.maketraceone(ρ[key])
    end
    return ρ
end

## Examples

""" http://arxiv.org/abs/1210.0111v2"""
function get_ρ_T1CD12(i)
    ρ = Dict()
    # |00><00| + |11><11|
    ρ[1] =  psep(1,1,2) + psep(2,2,2)
    # |00><00| + |11><11| + (|0> + |1>)⊗(|0> + |1>)(<0| + <1|)⊗(<0| + <1|)
    ρ[2] = psep(1,1,2) + psep(2,2,2) + sq(ψ(2))
    return ρ[i]
end

""" http://arxiv.org/abs/1210.0111v2"""
function get_ρ_T2CD12(i)
    ρ = Dict()
    # |00><00| + |11><11| + |12><12|
    ρ[1] = psep(1,1,3) + psep(2,2,3) + psep(2,3,3)
    # |00><00| + |01><01| + |11><11| + |12><12|
    ρ[2] = psep(1,1,3) + psep(1,2,3) + psep(2,2,3) + psep(2,3,3)
    # |00><00| + |01><01| + |02><02| + |11><11| + |12><12|
    ρ[5] = psep(1,1,3) + psep(1,2,3) +  psep(1,3,3) + psep(2,2,3) + psep(2,3,3)
    return ρ[i]
end

""" http://arxiv.org/abs/1210.0111v2"""
function get_ρ_Ex26CD12()
    ρ = [4 0 0 0 0 0
         0 4 2 0 0 2
         0 2 2 1 -1 0
         0 0 1 2 1 -1
         0 0 -1 1 5 1
         0 2 0 -1 1 2]
    return ρ
end

"""https://arxiv.org/abs/2011.08132v1"""
function get_ρ_DNY20()
    H_flat = [i₁*j₁ + i₂*j₂ for i₁ in 1:3 for i₂ in 1:3 for j₁ in 1:3 for j₂ in 1:3 ]
    return  reshape(H_flat,9,9)
end

"""https://arxiv.org/abs/quant-ph/9904005v2 page 7 eq. (11)"""
function get_ρ_DTT00a(n)
    #@assert 0 ≤ f ≤ 1/n  # to be separable
    f=1/(n)
    Ψn = (1/sqrt(n))*sum([kron(us.eᵢ(n,i),us.eᵢ(n,i)) for i in 1:n])
    ρᵥᵥf = f* us.sq(Ψn) + (1-f)/(n^2)*I(n^2)
    ρf =  us.take_Pᵀ(ρᵥᵥf,1,(n,n))
end

"""https://arxiv.org/abs/quant-ph/9904005v2 page 9 eq. (13)"""
function get_ρ_DTT00b()
    N = 2/(sqrt(5 + sqrt(5)))
    h = 0.5*(sqrt(1 + sqrt(5)))
    Ñ = sqrt(2/sqrt(5))

    v(i) = N .* [cos(2*π*i/5),sin(2*π*i/5),h]
    w(j) = Ñ .* sqrt(cos(π/5))*[cos(2*π*j/5),sin(2*π*j/5),cos(4*π*j/5),sin(4*π*j/5)]

    (1/140)*(I(12) - sum([kron(us.sq(v(i)),us.sq(w(i))) for i in 1:4 ]))
end

"""generates a matrix of the form ∑ʳaᵀa⊗bᵀb, with a,b ∈ uniform random entry-wise.  d ∈ {2,3,4} , r ∈ [9]"""
function gen_ρ_RAND(d, r::Integer, isℂ = false)
    Random.seed!(343) ; d₁,d₂ = d ; ρ = zeros(d₁*d₂,d₁*d₂)
    if isℂ
        a = [randn(1,d₁) + im*randn(1,d₁) for i in 1:r ]
        b = [randn(1,d₂) + im*randn(1,d₂) for i in 1:r ]
    else
        a = [randn(1,d₁) for i in 1:r ]
        b = [randn(1,d₂) for i in 1:r ]
    end
    for i in 1:r
        A = transpose(a[i])*conj(a[i])
        B = transpose(b[i])*conj(b[i])
        ρ = kron(A,B) + ρ
    end
    return ρ
end



# """ http://arxiv.org/abs/quant-ph/9605038v2 """
# function get_ρ_HHH1(p  = 0.5)
#     # Random.seed!(343)
#     a  = 0.5; b = 0.8; # arbitary possitive numbers
#     ψ₁ = a*Utils_states.ψ(2,2,2) + b*Utils_states.ψ(1,1,2)
#     ψ₂ = a*Utils_states.ψ(2,1,2) + b*Utils_states.ψ(1,2,2)
#     ρ_temp = p*Utils_states.sq(ψ₁) + (1 - p)*Utils_states.sq(ψ₂)
#     return ρ_temp
# end

end
