module Moments

export  eᵢ,
        var_kron,
        var_kron_C,
        make_mon_expo,
        make_xxᵀ⨂yyᵀ,
        get_ℝ_block_diag,
        make_xx̄ᵀ⨂yȳᵀ,
        get_γδ_dict,
        get_ℝℂ_coef_expo,
        get_coef_expoᴿ,
        get_ℂ_block_diag

## Utils
"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [Int(j==i) for j in 1:n]

""" input: A,B (arrays of integer tupples) output: A ⊗ B """
function var_kron(A,B)
   C_temp = [ a + b for  b in B, a ∈ A]
   rcs =  [(i,j) for i in 1:size(A)[1] for j in 1:size(B)[1]]
   return [C_temp[ij[2],kl[2],ij[1],kl[1]] for ij in rcs , kl in rcs ]
end

function var_kron_C(A,B,fun) # still a shitty fix
    C_temp = [isempty(𝐴) || isempty(𝐵) ? [] : fun(𝐴,𝐵)  for 𝐵 in B, 𝐴 in A]
    rcs =  [(i,j) for i in 1:size(A)[1] for j in 1:size(B)[1]]
    return  [C_temp[ij[2],kl[2],ij[1],kl[1]] for ij in rcs , kl in rcs ] # is the form correct ?
end

## Moments
"""[x]≦ₜ or [x]₌ₜ"""
function make_mon_expo(n::Int,t::Int; isle::Bool = true)
    @assert typeof(t) == Int64
    t == 0 ? (return [eᵢ(n,0)]) : 0
    tmp = make_mon_expo(n,t-1;isle=isle)
    M_vec = reshape([m + eᵢ(n,i) for i ∈ 1:n, m ∈ tmp],:,1)
    return unique(isle ? vcat(tmp,M_vec) : M_vec)
end

"""[x]≦ₜ[x]ᵀ≦ₜ or [x]₌ₜ[x]ᵀ₌ₜ"""
function make_mon_expo(n::Int,t::Tuple{Int,Int}; isle::Bool = true)
    M_vec1      = make_mon_expo(n,t[1]; isle=isle)
    M_vec2      = make_mon_expo(n,t[2]; isle=isle)
    M_mat_flat  = [mi+mj for mi in M_vec1, mj in M_vec2]
end

## Complex moments
"""[x,̄x]≦ₜ or [x,̄x]ᵀ₌ₜ"""
make_mon_expo(d::Tuple{Int,Int},t::Int; isle::Bool = true) = make_mon_expo(sum(2 .* d),t; isle=isle)

"""[x,̄x]≦ₜ[x,̄x]*≦ₜ or [x,̄x]₌ₜ[x,̄x]*₌ₜ but is actually [xᵣₑ,xᵢₘ]≦ₜ[xᵣₑ,xᵢₘ]ᵀ≦ₜ or [xᵣₑ,xᵢₘ]₌ₜ[xᵣₑ,xᵢₘ]ᵀ₌ₜ """
make_mon_expo(d::Tuple{Int,Int},t::Tuple{Int,Int}; isle::Bool = true) = make_mon_expo(sum(2 .* d),t; isle=isle)

## Real
"""xxᵀ⊗yyᵀ"""
function make_xxᵀ⨂yyᵀ(d)
    n = sum(d); p_ex = make_mon_expo(n,1;isle = false)
    x, y = p_ex[1:d[1]], p_ex[d[1]+1:n]
    return var_kron(x .+ reshape(x,1,:), y .+ reshape(y,1,:))
end

"""Returns diagonal block of the momoment matrix with partition:
    xᵃ⁺ᶜyᵇ⁺ᵈ where |a+c|,|b+d| ∈ 2N """
function get_ℝ_block_diag(d::Tuple{Int,Int},t::Tuple{Int,Int})
    mom_mat = Moments.make_mon_expo(sum(d),t)
    mom_vec = mom_mat[1,:]
    halfsum(arr) = [sum(arr[1:d[1]]),sum(arr[(d[1]+1):end])]
    tf(arr)      = map(x -> (isodd(x[1]),isodd(x[2])), halfsum.(arr))
    tf2(arr)     = map(p -> p[1], findall(arr))

    isevev(pair) = pair == (false,false)
    isodod(pair) = pair == (true,true)
    isevod(pair) = pair == (false,true)
    isodev(pair) = pair == (true,false)

    ee = tf2(isevev.(tf(mom_vec)))
    oo = tf2(isodod.(tf(mom_vec)))
    eo = tf2(isevod.(tf(mom_vec)))
    oe = tf2(isodev.(tf(mom_vec)))

    ℝ_block_diag = Dict("ee" => mom_mat[ee,ee],
                        "oo" => mom_mat[oo,oo],
                        "eo" => mom_mat[eo,eo],
                        "oe" => mom_mat[oe,oe])

    for key in ["ee","oo","eo","oe"]
         isempty(ℝ_block_diag[key]) ? delete!(ℝ_block_diag,key) : continue
    end
    return ℝ_block_diag
end

## Complex
""" ℜe(xx*⊗yy*), ℑm(xx*⊗yy*)
xx*⊗yy* = (xᵣₑ + i xᵢₘ)(xᵣₑ - i xᵢₘ)ᵀ ⊗ (yᵣₑ + i yᵢₘ)(yᵣₑ - i yᵢₘ)ᵀ
= (xᵣₑxᵣₑᵀ + xᵢₘxᵢₘᵀ + i (xᵢₘxᵣₑᵀ - xᵣₑxᵢₘᵀ)) ⊗ (yᵣₑyᵣₑᵀ + yᵢₘyᵢₘᵀ + i (yᵢₘyᵣₑᵀ - yᵣₑyᵢₘᵀ))
⟹ cℝxx̄ᵀ⨂yȳᵀ ∘ eℝxx̄ᵀ⨂yȳᵀ := ℜe(xx*⊗yy*) = (xᵣₑxᵣₑᵀ + xᵢₘxᵢₘᵀ ) ⊗ (yᵣₑyᵣₑᵀ + yᵢₘyᵢₘᵀ ) - (xᵢₘxᵣₑᵀ - xᵣₑxᵢₘᵀ) ⊗ (yᵢₘyᵣₑᵀ - yᵣₑyᵢₘᵀ)
⟹ c𝕀xx̄ᵀ⨂yȳᵀ ∘ e𝕀xx̄ᵀ⨂yȳᵀ := ℑm(xx*⊗yy*) = (xᵣₑxᵣₑᵀ + xᵢₘxᵢₘᵀ) ⊗ (yᵢₘyᵣₑᵀ - yᵣₑyᵢₘᵀ) + (xᵢₘxᵣₑᵀ - xᵣₑxᵢₘᵀ)⊗(yᵣₑyᵣₑᵀ + yᵢₘyᵢₘᵀ) """

function make_xx̄ᵀ⨂yȳᵀ(d,γδ_dict)
    p_ex = make_mon_expo(sum(2 .*d),1;isle = false)
    x, x̄, y, ȳ = split_expo(p_ex,d...)
    xx̄ᵀ⨂yȳᵀ = var_kron(x .+ reshape(x̄,1,:), y .+ reshape(ȳ,1,:))
end


## Block diagonalization
# The complex blocks
"""Splits a vector of exponents into the register components ααᶥββᶥ →  (α,αᶥ,β,βᶥ)"""
split_expo(ααᶥββᶥ,d₁,d₂) = (ααᶥββᶥ[1:d₁],ααᶥββᶥ[d₁+1:2*d₁],
                            ααᶥββᶥ[1+2*d₁:2*d₁+d₂],ααᶥββᶥ[1+2*d₁+d₂:end])

"""xᵅ̄xᵃyᵝ̄yᵇ -> xᵃx̄ᵅyᵇ̄yᵝ̄"""
function conj_expo(ααᶥββᶥ,d₁,d₂)
    α,αᶥ,β,βᶥ = split_expo(ααᶥββᶥ,d₁,d₂) ; return vcat(αᶥ,α,βᶥ,β)
end

"""Gets the principal block of the moment matrix after block diagonalization"""
function get_xx̄yȳMM_blocks(d,t)
    Iᵗ = make_mon_expo(d,t[1])
    Iᵗ_temp = map(x ->  sum.(split_expo(x,d...)), Iᵗ)
    r_list  = map(v -> v[1]-v[2],Iᵗ_temp) ; s_list  = map(v -> v[3]-v[4],Iᵗ_temp)
    P = Dict()
    for r ∈ -t[1]:t[1], s ∈ -t[1]:t[1]
        T = Iᵗ[(r_list .== r) .& (s_list .== s)]
        isempty(T) ? nothing : P[r,s] = T
    end
    build_block(k) = P[k] .+ reshape(map(x->conj_expo(x,d...),P[k]),1,:)
    Dict(zip(keys(P), [build_block(k) for k in keys(P)]))
end

## The real analogue
"""Returns all  (γ,γᶥ,ζ,ζᶥ),(δ,δᶥ,η,ηᶥ) ∈ (ℕᵈ)⁴ s.t.
(γ,γᶥ,ζ,ζᶥ)+(δ,δᶥ,η,ηᶥ)=(α,αᶥ,β,βᶥ)  ∀  ααᶥββᶥ ∈ Moment matrixₜ
"""
function get_γδ(m)
    # nz_coord = findall(.!(m .== 0)) ; n = length(m)
    # γ = [m] ; δ = [eᵢ(n,0)]
    # for k in nz_coord
    #     n_count = m[k]
    #     while n_count > 0
    #         mtemp = eᵢ(n,k)
    #         push!(γ, γ[end] - mtemp)
    #         push!(δ, δ[end] + mtemp)
    #         n_count -= 1
    #     end
    # end
    # return unique(vcat(γ,δ))
    n = length(m)
    nz_coord = findall(.!(m .== 0))
    isempty(nz_coord) ? (return [m]) : nothing
    c = pop!(nz_coord)
    Nar_old = [ j*eᵢ(n,c) for j in 0:m[c]]
    global Nar_old
    while  !isempty(nz_coord)
        c = pop!(nz_coord)
        Nar_new = [ j*eᵢ(n,c) for j in 0:m[c]]
        Nar_old = [ n + m for n in Nar_old for m in Nar_new ]
    end
    return Nar_old
end

function get_γδ_dict(d,t)
    MM_vec = make_mon_expo(d,t[1]*2)
    return Dict(zip(MM_vec,map(m ->([m] .- get_γδ(m),get_γδ(m)) ,MM_vec)))
end

# function get_γδ_dict(d,t)
#     M = Moments.make_mon_expo(d,2 .* t)
#     M_vec = M[:,1]
#     γδ_dict = Dict()
#     for k in Moments.make_mon_expo(d,2*t[1])
#         fd = hcat(map(x -> [x[1],x[2]],findall(M .== [k]))...)
#         γδ_dict[k] = [M_vec[fd[1,:]], M_vec[fd[2,:]]]
#     end
#     return γδ_dict
# end








"""The coefficients function
-1^|δᶥ+ηᶥ|
i^|δ+δᶥ+η+ηᶥ|
α!αᶥ!β!βᶥ! / γ!γᶥ!ζ!ζᶥ!δ!δᶥ!η!ηᶥ!
"""
f_temp(X) = prod(map(x-> prod(factorial.(x)),X))
function get_coef(p1,p2,d)
    γ,γᶥ,ζ,ζᶥ = Moments.split_expo(p1,d...) ; δ,δᶥ,η,ηᶥ = Moments.split_expo(p2,d...)
    a = ((-1.0)^sum(sum.([δᶥ,ηᶥ])))*((1.0im)^sum(sum.([δ,δᶥ,η,ηᶥ])))*f_temp([[γ,γᶥ,ζ,ζᶥ]+[δ,δᶥ,η,ηᶥ]...])
    return a/f_temp([γ,γᶥ,ζ,ζᶥ,δ,δᶥ,η,ηᶥ])
end
get_coef(p1p2_arr,d) = map((x,y)->get_coef(x,y,d),p1p2_arr...)

""" get_expo(γγᶥδδᶥ,ζζᶥηηᶥ,d) = [γ+γᶥ,δ+δᶥ,η+ηᶥ,ζ+ζᶥ]"""
function get_expo(p1,p2,d)
    γ,γᶥ,ζ,ζᶥ = split_expo(p1,d...) ; δ,δᶥ,η,ηᶥ = split_expo(p2,d...)
    return vcat(γ+γᶥ,δ+δᶥ,ζ+ζᶥ,η+ηᶥ)
end
get_expo(p1p2_arr,d) = map((x,y)->get_expo(x,y,d),p1p2_arr...)

"""f"""
function get_coef_expo(d,ααᶥββᶥ::Array{Int64,1},γδ_dict)
    γδ = γδ_dict[ααᶥββᶥ]
    return get_coef(γδ,d),get_expo(γδ,d)
end

function get_ℝℂ_coef_expo(coef,expo)
    ℝmask = isreal.(coef)
    ℝcoef = real(coef[ℝmask])   ; ℝexpo = expo[ℝmask]
    ℂcoef = imag(coef[.!ℝmask]) ; ℂexpo = expo[.!ℝmask]
    return (ℝcoef,ℝexpo),(ℂcoef,ℂexpo)
end

function get_ℝℂ_coef_expo(d,B,γδ_dict)
    s₁,s₂ = size(B)
    MM = map(x-> get_coef_expo(d,x,γδ_dict),B)
    MMℝℂ = map(x-> get_ℝℂ_coef_expo(x...),MM)

    ℝcoef = [MMℝℂ[i,j][1][1] for i in 1:s₁, j in 1:s₂]
    ℝexpo = [MMℝℂ[i,j][1][2] for i in 1:s₁, j in 1:s₂]
    ℂcoef = [MMℝℂ[i,j][2][1] for i in 1:s₁, j in 1:s₂]
    ℂexpo = [MMℝℂ[i,j][2][2] for i in 1:s₁, j in 1:s₂]
    return (ℝcoef,ℝexpo),(ℂcoef,ℂexpo)
end

"""(ℝcoef,ℝexpo),(ℂcoef,ℂexpo) ->
MMCoefᴿ = (ℝcoef -ℂcoef)        MMexᴿ :=     (ℝexpo ℂexpo)
          (ℂcoef  ℝcoef)    and              (ℂexpo ℝexpo)
 """
function get_coef_expoᴿ(d,B,γδ_dict)
    (ℝcoef,ℝexpo),(ℂcoef,ℂexpo) = get_ℝℂ_coef_expo(d,B,γδ_dict)
    MMCoefᴿ = hcat(vcat(ℝcoef,ℂcoef),vcat((-1.0)*ℂcoef,ℝcoef))
    MMexᴿ   = hcat(vcat(ℝexpo,ℂexpo),vcat(ℂexpo,ℝexpo))
    return MMCoefᴿ,MMexᴿ
end

""""""
function get_ℂ_block_diag(d,t;noBlock = false)
    if noBlock
        MMexᴿ   = Dict("Default" => map(x->[x],make_mon_expo(d,t) .+ reshape(make_mon_expo(d,t),1,:) ))
        MMCoefᴿ = Dict("Default" => fill([1.0],size(MMexᴿ["Default"])...))
        return MMCoefᴿ,MMexᴿ
    end
    xx̄yȳMM_blocks = get_xx̄yȳMM_blocks(d,t)
    γδ_dict = get_γδ_dict(d,t) ## This should be moved outside so that it can be precomputed.

    MMexᴿ= Dict() ; MMCoefᴿ = Dict()
    for B in keys(xx̄yȳMM_blocks)
        MMCoefᴿ[B],MMexᴿ[B] = get_coef_expoᴿ(d,xx̄yȳMM_blocks[B],γδ_dict)
    end
    return MMCoefᴿ,MMexᴿ
end


end


# U(vec) = vec .+ reshape(vec,1,:)
# W(vec1,vec2) = vec1 .+ reshape(vec2,1,:)
# function make_xx̄ᵀ⨂yȳᵀ(d)
#     d₁,d₂ = d
#     n = sum(2 .* d)
#     pre_expo = make_mon_expo(n,1;isle =false)
#     xᵣₑ, xᵢₘ = pre_expo[1:d₁], pre_expo[d₁+1:2*d₁]
#     yᵣₑ, yᵢₘ = pre_expo[2*d₁+1:2*d₁+d₂], pre_expo[2*d₁+d₂+1:2*d₁+2*d₂]
#
#     xᵣₑxᵣₑᵀ, xᵢₘxᵢₘᵀ, yᵣₑyᵣₑᵀ, yᵢₘyᵢₘᵀ = U(xᵣₑ), U(xᵢₘ), U(yᵣₑ), U(yᵢₘ)
#     xᵢₘxᵣₑᵀ, xᵣₑxᵢₘᵀ, yᵢₘyᵣₑᵀ, yᵣₑyᵢₘᵀ = W(xᵢₘ,xᵣₑ), W(xᵣₑ,xᵢₘ), W(yᵢₘ,yᵣₑ), W(yᵣₑ,yᵢₘ)
#
#     Real_dict =  Dict(("+ 1")  => var_kron(xᵣₑxᵣₑᵀ,yᵣₑyᵣₑᵀ),
#                       ("+ 2")  => var_kron(xᵣₑxᵣₑᵀ,yᵢₘyᵢₘᵀ),
#                       ("+ 3")  => var_kron(xᵢₘxᵢₘᵀ,yᵣₑyᵣₑᵀ),
#                       ("+ 4")  => var_kron(xᵢₘxᵢₘᵀ,yᵢₘyᵢₘᵀ),
#                       ("- 5")  => var_kron(xᵢₘxᵣₑᵀ,yᵢₘyᵣₑᵀ),
#                       ("+ 6")  => var_kron(xᵢₘxᵣₑᵀ,yᵣₑyᵢₘᵀ),
#                       ("+ 7")  => var_kron(xᵣₑxᵢₘᵀ,yᵢₘyᵣₑᵀ),
#                       ("- 8")  => var_kron(xᵣₑxᵢₘᵀ,yᵣₑyᵢₘᵀ))
#
#     Imag_dict =  Dict(("+ 1")  => var_kron(xᵣₑxᵣₑᵀ,yᵢₘyᵣₑᵀ),
#                       ("- 2") => var_kron(xᵣₑxᵣₑᵀ,yᵣₑyᵢₘᵀ),
#                       ("+ 3")  => var_kron(xᵢₘxᵢₘᵀ,yᵢₘyᵣₑᵀ),
#                       ("- 4") => var_kron(xᵢₘxᵢₘᵀ,yᵣₑyᵢₘᵀ),
#                       ("+ 5")  => var_kron(xᵢₘxᵣₑᵀ,yᵣₑyᵣₑᵀ),
#                       ("+ 6")  => var_kron(xᵢₘxᵣₑᵀ,yᵢₘyᵢₘᵀ),
#                       ("- 7") => var_kron(xᵣₑxᵢₘᵀ,yᵣₑyᵣₑᵀ),
#                       ("- 8") => var_kron(xᵣₑxᵢₘᵀ,yᵢₘyᵢₘᵀ))
#
#     return Dict("real" => Real_dict,
#                 "imag" => Imag_dict)
#
#
# end
