module C_constraints
    include(dirname(@__FILE__)*"\\Utils_cons.jl")
    include(dirname(dirname(@__FILE__))*"\\Moments.jl")
    using .Utils_cons
    using .Moments
    using LinearAlgebra
    const la = LinearAlgebra
    const uc = Utils_cons
    const mom = Moments

    export make_mon_expo_keys,
           make_PSD_con,
           make_ord4_con,
           make_loc_cons_S_inf,
           make_loc_cons_S₂,
           make_loc_cons_S₂₁,
           make_Gᴿ_con

    make_mon_expo_keys(d,t::Int) = mom.make_mon_expo(d,t*2)

    """L([xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜᵀ) ⪰ 0"""
    function make_PSD_con(d,t,Lx;noBlock=false)
        MMCoefᴿ,MMexᴿ = mom.get_ℂ_block_diag(d,t;noBlock=noBlock)
        return  Utils_cons.idx2var_dic(Lx,MMexᴿ,MMCoefᴿ)
    end

    """L(xx*⊗yy*) = ρ  ⟺ Lᴿ_ℝxx̄ᵀ⨂yȳᵀ = ℜe(ρ), Lᴿ_𝕀xx̄ᵀ⨂yȳᵀ = ℑm(ρ)  """
    function make_ord4_con(d,Lx)
        γδ_dict = mom.get_γδ_dict(d,(2,2))
        xx̄ᵀ⨂yȳᵀ= mom.make_xx̄ᵀ⨂yȳᵀ(d,γδ_dict)
        (cℝxx̄ᵀ⨂yȳᵀ,eℝxx̄ᵀ⨂yȳᵀ),(c𝕀xx̄ᵀ⨂yȳᵀ,e𝕀xx̄ᵀ⨂yȳᵀ) = mom.get_ℝℂ_coef_expo(d,xx̄ᵀ⨂yȳᵀ,γδ_dict)
        Lᴿ_ℝxx̄ᵀ⨂yȳᵀ = uc.idx2var_a(Lx,eℝxx̄ᵀ⨂yȳᵀ,cℝxx̄ᵀ⨂yȳᵀ)
        Lᴿ_𝕀xx̄ᵀ⨂yȳᵀ = uc.idx2var_a(Lx,e𝕀xx̄ᵀ⨂yȳᵀ,c𝕀xx̄ᵀ⨂yȳᵀ)
        return Lᴿ_ℝxx̄ᵀ⨂yȳᵀ,Lᴿ_𝕀xx̄ᵀ⨂yȳᵀ
    end

    """ Lᴿ(gᴿ⋅[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁ᵀ) ⪰ 0
        gᴿ = √(ρₘₐₓ) - ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ²) , i ∈ [d_1],
        or   √(ρₘₐₓ) - ((yᵣₑ)ⱼ² + (yᵢₘ)ⱼ²) , j ∈ [d_2] """
    make_loc_con(Lx,sqρ,eᵤ,eᵥ,B,C) = sqρ*uc.idx2var_a(Lx,B,C) - (uc.idx2var_a(Lx,B,C,[2*eᵤ]) + uc.idx2var_a(Lx,B,C,[2*eᵥ]))  # L((√ρₘₐₓ - ((xᵣₑ)ᵢ²+(xᵢₘ)ᵢ²))⋅ηₜ₋₁)
    function make_loc_cons_S_inf(ρ,d,t,Lx;noBlock=false)
        d₁,d₂ = d ; n = sum(2 .*d) ; sqρ   = sqrt(maximum(norm.(ρ)))
        MMCoefᴿ,MMexᴿ = mom.get_ℂ_block_diag(d,t .- 1;noBlock=noBlock)
        loc_con = Dict()
        for b in keys(MMCoefᴿ)
            for k in 1:d₁ # Constraint:  Lᴿ( (√ρₘₐₓ-((xᵣₑ)ᵢ²-(xᵢₘ)ᵢ²))⋅ηₜ₋₁)) ⪰ 0 for k ∈ [d₁]
                loc_con[(b,"x²_$k")] = make_loc_con(Lx,sqρ,uc.eᵢ(n,k),uc.eᵢ(n,k+d₁),MMexᴿ[b],MMCoefᴿ[b])
            end
            for k in 1:d₂ # Constraint:  Lᴿ( (√ρₘₐₓ-((yᵣₑ)ᵢ²+(yᵢₘ)ᵢ²))⋅ηₜ₋₁) ⪰ 0 for k ∈ [d₂]
                loc_con[(b,"y²_$k")] = make_loc_con(Lx,sqρ,uc.eᵢ(n,k+2*d₁),uc.eᵢ(n,k+2*d₁+d₂),MMexᴿ[b],MMCoefᴿ[b])
            end
        end
        return loc_con
    end

    """Lᴿ(gᴿ⋅[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁ᵀ) ⪰ 0
        gᴿ   = √Tr(ρ) - ∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ²),
                ∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ²) - ∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ²)"""
    tmp(Lx,ind,coef,n,s,l) = sum([uc.idx2var_a(Lx,ind,coef,[2* uc.eᵢ(n,k)]) for k in s:l]) #### SOmething is wrong
    function make_loc_cons_S₂(ρ,d,t,Lx;noBlock=false)
        d₁,d₂ = d ; n = sum(2 .*d) ; stρ = sqrt(real(tr(ρ)))
        MMCoefᴿ,MMexᴿ = mom.get_ℂ_block_diag(d, t.- 1,noBlock=noBlock)

        loc_con = Dict(); g₂ = Dict()
        for b in keys(MMCoefᴿ)
            e_η = MMexᴿ[b] ; c_η = MMCoefᴿ[b]
            xRterm = tmp(Lx,e_η,c_η,n,1,d₁)    + tmp(Lx,e_η,c_η,n,d₁+1,2*d₁)       #  Lᴿ(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            yRterm = tmp(Lx,e_η,c_η,n,2*d₁+1,2*d₁+d₂) + tmp(Lx,e_η,c_η,n,2*d₁+d₂+1,2*d₁+2*d₂)  #  Lᴿ(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            loc_con[b,"x"] = stρ*uc.idx2var_a(Lx,e_η,c_η) - xRterm #  √Tr(ρ) ⋅ Lᴿ(ηₜ₋₁) - Lᴿ(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            #loc_con[b,"y"] = stρ*uc.idx2var_a(Lx,MMexᴿ[b],MMCoefᴿ[b]) - yRterm #  √Tr(ρ) ⋅ Lᴿ(ηₜ₋₁) - Lᴿ(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            g₂[b] = xRterm - yRterm #
        end
        return loc_con, g₂
    end

    """Lᴿ(gᴿ⋅[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁ᵀ) ⪰ 0
        gᴿ = Tr(ρ) - ∑ᵈᵢ(uₓᵢ² + vₓᵢ²),
                ∑ᵈᵢ(u_yⱼ² + v_yⱼ²)  - 1"""
    function make_loc_cons_S₂₁(ρ,d,t,Lx;noBlock=false)
        d₁,d₂ = d ; n = sum(2 .*d) ; tr_ρ  = real(tr(ρ))
        MMCoefᴿ,MMexᴿ = mom.get_ℂ_block_diag(d, t.- 1,noBlock=noBlock)
        loc_con    = Dict() ; loc_con_eq = Dict()
        for b in keys(MMCoefᴿ)
            e_η = MMexᴿ[b] ; c_η = MMCoefᴿ[b]
            xRterm = tmp(Lx,e_η,c_η,n,1,d₁)    + tmp(Lx,e_η,c_η,n,d₁+1,2*d₁) #  Lᴿ(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            yRterm = tmp(Lx,e_η,c_η,n,2*d₁+1,2*d₁+d₂) + tmp(Lx,e_η,c_η,n,2*d₁+d₂+1,2*d₁+2*d₂) #  Lᴿ(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            loc_con[b]    = tr_ρ*uc.idx2var_a(Lx,MMexᴿ[b],MMCoefᴿ[b]) - xRterm #  Tr(ρ)⋅L(ηₜ₋₁) - L(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ ) ⪰ 0
            loc_con_eq[b] = uc.idx2var_a(Lx, MMexᴿ[b],MMCoefᴿ[b]) - yRterm #  1⋅L(ηₜ₋₁) -  L(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ ) = 0
        end
        return loc_con, loc_con_eq
    end

    """ρ⊗L(η) - L( (xx*⊗yy*) ⊗ η )⪰ 0 ∀ η ∈ blocks of ([x,̄x,y,̄y]ₜ₋₂[x,̄x,y,̄y]*ₜ₋₂)
     L( ρℝ ⊗ η)  - L( Gℝ[k] ⊗ η)⪰ 0 ∀ η ∈ blocks of [uₓ,vₓ,u_y,v_y]ₜ₋₂[uₓ,vₓ,u_y,v_y]^Tₜ₋₂  """
    function make_Gᴿ_con(ρ,d,t,Lx;noBlock=false)
        ρᴿ = [real(ρ) -imag(ρ); imag(ρ) real(ρ)]
        γδ_dict = mom.get_γδ_dict(d,(2,2))
        xx̄ᵀ⨂yȳᵀ = Moments.make_xx̄ᵀ⨂yȳᵀ(d,γδ_dict)
        c_xx̄ᵀ⨂yȳᵀᴿ,e_xx̄ᵀ⨂yȳᵀᴿ = mom.get_coef_expoᴿ(d,xx̄ᵀ⨂yȳᵀ,γδ_dict)
        MMᴿcoef,MMᴿexᴿ = mom.get_ℂ_block_diag(d, t.- 2,noBlock=noBlock)

        spec_sum(𝐴,𝐵) = [a + b  for a in 𝐴 for b in 𝐵]
        spec_mult(𝐶,𝐷) = [c * d  for c in 𝐶 for d in 𝐷]
        LGℝη = Dict()
        for b in keys(MMᴿcoef)
            Ex_temp   = mom.var_kron_C(e_xx̄ᵀ⨂yȳᵀᴿ,MMᴿexᴿ[b],spec_sum)
            Coef_temp = mom.var_kron_C(c_xx̄ᵀ⨂yȳᵀᴿ,MMᴿcoef[b],spec_mult)

            TEMP1 = uc.idx2var_a(Lx,Ex_temp,Coef_temp)  # L(Gᴿ[k]⋅ηₜ₋₂)
            TEMP2 = la.kron(ρᴿ,uc.idx2var_a(Lx,MMᴿexᴿ[b],MMᴿcoef[b])) # L(ρᴿ ⊗ ηₜ₋₂)
            LGℝη[b]   = TEMP2 - TEMP1
        end
        return  LGℝη
    end
end



#
# for b in keys(MMᴿcoef)
#     TEMP1 = sum([tmp2(MMexᴿ[b],MMCoefᴿ[b],k) for k in 1:8]) # ∑ₖL(Gᴿ[k]⋅ηₜ₋₂)
#     TEMP2 = la.kron(ρᴿ,uc.idx2var_arr(Lx,MMexᴿ[b],MMCoefᴿ[b])) # L(ρᴿ ⊗ ηₜ₋₂)
#     LGℝη[b]   = TEMP2 - TEMP1 # L(ρᴿ ⊗ ηₜ₋₂)  - ∑ₖL(Gᴿ[k] ⊗ ηₜ₋₂)
# end


        # LGℝη         = Dict()
        # tmp2(B,C,k) = la.kron(sm[k],ones(D .*size(B))) .*
        #     uc.idx2var_arr(
        #         Lx,
        #         mom.var_kron_C(Gᴿ[k],B),
        #         mom.var_kron_C(fill(0.0,size(Gᴿ[k])...),C))


# make_ord4_con(d,Lx) = uc.idx2varxx̄ᵀtyȳᵀ(Lx,mom.make_xx̄ᵀ⨂yȳᵀ(d))
