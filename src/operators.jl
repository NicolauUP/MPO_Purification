

#Define σ+ and σ- operators
ITensors.op(::OpName"σ+", ::SiteType"Qubit") = [0 1; 0 0]
ITensors.op(::OpName"σ-", ::SiteType"Qubit") = [0 0; 1 0]

#Define the Identity MPO

function Identity_MPO(Sites::Vector{<:Index})

    NSites = length(Sites)

    Id_op_sum = OpSum()

    for i in 1:NSites
        Id_op_sum += 1, "Id", i
    end
    Id_MPO = MPO(Id_op_sum, Sites)
    return Id_MPO / NSites
end

#Define the Translation of 1 site MPO

function translation_MPO(sites::Vector{<:Index})

    NSites = length(sites)

    T_R_op_sum = OpSum()
    T_L_op_sum = OpSum()

    for l in 1:NSites  #Sum over each site

        opsum_R_Temp = OpSum()
        opsum_L_Temp = OpSum()

        
        opsum_R_Temp += 1, "σ+", l
        opsum_L_Temp += 1, "σ-", l

        for m in l+1:NSites
           opsum_R_Temp *= 1, "σ-", m 
           opsum_L_Temp *= 1, "σ+", m
        end

        for i in 1:l-1
            opsum_R_Temp += 1, "Id", i # Identity
            opsum_L_Temp += 1, "Id", i # Identity
        end
        
        T_R_op_sum += opsum_R_Temp
        T_L_op_sum += opsum_L_Temp
    end

    T_R_MPO = MPO(T_R_op_sum, sites) #Right translation MPO
    T_L_MPO = MPO(T_L_op_sum, sites) #Left translation MPO  

    return T_R_MPO, T_L_MPO
end



function Build_H_MF_1D(T::Function, eltype::Type{<:Number}, sites::Vector{<:Index}, ϵ::Float64,χ::Int64)

    _, HopMPO, _ = QuanticsTCI(T, eltype, sites, ϵ) #MPO for the hopping function 

    T_R, T_L = translation_MPO(sites) #Translation MPOs

    H_T_R = apply(HopMPO, T_R; cutoff = ϵ, maxdim = χ) #Apply HoppingMPO to the translation MPO
    H_T_L = apply(T_L, ITensors.dag(HopMPO); cutoff = ϵ, maxdim = χ) #Hermitian COnjugate

    return +(H_T_R, H_T_L; cutoff = ϵ, maxdim = χ)
end
