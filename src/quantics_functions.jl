"""

Functions for the quantics (Q) tensor cross interpolation (TCI).

"""



function Convert_To_Binary(n::Int64, NSites::Int64)
    bits = last(bitstring(n), NSites)
    VectorString = []
    for i in bits
        push!(VectorString, string(i))
    end
    return VectorString
end


function BasisStateMPS(n::Int64, sites::Vector{<:Index})
    # Create a random MPS with the given sites and binary vector
    return random_mps(sites, Convert_To_Binary(n, length(sites)))
end

"""
     MatrixChecker(mpo::MPO, sites::Vector{<:Index}, i::Int64, j::Int64)
Check the matrix element of an MPO between two basis states specified by indices i and j
"""
function MatrixChecker(mpo::MPO, sites::Vector{<:Index}, i::Int64, j::Int64)
    # Convert indices to MPS
    Ψ_i = BasisStateMPS(i, sites)
    Ψ_j = BasisStateMPS(j, sites)
    
    # Compute the inner product
    return inner(Ψ_i', mpo, Ψ_j)
end


"""
        Quantics_TCI(f::Function, eltype::Type{<:Number}, sites::Vector{<:Index}, ϵ::Float64)
Obtain the quantics tensor cross interpolation of a function f over the basis states defined by the sites.
Returns the QTT, the corresponding MPO, and the MPS.
"""
function Quantics_TCI(f::Function, eltype::Type{<:Number}, sites::Vector{<:Index}, ϵ::Float64)
   # Get the number of sites
   NSites = length(sites)

   #Define the support/domain 
   
   XVals = collect(range(0, 2^NSites-1 ;length=2^NSites))

   #Obtain quantics tensor train, (ignoring the error and the number of iterations)
   QTT, _, _ = QuanticsTCI.quanticscrossinterpolate(eltype, f, XVals; tolerance = ϵ)

   #Convert the quantics tensor train to a tensor train
    TT = TCI.tensortrain(QTT.tci)


    mps = MPS(TT; sites)
    mpo = outer(mps', mps)
    for i in 1:NSites
        mpo.data[i] = Quantics._asdiagonal(mps.data[i], sites[i])
    end
    return QTT, mpo, mps

end



"""
        Get_ChargeDensity(ρ::MPO, sites::Vector{<:Index}, i::Int64,ϵ::Float64)
Get the charge density at site i from the density MPO ρ using quantics tensor cross interpolation.
Returns the QTT, the corresponding MPO, and the MPS.
"""
function Get_ChargeDensity(ρ::MPO, sites::Vector{<:Index}, i::Int64,ϵ::Float64)
    NSites = length(sites)

    f(α) = inner(BasisStateMPS(mod(Int(α), 2^NSites), sites)', ρ, BasisStateMPS(mod(Int(α), 2^NSites), sites) )
    # Get the quantics tensor cross interpolation
    QTT, mpo, mps = QuanticsTCI(f, eltype(ρ[1]), sites, ϵ)
    return QTT, mpo, mps
end


"""
        Get_HartreeTerm(ρ::MPO, sites::Vector{<:Index},ϵ::Float64)
Get the Hartree term from the density MPO ρ using quantics tensor cross interpolation.
Returns the QTT, the corresponding MPO, and the MPS.   
Specific for nearest-neighbor interactions.
"""
function Get_HartreeTerm(ρ::MPO, sites::Vector{<:Index},ϵ::Float64)
    NSites = length(sites)
    #Verify if this is the correct function or it is the other way around.  
    factor_left(α) = Int( (α - 1 ) >= 0) #PBC to the left
    factor_right(α) = Int( (α + 1 ) < 2^NSites) #PBC to the right 
    f(α) = factor_left(α) *  inner(BasisStateMPS(mod(Int(α-1), 2^NSites), sites)', ρ, BasisStateMPS(mod(Int(α-1), 2^NSites), sites) ) 
           +factor_right(α) *  inner(BasisStateMPS(mod(Int(α+1), 2^NSites), sites)', ρ, BasisStateMPS(mod(Int(α+1), 2^NSites), sites) ) 

    # Get the quantics tensor cross interpolation
    QTT, mpo, mps = Quantics_TCI(f, eltype(ρ[1]), sites, ϵ)
    return QTT, mpo, mps
end 

"""
        Get_FockTerm(ρ::MPO, sites::Vector{<:Index},ϵ::Float64)
Get the Fock term from the density MPO ρ using quantics tensor cross interpolation.
Returns the QTT, the corresponding MPO, and the MPS.   
Specific for nearest-neighbor interactions.
"""

function Get_FockTerm(ρ::MPO, sites::Vector{<:Index},ϵ::Float64)
    NSites = length(sites)
    #Verify if this is the correct function or it is the other way around.  
    f(α) = inner(BasisStateMPS(mod(Int(α+1), 2^NSites), sites)', ρ, BasisStateMPS(mod(Int(α), 2^NSites), sites) )
    # Get the quantics tensor cross interpolation
    QTT, mpo, mps = Quantics_TCI(f, eltype(ρ[1]), sites, ϵ)
    return QTT, mpo, mps
end




