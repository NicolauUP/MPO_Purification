begin 
using ITensors, ITensorMPS
using CairoMakie
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using LinearAlgebra
end

begin
include("operators.jl")
include("quantics_functions.jl")
end

function construct_rho_0(N::Int, H::MPO, H_max::Float64, H_min::Float64, Ne::Int, sites::Vector{<:Index})
    # Calculate mean energy
    μ = real(tr(H) / N)
    
    # Calculate the slope lambda to ensure eigenvalues stay in [0, 1]
    # This logic remains correct
    λ = minimum((Ne /(H_max - μ), (N - Ne)/(μ - H_min)))

    # Correct Formula: ρ0 = (Ne/N)I - (λ/N)(H - μI)
    # Rearranging terms for MPO addition: 
    # ρ0 = [ (Ne + λμ)/N ] * I  -  [ λ/N ] * H
    
    coeff_I = (Ne + λ * μ) / N
    coeff_H = - (λ / N)

    rho_0 = +(coeff_I * Identity_MPO(sites), coeff_H * H; cutoff = 1e-10, maxdim = 100)
    
    return rho_0
end

begin
t(x) = 1.0
v(x) =  1.0*(-1.0)^(isodd(x)) + 0.5*cos(2pi * (sqrt(5)-1)/2.0 * (x.-1)) + 1.0*cos(2pi * (sqrt(5)-1)/2.0 /(2^13)* (x.-1))
L = 15
N = 2^L

sites = siteinds("Qubit", L)

end


H0 = Build_H_MF_1D(v,t, ComplexF64, sites, 1e-8, 300)
begin
ϵ = 1e-8
maxχ = 300


H_min = -3.1
H_max = 3.1


Ne = div(N,2)
end

begin 
ρ0 = construct_rho_0(N, H0, H_max, H_min, Ne,sites)
println("Initial trace: ", real(tr(ρ0)))
println("Initial MPO constructed.\n")
denom_p = 200
for i in 1:40
    println("--- Step $i ---")
    T1_p = real(tr(ρ0))

    use_mcweeny = false

    P2 = apply(ρ0, ρ0; cutoff = ϵ, maxdim = maxχ)


    P3 = apply(ρ0, P2;  cutoff = ϵ, maxdim = maxχ)
    println("χ_1 : ", maximum(linkdims(ρ0)), " χ_2 : " , maximum(linkdims(P2)), " χ_3 : ", maximum(linkdims(P3)))

    T1 = real(tr(ρ0))
    T2 = real(tr(P2))
    T3 = real(tr(P3))

    denom = T1 - T2

    cn = (T2 - T3) / denom

    if abs(denom) < 1e-4 || (i > 5 && abs(cn - 0.5) < 0.05) && !use_mcweeny
        use_mcweeny = true
        println("Step $i: Switching to McWeeny (Fixed c_n = 0.5) for stability.")
    end
    
    
    
    if use_mcweeny
    # STABLE MODE: McWeeny Purification
    # Formula: 3ρ^2 - 2ρ^3
    # We calculate this efficiently: 3*P2 - 2*P3
    
    ρ0 = +(3.0 * P2, -2.0 * P3; cutoff=1e-10, maxdim=100)
    else

    c_ρ, c_p2, c_p3 = 0.0, 0.0 ,0.0



    if cn < 0.5
        inv_1_minus_cn = 1.0 / (1.0 - cn)
        c_ρ = (1 - 2*cn) * inv_1_minus_cn
        c_p2 = (1.0 + cn) * inv_1_minus_cn
        c_p3 = -1.0 * inv_1_minus_cn
    else
        inv_cn = 1.0 / cn
        c_ρ = 0.0
        c_p2 = (1.0 + cn) * inv_cn
        c_p3 = -1.0 * inv_cn
    end

        ρ0 = +(c_ρ * ρ0, c_p2 * P2, c_p3 * P3; cutoff = ϵ, maxdim = maxχ)
    
end
    rel_error = (T1 - T2) / T1


    println("  Trace (Ne): $T1") 
    println("  Abs Error: $denom")
    println("  Rel Error: $rel_error") # This should be < 1.0 and drop to 1e-6

    if rel_error < 1e-5 || (abs(T1 - T1_p) < 1e-6 && i>5)
        println("\nConverged!\n")
        break
    end
    denom_p = denom


end

println("Final trace: ", real(tr(ρ0))," Expected Ne: ", Ne)

end
begin
sites_lattice = LinRange(1, N, 512)
diagonal = zeros(Float64, length(sites_lattice))
for i in 1:length(sites_lattice)
    diagonal[i] = MatrixChecker(ρ0, sites, i,i)
end
fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], xlabel = "Site", ylabel = "Charge Density", title = "Charge Density Profile")
lines!(ax, sites_lattice, diagonal)
fig |> display
end 