
using ITensors, ITensorMPS
using CairoMakie
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using LinearAlgebra


include("operators.jl")
include("quantics_functions.jl")


t(x) = 1.0
v(x) = 1e-15
L = 10
N = 2^L

sites = siteinds("Qubit", L)




H0 = Build_H_MF_1D(v,t, ComplexF64, sites, 1e-10, 100)



H_min = -2.1
H_max = 2.1


μ = tr(H0)/N

ρ = 

