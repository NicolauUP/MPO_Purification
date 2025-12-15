
using ITensors, ITensorMPS
using CairoMakie
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using SparseArrays, LinearAlgebra, FFTW, Statistics, LsqFit


include("operators.jl")
include("quantics_functions.jl")


t(x) = 1.0

L = 10
N = 2^L

sites = siteinds("Qubit", L)

