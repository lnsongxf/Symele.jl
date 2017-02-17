type Meta # Because I want to be able to scale solution order
    numstates::Int
    numvars::Int
    numshocks::Int
    numparameters::Int
    solutionorder::Int
    solutionmethod::String
end

immutable Parameters
    names::Array{String,1}
    values::Array{Float64,1}
    dictionary::Dict{String,Int}
    syms::Array{SymPy.Sym,1}
end

immutable Variables
    names::Array{String,1}
    dictionary::Dict{String,Int}
    isstate::Array{Bool,1}
    state_syms::Array{SymPy.Sym,1}
    lead_syms::Array{SymPy.Sym,1}
    contemp_syms::Array{SymPy.Sym,1}
    forwardlooking::Array{Bool,1}
end

immutable Shocks
    names::Array{String,1}
    dictionary::Dict{String,Int}
    syms::Array{SymPy.Sym,1}
end

immutable Equations
    stringversions::Array{String,1}
    syms::Array{SymPy.Sym,1}
end

immutable SteadyState
    system_sym::Array{SymPy.Sym,1}
    system_lambda::Array{Function,1}
    values::Array{Float64,1}
    userdefined::Bool
    ssfunction::Function
end

immutable FirstOrderPerturbation{T}
    FYp::Array{T,2}
    FXp::Array{T,2}
    FY::Array{T,2}
    FX::Array{T,2}
    FXm::Array{T,2}
    FE::Array{T,2}
end
 

immutable NthOrderPerturbation{T}
    FYp::SparseMatrixCSC{T}
    FXp::SparseMatrixCSC{T}
    FY::SparseMatrixCSC{T}
    FX::SparseMatrixCSC{T}
    FXm::SparseMatrixCSC{T}
    FE::SparseMatrixCSC{T}
end

type PerturbationSystem{T} # Not immutable, so you can add higher orders of perturbation afterwards
    firstorder::FirstOrderPerturbation{T}
    higherorders::Array{NthOrderPerturbation{T},1}
end

type DSGEModel
    meta::Meta
    parameters::Parameters
    variables::Variables
    shocks::Shocks
    equations::Equations
    steadystate::SteadyState
    perturbationsystem_sym::PerturbationSystem{SymPy.Sym}
    perturbationsystem_lambda::PerturbationSystem{Function}
end