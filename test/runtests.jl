using Symele
using JLD
using Base.Test

# Load up the two example models, because I'm going to want to keep re using them (loading does get tested, just later)
threeEQNK = DSGEModel("../examples/threeEQNK.txt")
RBC = DSGEModel("../examples/RBC.txt")

# Run the tests of the model constructor
@testset "Model setup" begin
    include("modelsetup.jl")
end

# Check model utils work properly
@testset "Model utilities" begin
    include("modelutils.jl")
end

# Run the tests to solve the models
@testset "Model solutions (first-order)" begin
    include("solve.jl")
end

@testset "IRFs" begin
    include("irfs.jl")
end

@testset "Estimation" begin
    include("estimation.jl")
end
NKdata = load("data/NKdata.jld")
NKdata = NKdata["df"]
RBCdata = load("data/RBCdata.jld")
RBCdata = RBCdata["df"]
@testset "Simulation" begin
    include("simulation.jl")
end

@testset ("Shock decomposition") begin
    include("shockdecomposition.jl")
end