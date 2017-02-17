NKdata = load("data/NKdata.jld")
NKdata = NKdata["df"]
RBCdata = load("data/RBCdata.jld")
RBCdata = RBCdata["df"]

@testset "MLE" begin
    # Single parameter estimation
    parameters = ["theta"]
    startpoint = [0.8]
    lbd = [0.0+eps(0.0)]
    ubd = [1.0 - eps(1.0)]
    estimation = MLEEstimation(threeEQNK,parameters,startpoint,lbd,ubd,NKdata)
    MLEresults = estimate(estimation,threeEQNK)
    @test MLEresults.mode["theta"] ≈ 0.7652012392878662
    
    # Two parameter estimation
    parameters = ["rho","betta"]
    startpoint = [0.0,0.99]
    lbd = [0.0,0.0]
    ubd = [1.0,1.0]
    estimation = MLEEstimation(RBC,parameters,startpoint,lbd,ubd,RBCdata)
    MLEresults = estimate(estimation,RBC)
    @test MLEresults.mode["betta"] ≈ 0.9602637168127645 
    @test MLEresults.mode["rho"] ≈ 0.9035878776313421
end

@testset "RWMH" begin
    using Distributions
    parameters = ["theta"]
    priors = Dict("theta"=>Uniform(0,1))

    estimation = RWMHEstimation(threeEQNK,parameters,priors,NKdata)
    RWMHresults = estimate(estimation,threeEQNK)
    @test RWMHresults.median["theta"] ≈ 0.7567788408831552 
    @test RWMHresults.mean["theta"] ≈ 0.754462791022537
    
    # Two parameter estimation
    parameters = ["rho","betta"]
    priors = Dict("rho"=>Uniform(0,1),"betta"=>Uniform(0.85,1.0))
    estimation = RWMHEstimation(RBC,parameters,priors,RBCdata)
    RWMHresults = estimate(estimation,RBC)
    @test RWMHresults.median["rho"] ≈ 0.99
    @test RWMHresults.mean["rho"] ≈ 0.99
    @test RWMHresults.median["betta"] ≈ 0.0
    @test RWMHresults.mean["betta"] ≈ 0.0
end