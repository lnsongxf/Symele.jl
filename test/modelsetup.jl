function testloading(modelfile;kwargs...)
    DSGEModel(modelfile;kwargs...)
    return true
end

function getRBCSS(model,parameters)
    alpha = model.parameters.values[model.parameters.dictionary["alpha"]]
    beta = model.parameters.values[model.parameters.dictionary["betta"]]
    delta = model.parameters.values[model.parameters.dictionary["delta"]]
    sigma = model.parameters.values[model.parameters.dictionary["sigma"]]
    RBCSS = Array(Float64,RBC.meta.numvars)
    RBCSS[RBC.variables.dictionary["A"]] = 1
    RBCSS[RBC.variables.dictionary["k"]]  = log(((1/beta - 1 + delta)/alpha)^(1/(alpha-1)))
    RBCSS[RBC.variables.dictionary["Y"]]  = log(((1/beta - 1 + delta)/alpha)^(alpha/(alpha-1)))
    RBCSS[RBC.variables.dictionary["Inv"]]  = log(delta*((1/beta - 1 + delta)/alpha)^(1/(alpha-1)))
    RBCSS[RBC.variables.dictionary["r"]]  = log(1/beta - 1 + delta)
    RBCSS[RBC.variables.dictionary["c"]]  = log(((1/beta - 1 + delta)/alpha)^(alpha/(alpha-1)) - delta*((1/beta - 1 + delta)/alpha)^(1/(alpha-1)))
    RBCSS[RBC.variables.dictionary["lammbda"]]  = log((((1/beta - 1 + delta)/alpha)^(alpha/(alpha-1)) - delta*((1/beta - 1 + delta)/alpha)^(1/(alpha-1)))^(-sigma))
    return RBCSS
end

# Test setting up some basic examples

@test testloading("../examples/threeEQNK.txt")
# Check that the steady state is right
@test threeEQNK.steadystate.values ≈ zeros(Float64,threeEQNK.meta.numvars)

@test testloading("../examples/RBC.txt")
# RBC model has an analytical steady state. Check that the numerically computed SS is (within numerical tolerance) correct.
wrappedRBCSS = x->getRBCSS(RBC,x)
@test RBC.steadystate.values ≈ wrappedRBCSS(RBC.parameters.values)

# Try one that supplies a steady state function to check that works too
@test testloading("../examples/RBC.txt";ssfunction=wrappedRBCSS,ssfunctionordered=true)
ssfunctionRBC = DSGEModel("../examples/RBC.txt";ssfunction=wrappedRBCSS,ssfunctionordered=true)
@test RBC.steadystate.values ≈ wrappedRBCSS(RBC.parameters.values)