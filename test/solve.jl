# Check that solve works
function checksolves(model, solutionfile)
    sol = solve(model,true)
    return true
end

# Solve and check the solutions for the example models
function returncomparisons(model, solutionfile)
    sol = solve(model,true)
    d = load(solutionfile)
    preloadedsolution = d["solution"]
    return preloadedsolution, sol
end

# Run the tests

@test checksolves(threeEQNK,"data/threeEQNKsolution.jld")
preloadedsolution,sol = returncomparisons(threeEQNK,"data/threeEQNKsolution.jld")
@test preloadedsolution.MX ≈ sol.MX
@test preloadedsolution.ME ≈ sol.ME
@test preloadedsolution.MC ≈ sol.MC

@test checksolves(RBC,"data/RBCsolution.jld")
preloadedsolution,sol = returncomparisons(RBC,"data/RBCsolution.jld")
@test preloadedsolution.MX ≈ sol.MX
@test preloadedsolution.ME ≈ sol.ME
@test preloadedsolution.MC ≈ sol.MC