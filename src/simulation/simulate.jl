"""
    simulate(model::DSGEModel,shocks::Array{Float64,2})
Simulates the model for T periods using the T by k set of shocks (where k is the number of shocks in the model).
Returns a T by n Array, where n is the number of variables in the model.
"""
function simulate(model::DSGEModel,shocks::Array{Float64,2})
    T = size(shocks,1)
    shocks = permutedims(shocks,[2 1]) # Speed up memory access during simulation. This may not be quicker. I should profile.
    solution = solve(model,true,false)
    out = Array(Float64,model.meta.numvars,T)
    out[:,1] = solution.ME*shocks[:,1]
    for t = 2:T
        out[:,t] = solution.MX*(out[:,t-1] - model.steadystate.values) + solution.ME*shocks[:,t] + model.steadystate.values
    end
    out = permutedims(out,[2 1])
    return out
end

"""
    simulate(model::DSGEModel,T::Int
Simulates the model for T periods by randomly simulating shocks from a normal distribution. Returns a T by n Array, where n is the number of variables in the model.
"""
function simulate(model::DSGEModel,T::Int)
    shocks = randn(T,model.meta.numshocks)
    return simulate(model,shocks)
end
