type MLEEstimation
    parameters::Array{String,1}
    startpoint::Array{Float64,1}
    lbd::Array{Float64,1}
    ubd::Array{Float64,1}
    observed::Array{String,1}
    observedseries::Array{Float64}
    parameterindices::Array{Int,1}
    algorithm::Symbol
end

type MLEEstimationResults
    mode::Dict{String,Float64}
end

"""
   MLEEstimation(model::DSGEModel,parameters::Union{String,Array{String,1}},startpoint::Union{Number,Array{Float64,1}},lbd::Union{Number,Array{Float64,1}},ubd::Union{Number,Array{Float64,1}},observed::Union{String,Array{String,1}},observedseries::Array{Float64};algorithm::Symbol=:LN_COBYLA)
Define the parameters for a MLE estimation.
parameters are the parameters you wish to estimation. startpoint corresponds to an initial guess for the numerical optimiser. lbd and ubd are lower and upper bounds for the parameter; set to -Inf and Inf for unbounded (can be one sided). observed is an array of strings corresponding to which variables in your model are observed. observedseries is a T*k Array of data, where T is the length of your time series and the k rows correspond (in order) to the variables in observed. algorithm allows you to specific which optimisation algorithm NLopt uses. See NLopt's documentation for more information.
Returns an instance of type MLEEstimation, which can then be used to run the actual estimation.
"""
function MLEEstimation(model::DSGEModel,parameters::Union{String,Array{String,1}},startpoint::Union{Number,Array{Float64,1}},lbd::Union{Number,Array{Float64,1}},ubd::Union{Number,Array{Float64,1}},observed::Union{String,Array{String,1}},observedseries::Array{Float64};algorithm::Symbol=:LN_COBYLA)
    if isa(parameters,String)
        parameters = [parameters]
    end
    if isa(observed,String)
        observed = [observed]
    end
    if isa(startpoint,Number)
        startpoint = [startpoint]
    end
    if isa(lbd,Number)
        lbd = [lbd]
    end
    if isa(ubd,Number)
        ubd = [ubd]
    end
    # This is the vector of parameters we're tying to estimate
    parameterindices = Array(Int,length(parameters))
    for (i,name) in enumerate(parameters)
        if name in keys(model.parameters.dictionary)
            parameterindices[i] = model.parameters.dictionary[name]
        else
            error("$(name) is not a parameter in your model.")
        end
    end
    if length(lbd)!=length(parameters)
        error("You gave $(length(lbd)) lower bounds for $(length(parameters)) parameters.")
    end
    if length(ubd)!=length(parameters)
        error("You gave $(length(ubd)) lower bounds for $(length(parameters)) parameters.")
    end
    if length(observed) != size(observedseries,2)
        error("You specified $(length(observed)) observed variables, but your data series has $(size(observedseries,2)) columns.")
    end
    
    if length(observed) > model.meta.numshocks
        error("Stochastic singularity. You cannot have more observed than the number of shocks.")
    end
    
    return MLEEstimation(parameters,startpoint,lbd,ubd,observed,observedseries,parameterindices,algorithm)
end

"""
   MLEEstimation(model::DSGEModel,parameters::Union{String,Array{String,1}},startpoint::Union{Number,Array{Float64,1}},lbd::Union{Number,Array{Float64,1}},ubd::Union{Number,Array{Float64,1}},observedseries::DataFrames.DataFrame;algorithm::Symbol=:LN_COBYLA)
Define the parameters for a MLE estimation.
parameters are the parameters you wish to estimation. startpoint corresponds to an initial guess for the numerical optimiser. lbd and ubd are lower and upper bounds for the parameter; set to -Inf and Inf for unbounded (can be one sided). observedseries is a DataFrame with column names equal to which variables in your model are observed. algorithm allows you to specific which optimisation algorithm NLopt uses. See NLopt's documentation for more information.
Returns an instance of type MLEEstimation, which can then be used to run the actual estimation.
"""
function MLEEstimation(model::DSGEModel,parameters::Union{String,Array{String,1}},startpoint::Union{Number,Array{Float64,1}},lbd::Union{Number,Array{Float64,1}},ubd::Union{Number,Array{Float64,1}},observedseries::DataFrames.DataFrame;kwargs...)
    # Handle the dataframe
    observedsymbols = DataFrames.names(observedseries)
    observed = Array(String,length(observedsymbols))
    for i in eachindex(observed)
        observed[i] = string(observedsymbols[i])
    end
    observeddata = Array(Float64,size(observedseries,1),length(observedsymbols))
    for i in eachindex(observedsymbols)
        observeddata[:,i] = observedseries[observedsymbols[i]]
    end
    
    return MLEEstimation(model,parameters,startpoint,lbd,ubd,observed,observeddata;kwargs...)
end

"""
    estimation(estimation::MLEEstimation,model::DSGEModel;laypunovmaxattempts=10000)
estimation is an MLEEstimation that you have already defined. Call this function on that estimation to run the estimation job.
Returns a type MLEResults, containing the results from the estimation.
"""
function estimate(estimation::MLEEstimation,model::DSGEModel;laypunovmaxattempts=10000)
    MM = observationmatrix(model,estimation)

    estimationoptions = KalmanOptions(MM,laypunovmaxattempts)

    # Define the target function, grad doesn't get used, because I use a derivative-free 
    ll = (x,_) -> begin
        # The 1 gets the ll and ignores the shocks and states that get returned
        return kalmanloglikelihood(estimationoptions,model,x,estimation.parameterindices,estimation.observedseries,false)[1]
    end
    # Check initial values
    initialeval = ll(estimation.startpoint,[])
    if isnan(initialeval) || initialeval == Inf || initialeval == -Inf
        error("Steady state guess does not evaluate (e.g. due to divide-by-zero). Supply a different initial guess.")
    end

    # Define the optimisation problem
    opt = NLopt.Opt(estimation.algorithm, length(estimation.parameters))
    NLopt.max_objective!(opt, ll)
    NLopt.lower_bounds!(opt, estimation.lbd)
    NLopt.upper_bounds!(opt, estimation.ubd)
    optf,optx,ret = NLopt.optimize(opt, estimation.startpoint)
    if ret == :SUCCESS
        return MLEEstimationResults(Dict{String,Float64}(zip(estimation.parameters,optx)))
    elseif ret == :ROUNDOFF_LIMITED
        warn("Optimisation was roundoff limited. May not be optimal.")
        return MLEEstimationResults(Dict{String,Float64}(zip(estimation.parameters,optx)))
    else
        error("Numerical optimiser could not converge to mode. Perhaps try a different initial guess. NLopt return code: $(ret).")
    end
end
