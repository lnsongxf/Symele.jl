immutable StatesAndShocksEstimation
    observed::Array{String,1}
    series::Array{Float64,2}
end

"""
    getstatesandshocks(model::DSGEModel,observed::Array{String},series::Array{Float64};laypunovmaxattempts=10000,graph=true)
Returns the Kalman-filter-estimated shocks and variables for the data in the DataFrame 'observedseries'. laypunovmaxattempts restricts how many iterations the Kalman filter will short circuit at if it cannot converge to the unconditional variance (see documentation for kalmanfilter).

Returns two dictionarys: the first is the estimated variables (keys are variable names); the second is the estimated shocks (keys are shock names).
"""
function getstatesandshocks(model::DSGEModel,observed::Array{String},series::Array{Float64};laypunovmaxattempts=10000)
    estimation = StatesAndShocksEstimation(observed,series)
    MM = observationmatrix(model,estimation)
    estimationoptions = KalmanOptions(MM,laypunovmaxattempts)
    
    solution = solve(model,true,false)
    _,shocks,states = kalmanloglikelihood(estimationoptions,model,solution,estimation.series,true)
    
    # Format shocks and estimated variables into a nice dictionary
    outputshocks = Dict{String,Array{Float64}}()
    for (i,shock) in enumerate(model.shocks.names)
        outputshocks[shock] = shocks[i,:]
    end
    outputvars = Dict{String,Array{Float64}}()
    for (i,var) in enumerate(model.variables.names)
        outputvars[var] = states[i,:]
    end

    return outputvars,outputshocks
end

"""
    getstatesandshocks(model::DSGEModel,observed::Array{String},series::Array{Float64};laypunovmaxattempts=10000,graph=true)
Returns the Kalman-filter-estimated shocks and variables for the data in the DataFrame 'observedseries'. laypunovmaxattempts restricts how many iterations the Kalman filter will short circuit at if it cannot converge to the unconditional variance (see documentation for kalmanfilter).

Returns two dictionarys: the first is the estimated variables (keys are variable names); the second is the estimated shocks (keys are shock names).
"""
function getstatesandshocks(model::DSGEModel,observedseries::DataFrames.DataFrame;kwargs...)
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
    return getstatesandshocks(model,observed,observeddata;kwargs...)
end
