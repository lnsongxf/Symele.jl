"""
    shockdecomposition(model::DSGEModel,observed::Array{String},series::Array{Float64};laypunovmaxattempts=10000,graph=true)
Returns the Kalman-filter-estimated contribution of each shock to the observed variables. Specify which of the k series are observed by passing in a k-element array of strings corresponding to the variable names. The data are passed in as a T by k array 'series', with the k columns corresponding, in order, to the k observed series. laypunovmaxattempts restricts how many iterations the Kalman filter will short circuit at if it cannot converge to the unconditional variance (see documentation for kalmanfilter).

Returns a double dictionary of the decomposition. The first key indexes the observed series; the second key indexes the contribution to that series from each shock.
"""
function shockdecomposition(model::DSGEModel,observed::Array{String},series::Array{Float64};laypunovmaxattempts=10000,graph=true)
    vars,shocks = getstatesandshocks(model::DSGEModel,observed::Array{String},series::Array{Float64};laypunovmaxattempts=10000)
    # calculate contributions
    contributions = Dict{String,Array{Float64,2}}()
    for (name,value) in shocks
        feedin = zeros(Float64,length(value),model.meta.numshocks)
        feedin[:,model.shocks.dictionary[name]] = value
        contributions[name] = simulate(model,feedin)
    end
    
    # Now do the decomposition
    decomposition = Dict{String,Dict{String,Array{Float64,1}}}()
    for var in model.variables.names
        thisvardict = Dict{String,Array{Float64,1}}()
        initial = vars[var]
        for shock in model.shocks.names
            thisvardict[shock] = contributions[shock][:,model.variables.dictionary[var]]
            initial -= thisvardict[shock]
        end
        thisvardict["initialcondition"] = initial
        thisvardict[var] = vars[var]
        decomposition[var] = thisvardict
    end
    
    if graph
        plotdecomposition(decomposition)
    end
    return decomposition
end

"""
    shockdecomposition(model::DSGEModel,observedseries::DataFrames.DataFrame,vars::Array{String}=model.variables.names;laypunovmaxattempts=10000,graph=true)
Returns the Kalman-filter-estimated contribution of each shock to the observed variables in the DataFrame 'observedseries'. laypunovmaxattempts restricts how many iterations the Kalman filter will short circuit at if it cannot converge to the unconditional variance (see documentation for kalmanfilter).

Returns a double dictionary of the decomposition. The first key indexes the observed series; the second key indexes the contribution to that series from each shock.
"""
function shockdecomposition(model::DSGEModel,observedseries::DataFrames.DataFrame;kwargs...)
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
    return shockdecomposition(model,observed,observeddata;kwargs...)
end

"""
    plotdecomposition(decomposition)
Plots graphs showing the contribution of each shock (and initial conditions) to the variables. Takes the dictionary returned by shockdecomposition as input. No return value. Mostly an internal function; called by setting the keyword argument graph=true (which is the default) when calling shockdecomposition.
"""
function plotdecomposition(decomposition)
    for (var,decomp) in decomposition
        PyPlot.figure()
        PyPlot.title(var)
        legend = ()
        for (shock,values) in decomp
            if shock != var
                PyPlot.bar(1:length(values), values)
                legend = (legend...,shock)
            else
                PyPlot.plot(1:length(values),values)
                legend = (legend...,var)
            end
        end
        PyPlot.legend(legend)
    end
    return
end
