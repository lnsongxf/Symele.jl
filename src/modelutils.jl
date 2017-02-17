"""
   solveSS!(model,quiet=false)
Solves the steady state for the model as currently parameterised and updates model.steadystate.values. Setting quiet=true disable output; this is useful when repeatedly solving, such as during estimation)
"""
function solveSS!(model,quiet=false)
    if !model.steadystate.userdefined
        sswrapper = x->evaluateSS(model.steadystate.system_lambda,x,model.parameters.values)

        initialeval = sswrapper(model.steadystate.values)
        if any(isnan(initialeval)) || any(initialeval.== Inf) || any(initialeval.== -Inf)
            error("Steady state guess does not evaluate (e.g. due to divide-by-zero). Supply a different initial guess.")
        end
        nlout = NLsolve.nlsolve(NLsolve.not_in_place(sswrapper), model.steadystate.values)
        if !NLsolve.converged(nlout)
            if !quiet
                println("FAILED.")
                println("Residuals at initial guess:")
                println(sswrapper(model.steadystate.values))
                println("Parameters:")
                println(Dict(zip(model.parameters.names,model.parameters.values)))
            end
            error("Could not find steady state.")
        end
        model.steadystate.values[:] = nlout.zero
        return
    else
        model.steadystate.values[:] = model.steadystate.ssfunction(model.parameters.values)
    end
end

"""
   function updateparameters!(model::DSGEModel,newvalues::Array{Float64,1},indices::Array{Int,1},quiet=false)
Update model.parameters.values by setting the parameters given by indices to the values given by newvalues. Typically this is only used internally.
Function recalculates the steadystate; quiet disables output from the steady state solver. See the help file for solveSS!
"""
function updateparameters!(model::DSGEModel,newvalues::Array{Float64,1},indices::Array{Int,1},quiet=false)
    model.parameters.values[indices] = newvalues[:]
    solveSS!(model,quiet)
    return
end

"""
    function updateparameters!(model::DSGEModel,newvalues::Array{Float64,1},quiet=false)
Update model.parameters.values by setting it to the values given by newvalues. Function recalculates the steadystate; quiet disables output from the steady state solver. See the help file for solveSS!
"""
function updateparameters!(model::DSGEModel,newvalues::Array{Float64,1},quiet=false)
    model.parameters.values[:] = newvalues[:]
    solveSS!(model,quiet)
    return
end

"""
    function updateparameters!(model::DSGEModel,parameterlist::Array{String,1},newvalues::Array{Float64,1})
Update model.parameters.values by setting the parameters given by parameterlist to the corresponding value in newvalues. Function recalculates the steadystate; quiet disables output from the steady state solver. See the help file for solveSS!
"""
function updateparameters!(model::DSGEModel,parameterlist::Array{String,1},newvalues::Array{Float64,1})
    updateparameters!(model,newvalues,getparameterindices(model,parameterlist))
    return
end

updateparameters!(model::DSGEModel,parameterlist::String,newvalues::Number) = updateparameters!(model,[parameterlist],[newvalues])
    

function getparameterindices(model::DSGEModel,parameters::Array{String,1})
    parameterindices = Array(Int,length(parameters))
    for (i,name) in enumerate(parameters)
        if name in keys(model.parameters.dictionary)
            parameterindices[i] = model.parameters.dictionary[name]
        else
            error("$(name) is not a parameter in your model.")
        end
    end
    return parameterindices
end

function getvariableindices(model::DSGEModel,vars::Array{String,1})
    indices = Array(Int,length(vars))
    for (i,name) in enumerate(vars)
        if name in keys(model.variables.dictionary)
            indices[i] = model.variables.dictionary[name]
        else
            error("$(name) is not a variable in your model.")
        end
    end
    return indices
end

function getshockindices(model::DSGEModel,shocks::Array{String,1})
    indices = Array(Int,length(shocks))
    for (i,name) in enumerate(shocks)
        if name in keys(model.shocks.dictionary)
            indices[i] = model.shocks.dictionary[name]
        else
            error("$(name) is not a shocks in your model.")
        end
    end
    return indices
end
