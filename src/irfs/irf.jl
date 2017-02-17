"""
    function irf(model,responses,impulsess;periods=20,graph=true)
Construct the IRFs for the variables given by responses to the shocks given by impulses (both are arrays of strings).
Returns a double dictionary with first key corresponding to the shock and second to the response variable.
"""
function irf(model,responses,impulses;kwargs...)
    solution = solve(model,true)
    return irf(model,solution,responses,impulses;kwargs...)
end

"""
    function irf(models;periods=20,graph=true)
Construct IRFs for all shocks and variables in the model. Solves the model at the current parameter set.
Returns a double dictionary with first key corresponding to the shock and second to the response variable.
"""
irf(model;kwargs...) = irf(model,model.variables.names,model.shocks.names;kwargs...)

"""
    irf(model,solutions;periods=20,graph=true)
Construct IRFs for all shocks and variables in the model. Uses the supplied solution rather than solving the model again.
Returns a double dictionary with first key corresponding to the shock and second to the response variable.
"""
irf(model,solution;kwargs...) = irf(model,solution,model.variables.names,model.shocks.names;kwargs...)

"""
    irf(model,solution,responses,impulses;periods=20,graph=true)
Construct the IRFs for the variables given by responses to the shocks given by impulses (both are arrays of strings). Uses the supplied solution rather than solving the model again.
Returns a double dictionary with first key corresponding to the shock and second to the response variable.
"""
function irf(model,solution,responses,impulses;periods=20,graph=true)
    if isa(responses,String)
        responses = [responses]
    end
    if isa(impulses,String)
        impulses = [impulses]
    end
    # Check the responses and impulses
    responseindices = getvariableindices(model,responses)
    impulseindices = getshockindices(model,impulses)

    IRFs = Dict{String,Dict{String,Array{Float64,1}}}()

    for (impulse,impulsename) in zip(impulseindices,impulses)
        E = zeros(Float64,model.meta.numshocks,1)
        E[impulse,1] = 1.0

        Y = Array(Float64,model.meta.numvars,periods)
        Y[:,1] = solution.ME*E
        for t in 2:periods
            Y[:,t] = solution.MX*Y[:,t-1]
        end

        IRFs[impulsename] = Dict{String,Array{Float64,1}}()
        for (index,response) in zip(responseindices,responses)
            IRFs[impulsename][response] = Y[index,:]
        end
    end
    if graph
        graphirf(model,IRFs)
    end
    return IRFs
end

"""
    function graphirf(model,irf)
Graphs all the IRFs contained in irf (which must be in the same format as the returned output from a call to irf) for the given model. Called by irf if graph=true is set (as it is by default).
"""
function graphirf(model,irf)
    for (shockname,irfs) in irf
        numplots = length(irfs)
        numfigures = ceil(Int,numplots/9)
        PyPlot.figure()
        plotsize = [0 0]
        plotct = 1
        j = 0
        for (varname,response) in irfs
            j += 1
            if plotct == 1
                if numplots - j + 1 > 6
                    plotsize = [3 3]
                    elseif numplots - j + 1 > 4
                    plotsize = [3 2]
                    elseif numplots - j + 1 > 2
                    plotsize = [2 2]
                    elseif numplots - j + 1 > 1
                    plotsize = [2 1]
                else
                    plotsize = [1 1]
                end
                if numfigures > 1
                    temptitle = "IRFs to $(shockname) - $(figurect) of $(numfigures)"
                    else 
                    temptitle = "IRFs to $(shockname)"
                end
                PyPlot.suptitle(temptitle)
            end
            PyPlot.subplot(plotsize[1],plotsize[2],plotct)
            x = 1:size(response,1)
            PyPlot.plot(x,response)
            ax = PyPlot.gca()
            ax[:set_xlim]([0,size(response,1)])
            PyPlot.title(varname)
            plotct += 1
            if plotct > plotsize[1]*plotsize[2]
                plotct = 1
                PyPlot.figure()
            end
        end
    end
    return
end
