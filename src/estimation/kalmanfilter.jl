type KalmanOptions
    MM::Array{Float64,2}
    maxattempts::Int # When to short circuit out on the laypunov iteration
end

function kalmanloglikelihood(estimationoptions::KalmanOptions,model::DSGEModel,solution::DSGESolution,observedseries::Array{Float64},returnstates::Bool=false)
    T = size(observedseries,1)
    sbar = copy(model.steadystate.values)
    P = solution.ME*eye(Float64,size(solution.ME,2))*solution.ME'
    

    # Iterate to convergence to get the unconditional variance
    crit = 1.0
    tol = 1e-8
    i = 0
    while crit > tol && i < estimationoptions.maxattempts
        Pn = solution.MX*P*solution.MX' + solution.ME*solution.ME'
        crit = maxabs(Pn-P)
        P,Pn = Pn,P
        i+=1
    end

    logpdf = (x,mean,sigma) -> (-size(x,1)/2)*(log(2) + log(pi)) - 0.5*log(det(sigma)) - 0.5*dot((x-mean),(sigma\(x-mean)))

    ll = 0
    if returnstates
        shocks = Array(Float64,model.meta.numshocks,T)
        states = Array(Float64,model.meta.numvars,T)
    else
        shocks = Array(Float64,0)
        states = Array(Float64,0)
    end
    for t= 1:T
        # Estimated state, conditional on last period
        sbarm = solution.MX*(sbar-model.steadystate.values) + model.steadystate.values
        Pm = solution.MX*P*solution.MX' + solution.ME*solution.ME'
        
        # Expected observed, given estimated state
        ym = estimationoptions.MM*sbarm
        Fm = estimationoptions.MM*Pm*estimationoptions.MM'
        #return observedseries[t,:],ym,Fm
        ll += logpdf(observedseries[t,:],ym,Fm)
        
        # Now update estimated state
        sbar = sbarm + Pm*estimationoptions.MM'*(Fm\(observedseries[t,:]-ym))
        P = Pm - Pm*estimationoptions.MM'*(Fm\estimationoptions.MM)*Pm
        
        if returnstates
            # Don't do it otherwise, to avoid unnecessary allocations
            shocks[:,t] = solution.ME\(sbar - solution.MX*sbarm)
            states[:,t] = sbar
        end
    end
    return ll,shocks,states
end

function kalmanloglikelihood(estimationoptions::KalmanOptions,model::DSGEModel,newparameters::Array{Float64},parameterindices::Array{Int},observedseries::Array{Float64},returnstates::Bool=false)
    updateparameters!(model,newparameters,parameterindices,true)
    local solution
    try
        solution = solve(model,true,false)
    catch
        # Probably failed to satisfy BK at these parameters, that or SS didn't converge. Return -Inf for the likelihood.
        return -Inf,Array(Float64,0),Array(Float64,0)
    end
    return kalmanloglikelihood(estimationoptions,model,solution,observedseries,returnstates)
end

function kalmanloglikelihood(estimationoptions::KalmanOptions,model::DSGEModel,newparameters::Array{Float64},observedseries::Array{Float64})
    if length(newparameters)!=model.meta.numparameters
        error("Parameter vector is wrong size: $(length(newparameters)) instead of $(model.meta.numparameters).")
    end
    updateparameters!(model,newparameters,true)
    local solution
    try 
        solution = solve(model,true,false)
    catch
        # Failed to satisfy BK at these parameters. Return -Inf for the likelihood. (And empty arrays for the shocks and states)
        return -Inf,Array(Float64,0),Array(Float64,0)
    end
    return kalmanloglikelihood(estimationoptions,model,solution,observedseries,returnstates)
end  

kalmanloglikelihood(estimationoptions::KalmanOptions,model::DSGEModel,observedseries::Array{Float64},returnstates::Bool=false) = kalmanloglikelihood(estimationoptions,model,model.parameters.values,collect(1:model.meta.numvars),observedseries,returnstates)
