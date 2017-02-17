type RWMHEstimation
    parameters::Array{String,1}
    priors::Array{Any,1}
    priors_lambda::Array{Function,1}
    parameterindices::Array{Int,1}
    observed::Array{String,1}
    observedseries::Array{Float64}
    nchains::Int
    chainlength::Int
    burnin::Int
    proposalcovariance::Array{Float64,2}
    jumpscale::Float64
end

immutable RWMHChain
    mean::Array{Float64}
    median::Array{Float64}
    draws::Array{Float64,2}
    acceptancerate::Float64
end

immutable RWMHEstimationResults
    mean::Dict{String,Float64}
    median::Dict{String,Float64}
    posterior::Array{Float64,2}
    chains::Array{RWMHChain,1}
end

"""
   RWMHEstimation(model::DSGEModel,parameters::Union{String,Array{String,1}},priors::Dict,observed::Union{String,Array{String,1}},observedseries::Array{Float64};nchains::Int=2,chainlength::Int=10000,burnin::Int=1000,proposalcovariance::Array{Float64,2}=zeros(Float64,0,0),jumpscale::Float64=0.2)
Define the settings for a random walk metropolis hastings estimation.
parameters are the parameters you wish to estimation. priors is a dictionary with keys of the parameters you are estimating (as strings) and values as the prior distribution. Using distributions from the Distributions package is supported (and strongly recommended). You can, alternatively, supply your own anonymous function of your desired PDF. observed is an array of strings corresponding to which variables in your model are observed. observedseries is a T*k Array of data, where T is the length of your time series and the k rows correspond (in order) to the variables in observed.
Optional keyword arguments are:
nchains: the number of chains used to sample from the prior
chainlength: the length of each chain
burnin: how much to drop from the start of the chain as burn in
proposalcovariance: the covariance matrix of the proposal distribution. If unsupplied, the inverse of the hessian at the mode is used as default.
jumpscale: how much to scale the jumps in the proposal distribution. This can be used to tune the acceptance rate. Setting jumpscale higher will tend to lower the acceptance rate.
Returns an instance of type RWMHEstimation, which can then be used to run the actual estimation.
"""
function RWMHEstimation(model::DSGEModel,parameters::Union{String,Array{String,1}},priors::Dict,observed::Union{String,Array{String,1}},observedseries::Array{Float64};nchains::Int=2,chainlength::Int=10000,burnin::Int=1000,proposalcovariance::Array{Float64,2}=zeros(Float64,0,0),jumpscale::Float64=0.2)
    if isa(parameters,String)
        parameters = [parameters]
    end
    if isa(observed,String)
        observed = [observed]
    end
    # This is the vector of parameters we're tying to estimate
    parameterindices = getparameterindices(model,parameters)
    
    if length(priors)!=length(parameters)
        error("Number of priors ($(length(priors))) does not match number of parameters to be estimated ($(length(parameters))).")
    end
    priors_ordered,priors_lambda = lambdifypriors(parameters,priors)
    
    if length(observed) > model.meta.numshocks
        error("Stochastic singularity. You cannot have more observed than the number of shocks.")
    end
    
    return RWMHEstimation(parameters,priors_ordered,priors_lambda,parameterindices,observed,observedseries,nchains,chainlength,burnin,proposalcovariance,jumpscale)
end

"""
   RWMHEstimation(model::DSGEModel,parameters::Union{String,Array{String,1}},priors::Dict,observedseries::DataFrames.DataFrame;nchains::Int=2,chainlength::Int=10000,burnin::Int=1000,proposalcovariance::Array{Float64,2},jumpscale::Float64=0.2)
Define the settings for a random walk metropolis hastings estimation.
parameters are the parameters you wish to estimation. priors is a dictionary with keys of the parameters you are estimating (as strings) and values as the prior distribution. Using distributions from the Distributions package is supported (and strongly recommended). You can, alternatively, supply your own anonymous function of your desired PDF. observedseries is a DataFrame with column names equal to which variables in your model are observed.
Optional keyword arguments are:
nchains: the number of chains used to sample from the prior
chainlength: the length of each chain
burnin: how much to drop from the start of the chain as burn in
proposalcovariance: the covariance matrix of the proposal distribution. If unsupplied, the inverse of the hessian at the mode is used as default.
jumpscale: how much to scale the jumps in the proposal distribution. This can be used to tune the acceptance rate. Setting jumpscale higher will tend to lower the acceptance rate.
Returns an instance of type RWMHEstimation, which can then be used to run the actual estimation.
"""
function RWMHEstimation(model::DSGEModel,parameters::Union{String,Array{String,1}},priors::Dict,observedseries::DataFrames.DataFrame;kwargs...)
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
    
    return RWMHEstimation(model,parameters,priors,observed,observeddata;kwargs...)
end

function RWMHChain(estimation::RWMHEstimation,ll::Function,H::Array{Float64,2},mode::Array{Float64,1},c::Int)
    # This function runs the actual posterior sampler
    draws = Array(Float64,length(estimation.parameters),estimation.chainlength-estimation.burnin)
    # Preallocate the random acceptance thresholds
    acceptthresholds = rand(Float64,estimation.chainlength)
    # Preallocate the jumps
    jumpD = Distributions.MvNormal(H)
    jumps = Distributions.rand(jumpD,estimation.chainlength)
    current = mode
    proposal = similar(current)
    llcurrent = ll(current)
    acceptances = 0
    ProgressMeter.@showprogress 1 "Sampling from posterior, chain $(c) of $(estimation.nchains)..." for i in 1:estimation.chainlength
        proposal[:] = current[:] + estimation.jumpscale*jumps[:,i]
        llproposal = ll(proposal)
        alpha = exp(llproposal-llcurrent)
        if alpha > acceptthresholds[i]
            current,proposal = proposal,current
            llcurrent = llproposal
            acceptances += 1
        end
        if i > estimation.burnin
            draws[:,i-estimation.burnin] = current
        end
    end
    return RWMHChain(mean(draws,2),median(draws,2),draws,acceptances/estimation.chainlength)
end

"""
    function estimation(estimation::RWMHEstimation,model::DSGEModel;laypunovmaxattempts=10000)
Runs the RWMH estimation for a previously defined estimation (which sets the estimation parameters).
Returns an instance of type RWMHResults which contains the full posterior, and summary statistics of.
"""
function estimate(estimation::RWMHEstimation,model::DSGEModel;laypunovmaxattempts=10000)
    MM = observationmatrix(model,estimation)

    estimationoptions = KalmanOptions(MM,laypunovmaxattempts)
    
    # Construct the likelihood function
    ll = x -> begin
        priors = evaluatepriors(estimation.priors_lambda,x)
        if isfinite(priors)
            l,_,_ = kalmanloglikelihood(estimationoptions,model,x,estimation.parameterindices,estimation.observedseries,false)
            return priors + l
        else
            return priors
        end
    end
    
    H,mode = findmode_MH(estimation,ll)

    chains = Array(RWMHChain,estimation.nchains)
    for c in 1:estimation.nchains
        chains[c] = RWMHChain(estimation,ll,H,mode,c)
    end
    # Now stitch the chains together to create aggregate posterior
    posterior = Array(Float64,length(estimation.parameters),estimation.nchains*(estimation.chainlength-estimation.burnin))
    for c in 1:estimation.nchains
        posterior[:,((c-1)*(estimation.chainlength-estimation.burnin)+1):c*(estimation.chainlength-estimation.burnin)] = chains[c].draws
    end
    
    return RWMHEstimationResults(Dict(zip(estimation.parameters,mean(posterior,2))),Dict(zip(estimation.parameters,median(posterior,2))),posterior,chains)
end


function findmode_MH(estimation,ll;NLopt_xtol_rel=1e-6)
    ll_wgrad = (x,_)->ll(x)
    # Find the mode
    opt = NLopt.Opt(:LN_COBYLA, length(estimation.parameters))
    NLopt.max_objective!(opt, ll_wgrad)
    NLopt.xtol_rel!(opt,NLopt_xtol_rel)
    # set startpoint to the mean of the priors
    startpoint = Array(Float64,length(estimation.parameters))
    lbd = similar(startpoint)
    ubd = similar(startpoint)
    for i in eachindex(startpoint)
        if isa(estimation.priors[i],Distributions.Distribution)
            startpoint[i] = Distributions.mean(estimation.priors[i])
            lbd[i] = Distributions.minimum(estimation.priors[i]) + eps(Distributions.minimum(estimation.priors[i]))
            ubd[i] = Distributions.maximum(estimation.priors[i]) - eps(Distributions.maximum(estimation.priors[i]))
        else
            warn("Finding the mean and support of user-supplied log PDF is not supported. Using 1 as start point for mode finding and allowing it to be unbounded.")
            startpoint[i] = 1.0
            lbd[i] = -Inf
            ubd[i] = Inf
        end
    end
    NLopt.lower_bounds!(opt, lbd)
    NLopt.upper_bounds!(opt, ubd)
    _,mode,_ = NLopt.optimize(opt, startpoint)
    if isempty(estimation.proposalcovariance)
        # Take the hessian as the covariance matrix of the proposal distributon
        H = -inv(Calculus.hessian(ll, mode))
        if any(H.==Inf) || any(H.==-Inf) || any(isnan(H))
            println("Could not invert hessian.")
            println(Calculus.hessian(ll, mode))
            println(mode)
            error("Cannot continue.")
            #warn("Hessian could not be inverted at mode. Using identity covariance matrix instead.")
            #H = eye(Float64,size(H)...)
        end
        if !issymmetric(H)
            warn("Inverse of hessian at mode is not symmetric.\nInverse: $H \nHessian: $(Calculus.hessian(ll, mode))")
            warn("This is due to a known, but unresolved, issue in Calculus (#91).")
            warn("Proceeding by splitting the difference (H = (H + H')/2)...")
            H = (H + H')/2
        end
    else
        H = estimation.proposalcovariance
    end
    return H,mode
end

function evaluatepriors(priors_lambda,x)
    ll = 0
    for i in eachindex(x)
        ll += priors_lambda[i](x[i])
    end
    return ll
end

function lambdifypriors(parameters,priors)
    priors_lambda = Array(Function,length(parameters))
    priors_ordered = Array(Any,length(parameters))
    for (i,name) in enumerate(parameters)
        if name in keys(priors)
            if isa(priors[name],Distributions.Distribution)
                priors_lambda[i] = x->Distributions.logpdf(priors[name],x)
            elseif isa(priors[name],Function)
                priors_lambda[i] = priors[name]
            end
            priors_ordered[i] = priors[name]
        else
            error("You have not set a prior for $(name).")
        end
    end
    return priors_ordered,priors_lambda
end
