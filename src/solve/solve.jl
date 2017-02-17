"""
    function solve(model::DSGEModel,quiet=false,permitsunspots=false;method="default")
Find the solution to the model (at the parameters as currently set). quiet disables any form of output from the solver (other than errors); this is useful when you are repeatedly solving the model, such as during estimation. permitsunspots allows sunspots (only in qzsolve). THIS SHOULD BE DISABLED. method allows you to choose how to solve the model. Possible options: "default" uses the method chosen by the model; "qz" uses the qz decomposition method.
Returns an instance of type DSGESolution.
"""
function solve(model::DSGEModel,quiet=false,permitsunspots=false;method="default")
    if !quiet
        println("Solving model to order $(model.meta.solutionorder)")
    end
    if method == "default"
        method = model.meta.solutionmethod
    end
    # Farm out to the different possible methods here I guess?
    if method == "qz" && model.meta.solutionorder == 1
        return qzsolve(model,permitsunspots)
    else
        error("Unknown solution method $(method).")
    end
end

function createarglist(model)
    arglist = vcat(model.steadystate.values,model.parameters.values)
    return arglist
end

function evaluatematrix(model,arglist,functionmatrix)
    if isa(functionmatrix,SparseMatrixCSC)
        out = spzeros(Float64,size(functionmatrix)...)
    else
        out = Array(Float64,size(functionmatrix)...)
    end
    for i in eachindex(functionmatrix,out)
        out[i] = functionmatrix[i](arglist...)
    end
    return out
end

function qzsolve(model,permitsunspots)
    arglist = createarglist(model)
    FYp = evaluatematrix(model,arglist,model.perturbationsystem_lambda.firstorder.FYp)
    FXp = evaluatematrix(model,arglist,model.perturbationsystem_lambda.firstorder.FXp)
    FY = evaluatematrix(model,arglist,model.perturbationsystem_lambda.firstorder.FY)
    FX = evaluatematrix(model,arglist,model.perturbationsystem_lambda.firstorder.FX)
    FXm = evaluatematrix(model,arglist,model.perturbationsystem_lambda.firstorder.FXm)
    FE = evaluatematrix(model,arglist,model.perturbationsystem_lambda.firstorder.FE)
    
    nfl = countnz(model.variables.forwardlooking)
    ny = model.meta.numvars - model.meta.numstates
    nt = model.meta.numvars + nfl
    nfly = countnz(!model.variables.isstate & model.variables.forwardlooking)
    
    LHS = zeros(Float64,nt,nt)
    RHS = zeros(Float64,nt,nt)
    ME = zeros(Float64,nt,model.meta.numshocks)
    MD = vcat(zeros(Float64,model.meta.numvars,nfl),eye(Float64,nfl,nfl))
    
    LHS[1:model.meta.numvars,1:ny] = FY[:,:]
    LHS[1:model.meta.numvars,(ny+1):(model.meta.numvars)] = FX[:,:]
    LHS[1:model.meta.numvars,(model.meta.numvars+1):(model.meta.numvars+nfly)] = FYp[:,model.variables.forwardlooking[!model.variables.isstate]]
    LHS[1:model.meta.numvars,(model.meta.numvars+nfly+1):end] = FXp[:,model.variables.forwardlooking[model.variables.isstate]]
    LHS[(model.meta.numvars+1):end,1:nfl] = eye(Float64,nfl,nfl)
    
    RHS[1:model.meta.numvars,(1:model.meta.numvars)[model.variables.isstate]] = -FXm[:,:]
    RHS[(model.meta.numvars+1):end,(model.meta.numvars+1):end] = eye(Float64,nfl,nfl)
    
    ME[1:model.meta.numvars,:] = -FE[:,:]
    
    F = schurfact(RHS,LHS)
    # schurfact gives: LHS = F[:Q]*F[:T]*F[:Z]'; RHS = F[:Q]*F[:S]*F[:Z]'
    # This is slightly different (transposes differ) to matlab's way of doing it
    unstableroots = abs(F[:alpha]./F[:beta]) .> 1
    no = countnz(unstableroots)
    ni = nt - no # number of stable roots
    if no == nfl
        ordschur!(F, !unstableroots)
        T11 = F[:T][1:ni,1:ni]
        S11 = F[:S][1:ni,1:ni]
        S12 = F[:S][1:ni,(ni+1):end]
        S22 = F[:S][(ni+1):end,(ni+1):end]
        Q1 = F[:Q][:,1:ni]'
        Q2 = F[:Q][:,(ni+1):end]'
        Z1 = F[:Z][1:ni,:]
        U1,S,V1 = svd(Q2*MD)
        D11 = diagm(S)
        capPsi = Q1*MD*V1*(D11\(U1'))
        MX = Z1*vcat(hcat(T11\S11,T11\(S12 - capPsi*S22)), zeros(Float64,no,nt))*Z1'
        ME = Z1*vcat(T11\(Q1 - capPsi*Q2),zeros(no,nt))*ME # * covariance matrix
        MC = reshape(copy(model.steadystate.values),model.meta.numvars,1)
    elseif no > nfl
        error("Model is explosive ($(no) explosive eigenvalues for $(nfl) forward looking variables).")
    else
        if permitsunspots
            error("Sunspots not yet implemented.")
        else
            error("Model is indeterminate ($(no) explosive eigenvalues for $(nfl) forward looking variables).")
        end
    end
        
    
    return DSGESolution(1,MX,ME,MC)
end
