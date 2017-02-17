function subSSsym(vars,expr)
    subsargs = pairedtuplearray(vcat(vars.lead_syms,vars.state_syms),vcat(vars.contemp_syms,vars.contemp_syms[vars.isstate]))
    return SymPy.subs(expr,subsargs...)
end

function lambdify(vars,parameters,a::Array{SymPy.Sym})
    out = Array(Function,size(a)...)
    for i in eachindex(a,out)
        # First sub out leads and lags for contemps
        expr = subSSsym(vars,a[i])
        out[i] = SymPy.lambdify(expr,[vars.contemp_syms;parameters.syms])
    end
    return out
end

function lambidfy(vars,parameters,pert::NthOrderPerturbation{SymPy.Sym})
    FYp = lambdify(vars,parameters,pert.FYp)
    FXp = lambdify(vars,parameters,pert.FXp)
    FY = lambdify(vars,parameters,pert.FY)
    FX = lambdify(vars,parameters,pert.FX)
    FXm = lambdify(vars,parameters,pert.FXm)
    FE = lambdify(vars,parameters,pert.FE)
    return NthOrderPerturbation(FYp,FXp,FX,FY,FE)
end

function lambdify(vars,parameters,meta,pert::PerturbationSystem{SymPy.Sym})
    FYp = lambdify(vars,parameters,pert.firstorder.FYp)
    FXp = lambdify(vars,parameters,pert.firstorder.FXp)
    FY = lambdify(vars,parameters,pert.firstorder.FY)
    FX = lambdify(vars,parameters,pert.firstorder.FX)
    FXm = lambdify(vars,parameters,pert.firstorder.FXm)
    FE = lambdify(vars,parameters,pert.firstorder.FE)
    firstorder = FirstOrderPerturbation(FYp,FXp,FX,FY,FXm,FE)
    higherorders = Array(NthOrderPerturbation{Function},meta.solutionorder-1)
    for i in eachindex(higherorders)
        higherorders[i] = lambidfy(vars,parameters,pert.higherorders[i])
    end
    return PerturbationSystem(firstorder,higherorders)
end

function PerturbationSystem(meta,parameters,vars,shocks,equations,steadystate)
    firstorder = FirstOrderPerturbation(meta,equations,vars,shocks)
    higherorders = Array(NthOrderPerturbation{SymPy.Sym},meta.solutionorder)
    
    pert_sym = PerturbationSystem(firstorder,higherorders)
    pert_lambda = lambdify(vars,parameters,meta,pert_sym)
    return pert_sym,pert_lambda
end

function FirstOrderPerturbation(meta::Meta,equations::Equations,vars::Variables,shocks::Shocks)
    FYp = jacobian(equations.syms,vars.lead_syms[!vars.isstate])
    FXp = jacobian(equations.syms,vars.lead_syms[vars.isstate])
    FX = jacobian(equations.syms,vars.contemp_syms[vars.isstate])
    FY = jacobian(equations.syms,vars.contemp_syms[!vars.isstate])
    FXm = jacobian(equations.syms,vars.state_syms)
    FE = jacobian(equations.syms,shocks.syms)
    
    return FirstOrderPerturbation(FYp,FXp,FX,FY,FXm,FE)
end