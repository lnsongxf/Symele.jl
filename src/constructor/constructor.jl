function sympify_catch(exp,eqno)
    try
        return Expr = Sym(exp)
    catch err
        if string(err.val) == "PyObject SympifyError()"
            # This is a terrible way to do this, but I can't figure a better way...
            println("Syntax in equation $(eqno) is invalid.")
            rethrow(err)
        else
            println("Error sympifying equation $(eqno):")
            rethrow(err)
        end
    end
end

function pairedtuplearray{T}(a::Array{T,1},b::Array{T,1})
    if length(a)!=length(b)
        error("Arrays must have same length. ($(length(a)) versus $(length(a)) elements)")
    end
    out = Array(Tuple{T,T},length(a))
    for i in eachindex(out,a,b)
        out[i] = (a[i],b[i])
    end
    return out
end
    
function checksymbols(equations,vars,shocks,parameters)
    # This function checks all the expressions for symbols that are unrecognised.
    # Catching this early is helpful for avoiding errors later.
    validsymbols = vcat(vars.lead_syms,shocks.syms,vars.contemp_syms,vars.state_syms,parameters.syms[:])
    for (i,expression) in enumerate(equations.syms)
        symbollist = SymPy.free_symbols(expression)
        for symbol in symbollist
            if !(symbol in validsymbols)
                error("Unknown symbol $(symbol) in equation $(i).")
            end
        end
    end
    return
end

function createsymversions(vars,suffix)
    out = Array(SymPy.Sym,length(vars))
    for (i,var) in enumerate(vars)
        try
            out[i] = Sym(string(var,suffix))
        catch
            error("Unable to convert $(string(var,suffix)) to a symbolic expression.")
        end
    end
    return out
end    
# Create method to allow call without optional suffix
createsymversions(vars) = createsymversions(vars,"")

function checknames(names)
    # This function throws an error for names that are unacceptable/reserved in sympy (i.e. pi, inf etc)
    badnames = Array(String,0)
    push!(badnames,"Inf")
    push!(badnames,"inf")
    push!(badnames,"exp")
    push!(badnames,"pi")
    push!(badnames,"log")
    push!(badnames,"beta")
    push!(badnames,"lambda")
    push!(badnames,"zeta")

    for name in names
        if any(badnames.==name)
            error("Name: $(name) is a reserved name in SymPy. Please change it.")
        end
    end
    return 
end

function sympifyallexpressions(expressions,leads_sym,vars)
    warn("I am ugly. Make me better")
    expressions_sym = Array(SymPy.Sym,length(expressions))

    lags = createsymversions(vars,"_lag") # Create a list of all lagged variable - these are _potential_ states
    writtenlags = createsymversions(vars,"(-1)") # Preallocating this saves time, by not having to repeatedly call Sym
    writtenleads = createsymversions(vars,"(+1)") # As above

    for (eqno,exp) in enumerate(expressions)
        Expr = sympify_catch(exp,eqno)
        expressions_sym[eqno] = Expr
    end
    
    subsargs = pairedtuplearray(vcat(writtenlags,writtenleads),vcat(lags,leads_sym))
    if !isempty(subsargs)
        expressions_sym = SymPy.subs(expressions_sym,subsargs...)
    end
    
    isstate = falses(length(lags))
    
    for expression in expressions_sym
        isstate = isstate | expressioncontains(expression,lags)
    end
    
    # Now create the states array
    numstates = countnz(isstate)
    states_sym = Array(SymPy.Sym,numstates)
    counter = 0
    for (i,state) in enumerate(lags)
        if isstate[i]
            counter += 1
            states_sym[counter] = state
        end
    end

    return expressions_sym,isstate,states_sym
end

function evaluateSS(ssfunction,varvalues,parametervalues)
    out = Array(Float64,length(ssfunction))
    for (i,f) in enumerate(ssfunction)
        out[i] = f(varvalues...,parametervalues...)
    end
    return out
end

function lambdifyss(ss,vars,parameters)
    # Takes in sympy expressions, and symbolic object versions of the variables and parameters
    # Converts to a function that takes in arrays of variable values and parameter values
    # and returns equation errors
    functionarray = Array(Function,length(ss))
    for (i,expression) in enumerate(ss)
        functionarray[i] = SymPy.lambdify(expression,[vars;parameters])
    end
    return functionarray
end

function expressioncontains(expression::SymPy.Sym,symbol::SymPy.Sym)
    if symbol in SymPy.free_symbols(expression)
        return true
    end
    return false
end

function expressioncontains(expression::SymPy.Sym,symbols::Array{SymPy.Sym})
    count = 0
    isin = falses(length(symbols))
    for (i,symbol) in enumerate(symbols)
        if symbol in SymPy.free_symbols(expression)
            isin[i] = true
        end
    end
    return isin
end

function Parameters(names,values)
    if length(names)!=length(values)
        error("You have $(length(names)) parameters, but $(length(values)) parameter values.")
    end
    return Parameters(names,values,Dict(zip(names,1:length(names))),createsymversions(names))
end

function Shocks(names)
    dictionary = Dict(zip(names,1:length(names)))
    syms = createsymversions(names)
    return Shocks(names,dictionary,syms)
end

function rearrangeequations(equations)
    out = similar(equations)
    for (i,eq) in enumerate(equations)
        numequals = length(matchall(r"=", eq))
        if numequals > 1
            error("You have more than 1 equals sign in equation $(i).")
        elseif numequals == 1
            index = search(eq,'=')
            out[i] = string(eq[1:(index-1)]," - (", eq[(index+1):end] , ")")
        else
            out[i] = eq
        end
    end
    return out
end

function findforwardlooking(eqsyms,leads_sym)
    isfl = falses(length(leads_sym))
    for (i,var) in enumerate(leads_sym)
        for expression in eqsyms
            if !isfl[i] && expressioncontains(expression,var)
                isfl[i] = true
            end
        end
    end
    return isfl
end    

function parseequationsandvars(vars,equations)
    if length(equations)!=length(vars)
        error("You have $(length(equations)) equations for $(length(vars)) endogenous variables.")
    end
    
    # Define arrays for all the symbolic objects, except states, because we don't know them yet
    leads_sym = createsymversions(vars,"_lead")
    contemp_sym = createsymversions(vars)
    
    # Check the equations have zero or one equals signs. If one, convert to a one-sided expressions
    zeroedequations = rearrangeequations(equations)
    
    # Sympify the expressions and detect state variables while doing so
    eqsyms,var_isstate,states_sym = sympifyallexpressions(zeroedequations,leads_sym,vars)
    # Find which vars are forwardlooking
    var_isforward = findforwardlooking(eqsyms,leads_sym)
    
    equations = Equations(equations,eqsyms)
    # Re order to be [Y X]
    vars = vcat(vars[!var_isstate],vars[var_isstate])
    leads_sym = vcat(leads_sym[!var_isstate],leads_sym[var_isstate])
    contemp_sym = vcat(contemp_sym[!var_isstate],contemp_sym[var_isstate])
    var_isforward = vcat(var_isforward[!var_isstate],var_isforward[var_isstate])
    var_isstate = vcat(var_isstate[!var_isstate],var_isstate[var_isstate])
    
    var_dictionary = Dict(zip(vars,1:length(vars)))
    
    vars = Variables(vars,var_dictionary,var_isstate,states_sym,leads_sym,contemp_sym,var_isforward)
    return equations,vars
end

function SteadyState(meta::Meta,equations::Equations,vars::Variables,parameters::Parameters,shocks::Shocks,ss_guess::Array{Float64,1})
    if length(ss_guess) != meta.numvars
        error("Your supplied guess for the steady state has $(length(ss_guess)) entries, but you have $(meta.numvars) variables.")
    end
    # Need to reorder the steady state guess, becuase the variables get reordered
    ss_guess = vcat(ss_guess[!vars.isstate],ss_guess[vars.isstate])
    
    subsargs = pairedtuplearray(vcat(vars.lead_syms,vars.state_syms,shocks.syms),vcat(vars.contemp_syms,vars.contemp_syms[vars.isstate],zeros(Float64,meta.numshocks)))
    expressionsSS = SymPy.subs(equations.syms,subsargs...)
    lambdifiedSS = lambdifyss(expressionsSS,vars.contemp_syms,parameters.syms)
    
    return SteadyState(expressionsSS,lambdifiedSS,ss_guess,false,x->x)
end

function reorderss!(x,isstate)
    x[:] = vcat(x[!isstate],x[isstate])
    return x
end

function SteadyState(variables,ssfunction,ordered,ss_guess)
    if !ordered
        orderedssfunction = x->begin
            ss = ssfunction(x)
            reorderss!(ss,variables.isstate)
            return ss
        end
        return SteadyState(Array(SymPy.Sym,0),Array(Function,0),ss_guess,true,orderedssfunction)
    else
        return SteadyState(Array(SymPy.Sym,0),Array(Function,0),ss_guess,true,ssfunction)
    end
end

"""
    function DSGEModel(var_names,shock_names,parameter_names,equations,parameter_values,ss_guess=zeros(length(var_names));solutionorder::Int=1,solutionmethod::String="qz",ssfunction=false,ssfunctionordered::Bool=false)
Constructor for DSGE model. var_names, shocks_names, parameter_names and equations are arrays of strings. Parameter values are a numeric array. ss_guess is an optional guess at the model's steady state; having a good guess can help the numerical solver converge. If unspecified it defaults to all zeros. The guess should be in the same order as the variables in var_names. ssfunction is an optional alternative to using the built in steady state solver (this is useful in cases where the steady state has an analytical solution for instance). It should take as arguments a vector of parameter values and return a vector corresponding to the steady state (in the same order as var_names). Alternatively, it can return the steady state variables in the order that the model object re-orders the variables to by setting ssfunctionordered=true. (This requires you to know the variable order in advance, and is for advanced use only.
Returns your model as an instance of type DSGEModel. 
"""
function DSGEModel(var_names,shock_names,parameter_names,equations,parameter_values,ss_guess=zeros(length(var_names));solutionorder::Int=1,solutionmethod::String="qz",ssfunction=false,ssfunctionordered::Bool=false)
    # Some useful constants first
    numvars = length(var_names)
    numshocks = length(shock_names)
    numparameters = length(parameter_names)
    
    # Check all the names are ok
    checknames(var_names)
    checknames(shock_names)
    checknames(parameter_names)
    
    shocks = Shocks(shock_names)
    equations,vars = parseequationsandvars(var_names,equations)
    meta = Meta(countnz(vars.isstate),numvars,numshocks,numparameters,solutionorder,solutionmethod)
    parameters = Parameters(parameter_names,parameter_values)
    # Check for unknown symbols in any of the equations
    checksymbols(equations,vars,shocks,parameters)   

    # Create the steady state system
    if isa(ssfunction,Bool)
        steadystate = SteadyState(meta,equations,vars,parameters,shocks,ss_guess)
    elseif isa(ssfunction,Function)
        # user supplied steady state function
        steadystate = SteadyState(vars,ssfunction,ssfunctionordered,ss_guess)
    end
    # Create the perturbation system
    perturbationsystem_sym,perturbationsystem_lambda = PerturbationSystem(meta,parameters,vars,shocks,equations,steadystate)
    
    model = DSGEModel(meta,parameters,vars,shocks,equations,steadystate,perturbationsystem_sym,perturbationsystem_lambda)
    solveSS!(model)
    return model
end
