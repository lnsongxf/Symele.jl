function observationmatrix(model,estimation)
    MM = zeros(Float64,length(estimation.observed),model.meta.numvars)
    for (i,name) in enumerate(estimation.observed)
        if name in keys(model.variables.dictionary)
            MM[i,model.variables.dictionary[name]] = 1.0
        else
            error("$(name) cannot be observed variable because it is not a variable name in your model.")
        end
    end
    return MM
end
