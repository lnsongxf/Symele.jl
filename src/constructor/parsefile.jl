"""
   DSGEModel(file;kwargs...)
Construct a DSGE model from a model file. Pass in the path as a string to the file location. All of the kwargs that are valid for the DSGEModel constructor are valid here. See that help file.
For information on how to write a model file, see the model file documentation, or the examples contained in examples/
Returns your model as an instance of type DSGEModel.
"""
function DSGEModel(file::String;kwargs...)

    f = open(file)
	modelfilestring = readlines(f)
    close(f)

	# Find the indexes of important things
	ends = Bool[strip(line) == "end" for line in modelfilestring]
    ends = (1:length(modelfilestring))[ends]
	parameters = findline("parameters:",modelfilestring)
	parametervalues = findline("parametervalues:",modelfilestring)
	shocks = findline("shocks:",modelfilestring)
	vars = findline("vars:",modelfilestring)
	equations = findline("equations:",modelfilestring)
	ssguess = findline("ssguess:",modelfilestring)
	
    equationend = minimum(ends[ends.>equations])
    ssguessend = minimum(ends[ends.>ssguess])
    valuesend = minimum(ends[ends.>parametervalues])
	
	vars = parsevars(modelfilestring[vars])
	shocks = parseshocks(modelfilestring[shocks])
	parameters = parseparameters(modelfilestring[parameters])
	parametervalues_dict = splitonequals(modelfilestring[(parametervalues+1):(valuesend-1)])
	equations = parseequations(modelfilestring[(equations+1):(equationend-1)])
	ssguess_dict = splitonequals(modelfilestring[(ssguess+1):(ssguessend-1)])

    # Now order the parameter values and ssguess
	parametervalues = Array(Float64,length(parameters))
	unset = trues(length(parameters))
	for (key,val) in parametervalues_dict
		if countnz(parameters.==key) == 0
			error("You have set a parameter value for $key, but it is not a declared parameter.")
        else
            parametervalues[parameters.==key] = val
			unset[parameters.==key] = false
        end
	end
	if any(unset)
		error("You have not set a value for parameter(s): $(parameters[unset]).")
	end
	
	ssguess = zeros(Float64,length(vars)) # Default to zero if unset
	for (key,val) in ssguess_dict
		if countnz(vars.==key) == 0
			error("You have set a steady state guess for $key, but it is not a declared variable.")
        else
            ssguess[vars.==key] = val
        end
	end
	
    return DSGEModel(vars,shocks,parameters,equations,parametervalues,ssguess;kwargs...)
end

function findline(keyword,modelfilestring)
	indices = Bool[startswith(strip(line),keyword) for line in modelfilestring]
	if countnz(indices) != 1
		error("Model file contains $(countnz(indices)) occurences of the $keyword keyword. It must appear only once.")
	end
	
	return (1:length(modelfilestring))[indices][1] # The extra one gets the scalar, rather than a 1 element array
end

function parsevars(vars)
	vars = strip(vars[6:end]) # Chop off the keyword and whitespace
	return splitoncomma(vars)
end

function parseshocks(shocks)
	shocks = strip(shocks[8:end]) # Chop off the keyword and whitespace
	return splitoncomma(shocks)
end

function parseparameters(parameters)
	parameters = strip(parameters[12:end]) # Chop off the keyword and whitespace
	return splitoncomma(parameters)
end

function splitonequals(input)
	out = Dict{String,Float64}()
	for i in eachindex(input)
		pair = split(input[i],"=")
		if length(pair)!=2
			error("Error reading line '$(in[i])'; wrong number of equals signs.")
        else
            out[strip(pair[1])] = parse(Float64,strip(pair[2]))
        end
	end
	return out
end

function parseequations(equations)
	for i in eachindex(equations)
		equations[i] = strip(equations[i])
	end
	return equations
end


function splitoncomma(in)
	list = split(in,",")
	for i in eachindex(list)
		list[i] = strip(list[i])
	end
	return list
end
