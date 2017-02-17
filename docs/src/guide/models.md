# Dealing with Models

## The DSGEModel object

The core of the Symele package is the DSGEModel object. This is a container for all of the information about your model: it's parameters, variables, shocks, equations and the perturbation approximation.

Once constructed a model is a variable like any other. You can copy it; have many in memory, corresponding to different things; update parameters on the fly; etc. This permits you to experiment with your model, interactively at the Julia REPL.

## Constructing a model

There are two ways to construct a model in Symele. The first uses a model file, which will be familiar to Dynare users. The file specifies parameters, parameter values, variables, shocks, the model's first order conditions and, optionally, a guess at the deterministic steady state. Symele take its from there. This is a natural and easy way to write down your model and have the computer understand it.

The second way is to use the more direct constructor. Here you instead pass Symele a series of arrays corresponding to parameters, parameter values etc. While somewhat less elegant than passing in a model file, it can be useful (for instance, constructing the model programatically where you have many sectors with the same first order conditions but different variable names).

The examples/ folder contains example model files. Model files are plain text files, and can have extension. To construct a model using a model file simply call (running this command requires your present working directory to be the Symele package folder, try cd(Pkg.dir("Symele"))):
```julia
RBC = DSGEModel("examples/RBC.txt")
'''



## Working with models



## Solving models

Solving a model is simple. Just call `solve(model)'. Currently, only the QZ decomposition solution method is implemented, so that will be used. In future, other solution methods might be made available, accessible via keyword argument.

`solve' returns a DSGESolution. This is a collection of matrices such that:
```julia
X_t = solution.MX*(X_{t-1} - SS) + solution.ME*E + solution.MC
'''
where `SS' is the deterministic steady state and `E' is the column vector of shocks in the model.

Higher-order solutions will be implemented in later versions.

As with model, you can have many solutions in memory. Once created, a solution corresponds to the model _as it was when you solved it_. If you subsequently update the model parameters, the solution does not change (you instead need to resolve the model).

This is useful for comparing solutions under different parameters, e.g.:
```julia
updateparameters!(model,"someparameter",0.9)
highparametersolution = solve(model)
updateparameters!(model,"someparameter",0.1)
lowparametersolution = solve(model)
'''

## IRFs and simulations

