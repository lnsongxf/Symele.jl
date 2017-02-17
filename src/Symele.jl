module Symele

# Required packages
import Base.copy
import SymPy
import SymPy: Sym, jacobian # Used a lot, and doesn't lead to conflict
import Calculus
import DataFrames
import Distributions
import NLopt
import NLsolve
import ProgressMeter
import PyPlot

include("constructor/modeltypes.jl")
include("constructor/constructor.jl")
include("constructor/parsefile.jl")
include("constructor/perturbation.jl")
include("modelutils.jl")
include("solve/solutiontypes.jl")
include("solve/solve.jl")
include("simulation/simulate.jl")
include("estimation/utils.jl")
include("estimation/kalmanfilter.jl")
include("estimation/mle.jl")
include("estimation/rwmh.jl")
include("estimation/getstatesandshocks.jl")
include("estimation/shockdecomposition.jl")
include("irfs/irf.jl")

include("estimation/hpfilter.jl")

# Functions to be exported
export DSGEModel, solve, estimate, MLEEstimation, RWMHEstimation, hpfilter, irf, graphirf, updateparameters!, simulate, shockdecomposition, getstatesandshocks

end # module
