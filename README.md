# Symele

[![Build Status](https://travis-ci.org/angusmoore/Symele.jl.svg?branch=master)](https://travis-ci.org/angusmoore/Symele.jl)

[![Coverage Status](https://coveralls.io/repos/angusmoore/Symele.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/angusmoore/Symele.jl?branch=master)

[![codecov.io](http://codecov.io/github/angusmoore/Symele.jl/coverage.svg?branch=master)](http://codecov.io/github/angusmoore/Symele.jl?branch=master)

Symele is a symbolic modelling engine for DSGE models. It aims to automate the most tedious and error-prone part of solving DSGE models by perturbation methods: taking the perturbation. In that sense, it is very similar to existing tools such as Dynare. However, it is built with a more interactive workflow in mind than Dynare: models are ojbects, many of which can be kept in memory at once. This design allows you to experiment and test the parameters of your model interactively, at Julia's REPL, rather than executing a static 'model file'.

