function testdecomp(model,data)
    shockdecomposition(model,data)
    return true
end

@test testdecomp(threeEQNK,NKdata)
@test testdecomp(RBC,RBCdata)


NKshocks = load("data/NKdata.jld")
NKactualshocks = NKshocks["actualshocks"]
RBCshocks = load("data/RBCdata.jld")
RBCactualshocks = RBCshocks["actualshocks"]

_,estshocks = getstatesandshocks(threeEQNK,NKdata)
for (i,name) in enumerate(threeEQNK.shocks.names)
    @test estshocks[name] == NKactualshocks[:,i]
end

_,estshocks = getstatesandshocks(RBC,RBCdata)
for (i,name) in enumerate(RBC.shocks.names)
    @test estshocks[name] == RBCactualshocks[:,i]
end
