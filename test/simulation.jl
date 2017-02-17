function checksimulation(model)
    simulate(model,100)
    return true
end

@test checksimulation(threeEQNK)
@test checksimulation(RBC)