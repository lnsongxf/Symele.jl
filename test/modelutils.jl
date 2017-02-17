function testupdateparameters(model,p1,p2)
    updateparameters!(model,p1,p2)
    return true
end

@test testupdateparameters(threeEQNK,"betta",0.98)
@test testupdateparameters(threeEQNK,["betta","rhoR"],[0.99,0.9])
#@test Symele.getshockindices(threeEQNK,["e_A"]) == [2]
#@test Symele.getvariableindices(threeEQNK,["Y","inflation"]) == [1,3]