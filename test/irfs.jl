function checkirfs(model)
    irf(model;periods=20,graph=true)
    return true
end

@test checkirfs(RBC)
AeA = irf(RBC,["A"],["e_A"];periods=20)
correct = Array(Float64,20,1)
correct[1] = 1.
for i = 2:20
    correct[i] = RBC.parameters.values[RBC.parameters.dictionary["rho"]]*correct[i-1]
end
@test AeA["e_A"]["A"] ≈ correct