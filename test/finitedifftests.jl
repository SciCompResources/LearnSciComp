@info "Starting finite difference tests"
using LearnSciComp

@testset "weights from fornberg(1998)" begin
    order = 2; z = 0; x = [-1, 0, 1.0];
    @test fornberg(order, z, x) == [1,-2,1]  # central difference of second-derivative with unit-step

    order = 1; z = 0; x = [-1., 1.0];
    @test fornberg(order, z, x) == [-0.5,0.5] # central difference of first-derivative with unit step

    order = 1; z = 0; x = [0, 1];
    @test fornberg(order, z, x) == [-1, 1] # forward difference

    order = 1; z = 1; x = [0, 1];
    @test fornberg(order, z, x) == [-1, 1] # backward difference

    order = 3; z = 0; x = [-2,-1,0,1,2]   # third derivative central-diff
    @test fornberg(order, z, x) == [-1/2, 1, 0, -1, 1/2]

    # forward-diff of third derivative with order of accuracy == 3
    order = 3; z = 0; x = [0,1,2,3,4,5]
    @test fornberg(order, z, x) == [-17/4,	71/4	,−59/2,	49/2,	−41/4,	7/4]
end
