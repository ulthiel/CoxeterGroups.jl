using Test
using CoxeterGroups

@testset "Basic Coxeter system operations" begin
    # Create a system of type H3 x A2
    coxeter_mat = [
        1 5 2 2 2
        5 1 3 2 2
        2 3 1 2 2
        2 2 2 1 3
        2 2 2 3 1
    ]
    H3A2 = coxeter_system(coxeter_mat)

    # Test usual Coxeter system accessors
    @test coxeter_name(H3A2) == "H3 x A2"
    @test coxeter_matrix(H3A2) == coxeter_mat
    @test rank(H3A2) == 5
    @test is_irreducible(H3A2) == false
    @test is_finite_type(H3A2) == true
    @test is_affine_type(H3A2) == false
    @test order(H3A2) == 120 * 6
    @test number_of_reflections(H3A2) == 18
    @test degrees(H3A2) == [2, 2, 3, 6, 10]     # Should be sorted
    @test exponents(H3A2) == [1, 1, 2, 5, 9]    # Should be sorted

    # Shouldn't be able to take the Coxeter number of a reducible system.
    @test_throws ErrorException coxeter_number(H3A2)
end
