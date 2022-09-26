using Test
using CoxeterGroups

# Check our predicates are working.
@test is_coxeter_matrix(Array{Int64}(undef, 0, 0))
@test is_coxeter_matrix([1;;])
@test !is_coxeter_matrix([2;;])


# Construct every element of the group by right-multiplying by generators.
# Used for some basic tests to see if the group is being constructed correctly.
function enumerate_whole_group(G::CoxGrp)
    id = one(G)
    seen = Set([id])
    queue = [id]
    gens = generators(G)

    i = 1
    while i <= length(queue)
        for s = gens
            ws = queue[i] * s
            if ws ∉ seen
                push!(seen, ws)
                push!(queue, ws)
            end
        end

        i += 1
    end

    return queue
end

# Test each of the different implementations
group_constructors = [
    coxeter_group_min,
    coxeter_group_recursive,
]

# For each implementation, check they generate the correct number of elements in the symmetric group.
for group_constructor = group_constructors
    for rank = 0:5
        type = coxeter_type("A", rank)
        W, gens = group_constructor(type.gcm)
        Welts = enumerate_whole_group(W)
        @test length(Welts) == type.weyl_order
    end
end

# Test basic multiplication in symmetric groups.
for n in 0:5
    W, gens = symmetric_group(n)
    Welts = enumerate_whole_group(W)
    @test length(Welts) == factorial(n)
end

# Check an implementation supports all operations appearing in the documentation (not necessarily that
# the implementation of these operations is correct).
groups = [
    coxeter_group_min(Array{Int64}(undef, 0, 0)),
    coxeter_group_min([2;;]),                       # A1, GCM
    coxeter_group_min([2 -1; -1 2]),                # A2, GCM
    coxeter_group_min([2 -1 0; -1 2 -1; 0 -1 2]),   # A3, GCM
    coxeter_group_min([1 3 2; 3 1 3; 2 3 1]),       # A3, Coxeter matrix
    coxeter_group_min([1 5 2; 5 1 3; 2 3 1]),       # H3, Coxeter matrix
    symmetric_group(1),
    symmetric_group(2),
    symmetric_group(3),
    symmetric_group(4),
]

@testset "Supported operations $W" for (W, gens) in groups
    @test gens == generators(W)
    @test rank(W) == length(gens)
    @test isone(one(W))
    if is_finite(W)
        w0 = longest_element(W)
        @test all([length(w0 * gen) < length(w0) for gen in gens])
    end

    id = one(W)
    @test W == parent(id)
    @test inv(id) == id
    @test length(id) == 0
    @test sign(id) == 1
    @test length(short_lex(id)) == 0
    @test length(inverse_short_lex(id)) == 0

    if rank(W) > 0
        coxelt = *(gens...)
        @test length(coxelt) == rank(W)
        @test is_left_descent(1, coxelt)
        @test is_right_descent(coxelt, rank(W))
        @test left_multiply(1, coxelt) == gens[1] * coxelt
        @test right_multiply(coxelt, 1) == coxelt * gens[1]
        @test short_lex(coxelt) == 1:rank(W)
        @test inverse_short_lex(coxelt) == reverse(short_lex(coxelt^-1))
        @test coxelt^5 == coxelt^2 * coxelt^3
        @test coxelt^-2 == inv(coxelt)^2
    end
end

matrices = Dict{String,Matrix{Int64}}()
matrices["A0"] = Array{Int64}(undef, 0, 0)
matrices["A1"] = [1;;]
matrices["A2"] = [1 3;3 1]
matrices["A3"] = [1 3 2;3 1 3;2 3 1]
matrices["A4"] = [1 3 2 2;3 1 3 2;2 3 1 3;2 2 3 1]
matrices["A~1"] = [1 0;0 1]
matrices["A~2"] = [1 3 3;3 1 3;3 3 1]
matrices["A~3"] = [1  3  2  3;3  1  3  2;2  3  1  3;3  2  3  1]
matrices["A~4"] = [1  3  2  2  3;3  1  3  2  2;2  3  1  3  2;2  2  3  1  3;3  2  2  3  1]
matrices["B2"] = [1  4;4  1]
matrices["B3"] = [1  3  2;3  1  4;2  4  1]
matrices["B4"] = [ 1  3  2  2;3  1  3  2;2  3  1  4;2  2  4  1]
matrices["B~1"] = [1  4;4  1]
matrices["B~2"] = [1  4  2;4  1  4;2  4  1]
matrices["B~3"] = [1  2  3  2;2  1  3  2;3  3  1  4;2  2  4  1]
matrices["B~4"] =[1  2  3  2  2;2  1  3  2  2;3  3  1  3  2;2  2  3  1  4;2  2  2  4  1]
matrices["C2"] = [1  4;4  1]
matrices["C3"] = [1  3  2;3  1  4;2  4  1]
matrices["C4"] = [1  3  2  2;3  1  3  2;2  3  1  4;2  2  4  1]
matrices["C~1"] = [1  4;4  1]
matrices["C~2"] = [1  4  2;4  1  4;2  4  1]
matrices["C~3"] = [1  4  2  2;4  1  3  2;2  3  1  4;2  2  4  1]
matrices["C~4"] = [1  4  2  2  2;4  1  3  2  2;2  3  1  3  2;2  2  3  1  4;2  2  2  4  1]
matrices["D2"] = [1  2;2  1]
matrices["D3"] = [1  3  3;3  1  2;3  2  1]
matrices["D4"] = [1  3  2  2;3  1  3  3;2  3  1  2;2  3  2  1]
matrices["D~3"] = [1  3  2  3;3  1  3  2;2  3  1  3;3  2  3  1]
matrices["D~4"] = [1  2  3  2  2;2  1  3  2  2;3  3  1  3  3;2  2  3  1  2;2  2  3  2  1]
matrices["E6"]=[1  3  2  2  2  2;3  1  3  2  2  2;2  3  1  3  2  3;2  2  3  1  3  2;2  2  2  3  1  2;2  2  3  2  2  1]
matrices["E7"]=[1  3  2  2  2  2  2;3  1  3  2  2  2  2;2  3  1  3  2  2  3;2  2  3  1  3  2  2;2  2  2  3  1  3  2;2  2  2  2  3  1  2;2  2  3  2  2  2  1]
matrices["E8"]=[1  3  2  2  2  2  2  2;3  1  3  2  2  2  2  2;2  3  1  3  2  2  2  2;2  2  3  1  3  2  2  2;2  2  2  3  1  3  2  3;2  2  2  2  3  1  3  2;2  2  2  2  2  3  1  2;2  2  2  2  3  2  2  1]
matrices["E~6"]=[1  3  2  2  2  2  2;3  1  3  2  2  2  2;2  3  1  3  2  3  2;2  2  3  1  3  2  2;2  2  2  3  1  2  2;2  2  3  2  2  1  3;2  2  2  2  2  3  1]
matrices["E~7"]=[1  3  2  2  2  2  2  2;3  1  3  2  2  2  2  2;2  3  1  3  2  2  2  2;2  2  3  1  3  2  2  3;2  2  2  3  1  3  2  2;2  2  2  2  3  1  3  2;2  2  2  2  2  3  1  2;2  2  2  3  2  2  2  1]
matrices["E~8"]=[1  3  2  2  2  2  2  2  2;3  1  3  2  2  2  2  2  2;2  3  1  3  2  2  2  2  2;2  2  3  1  3  2  2  2  2;2  2  2  3  1  3  2  2  2;2  2  2  2  3  1  3  2  3;2  2  2  2  2  3  1  3  2;2  2  2  2  2  2  3  1  2;2  2  2  2  2  3  2  2  1]
matrices["F4"] = [1  3  2  2;3  1  4  2;2  4  1  3;2  2  3  1]
matrices["F~4"] = [1  3  2  2  2;3  1  4  2  2;2  4  1  3  2;2  2  3  1  3;2  2  2  3  1]
matrices["G2"]=[1  6;6  1]
matrices["G~2"] = [1  6  2;6  1  3;2  3  1]
matrices["H2"]=[1  5;5  1]
matrices["H3"] =[1  5  2;5  1  3;2  3  1]
matrices["H4"] = [1  5  2  2;5  1  3  2;2  3  1  3;2  2  3  1]
matrices["I2(0)"] = [1  0;0  1]
matrices["I2(2)"] = [1  2;2  1]
matrices["I2(3)"] = [1  3;3  1]

@testset "Coxeter matrix generation of $group_type" for (group_type, coxeter_matrix) in matrices
    CG, _ = coxeter_group_recursive(group_type)
    @test all(coxeter_matrix .== CG.M)
end

include("Test_CoxGrpMin_ReflTable.jl")
include("Test_CoxeterSystems.jl")
