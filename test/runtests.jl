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
    CoxeterGroup,
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
