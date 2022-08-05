using Test
using CoxeterGroups

# Construct every element of the group by right-multiplying by generators.
# Used for some basic tests to see if the group is being constructed correctly.
function enumerate_whole_group(G::CoxeterGroup)
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


# We should be able to construct the trivial group with a 0x0 Coxeter matrix.
S1, () = CoxeterGroup(Array{Int64}(undef, 0, 0), String[])
S1Elts = enumerate_whole_group(S1)
@test length(S1Elts) == 1


S2, (s,) = CoxeterGroup([1;;], ["s"])
S2Elts = enumerate_whole_group(S2)
@test length(S2Elts) == 2


S3, (s, t) = CoxeterGroup([1 3; 3 1], ["s", "t"])
S3Elts = enumerate_whole_group(S3)
@test length(S3Elts) == 6
