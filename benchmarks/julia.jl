using CoxeterGroups

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

function multiply_all_pairs(Gelts)
    sum = 0
    for x = Gelts, y = Gelts
        sum += length(x * y)
    end
    return sum
end

# CoxeterGroup
# For A5, took about 3 us per element. Magma takes about 0.2 us per element.
# For A6, took about 6.8 us per element. Magma takes about 0.2 us per element.
# G, S = CoxeterGroup(coxeter_type("A", 4).coxeter_mat)

# CoxGrpMin
# For A5, took about 0.25 us per element. Magma takes about 0.20 us per element.
# For A6, took about 0.33 us per element. Magma takes about 0.26 us per element.
#G, S = coxeter_group_min(coxeter_type("A", 4).gcm)

# CoxGrpSym
# For A6 = S7, took about 0.075 us per element.
G, S = symmetric_group(7)

Gelts = enumerate_whole_group(G)
time = @elapsed lengthsum = multiply_all_pairs(Gelts)
us_per_elt = time * 1_000_000 / length(Gelts)^2
print("Took about $us_per_elt us per element\n")
print(lengthsum)
