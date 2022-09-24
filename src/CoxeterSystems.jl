using Match

export classify_coxeter_matrix, classify_gcm
export is_irreducible, is_affine_type, is_finite_type
export coxeter_system, coxeter_name, degrees, order, coxeter_number, number_of_reflections, exponents

#=
    Classification of finite and affine-type Dynkin diagrams, and finite-type Coxeter systems.

The functions classify_coxeter_matrix(...) and classify_gcm(...) take a Coxeter matrix or GCM, and return a list their
connected components. Each connected component is a pair, with the first element giving the type, like (:A, 4), and the
second element of the pair giving the vertex reading order which takes the given matrix to Kac' convention. These are
meant to be internal library functions, used when constructing Coxeter systems or Cartan systems.

For example, classification of a Coxeter matrix (1's and 2's are omitted)
    .  3  .  .  .  .  .  .  .
    3  .  .  .  .  .  .  3  .
    .  .  .  3  .  .  .  .  .
    .  .  3  .  3  .  .  .  .
    .  .  .  3  .  4  .  .  .
    .  .  .  .  4  .  .  .  .
    .  3  .  .  .  .  .  .  .
    .  .  .  .  .  .  .  .  .
whose underlying Coxeter graph looks like
    1 -- 2 -- 7    8     6 === 5 -- 4 -- 3
would return the list of components
    [
        ((:A, 3), [1, 2, 7]),
        ((:C, 4), [3, 4, 5, 6]),
        ((:A, 1), [8]),
    ]
The classify_coxeter_matrix(...) function can output finite types like (:A, 3), (:C, 4), (:I, 2, m), or affine types
such as (:A, 3, :aff). It will output type C rather than B, and only uses (:I, 2, m) for m = 0 or m ≥ 7. The other
classify_gcm(...) function can output a larger range of types, in particular it can output the "dual affine" types,
such as (:B, 3, :dualaff), and the interesting type A affine types, which are (:BC, n, :aff). Connected components of
indefinite type are simply classified as (:Unknown, rank).

The following table is a list of all finite and affine type Dynkin diagrams and finite type Coxeter systems. The I(GCM)
and I(Cox) columns specify which types should be taken to make an irredundant list of finite and affine type GCMs or
Coxeter systems - when a Coxeter type should be omitted, an isomorphic type is placed in there instead. The Shape column
classifies the underlying undirected graph, which is either a path, cycle, or tree. When the underlying graph is a tree,
the vertex degrees which are ≥ 3 are listed. The Bond column is the multiset of Coxeter bonds, excluding m_st = 2 or 3.
The ~ symbol means affinisation, and the @ symbol means dual affinisation.

Type         I(GCM)    I(Cox)   Shape       Bond
----         ------    ------   -----       ----
An           n ≥ 1     n ≥ 1    Path        -
Bn           n ≥ 3     Cn       Path        4
Cn           n ≥ 2     n ≥ 2    Path        4
Dn           n ≥ 4     n ≥ 4    Tree(3)     -
E6           Yes       Yes      Tree(3)     -
E7           Yes       Yes      Tree(3)     -
E8           Yes       Yes      Tree(3)     -
F4           Yes       Yes      Path        4
G2           Yes       Yes      Path        6

Hn           n=2,3,4   n=2,3,4  Path        5
I2(m)        m ≥ 7     m ≥ 7    Path        m

A~1          Yes       Yes      Path        ∞
A~n          n ≥ 2     n ≥ 2    Cycle       -
B~n          n ≥ 3     n ≥ 3    Tree(3)     4
C~n          n ≥ 2     n ≥ 2    Path        4,4                             Bonds point in.
D~4          Yes       Yes      Tree(4)     -
D~n          n ≥ 5     n ≥ 5    Tree(3,3)   -
E~6          Yes       Yes      Tree(3)     -
E~7          Yes       Yes      Tree(3)     -
E~8          Yes       Yes      Tree(3)     -
F~4          Yes       Yes      Path        4
G~2          Yes       Yes      Path        6

BC~1         Yes       A~1      Path        ∞       Kac: A2(2)
BC~n         n ≥ 2     C~n      Path        4,4     Kac: A_{2(n - 1)}(2)    Bonds point same direction.
B@n          n ≥ 3     B~n      Tree(3)     4       Kac: A_{2n - 1}(2)
C@n          n ≥ 2     C~n      Path        4,4     Kac: D_{n-1}(2)         Bonds point out.
F@4          Yes       F~4      Path        4,4     Kac: E6(2)
G@2          Yes       G~2      Path        6       Kac: D4(3)
=#

"""
    classify_coxeter_matrix(mat)

Classify a Coxeter matrix as one of the finite or affine type Coxeter systems. A list of pairs is returned, with each
pair giving a classification such as (:A, 4) or (:D, 6, :aff), together with an ordered list of vertices making up that
component in the Kac convention.
"""
function classify_coxeter_matrix(coxmat::AbstractMatrix{T}) where {T <: Integer}
    is_coxeter_matrix(coxmat) || error("Argument was not a Coxeter matrix")

    # Rank of the large Coxeter matrix.
    coxrank = size(coxmat)[1]

    # Adjacency list and degrees for the underlying undirected graph.
    adj = Dict(s => Int[t for t in 1:coxrank if s != t && coxmat[s, t] != 2] for s in 1:coxrank)
    degrees = [length(adj[s]) for s in 1:coxrank]


    # Create a list of the connected components of the underlying undirected graph.
    components = Vector{Vector{Int64}}()
    seen = zeros(Bool, coxrank)
    for start in 1:coxrank
        seen[start] && continue

        component = [start]
        seen[start] = true
        pos = 1
        while pos <= length(component)
            node = component[pos]
            pos += 1
            for t in 1:coxrank
                if !seen[t] && coxmat[node, t] != 2 && coxmat[node, t] != 1
                    push!(component, t)
                    seen[t] = true
                end
            end
        end
        push!(components, component)
    end

    # Perform a breadth-first search from a point, returning a triple (vertices, dist, pred) of the vertices visited
    # in BFS order, their distances from the source, and a dictionary mapping each vertex to its predecessor on a
    # shortest path to the source.
    function bfs(start)
        dist = Dict(start => 0)
        pred = Dict(start => start)
        order = [start]
        pos = 1
        while pos <= length(order)
            s = order[pos]
            pos += 1
            for t in 1:coxrank
                if s != t && coxmat[s, t] != 2 && !haskey(dist, t)
                    dist[t] = dist[s] + 1
                    pred[t] = s
                    push!(order, t)
                end
            end
        end
        return (order, dist, pred)
    end

    # Given a pred dictionary output by BFS, and a starting vertex v, return the shortest path [v, ..., BFS source.]
    function getpath(pred, v)
        path = [v]
        while true
            v = path[end]
            pred[v] == v && return path
            push!(path, pred[v])
        end
    end

    # Given a vertex in a cycle, return a list of vertices in that cycle, starting from start, in cycle order.
    function getcycle(start)
        degrees[start] == 2 || error("Degree of start should be 2 for getcycle")
        order = [start]

        prev = start
        next = minimum(adj[start])
        while next != start
            push!(order, next)
            degrees[next] == 2 || error("Not a cycle")
            (prev, next) = (next, [s for s in adj[next] if s != prev][1])
        end

        return order
    end


    # Classify a component, returning something like (:A, [4, 1]) or (:Unknown, []).
    function classify_component(comp)
        # The rank of this irreducible component.
        rank = length(comp)

        # List, with multiplicity, the bonds m ∈ {0, 4, 5, ...}
        multbonds = sort(Int[m for m in [coxmat[comp[s], comp[t]] for s in 1:rank for t in s+1:rank] if m != 2 && m != 3])

        # A connected graph is a tree iff |E| = |V| - 1, and a cycle iff deg(v) = 2 for all v ∈ V.
        nedges = sum(degrees[s] for s in comp) ÷ 2
        is_tree = rank - 1 == nedges
        is_cycle = all(degrees[s] == 2 for s in comp)

        # The only non-tree graph we accept is A~n for n ≥ 2, which is a cycle.
        if !is_tree
            # If we have a cycle where every edge has bond multiplicity m = 3, then we have affine type A.
            is_cycle && multbonds == [] && return ((:A, rank-1, :aff), getcycle(minimum(comp)))

            # Otherwise, unknown type.
            return ((:Unknown, rank), comp)
        end

        # If the tree is made up of one vertex, then we have type A1.
        nedges == 0 && return ((:A, 1), comp)


        # List the leaves of the tree, and for each leaf, the bond multiplicity incident on it.
        leaves = sort([s for s in comp if degrees[s] == 1])
        leaf_bond = Dict(s => [coxmat[s, t] for t in comp if s != t && coxmat[s, t] != 2][1] for s in leaves)

        # The degrees of the component which are ≥ 3.
        bigdegrees = sort(Int[degrees[s] for s in comp if degrees[s] >= 3])

        # A path is a tree where all degrees are ≤ 2, so bigdegrees will be empty.
        if bigdegrees == []
            # Reorder the leaves so that the one with the lower bond multiplicity comes first.
            sort!(leaves; by=leaf -> leaf_bond[leaf])

            #                       mleft               mright
            # Let the path be left ------- (some path) -------- right, where mleft ≤ mright. Reorder comp so that
            # vertices are in path order [left, ..., right].
            (left, right) = leaves
            (comp, _, _) = bfs(left)
            (comp[1] == left && comp[rank] == right) || error("Comp should go from left to right now.")

            # The bond multiplicity incident on the left and right leaves.
            mleft = leaf_bond[left]
            mright = leaf_bond[right]


            ## Finite-type Coxeter systems which are paths.

            # An: No bonds of multiplicity ≂̸ 3.
            multbonds == [] && return ((:A, rank), comp)

            # Cn: Extra bond multiplicities {4}, and that bond is incident on the right vertex.
            multbonds == [4] && mright == 4 && return ((:C, rank), comp)

            # F4: Extra bond multiplicities {4}, rank 4, and both leaves have simple bonds.
            multbonds == [4] && rank == 4 && mleft == mright == 3 && return ((:F, 4), comp)

            # G2: Extra bond multiplicities {6}, rank 2.
            multbonds == [6] && rank == 2 && return ((:G, 2), comp)

            # Hn: Extra bond multiplicities {5}, that 5 is incident on the right vertex, and rank is 2, 3, or 4.
            multbonds == [5] && mright == 5 && 2 <= rank <= 4 && return ((:H, rank), comp)

            # I2(m): Extra bond multiplicities {m} for some m ≥ 7, rank 2.
            m = (length(multbonds) == 1) ? multbonds[1] : -1
            rank == 2 && m >= 7 && return ((:I, 2, m), comp)


            ## Affine-type Coxeter systems which are paths.

            # A~1: Extra bond multiplicities {∞}, rank 2.
            rank == 2 && multbonds == [0] && return ((:A, 1, :aff), comp)

            # C~n: Extra bond multiplicities {4, 4}, each incident on a leaf.
            # Coxeter graph: (n+1) == 1 -- 2 ... -- (n-2) -- (n-1) == n
            multbonds == [4, 4] && mright == mleft == 4 && return ((:C, rank-1, :aff), [comp[2:n+1]; comp[1]])

            # F~4: Extra bond multiplicities {4}, not incident on a leaf, and rank 5. In Kac, the bond in F4 points to
            # the right, and so in the affinisation we should have [3, 3, 4, 3] as the path bonds.
            if multbonds == [4] && mleft == 3 && mright == 3
                # Kac labelling should have 5 -- 1 -- 2 == 3 -- 4
                comp = coxmat[comp[3], comp[4]] == 4 ? comp : reverse(comp)
                return ((:F, 4, :aff), [comp[2:5]; comp[1]])
            end

            # G~2: Extra bond multiplicities {6}, incident on the right leaf, with rank 3.
            # Currently our path looks like 1 -- 2 ≡≡ 3, we need to return 3 -- 1 ≡≡ 2
            multbonds == [6] && rank == 3 && mright == 6 && return ((:G, 2, :aff), [comp[2:3]; comp[1]])

            ## No more finite or affine-type Coxeter systems which are paths.
            return ((:Unknown, rank), comp)
        end

        # Now for those trees with a unique vertex of degree 3, and no other vertices of higher degree.
        if bigdegrees == [3]
            # A tree with a unique vertex of degree 3 is essentially three paths joined to a special "star" point.
            # It has three leaves: let these leaves be p, q, r ordered by distance from the star point, with ties
            # broken according to their incident bond.

            # Perform a BFS to get distances from the star vertex.
            star = [s for s in comp if degrees[s] == 3][1]
            (_, dist, pred) = bfs(star)

            # Reorder the leaves by the criterion above, and take vectors of the resulting distances and bond mults.
            sort!(leaves; by=leaf -> (dist[leaf], leaf_bond[leaf]))
            leaf_dists = [dist[leaf] for leaf in leaves]
            leaf_bonds = [leaf_bond[leaf] for leaf in leaves]

            ## Finite-type Coxeter systems which are trees.

            # Dn: All simple edges, with at least two leaves having distance 1 to the star vertex.
            if multbonds == [] && leaf_dists[1:2] == [1, 1]
                # In type D4 in Kac, the star vertex is the second vertex.
                rank == 4 && return ((:D, 4), [star ; leaves])

                # Otherwise, BFS from that unique leaf at distance ≥ 2 from the star vertex to get the ordering.
                (comp, _, _) = bfs(leaves[3])
                return ((:D, rank), comp)
            end

            # En: All simple edges, ranks 6, 7, 8, with leaves at distances [1, 2, ?] from the star vertex.
            if multbonds == [] && 6 <= rank <= 8 && leaf_dists[1:2] == [1, 2]
                # For E6 and E7, we want a path 1 -- 2 -- 3star -- 4 -- 5 ( -- 6 ), with the star leaf attached to 3.
                # For E8 however, we want 1 -- 2 -- 3 -- 4 -- 5star -- 6 -- 7 to fit with Kac' labelling.
                order = (
                    (rank == 6 || rank == 7)
                    ? [getpath(pred, leaves[2]) ; getpath(pred, leaves[3])[end-1:-1:1] ; leaves[1]]
                    : [getpath(pred, leaves[3]) ; getpath(pred, leaves[2])[end-1:-1:1] ; leaves[1]]
                )
                return ((:E, rank), order)
            end


            ## Affine-type Coxeter systems which are trees.

            # B~n: Extra bond multiplicities {4}, two leaves incident on simple bonds at distance 1, and the last
            # leaf incident on the bond with multiplicity 4.
            multbonds == [4] && leaf_dists[1:2] == [1, 1] && leaf_bonds[3] == 4 && (
                return ((:B, rank - 1, :aff), [leaves[1] ; getpath(pred, leaves[3])[end:-1:1] ; leaves[2]])
            )

            # E~6: All bonds m = 3, and each leaf at distance 2.
            if multbonds == [] && leaf_dists == [2, 2, 2]
                # Extended diagram is
                #             aff
                #              |
                #              6
                #              |
                #  1 --- 2 --- 3 --- 4 --- 5
                nodes = [
                    getpath(pred, leaves[1]);
                    getpath(pred, leaves[2])[end:-1:2];
                    getpath(pred, leaves[3])[end:-1:2]
                ]
                return ((:E, 6, :aff), nodes)
            end

            # E~7: All bonds m=3, leaf distances 1, 3, and 3.
            if multbonds == [] && leaf_dists == [1, 3, 3]
                # Extended diagram is
                #                     7
                #                     |
                # aff --- 1 --- 2 --- 3 --- 4 --- 5 --- 6
                (aff, left...) = getpath(pred, leaves[2])
                right = getpath(pred, leaves[3])[end:-1:2]
                return ((:E, 7, :aff), [left ; right ; leaves[1] ; aff])
            end

            # E~8: All bonds m=3, leaf distances 1, 2, 5.
            if multbonds == [] && leaf_dists == [1, 2, 5]
                # Extended diagram is
                #                                 8
                #                                 |
                # aff --- 1 --- 2 --- 3 --- 4 --- 5 --- 6 --- 7
                (aff, left...) = getpath(pred, leaves[3])
                right = getpath(pred, leaves[2])[end:-1:2]
                return ((:E, 8, :aff), [left ; right ; leaves[1] ; aff])
            end


            ## No more finite or affine type Dynkin diagrams of this kind.
            return ((:Unknown, rank), comp)
        end

        # D~4: Branching degrees {4}, simple bonds throughout, rank 5.
        if bigdegrees == [4] && multbonds == [] && rank == 5
            # The 4-valent vertex comes second.
            star = [s for s in comp if degrees[s] == 4][1]
            return ((:D, 4, :aff), [leaves[1] ; star ; leaves[2:end]])
        end

        # D~n for n ≥ 5: Branching degrees {3, 3}, simple bonds, each leaf incident on a branching vertex.
        leaf_neighbour_degrees = [degrees[adj[leaf][1]] for leaf in leaves]
        if bigdegrees == [3, 3] && multbonds == [] && leaf_neighbour_degrees == [3, 3, 3, 3]
            # Extended diagram is
            #       aff             n
            #        |              |
            #  1 --- 2 --- ... --- n-2 --- n-1
            (order, _, _) = bfs(leaves[1])
            affleaf = [leaf for leaf in leaves if leaf in adj[leaves[1]]]
            return ((:D, rank - 1, :aff), [[s for s in order if s != affleaf] ; affleaf])
        end

        # This exhausts all finite and affine type Coxeter systems.
        return ((:Unknown, rank), comp)
    end

    return [classify_component(component) for component in components]
end


function classify_gcm(gcm::AbstractMatrix{T}) where {T <: Integer}
    is_gcm(gcm) || error("Was not given a GCM")

    # First run the classification of the underlying Coxeter graph, then check bond directions to distinguish
    # between dual cases. This relies on the convention we know which is being returned in the classification: eg if
    # a Coxeter graph is classified as type C, then the double bond is at the end.
    cox_components = classify_coxeter_matrix(gcm_to_coxeter_matrix(gcm))

    convert_component(type, comp) = @match type begin
        # Cases where we can just forward the classification of the Coxeter matrix.
        (:Unknown, n)                   => (type, comp)
        (:A, n)                         => (type, comp)
        (:A, n, :aff), if n >= 2 end    => (type, comp)
        (:D, n)                         => (type, comp)
        (:D, n, :aff)                   => (type, comp)
        (:E, n)                         => (type, comp)
        (:E, n, :aff)                   => (type, comp)

        # Non-simply laced finite types. Note the Coxeter classifier will never return type B, only type C.
        (:C, n) => gcm[comp[n], comp[n-1]] == -2 ? ((:B, n), comp) : ((:C, n), comp)
        (:F, 4) => gcm[comp[3], comp[2]] == -2 ? ((:F, 4), comp) : ((:F, 4), reverse(comp))
        (:G, 2) => gcm[comp[2], comp[1]] == -3 ? ((:G, 2), comp) : ((:G, 2), reverse(comp))

        # Non-simply laced affine types
        (:A, 1, :aff) => begin
            # If we have the (-2, -2) Cartan matrix, then it is symmetric so nothing more to do.
            gcm[comp[1], comp[2]] == -2 && return ((:A, 1, :aff), comp)

            # If we have the (-1, -4) Cartan matrix, the affine vertex should be the one with the -1 in its row.
            return (gcm[comp[2], comp[1]] == -1) ? ((:BC, 1, :aff), comp) : ((:BC, 1, :aff), reverse(comp))
        end

        (:B, n, :aff) => begin
            # We know n ≥ 3. We need to figure out whether the bond between n-1 and n is pointing to the
            # right (type B~n) or left (type B@n).
            return (gcm[comp[n], comp[n-1]] == -2) ? ((:B, n, :aff), comp) : ((:B, n, :dualaff), comp)
        end

        (:C, n, :aff) => begin
            # We know n ≥ 2, and Coxeter graph:   (n+1) == 1 -- 2 -- ... -- (n-1) == n.

            # Both arrows pointing in is type C~n.
            gcm[comp[1], comp[n+1]] == -2 && gcm[comp[n-1], comp[n]] == -2 && return ((:C, n, :aff), comp)

            # Both arrows pointing out is type C@n.
            gcm[comp[1], comp[n+1]] == -1 && gcm[comp[n-1], comp[n]] == -1 && return ((:C, n, :dualaff), comp)

            # Both arrows pointing the same way is type BC~n. The two arrows both point away from the affine vertex.
            # (n+1) =>= 1 -- 2 ... -- (n-1) =>= n.
            gcm[comp[1], comp[n+1]] == -2 && return ((:BC, n, :aff), comp)

            # At this point we have the right diagram but in the wrong direction:
            # (n+1) =<= 1 -- 2 ... -- (n-1) =<= n.
            # The correct order is (n-1), (n-2), ..., 1, n+1, n.
            return ((:BC, n, :aff), [comp[n-1:-1:1]; comp[n+1]; comp[n]])
        end

        (:F, n, :aff) => begin
            # Arrow between 2 and 3 pointing right is usual affinisation, otherwise dual affinisation.
            return (gcm[comp[3], comp[2]] == -2) ? ((:F, 4, :aff), comp) : ((:F, 4, :dualaff), comp)
        end

        (:G, n, :aff) => begin
            # Arrow between 1 and 2 pointing right is usual affinisation, otherwise dual.
            return (gcm[comp[2], comp[1]] == -3) ? ((:G, 2, :aff), comp) : ((:G, 2, :dualaff), comp)
        end
    end

    return [convert_component(comp...) for comp in cox_components]
end

# Convert a type (:A, 4) to a string like "A4".
type_to_string(type::Tuple) = @match type begin
    (letter, rank, :aff) => "$letter~$rank"
    (letter, rank, :dualaff) => "$letter@$rank"
    (:I, 2, m) => "I2(m)"
    _ => join(map(string, type))
end

# Convert a list of types to a string like "A4 x H3"
type_to_string(types::Vector) = join([type_to_string(type) for (type, comp) in types], " x ")

# Check if a particular type like (:A, 5) or (:C, 5, :aff) is finite type.
is_finite_type(type::Tuple) = @match type begin
    (:Unknown, n) => false
    (:A, n) => true
    (:B, n) => true
    (:C, n) => true
    (:D, n) => true
    (:E, n) => 6 <= n <= 8
    (:F, 4) => true
    (:G, 2) => true
    (:H, n) => 2 <= n <= 4
    (:I, 2, m) => m >= 2
    _ => false
end

# Check if a list of components like [((:A, 5), [1, 2, 3, 4, 5]), ((:B, 3, :aff), [6, 7, 8])] is finite type.
is_finite_type(composite_type::Vector) = all(is_finite_type(type) for (type, comp) in composite_type)

# Check if a particular type (:A, 5, :aff) is affine type.
is_affine_type(type::Tuple) = @match type begin
    (_, _, :aff) => true
    (_, _, :dualaff) => true
    _ => false
end

# Check if a list of types is affine type. This only occurs when there is a single type, and it is affine.
is_affine_type(composite_type::Vector) = length(composite_type) == 1 && is_affine_type(composite_type[1][1])

# The degrees of a finite type irreducible Coxeter system.
# From Table 1 in §3.7 of Humphrey's "Reflection groups and Coxeter groups".
degrees(type::Tuple) = @match type begin
    (:A, n), if n >= 1 end => Vector(2:n+1)
    (:B, n), if n >= 2 end => Vector(2:2:2*n)
    (:C, n), if n >= 2 end => Vector(2:2:2*n)
    (:D, n), if n >= 2 end => sort([Vector(2:2:2*n-2); n])
    (:E, 6) => [2, 5, 6, 8, 9, 12]
    (:E, 7) => [2, 6, 8, 10, 12, 14, 18]
    (:E, 8) => [2, 8, 12, 14, 18, 20, 24, 30]
    (:F, 4) => [2, 6, 8, 12]
    (:G, 2) => [2, 6]
    (:H, 3) => [2, 6, 10]
    (:H, 4) => [2, 12, 20, 30]
    (:I, 2, m), if m >= 2 end => [2, m]
    _ => error("$(type) is not finite type")
end


# A Coxeter system is a classified Coxeter matrix.
struct CoxeterSystem
    coxeter_matrix::Matrix{Int}
    components::Vector{Tuple{Tuple, Vector{Int}}}
end

"""
    coxeter_system(mat)

Create a Coxeter system from a Coxeter matrix or Cartan matrix.
"""
function coxeter_system(mat::AbstractMatrix{T}) where {T <: Integer}
    if is_gcm(mat)
        mat = gcm_to_coxeter_matrix(mat)
    end
    is_coxeter_matrix(mat) || error("Given matrix was not a Coxeter or Cartan matrix.")
    components = classify_coxeter_matrix(mat)
    return CoxeterSystem(mat, components)
end

# Pretty-print some data about a Coxeter system.
function Base.show(io::IO, cox::CoxeterSystem)
    adjectives = [
        "rank $(rank(cox))",
        length(cox.components) == 1 ? "irreducible" : "reducible",
        is_finite_type(cox.components) ? "finite type" : is_affine_type(cox.components) ? "affine type" : "indefinite type",
    ]
    write(io, "Coxeter system ($(join(adjectives, ", "))) of type $(type_to_string(cox.components))")
end

@doc raw"""
    coxeter_matrix(cox::CoxeterSystem)

The Coxeter matrix of a Coxeter system ``(W, S)`` is the matrix ``m_{s, t} = \operatorname{ord}(st)`` for ``s, t ∈ S``.
"""
coxeter_matrix(cox::CoxeterSystem) = copy(cox.coxeter_matrix)

"""
    rank(cox::CoxeterSystem)

The rank of a Coxeter system ``(W, S)`` is the number ``|S|`` of simple generators.
"""
rank(cox::CoxeterSystem) = size(cox.coxeter_matrix)[1]

"""
    coxeter_name(cox::CoxeterSystem)

A string representing the type of the system, such as "H3 x A2".
"""
coxeter_name(cox::CoxeterSystem) = type_to_string(cox.components)

"""
    is_irreducible(cox::CoxeterSystem)

A Coxeter system ``(W, S)`` is irreducible if its associated underyling graph consists of a single connected component.
In particular, ``|S| ≥ 1``.
"""
is_irreducible(cox::CoxeterSystem) = length(cox.components) == 1

"""
    is_finite_type(cox::CoxeterSystem)

A Coxeter system ``(W, S)`` is finite type if the group ``W`` is finite.
"""
is_finite_type(cox::CoxeterSystem) = is_finite_type(cox.components)

"""
    is_affine_type(cox::CoxeterSystem)

A Coxeter system ``(W, S)`` is affine type if is irreducible, and the unique component is of affine type.
"""
is_affine_type(cox::CoxeterSystem) = is_affine_type(cox.components)

@doc raw"""
    degrees(cox::CoxeterSystem)

The degrees ``d_1 ≤ ⋯ ≤ d_r`` of a Coxeter system ``(W, S)`` are the degrees of the fundamental invariants for the
action of ``W`` on ``\operatorname{Sym}(V^*)``, where ``V`` is an irreducible representation of ``W`` with ``S`` acting
by reflections.
Throws an error if the system is not finite type.
"""
function degrees(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    return sort([degree for (type, comp) in cox.components for degree in degrees(type)])
end

@doc raw"""
    exponents(cox::CoxeterSystem)

The exponents ``m_1 ≤ ⋯ ≤ m_r`` of a Coxeter system ``(W, S)`` describe the eigenvalues (with multiplicity) of a Coxeter
element acting on the Tits representation ``V``, which are roots of unity ``\exp(2 \pi i m_i / h)`` where ``h`` is the
Coxeter number.
Throws an error if the system is not finite type.
"""
function exponents(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    return (degrees(cox) .- 1)
end

"""
    order(cox::CoxeterSystem)

The order of a Coxeter system ``(W, S)`` is the order of the group ``W``.
Throws an error if the system is not finite type.
"""
function order(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    return reduce(*, degrees(cox); init=BigInt(1))
end

"""
    number_of_reflections(cox::CoxeterSystem)

The number of reflections in a Coxeter system ``(W, S)`` is the number of conjugates of the simple generators ``S``, or
equivalently the number of roots in a root system for ``W``.
Throws an error if the system is not finite type.
"""
function number_of_reflections(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    return sum(degrees(cox)) - rank(cox)
end

"""
    coxeter_number(cox::CoxeterSystem)

The Coxeter number ``h`` of an irreducible Coxeter system ``(W, S)`` is the order of any Coxeter element (product of all
simple reflections), or ``2 |Φ^+| / |S|``, or the largest degree ``d_{|S|}``.
If the Coxeter system is not irreducible then these definitions all diverge, so this function will throw an error on a
reducible system.

Throws an error if the system is reducible or not of finite type.
"""
function coxeter_number(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    is_irreducible(cox) || error("Coxeter system is reducible")

    return degrees(cox)[end]
end
