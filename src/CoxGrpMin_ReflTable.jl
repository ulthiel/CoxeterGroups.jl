#=
This purpose of this file is to implement the function create_reflection_table_coxeter(), which takes an arbitrary
Coxeter matrix and outputs the minimal root reflection table. It follows the general approach outlined by Casselman
in "Multiplication in Coxeter Groups II", and one certainly needs to read that paper to understand the code below. I
would also recommend reading the more straightforward create_reflection_table_gcm(), which can be used when the root
system is specified as a GCM (and so in particular is crystallographic).

A key supporting data structure is a QInt, which knows just enough about quantum integers at roots of unity to be useful
to the algorithm. The QInt structure should not be used outside of this file.

Let q ≥ 1 be an integer, and let [n]_q denote the quantum integer (z^n - z^-n) / (z - z^-1) where z = exp(2 π / q), so:
  [0]_q = 0
  [1]_q = 1
  [2]_q = z + z^-1        = 2 cos (2 π / q)
  [3]_q = z^2 + 1 + z^-2
and so on. A Coxeter matrix can be converted to a symmetric Cartan matrix containing roots of unity by turning a bond
multiplicity m into the quantum integer -[2]_2m (with special treatment for the diagonal and infinite entries). Note
in the example below that [2]_4 = 0 and [2]_6 = 1.

  [ 1  3  5  2 ]      [   2      -[2]_6   -[2]_10   0     ]
  [ 3  1  7  3 ]  =>  [ -[2]_6     2      -[2]_14  -[2]_6 ]
  [ 5  7  1  4 ]  =>  [ -[2]_10  -[2]_14    2      -[2]_8 ]
  [ 2  3  4  1 ]      [   0      -[2]_6   -[2]_8     2    ]

In the usual reflection representation (the Tits representation) of a Coxeter group, the coefficient of any minimal root
has one of the following forms:
- An integer, i.e. a QInt with q = 1,
- An integer multiple of [2]_2m for some m ≥ 3 in the Cartan matrix,
- Exactly [n]_2m for n > 2 (only occurs for dihedral roots [n] α_s + [n+1] α_t, where m = m_st),
- a + b[2] for integers a, b: only occurs when m_st = 5.

A QInt can faithfully represent a linear combination of quantum integers over a certain q provided that q=1 (where a
QInt is just an integer), or provided that 2 ≤ q = 2m ≤ 12, i.e. for m_st in {2, 3, 4, 5, 6}. For these kind of QInts,
equality of QInts works the same as equality for the underlying numbers. For more general q, the QInt data structure
does not know how to simplify numbers (it doesn't have any cyclotomic polynomial knowledge), but we are relying on the
fact that this should never be a problem, by the classification of the special forms of a minimal root listed above.

As the algorithm proceeds, it will happily take linear combinations of QInts at q=1 (i.e. regular integers) with QInts
at q = 2m ≥ 2, but it will throw an error if mixed cyclotomy (adding QInts with different q ≂̸ 1) would ever occur. There
are enough special cases in the algorithm such that this should never happen.

The other operation we need on QInts (and for a slightly larger class of QInts than those listed above) is to check
whether they lie in the range (-2, 2), since whenever the product ⟨α, β⟩ lies outside this range, then α and β generate
an infinite group (and hence one cannot be a minimal root). Again, we try to implement this operation for a class of
QInts just large enough for what we need.
=#

# Make it a subtype of Number, so that we can use a Vector of QInt in the actual algorithm.
struct QInt <: Number
    q::Int
    coeffs::Vector{Tuple{Int, Int}}

    """Create a QInt from a list of arbitrary (n, coeff) pairs."""
    function QInt(q::Int, pairs::Vector{Tuple{Int, Int}})
        q >= 1 || error("q = $q should be positive")
        coeffs = spzeros(Int, q)
        for (n, c) in pairs
            (n, c) = qint_reduce(q, n, c)

            # Since [0]_q = 0, ignore the coefficient in front of zero.
            if n != 0
                coeffs[n] += c
            end
        end

        # Everything zero? Return an integer-qint.
        iszero(coeffs) && return new(1, [(1, 0)])

        # Only the 1-coefficient remains? Return an integer-qint.
        support = findall(!iszero, coeffs)
        support == [1] && return new(1, [(1, coeffs[1])])

        # General case.
        return new(q, [(n, coeffs[n]) for n in support])
    end

    """Create the QInt [n]_q."""
    QInt(q::Int, n::Int) = QInt(q, [(n, 1)])

    """Create a QInt from an integer, [n]_1."""
    QInt(n::Int) = QInt(1, n)
end

Base.:(==)(x::QInt, y::QInt) = x.q == y.q && x.coeffs == y.coeffs
Base.hash(x::QInt, h::UInt) = hash(x.coeffs, hash(x.q, h))

Base.:(-)(x::QInt) = QInt(x.q, [(n, -c) for (n, c) in x.coeffs])

function Base.:(+)(x::QInt, y::QInt)
    x.q == 1 || y.q == 1 || x.q == y.q || error("Cannot add quantum integers with different q.")
    return QInt(max(x.q, y.q), [x.coeffs; y.coeffs])
end

function Base.:(-)(x::QInt, y::QInt)
    x.q == 1 || y.q == 1 || x.q == y.q || error("Cannot add quantum integers with different q.")
    return QInt(max(x.q, y.q), [x.coeffs; [(n, -c) for (n, c) in y.coeffs]])
end

function Base.:(*)(x::QInt, y::QInt)
    x.q == 1 || y.q == 1 || x.q == y.q || error("Cannot multiply quantum integers with different q.")
    return QInt(max(x.q, y.q), [(nm, c*d) for (n, c) in x.coeffs for (m, d) in y.coeffs for nm in qint_prod(n, m)])
end

"""
Test if a QInt is in the range (-2, 2), or throw an error if this cannot be exactly determined.
It will always work (return true or false) for numbers of the form:
- Integers, i.e. integer multiples of [1]_1.
- Integer multiples of a quantum integer [n]_2m.
- ±(1 - [2]_10) = ±2 cos(2 pi / 5) for m = 5
- ±(1 - [3]_2m) = ±2 cos(2 pi / m) for m >= 7
- Any linear combination of two or more quantum integers, all having the same sign, eg -[2] - [4]
"""
function in_two_interval(x::QInt)
    # For integer QInts.
    x.q == 1 && return -2 < x.coeffs[1][2] < 2

    # At this point all we should be checking are q=2m for m ≥ 4.
    mod(x.q, 2) == 0 && x.q >= 8 || error("Cannot check if QInt $x in (-2, 2) with q = $q")

    # Extract coefficient list, and perhaps flip signs so that the first is positive.
    coeffs = x.coeffs[1][2] > 0 ? x.coeffs : [(n, -c) for (n, c) in x.coeffs]

    # For integer multiples of a quantum integer, we have [2] < 2 always, and [3] < 2 for q ≤ 10. Otherwise ≥ 2.
    if length(coeffs) == 1
        (n, c) = coeffs[1]
        return (n, c) == (2, 1) || (n, c) == (3, 1) && x.q <= 10
    end

    # If everything has the same sign, then we have at least 2[2] since len(coeffs) ≥ 2. So the result is ≥ 2.
    all(c > 0 for (n, c) in coeffs) && return false

    # [1] - [3] = -2 cos(2 π / m) is in the range (-2, 2)
    coeffs == [(1, 1), (3, -1)] && return true

    # [1]_10 - [2]_10 = -(golden ratio) is in the range (-2, 2)
    x.q == 10 && coeffs == [(1, 1), (2, -1)] && return true

    error("Cannot check if QInt $x is in the range (-2, 2)")
end



"""Transform a pair (n, c) representing c [n]_q into an equivalent pair, with n reduced as far as possible."""
function qint_reduce(q::Int, n::Int, c::Int)
    q >= 1 || error("q = $q should be ≥ 1")

    # [-n] = -[n] for all q.
    if n < 0
        (n, c) = (-n, -c)
    end

    # [n]_1 = n for all n.
    q == 1 && return 1, n*c

    # [n]_2 = (-1)^n * n, so [1]_2 = 1, [2]_2 = -2, [3]_2 = 3, etc.
    q == 2 && return n * (n%2 == 1 && c || -c)

    # Reduce n into the range [0, q) using periodicity [n + q]_q = [n]_q.
    if !(0 <= n < q)
        n = mod(n, q)
    end

    # Reduce into the range [0, q/2) using oddness: [q-n]_q = -[n]_q.
    if !(0 <= n < q/2)
        (n, c) = (q - n, -c)
    end

    # If q is even, reduce into the range [0, q/4] using symmetry about q/2.
    if q % 2 == 0 && !(0 <= n <= q / 4)
        n = q ÷ 2 - n
    end

    # Finally, if q happens to equal 12, we know that [3]_12 = 2.
    if q == 12 && n == 3
        (n, c) = (1, 2 * c)
    end

    return (n, c)
end

"""Return the coefficients appearing in the product [n][m] for n, m ≥ 0. For example, [n][2] = [n+1] + [n-1]."""
function qint_prod(n::Int, m::Int)
    n >= 0 && m >= 0 || error("n, m must be ≥ 0.")
    if n > m
        (n, m) = (m, n)
    end

    return [j+m for j in 1-n:2:n]
end

"""
Create the minimal root reflection table for an arbitrary Coxeter matrix.

This algorithm is indended to be read accompanied by *Computation in Coxeter groups. II: Constructing minimal roots* by
Casselman [Cas08].
"""
function create_reflection_table_coxeter(coxeter_mat::AbstractArray{T}) where {T <: Integer}
    is_coxeter_matrix(coxeter_mat) || error("Not a Coxeter matrix")

    # Convert the Coxeter matrix to a Cartan matrix with QInt entries, by replacing m with -[2]_2m for finite m, or
    # by the integer -2 for m = 0 (representing an infinite bond), or the integer 2 for m = 1 (diagonal entries).
    rank = size(coxeter_mat)[1]
    cartan_mat = map(m -> m == 0 ? QInt(-2) : m == 1 ? QInt(2) : QInt(2*m, -2), coxeter_mat)

    # Roots: A list of roots we encounter, in breadth-first search order. Initialised to the simple roots.
    # Depth: Map a root to its depth, with simple roots having depth 1. Its keys track roots as a set.
    # Refl: If refl[s, root] is zeroroot, then s(root) is negative or no longer minimal.
    #       If refl[s, root] is undefined, then s(root) has not yet been determined.
    #       Otherwise, refl[s, root] = s(root) where s(root) is minimal.
    # Locks: Take a root to the set of simple reflections which make s(root) positive and no longer minimal.
    zeroroot = [QInt(0) for _ in 1:rank]
    roots = [[QInt(i == j ? 1 : 0) for i in 1:rank] for j in 1:rank]
    depth = Dict(root => 1 for root in roots)
    refl = Dict((i, roots[i]) => zeroroot for i in 1:rank)
    locks = Dict(root => Set() for root in roots)

    # Apply the reflection s(root), using the usual formula s(λ) = λ - ⟨α_s^, λ⟩ α_s. Note that the cartan matrix
    # here is symmetric, so there is no difference between roots and coroots.
    reflect(s::Int, root::Vector{QInt}) = root - sum(cartan_mat[s, :] .* root) * roots[s]

    # Apply the reflection s(root), and install this information in the data structures.
    function install_reflection(s::Int, root::Vector{QInt})
        reflected = reflect(s, root)
        if !haskey(depth, reflected)
            push!(roots, reflected)
            depth[reflected] = depth[root] + 1
            locks[reflected] = Set()
        end
        refl[(s, root)] = reflected
        refl[(s, reflected)] = root
        union!(locks[reflected], locks[root])
    end

    # Record the fact that s(root) is no longer minimal.
    function install_lock(s::Int, root::Vector{QInt})
        refl[(s, root)] = zeroroot
        push!(locks[root], s)
    end

    # Decide what to do with the pair (s, root): one of "skip", "lock", or "reflect".
    function decide(s::Int, root::Vector{QInt})
        # If the reflection is already calculated (since the root is s(some_lower_root)), skip.
        haskey(refl, (s, root)) && return "skip"

        # supp: the set of simples over which the root is supported.
        # supp_link: the set of edges incident on s and contained within the support, indexed by the vertex != s.
        supp = [t for t in 1:rank if root[t] != QInt(0)]
        supp_link = [t for t in 1:rank if t != s && root[t] != QInt(0) && coxeter_mat[s, t] != 2]

        # If a root is dihedral (has a support of size 2) and s is one of the elements of this support, then go ahead
        # and install the reflection. Due to the way that QInts work, this will do the right thing when reflections loop
        # back around.
        length(supp) == 2 && s in supp && return "reflect"

        # When s is outside the support of the root...
        if s ∉ supp
            # If s is not in the support, nor linked to anything in the support, then s(root) = root. Install the reflection.
            length(supp_link) == 0 && return "reflect"

            # Corollary 8.3 of [Cas08]: The support of a minimal root is a tree containing no infinite bonds. Therefore
            # if we are about to add an infinite bond or a cycle, the result cannot be minimal.
            any(coxeter_mat[s, t] == 0 for t in supp_link) && return "lock"
            length(supp_link) >= 2 && return "lock"

            # Now the link consists of a single vertex t. If there is a simple link between s and t, then the pairing
            # <root, alpha_s^> = -root[t], so we can reflect iff 0 < root[t] < 2. If there is a multiple bond between
            # s and t, then we may reflect if and only if root[t] = 1, since otherwise both [2], root[t] >= sqrt(2),
            # and hence <root, alpha_s^> would be <= -2.
            t = supp_link[1]
            if coxeter_mat[s, t] == 3
                return in_two_interval(root[t]) ? "reflect" : "lock"
            else
                return root[t] == QInt(1) ? "reflect" : "lock"
            end
        end

        # At this point we know that s is within the support of the root, and the support size is not 2. It shouldn't be
        # 1 either, since that should only happen for simple roots, and we've already initialised the refl table for them.
        length(supp) >= 3 || error("Length of support is $(length(supp)), it should be ≥ 3.")

        # If root[s] = 1, then there is a lot to say about performing the reflection here. Since s in in the support of
        # the root, we also know that in its supp_link we will only encounter finite edges.
        if root[s] == QInt(1)
            # If there are four or more edges in the link, then the scalar product here will be ≤ -2.
            length(supp_link) >= 4 && return "lock"

            # If there are three edges in the link, then the only way to have a scalar product > -2 is for all of the
            # edges to be simple, and for the root to have coefficients 1 in those places.
            if length(supp_link) == 3
                return all(coxeter_mat[s, t] == 3 && root[t] == QInt(1) for t in supp_link) ? "reflect" : "lock"
            end

            # If there are two edges in the link, then we have three cases depending on whether these edges are both
            # multiple bonds, or one multiple and one simple, or both simple. This is 11.3 and 11.4 of [Cas08].
            if length(supp_link) == 2
                # Let t and u be the vertices in the link, with m_st ≤ m_su.
                (t, u) = sort(supp_link, by=x -> coxeter_mat[s, x])

                # If they are both multiple bonds, then since root[s] = 1 we should already have root[t], root[u] ≥ √2.
                if coxeter_mat[s, t] != 3
                    (root[t] != QInt(1) && root[u] != QInt(1)) || error("Unexpected root coeffs")
                    return "lock"
                end

                # If st is simple and su is multiple, by checking on a case-by-case basis we can only reflect here if
                # m_su is 4 or 5, root[t] = 1, and root[u] = [2].
                if coxeter_mat[s, t] == 3 && coxeter_mat[s, u] != 3
                    if root[t] == QInt(1) && 4 <= coxeter_mat[s, u] <= 5 && root[u] == QInt(2 * coxeter_mat[s, u], 2)
                        return "reflect"
                    else
                        return "lock"
                    end
                end

                # If both bonds are simple, then reflections can only take place if both ends are 1, or one is 1 and the
                # other 2, or one is 1 and the other is 1 + [2]_10 (specifically in m_su = 5).
                allowed_sets = [Set([QInt(1)]), Set([QInt(1), QInt(2)]), Set([QInt(1), QInt(1) + QInt(10, 2)])]
                actual = Set([root[t], root[u]])
                return any(actual ⊆ allowed for allowed in allowed_sets) ? "reflect" : "lock"
            end

            # If there is a single element in the support linked to s (which is not in the support), then the only two
            # cases to consider are a single bond or multiple bond.
            if length(supp_link) == 1
                t = supp_link[1]

                # We'll let the single bond case fall through to the end of the function.

                # In a multiple bond, the other coefficient shouldn't be one.
                if coxeter_mat[s, t] > 3
                    root[t] != QInt(1) || error("Encountered unexpected root with coefficient 1")

                    # If root[t] is [2], then 2 - [2][2] = 1 - [3] > -2 always.
                    # If root[t] is [3], then 2 - [2][3] = 2 - [2] - [4] > -2 only for m ≤ 6.
                    mst = coxeter_mat[s, t]
                    if root[t] == QInt(2*mst, 2) || root[t] == QInt(2*mst, 3) && mst <= 6
                        return "reflect"
                    end
                end
            end
        end

        pairing = sum(cartan_mat[s, :] .* root)
        return in_two_interval(pairing) ? "reflect" : "lock"
    end


    # Main loop. Pull roots off the queue one by one, in breadth-first reflection order. As each is pulled off its
    # inherited locks are calculated, then all simple reflections which are not inherited locks are applied.
    pos = 1
    while pos <= length(roots)
        root = roots[pos]
        pos += 1

        # Install the inherited locks
        for s in locks[root]
            install_lock(s, root)
        end

        # Calculate the new roots.
        for s in 1:rank
            if s in locks[root]
                continue
            end
            action = decide(s, root)
            if action == "skip"
            elseif action == "lock"
                install_lock(s, root)
            elseif action == "reflect"
                install_reflection(s, root)
            else
                error("Uknown action $action")
            end
        end
    end

    # Now build the minimal root reflection table.
    rootidx = Dict(root => i for (i, root) in enumerate(roots))
    rootidx[zeroroot] = 0
    refltab = zeros(Int64, rank, length(roots))
    for i = 1:length(roots), s = 1:rank
        refltab[s, i] = rootidx[get(refl, (s, roots[i]), zeroroot)]
    end

    return refltab
end
