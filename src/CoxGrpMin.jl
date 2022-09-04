# CoxMin contains the CoxGrpMin and CoxEltMin implementations which use the minimal root automaton.
import LinearAlgebra: I

export CoxGrpMin, CoxEltMin, coxeter_group_min

"A Coxeter group, implemented using the minimal root reflection table."
mutable struct CoxGrpMin <: CoxGrp
    "The Coxeter matrix. Zero entries stand for infinity."
    coxeter_matrix::Matrix{Int64}

    # refl_table[s, α] is the reflection s⋅α of the root. The simple roots are numbered from 1 up to the rank,
    # the other minimal roots are numbered arbitrarily. refl_table[s, s] = 0 for each simple root s, and also
    # refl_table[s, α] = 0 if s⋅α leaves the set of minimal roots.
    # Since Julia stores matrices in column-major order, this (rather than the transpose) seems like the most
    # sensible layout.
    refl_table::Matrix{Int64}

    # is_finite stores whether the Coxeter group is finite or not. This can also be determined by looking for
    # holes in the refl_table, but we will cache the result of this search here.
    is_finite::Bool
end


"""
A CoxEltMin is a representation of a Coxeter group element in ShortLex normal form: as a word in the generators,
this is the lexicographically least reduced expression. Fixing a particular normal form means that equality testing
and hashing of Coxeter group elements is straightforward, as each has a unique representation.
"""
struct CoxEltMin <: CoxElt
    # Parent group
    group::CoxGrpMin

    # ShortLex normal form of the word. (Supports Coxeter groups with up to 255 generators, numbered [1, 255]).
    word::Vector{UInt8}
end


"""
    coxeter_group_min(mat) -> CoxGrpMin, gens

Given a Coxeter matrix or GCM, create a Coxeter group implemented via the minimal roots automaton.
Currently only a GCM can be passed.
"""
function coxeter_group_min(mat::Matrix{T}) where T <: Integer
    if !is_gcm(mat)
        throw(DomainError(mat, "Is not a generalised Cartan matrix"))
    end

    # Rank of the system.
    n, _ = size(mat)

    # Generate the minimal root reflection table using the GCM algorithm.
    refl_table = create_reflection_table_gcm(mat)

    # There are n zeros in the reflection table always, coming from the simple reflections reflecting
    # their simple roots to negative. If there are more than n zeros, then there are some non-minimal
    # positive roots, and the resulting group is infinite.
    is_finite = count(x -> x == 0, refl_table) > n
    G = CoxGrpMin(gcm_to_coxeter_matrix(mat), refl_table, is_finite)
    return G, generators(G)
end


"""
    generators(W::CoxGrp)

Return the simple generators of ``W``.
"""
generators(grp::CoxGrpMin) = [CoxEltMin(grp, [s]) for s=1:rank(grp)]

"""
    rank(W::CoxGrp) -> Integer

Return the rank (number of generators) of the Coxeter group.
"""
rank(grp::CoxGrpMin) = size(grp.coxeter_matrix)[1]

"""
    one(W::CoxGrp) -> CoxElt

Return the identity in the Coxeter group ``W``.
"""
Base.one(grp::CoxGrpMin) = CoxEltMin(grp, [])

"""
    coxeter_matrix(W::CoxGrp) -> Matrix

Return the Coxeter matrix for this group.
"""
coxeter_matrix(grp::CoxGrpMin) = copy(grp.coxeter_matrix)

"""
    is_finite(W::CoxGrp) -> Bool

Return whether the group ``W`` is finite or not."""
is_finite(grp::CoxGrpMin) = grp.is_finite

"""
    parent(w::CoxElt) -> CoxGrp

Return the parent group of ``W``.
"""
Base.parent(w::CoxEltMin) = w.group

"""
    isone(w::CoxElt) -> Bool

Check if ``w`` is the identity of the group.
"""
Base.isone(w::CoxEltMin) = length(w.word) == 0

"""
    length(w::CoxElt) -> Integer

The length of the group element ``w``, i.e. the number of generators in any reduced expression for ``w``.
"""
Base.length(w::CoxEltMin) = length(w.word)

"Check equality between Coxeter group elements"
Base.:(==)(w::CoxEltMin, x::CoxEltMin) = w.word == x.word

Base.hash(w::CoxEltMin, h::UInt) = hash(w.word, h)

"""
    short_lex(w::CoxElt) -> Vector{Int8}

Return the `ShortLex` normal form of ``w``, meaning the lexicographically least reduced expression.
"""
short_lex(w::CoxEltMin) = copy(w.word)

"Print a Coxeter group element in ShortLex normal form"
#Base.show(io::IO, x::CoxEltMin) = join(io, x.word)

"""
    is_right_descent(w::CoxElt, t::Integer) -> Bool

Return true if ``t`` is a right descent of ``w``, i.e. if ``l(wt) < l(w)``.
"""
function is_right_descent(w::CoxEltMin, t::T) where {T <: Integer}
    if !(1 <= t <= rank(w.group))
        throw(DomainError(t, "$t is not a generator in the range [1, $(rank(w.group))]"))
    end

    # Let t be a reflection (not necessarily simple) with corresponding root α, and w any group element.
    # Then l(wt) < l(w) if and only if w⋅α < 0.
    # Hence given a reduced expression s(1) ... s(n), not necessarily in ShortLex normal form, we may test
    # if α becomes negative by running s(n)⋅α, s(n-1)⋅(s(n)⋅α), etc through the minimal root reflection table.
    root = Int64(t)

    # It seems that neither Base.reverse or Iterators.reverse can iterate over a vector in reverse without
    # causing an allocation?
    for i = length(w.word):-1:1
        s = w.word[i]

        # If we are about to reflect to a negative root, since the rest of the word is reduced we will never
        # make it back to a positive root. Hence root becomes negative in the result.
        if s == root
            return true

        # If we are about to reflect to a non-minimal root, then since the word is reduced we'll also never
        # make it back to the minimal roots.
        elseif w.group.refl_table[s, root] == 0
            return false
        end

        root = w.group.refl_table[s, root]
    end

    return false
end

"""
    is_left_descent(t::Integer, w::CoxElt) -> Bool

Return true if ``t`` is a left descent of ``w``, i.e. if ``l(tw) < l(w)``.
"""
function is_left_descent(t::T, w::CoxEltMin) where T <: Integer
    if !(1 <= t <= rank(w.group))
        throw(DomainError(t, "$t is not a generator in the range [1, $(rank(w.group))]"))
    end

    # See the comments in is_right_descent(), here we are just reversing the order in which we read the word.
    root = UInt8(t)
    for s = w.word
        if s == root
            return true
        elseif w.group.refl_table[s, root] == 0
            return false
        end

        root = w.group.refl_table[s, root]
    end

    return false
end

"""
Multiply a word on the right by a simple generator t, mutating the word.
The word must be reduced, and if it is a ShortLex normal form, then the resulting
word will also be a ShortLex normal form.
"""
function _mult_r!(refl::Matrix{Int64}, word::Array{UInt8}, t::UInt8)
    # We proceed along the word from right-to-left, looking for either an insertion point or a deletion point.
    # A priori, we don't know which will take place (insertion or deletion).
    # If a deletion eventually takes place, its location is unique and there is only one thing to do.
    # If an insertion is to take place there might be multiple valid places to insert a letter (the obvious place
    # is at the end), but we want to choose the place which makes the result as lexicographically least as possible.

    # Let word = s(1) ... s(n). We start with the root t, and incrementally calculate the reflections s(n)⋅t,
    # s(n-1) s(n) ⋅ t, and so on. For 1 <= i <= n, let α = s(i+1)...s(n)⋅t be the root so far, which we assume
    # inductively to be a minimal root (so in particular, still positive). We start out with the potential insertion
    # position at the back, and the potential insertion letter equal to t.
    # By considering s(i) and α there are several cases:
    #   1. s(i)⋅α < 0 iff α is the simple root corresponding to s(i). If this happens, then since the word is reduced
    #      the root will stay negative, and we are looking at a deletion. The reflection corresponding to the root
    #      α is s(i+1) ... s(n) t s(n) ... s(i+1) = s(i), and hence we have s(i+1) ... s(n) t = s(i) ... s(n). After
    #      replacing, we have an s(i)^2 which can be deleted. So in this case we are looking at a deletion at s(i).
    #
    #   2. s(i)⋅α is the simple root corresponding to some other generator r. In this case we get another insertion
    #      point: considering the reflection associated to α we have s(i) ... s(n) = r s(i) ... s(n), so the new
    #      insertion point is to insert r just before position i. If r < s(i) lexicographically then this gives a
    #      better insertion point than previously, and so we update our best insertion point. Otherwise, do nothing.
    #
    #   3. s(i)⋅α is still positive, but no longer minimal. In this case since the word is reduced we will forever
    #      stay in the set of non-minimal positive roots, so we have an insertion. Return the best insertion point
    #      so far.
    # If nothing special happens, then simply update the root and continue. If we reach the start of the word then
    # we have an insertion, so return the best insertion so far.

    insertion_point = length(word) + 1
    insertion_letter = t
    root = Int64(t)
    for i = length(word):-1:1
        # Deletion case
        if word[i] == root
            deleteat!(word, i)
            return
        end

        root = refl[word[i], root]
        if root == 0
            # No longer minimal: return the best insertion point.
            break
        end

        # Check if we have a better insertion point now. Since word[i] is a simple
        # root, if root < word[i] it must be simple.
        if root < word[i]
            insertion_point = i
            insertion_letter = UInt8(root)
        end
    end

    insert!(word, insertion_point, insertion_letter)
    return
end

"""
    *(w, x)

The product ``wx`` of Coxeter group elements."""
function Base.:(*)(w::CoxEltMin, x::CoxEltMin)
    if w.group != x.group
        throw(ArgumentError("The elements come from different Coxeter groups"))
    end

    refl = w.group.refl_table
    word = copy(w.word)
    for s = x.word
        _mult_r!(refl, word, s)
    end

    return CoxEltMin(w.group, word)
end


"""
    right_multiply(w::CoxElt, s::Integer) -> CoxElt

Return the product ``ws``, where ``s`` is the simple reflection indexed by ``s``.
"""
function right_multiply(w::CoxEltMin, s::T) where T <: Integer
    1 <= s <= rank(w.group) || throw(DomainError(s, "$s is not a generator in the range [1, $(rank(w.group))]"))
    word = copy(w.word)
    _mult_r!(w.group.refl_table, word, UInt8(s))
    return CoxEltMin(w.group, word)
end

"""
    left_multiply(s::Integer, w::CoxElt) -> CoxElt

Return the product ``sw``, where ``s`` is the simple reflection indexed by ``s``.
"""
function left_multiply(s::T, w::CoxEltMin) where T <: Integer
    1 <= s <= rank(w.group) || throw(DomainError(s, "$s is not a generator in the range [1, $(rank(w.group))]"))
    return CoxEltMin(w.group, [UInt8(s)]) * w
end


"""
    inv(w)

The inverse ``w^{-1}`` of a Coxeter group element."""
function Base.inv(w::CoxEltMin)
    refl = w.group.refl_table
    word = UInt8[]
    sizehint!(word, length(w.word))

    for s in reverse(w.word)
        _mult_r!(refl, word, s)
    end

    return CoxEltMin(w.group, word)
end

@doc raw"""
    create_reflection_table_gcm(gcm)

Create the minimal root reflection table for a Coxeter group defined by a GCM (generalised Cartan matrix).
This is a straightforward algorithm since we need only deal with integers, rather than algebraic integers.

The return value is a matrix `refl_table` of size R×T, where R is the rank of the group (the size of the square GCM),
and T is the number of minimal roots. The minimal roots are indexed from 1:T, with the first 1:R of them corresponding
to the simple roots, and the other indices arbitrary.

If ``β = s(α)`` is a minimal root, then `refl_table[s, α]` stores the index of β, and otherwise `refl_table[s, α] = 0`.
Note that `refl_table[s, s] = 0` for every simple root `s`.
"""
function create_reflection_table_gcm(gcm::Matrix{T}) where T <: Integer
    if !is_gcm(gcm)
        throw(DomainError(gcm, "Not a generalised Cartan matrix"))
    end

    # My convention for the Cartan matrix is that gcm[i, j] = ⟨α_i^, α_j⟩. Not that it particularly matters here,
    # since a transposed Cartan matrix yields the same Coxeter system.

    # We will keep track of root vectors in the simple root basis, and their associated coroots in the simple
    # coroot basis: both the simple roots and coroots start out as coordinate vectors. We then start generating
    # the root system by reflecting anywhere we can to get new roots. The pairing between the two bases is
    # using the Cartan matrix.
    # The condition ⟨β, α_s⟩ ≤ -2 which is used to detect whether s(β) is non-minimal in the symmetric case can
    # be replaced by the condition that ⟨β^, α_s⟩⟨α_s^, β⟩ ≥ 4. (Both these conditions are equivalent to the pair
    # {s, r_β} of reflections generating an infinite dihedral group).

    # We will keep the roots and coroots in parallel arrays, and keep an extra dictionary mapping a root
    # to its index in this array, which will also serve as a marker for whether we have seen a root before.
    n, _ = size(gcm)
    id = Matrix{T}(I, n, n)
    roots = [id[s, :] for s = 1:n]
    coroots = [id[s, :] for s = 1:n]
    rootidx = Dict(roots[s] => s for s = 1:n)

    # The reflection table maps (s, i) to the index of s(roots[i]). If s(roots[i]) is negative (this happens
    # iff roots[i] is the simple reflection corresponding to s, iff s = i) then refl[s, i] = 0, and if
    # s(roots[i]) is no longer a minimal root, then refl[s, i] = 0.
    refl = Dict((s, s) => 0 for s = 1:n)

    i = 1
    while i <= length(roots)
        for s = 1:n
            # If the reflection s(roots[i]) has already been calculated, skip. This will happen once
            # for each root, since it already knows what lower reflection first created it.
            if haskey(refl, (s, i))
                continue
            end

            # Calculate the pairing ⟨α_s^, roots[i]⟩ and the copairing ⟨coroots[i], α_s⟩.
            pairing = sum(gcm[s, :] .* roots[i])
            copairing = sum(coroots[i] .* gcm[:, s])
            if pairing * copairing >= 4
                refl[s, i] = 0
                continue
            end

            # Calculate the reflection s(roots[i]).
            new_root = copy(roots[i])
            new_root[s] -= pairing

            # If we've not seen this root before, then add it to our list.
            if !haskey(rootidx, new_root)
                push!(roots, new_root)

                new_coroot = copy(coroots[i])
                new_coroot[s] -= copairing
                push!(coroots, new_coroot)

                rootidx[new_root] = length(roots)
            end

            # Record the reflection data.
            si = rootidx[new_root]
            refl[s, i] = si
            refl[s, si] = i
        end
        i += 1
    end

    # Reassemble our reflection data into a matrix, where the rows are indexed by simple reflections
    # and the columns are indexed by roots. Iterating over all possible root indices and simple reflections
    # double-checks that we have defined all reflections.
    refltab = zeros(Int64, n, length(roots))
    for i = 1:length(roots), s = 1:n
        refltab[s, i] = refl[s, i]
    end

    return refltab
end
