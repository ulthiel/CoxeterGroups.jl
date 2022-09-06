export CoxGrpSym, CoxEltSym, symmetric_group

"""
CoxGrpSym implements Coxeter-theoretic operations on symmetric groups.
"""
struct CoxGrpSym <: CoxGrp
    n::Int

    CoxGrpSym(n) = (n >= 0) ? new(n) : error("Negative argument $n to CoxGrpSym.")
end

"""
A CoxEltSym is the element type for CoxGrpSym. A permutation σ: [n] → [n] is represented
as a pair of arrays, fwd = [σ(1), ..., σ(n)] and inv = [σ^-1(1), ..., σ^-1(n)]. Keeping both
arrays makes inversion and checking left/right descents of a particular generator into O(1) operations.

Multiplication σ * π is treated as composition (σ ∘ π) of permutations, therefore multiplying σ
on the right by the transposition (i, j) has the effect of swapping fwd[i] and fwd[j]. Right Multiplication
by (i, j) has the effect of swapping the positions of i and j in the inv array, i.e. if a and b are the
indices such that inv[a] = i and inv[b] = j, then a and b are swapped.

Currently the CoxGrpSym type is just an integer, so the group reference could be replaced by an integer,
but later we may want to add extra metadata to the group, like letters used for generators.
"""
struct CoxEltSym <: CoxElt
    group::CoxGrpSym
    fwd::Vector{Int}
    inv::Vector{Int}
end

@doc raw"""
    symmetric_group(n) -> CoxGrpSym, gens

Create the symmetric group on ``n`` elements (a Coxeter system of type ``A_{n-1}``), for ``n \geq 0``.
"""
function symmetric_group(n)
    n >= 0 || error("The argument $n to symmetric_group must be nonnegative.")
    G = CoxGrpSym(n)
    return G, generators(G)
end

# Note that both S0 and S1 are the trivial group, and so have Coxeter rank 0.
rank(grp::CoxGrpSym) = max(0, grp.n - 1)

is_finite(grp::CoxGrpSym) = true

Base.one(grp::CoxGrpSym) = CoxEltSym(grp, collect(1:grp.n), collect(1:grp.n))

function generators(grp::CoxGrpSym)
    fwds = [collect(1:grp.n) for _ in 1:grp.n-1]
    for (i, fwd) in enumerate(fwds)
        fwd[i], fwd[i+1] = fwd[i+1], fwd[i]
    end
    return [CoxEltSym(grp, fwd, fwd) for fwd in fwds]
end

Base.parent(w::CoxEltSym) = w.group

Base.isone(w::CoxEltSym) = w.fwd == 1:w.group.n

Base.:(==)(w::CoxEltSym, x::CoxEltSym) = w.group.n == x.group.n && w.fwd == x.fwd

Base.hash(w::CoxEltSym, h::UInt) = hash(w.fwd, h)

function Base.:(*)(w::CoxEltSym, x::CoxEltSym)
    w.group.n == x.group.n || error("Incompatible parent groups")
    return CoxEltSym(w.group, [w.fwd[i] for i in x.fwd], [x.inv[i] for i in w.inv])
end

Base.inv(w::CoxEltSym) = CoxEltSym(w.group, w.inv, w.fwd)

# Coxeter length. Currently the O(n^2) straightforward algorithm. If we decide to make
# CoxGrpSym handle very large n, this should be replaced with the O(n log n) inversion
# counting merge-sort algorithm, or something as fast.
Base.length(w::CoxEltSym) = sum(1 for i in 1:length(w.fwd) for j in i+1:length(w.fwd) if w.fwd[i] > w.fwd[j]; init = 0)

longest_element(grp::CoxGrpSym) = CoxEltSym(grp, collect(grp.n:-1:1), collect(grp.n:-1:1))

is_right_descent(w::CoxEltSym, t::T) where {T <: Integer} = w.fwd[t] > w.fwd[t+1]
is_left_descent(t::T, w::CoxEltSym) where {T <: Integer} = w.inv[t] > w.inv[t+1]

function right_multiply(w::CoxEltSym, t::T) where T <: Integer
    fwd = copy(w.fwd)
    inv = copy(w.inv)
    inv[fwd[t]], inv[fwd[t+1]] = inv[fwd[t+1]], inv[fwd[t]]
    fwd[t], fwd[t+1] = fwd[t+1], fwd[t]
    return CoxEltSym(w.group, fwd, inv)
end

function left_multiply(t::T, w::CoxEltSym) where T <: Integer
    fwd = copy(w.fwd)
    inv = copy(w.inv)
    fwd[inv[t]], fwd[inv[t+1]] = fwd[inv[t+1]], fwd[inv[t]]
    inv[t], inv[t+1] = inv[t+1], inv[t]
    return CoxEltSym(w.group, fwd, inv)
end


# Calculate the ShortLex word associated to an inverse permutation array in O(coxeter length) time.
# In other words, to get the ShortLex normal form of (fwd, inv), call short_lex_array(inv).
function short_lex_array(inv::Array{T}) where T <: Integer
    # We'll be mutating the array, so take a copy.
    inv_mut = copy(inv)

    shortlex = T[]
    i = 1
    while i <= length(inv_mut) - 1
        if inv_mut[i] > inv_mut[i+1]
            push!(shortlex, i)
            inv_mut[i], inv_mut[i+1] = inv_mut[i+1], inv_mut[i]
            i = max(i - 1, 1)
        else
            i += 1
        end
    end

    return shortlex
end

inverse_short_lex(w::CoxEltSym) = reverse!(short_lex_array(w.fwd))
short_lex(w::CoxEltSym) = short_lex_array(w.inv)
