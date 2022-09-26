module CoxeterGroups

# Abstract supertypes
export CoxGrp, CoxElt

# Group functions
export coxeter_matrix, rank, generators, longest_element, is_finite

# Element functions
export length, sign, is_right_descent, is_left_descent, right_multiply, left_multiply, short_lex, inverse_short_lex


"""The abstract supertype for Coxeter groups."""
abstract type CoxGrp end

"""The abstract supertype for Coxeter group elements."""
abstract type CoxElt end

include("CoxeterMatrices.jl")
include("CoxeterGroupData.jl")
include("CoxeterSystems.jl")
include("CoxGrpMin.jl")
include("CoxGrpMin_ReflTable.jl")
include("CoxGrpSym.jl")
include("CoxGrpRec.jl")

# Default implementations of some functions.

"""
    isone(w)

Return true if the group element is the identity."""
Base.isone(w::CoxElt) = w == one(parent(w))

"Return the numerically least left descent of w, or 0 if there is no left descent (iff w = id)."
function first_left_descent(w::CoxElt)
    for i in 1:rank(parent(w))
        if is_left_descent(i, w)
            return i
        end
    end
    return 0
end

"Return the numerically least right descent of w, or 0 if there is no right descent (iff w = id)."
function first_right_descent(w::CoxElt)
    for i in 1:rank(parent(w))
        if is_right_descent(w, i)
            return i
        end
    end
    return 0
end

"""
    length(w::CoxElt)

The Coxeter length of ``w``: the length of a shortest expression in the Coxeter generators for ``w``."""
function Base.length(w::CoxElt)
    length = 0
    while true
        s = first_left_descent(w)
        s == 0 && return length
        length += 1
        w = left_multiply(s, w)
    end
end

"""
    sign(w::CoxElt) -> Int

The sign homomorphism ``(-1)^{l(w)}``.
"""
Base.sign(w::CoxElt) = length(w) % 2 == 0 ? 1 : -1

"""
    short_lex(w)

Return the ShortLex normal form of w: the lexicographically least reduced expression for w.
"""
function short_lex(w::CoxElt)
    word = []
    while true
        s = first_left_descent(w)
        s == 0 && return word
        push!(word, s)
        w = left_multiply(s, w)
    end
end

"""
    inverse_short_lex(w)

Return the InverseShortLex normal form of `w`: the lexicographically least reduced expression for ``w``
when read backwards, or equivalently the reverse of the ShortLex normal form of ``w^{-1}``.
"""
function inverse_short_lex(w::CoxElt)
    word = []
    while true
        s = first_right_descent(w)
        s == 0 && return reverse!(word)
        push!(word, s)
        w = right_multiply(w, s)
    end
end

"""
    *(x, y)

Take the product of two Coxeter group elements."""
function Base.:(*)(x::CoxElt, y::CoxElt)
    parent(x) == parent(y) || error("Group elements $x and $y do not come from the same parent.")
    while true
        s = first_left_descent(y)
        s == 0 && return x
        x = right_multiply(x, s)
        y = left_multiply(s, y)
    end
end

# The product over one element sometimes comes up, eg when doing *(gens...) to get a Coxeter element.
Base.:(*)(x::CoxElt) = x

"""
    inv(x)

The inverse of a Coxeter group element."""
function Base.inv(x::CoxElt)
    inv = one(parent(x))
    while true
        s = first_left_descent(x)
        s == 0 && return inv
        x = left_multiply(s, x)
        inv = left_multiply(s, inv)
    end
end

"""
    ^(x::CoxElt, e::T) where T <: Integer

Return the power ``x^e``."""
function Base.:(^)(x::CoxElt, e::T) where T <: Integer
    e == 0 && return one(parent(x))
    if e < 0
        x = inv(x)
        e = -e
    end

    e == 1 && return x

    # Exponentiate by squaring
    acc = one(parent(x))
    pow2 = x
    while e != 0
        if e % 2 == 1
            acc *= pow2
        end
        e = div(e, 2)
        pow2 *= pow2
    end
    return acc
end

"""
    longest_element(W::CoxGrp)

Return the longest element if ``W`` is finite, otherwise error."""
function longest_element(grp::CoxGrp)
    is_finite(grp) || error("longest_element only works for finite Coxeter groups")
    w = one(grp)
    while true
        finished = true
        for i in 1:rank(grp)
            if !is_right_descent(w, s)
                w = right_multiply(w, s)
                finished = false
            end
        end

        finished && return w
    end
end




end # module
