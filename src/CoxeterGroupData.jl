import LinearAlgebra: diagind

export coxeter_type, is_coxeter_matrix

"""
    is_coxeter_matrix(mat)

Check if an integer matrix is a Coxeter matrix. A Coxeter matrix is a square symmetric matrix with integer entries,
where the diagonal elements are 1, and then off-diagonal elements are in {0} ∪ {2, 3, 4, ...}. Under this convention,
a zero represents ∞.
"""
function is_coxeter_matrix(mat::Matrix{T}) where T <: Integer
    n, m = size(mat)
    return (
        n == m
        && all(mat[i, i] == 1 for i = 1:n)
        && all(mat[i, j] == mat[j, i] for i = 1:n, j = 1:n if i != j)
        && all(mat[i, j] == 0 || mat[i, j] >= 2 for i = 1:n, j = 1:n if i != j)
    )
end


"""
    is_gcm(mat)

Check if an integer matrix is a generalised Cartan matrix. A generalised Cartan matrix (GCM) is a square
matrix with integer entries, such that:

1. All diagonal entries are 2,
2. All off-diagonal entries are 0 or negative, and
3. ``a_{ij} = 0`` if and only if ``a_{ji} = 0``.
"""
function is_gcm(mat::Matrix{T}) where T <: Integer
    n, m = size(mat)
    return (
        n == m
        && all(mat[i, i] == 2 for i = 1:n)
        && all(mat[i, j] <= 0 for i = 1:n, j = 1:n if i != j)
        && all(mat[i, j] == 0 for i = 1:n, j = 1:n if mat[j, i] == 0)
    )
end


"""
    gcm_to_coxeter_matrix(gcm)

Convert a generalised Cartan matrix to its corresponding Coxeter matrix."""
function gcm_to_coxeter_matrix(gcm::Matrix{T}) where T <: Integer
    if !is_gcm(gcm)
        throw(DomainError(gcm, "Is not a generalised Cartan matrix"))
    end

    # Replace each entry with a_{ij} a_{ji}, then map (0, 1, 2, 3, ≥ 4) to m_{ij} = (2, 3, 4, 6, 0 = ∞).
    # Finally, replace the diagonal with 1's and return.
    function a_to_m(a)
        a == 0 && return 2
        a == 1 && return 3
        a == 2 && return 4
        a == 3 && return 6
        return 0
    end

    cartan_mat = map(a_to_m, gcm .* gcm)
    cartan_mat[diagind(cartan_mat)] .= 1
    return cartan_mat
end


"""A CoxType holds metadata about a Coxeter group"""
struct CoxType
    """Rank of the Coxeter group."""
    rank::Int

    """Generalised Cartan matrix. Omitted for non-crystallographic types, like H3 and H4."""
    gcm::Matrix{Int}

    """Coxeter matrix."""
    coxeter_mat::Matrix{Int}

    """Whether the group is finite."""
    # is_finite::Bool

    """The number of reflections in the group (half the number of roots). Omitted for infinite groups."""
    # refl_count::Int

    """The number of elements in the group. Omitted for infinite groups."""
    weyl_order::BigInt

    """The exponents of the group, in increasing order. Omitted for infinite groups."""
    # exponents::Vector{Int}
end

"""Convenience method for constructing a CoxType from the GCM and weyl order only."""
function coxeter_type_from_gcm(gcm::Matrix{Int}, weyl_order::BigInt)
    @assert is_gcm(gcm)

    return CoxType(size(gcm, 1), gcm, gcm_to_coxeter_matrix(gcm), weyl_order)
end

"""Construct a CoxType from a letter and rank, like ("A", 4) for the symmetric group S5."""
function coxeter_type(letter::String, rank::Int)
    if letter == "A"
        gcm = zeros(Int64, rank, rank)
        gcm[diagind(gcm)] .= 2
        gcm[diagind(gcm, 1)] .= -1
        gcm[diagind(gcm, -1)] .= -1
        return coxeter_type_from_gcm(gcm, factorial(big(rank + 1)))
    end

    throw(DomainError(letter, "Unknown letter"))
end
