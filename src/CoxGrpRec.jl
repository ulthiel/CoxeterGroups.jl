# By T. Schmit (2021)

export CoxGrpRec, CoxEltRec, coxeter_group_recursive

"""
    CoxeterGroup

A **Coxeter group** is a group **G** together with a generating subset **S**, whose elements are subject to relations:
```math
\\begin{aligned}
s² &= 1 \\quad \\text{ for } s∈S \\\\

(st)^{m_{s,t}} = 1 \\quad \\text{ for } s≠t∈S,\\; m_{s,t}∈\\mathbb{N}_{≥2} \\\\
\\end{aligned}
```
A CoxeterGroup can be defined by its Coxeter matrix M=(m_{s,t})_{s,t∈S} and the generating subset **S**.
Elements in CoxeterGroup are represented by their expression in InverseShortLex form.
Once the CoxeterGroup has been established, you can calculate on it using "*".

# Example
```julia-repl
julia> G,(a,b,c) = CoxeterGroup([1 4 2; 4 1 3; 2 3 1], ["a", "b", "c"])
(CoxeterGroup([1 4 2; 4 1 3; 2 3 1], ["a", "b", "c"]), CoxeterElement[[1], [2], [3]])
julia> w = a*c
c*a
julia> w*c
a
```
"""
struct CoxGrpRec <: CoxGrp
    M::Matrix{Int64}
    S::Vector{String}
end

"""
    coxeter_group_recursive(M, [names])

Create a Coxeter group with elements represented in InverseShortLex normal form, with group operations
implemented with a recursive algorithm.
"""
function coxeter_group_recursive(M::Matrix{Int64}, S::Array{String, 1})
    if is_gcm(M)
        M = gcm_to_coxeter_matrix(M)
    end
    is_coxeter_matrix(M) || throw(ArgumentError("M has to be a Coxeter matrix"))
    length(S) == size(M,1) || throw(ArgumentError("For M of size $(size(M)), S should have length $(size(M,1))"))

    CG = CoxGrpRec(M,S)
    return CG, [CoxEltRec(CG, [i], true) for i=1:length(S)]
end

function coxeter_group_recursive(M::Matrix{Int64})
    return coxeter_group_recursive(M, ["<$i>" for i=1:size(M,1)])
end

"""
    CoxeterElement

A **Coxeter element** is an element of a Coxeter group. We use CoxeterElement in order to facilitate the calculation in a given Coxeter group.
A CoxeterElement is always represented by its InverseShortLex form.
"""
struct CoxEltRec <: CoxElt
    G::CoxGrpRec
    w::Array{Int,1}
    """
        CoxeterElement(G::CoxeterGroup, w::Array{Int,1}, isNormal=false::Bool)

    Returns the CoxeterElement corresponding to S[w] for S the generating set of G.
    If w is already in InverseShortLex form, you can set isNormal=true in order to save on computation.

    # Example
    ```julia-repl
    julia> G,S = CoxeterGroup([1 4 2; 4 1 3; 2 3 1], ["a", "b", "c"])
    (CoxeterGroup([1 4 2; 4 1 3; 2 3 1], ["a", "b", "c"]), CoxeterElement[[1], [2], [3]])
    julia> w = CoxeterElement(G, [1,2,3])
    abc
    julia> v = CoxeterElement(G, [3,2,3])
    bcb
    ```
    """
    function CoxEltRec(G::CoxGrpRec, w::Array{Int,1}, isNormal=false::Bool)
        if isNormal
            new(G, w)
        else
            new(G, NF(w,G))
        end
    end
end

#This function is only for convenience. It can only be called to produce the identity element
function CoxEltRec(G::CoxGrpRec, w::Array{Any,1}, isNormal=false::Bool)
    isempty(w) || throw(ArgumentError("w has to be of type Array{Int,1}; I received an element of type $(typeof(w))"))
    return CoxEltRec(G, Int[], true)
end


function Base.show(io::IO, mime::MIME"text/plain", w::CoxEltRec)
    output = ""
    for s in w
        output *= w.G.S[s]
    end
    if length(w) > 0
        print(io, output)
    else
        print(io, "<>")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", C::CoxGrpRec)
    println(io, "Coxeter Group with Coxeter Matrix:")
    show(io, mime, C.M)
    println(io, "\nand generating Set:")
    println(io, C.S)
end

"""The identity element of the Coxeter group."""
Base.one(C::CoxGrpRec) = CoxEltRec(C, Int[], true)

"""The Coxeter generators."""
generators(C::CoxGrpRec) = [CoxEltRec(C, [i], true) for i=1:rank(C)]


Base.length(w::CoxEltRec) = length(w.w)

Base.getindex(w::CoxEltRec, i) = getindex(w.w, i)


function Base.insert!(w::CoxEltRec, i::Int, s::Int)
    insert!(w.w, i, s)
end

function Base.deleteat!(w::CoxEltRec, i::Int)
    deleteat!(w.w, i)
end

Base.copy(w::CoxEltRec) = CoxEltRec(w.G, copy(w.w), true)

Base.:(==)(w::CoxEltRec, x::CoxEltRec) = w.G == x.G && w.w == x.w

Base.hash(w::CoxEltRec, h::UInt) = hash(w.w, h)

"""
    <(x::CoxeterElement, y::CoxeterElement)

Returns true if x<y according to the InverseShortLex ordering.
"""
function Base.:(<)(x::CoxEltRec, y::CoxEltRec)
    x.G.M==y.G.M || throw(ArgumentError("x and y do not belong to the same Coxeter group"))
    if length(x) < length(y)
        return true
    elseif length(x) > length(y)
        return false
    else
        i = length(x)
        while i>0
            if x[i] == y[i]
                i-=1
            else
                return (x[i]<y[i])
            end
        end
        return false
    end
end

"""
    >(x::CoxeterElement, y::CoxeterElement)

Returns true if x>y according to the InverseShortLex ordering.
"""
function Base.:(>)(x::CoxEltRec, y::CoxEltRec)
    return <(y,x)
end

"""
    *(x::CoxeterElement, y::CoxeterElement)

The operation of the Coxeter group. Returns the CoxeterElement in InverseShortLex form of x⋅y.
"""
function Base.:(*)(x::CoxEltRec, y::CoxEltRec)
    x.G.M==y.G.M || throw(ArgumentError("x and y do not belong to the same Coxeter group"))
    return NF(x,y)
end

"""
    coxeter_matrix(CG::CoxeterGroup)

Returns the associated Coxeter matrix for CG
"""
coxeter_matrix(CG::CoxGrpRec) = deepcopy(CG.M)

"""
    parent(w::CoxeterElement)

Returns the Coxeter group in which w currently resides.
"""
Base.parent(w::CoxEltRec) = w.G

"""
    rank(CG::CoxeterGroup)

Returns the size of the generating set for CG.
"""
rank(CG::CoxGrpRec) = length(CG.S)




#returns the normal form of S[x] in CG
function NF(x::Array{Int,1}, CG::CoxGrpRec)
    if length(x)>1
        return NF(x[1:end-1], CoxEltRec(CG, [x[end]], true))
    else
        return CoxEltRec(CG, x, true)
    end
end

#returns the normal form of x*y
function NF(x::CoxEltRec, y::CoxEltRec)
    result = copy(y)
    for i = length(x):-1:1
        result = NF!(x[i], result)
    end
    return result
end

#returns the normal form of x*y
function NF(x::Array{Int,1}, y::CoxEltRec)
    result = copy(y)
    for i = length(x):-1:1
        result = NF!(x[i], result)
    end
    return result
end


function NF(s::Int, w::CoxEltRec)
    return NF!(s, copy(w))
end


function NF!(s::Int, w::CoxEltRec)
    n = length(w)
    l = 1
    r = 1
    σ = s
    while r <= n
        c = w[r]
        r = r + 1
        if c >= σ
            σ = c
            t = Exchange(s, CoxEltRec(w.G, w.w[l:r-1], true))
            if t !== nothing
                if s == w[l]
                    deleteat!(w,l)
                    return(w)
                else
                    s = t
                    l = r
                    σ = t
                end
            end
        end
    end
    insert!(w, l, s)
    return(w)
end


function Exchange(s::Int, w::CoxEltRec)
    CG = w.G
    M = CG.M
    n = length(w)
    if n == 1
        if w[1] == s
            return s
        elseif w[1]<s || M[w[1],s]!=2
            return nothing
        else
            return s
        end
    end
    if n == M[w[1], s]-1
        i = 2
        while i<=n && (w[i]==s || w[i]==w[1])
            i += 1
        end
        if i == n + 1
            if w[n] > w[n-1]
                return w[n-1]
            else
                return nothing
            end
        end
    else
        sy = CoxEltRec(CG, pushfirst!(w[1:n-1], s), true)
        z = NF(w[1], sy) #s_1sy
        if length(z) < n
            m = 1
            while m<length(z) && z[m+1] == sy[m+1]
                m = m + 1
            end
            ω = NF(vcat(w[1], s, w[1:m-1]), CoxEltRec(CG, w[m+1:n], true))

            for t in ω.w
                if t > ω.w[end]
                    return ω.w[end]
                end
            end
        end
        return nothing
    end
end
