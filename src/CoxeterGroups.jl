# By T. Schmit (2021)

module CoxeterGroups

import AbstractAlgebra: MatElem#, MatrixAlgebra
import Nemo: matrix, ZZ, fmpz
import Base: *,<,>
export CoxeterGroup, CoxeterElement, setCoxeterGroup, setRelation!, getRelation, getCoxeterMatrix, getCoxeterGroup, rank, generators

"Abstract supertype for Coxeter groups"
abstract type CoxGrp end

"Abstract supertype for Coxeter group elements"
abstract type CoxElt end

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
struct CoxeterGroup <: CoxGrp
    M::MatElem
    S::Array{String,1}

    function CoxeterGroup(M::MatElem)
        is_coxeter_matrix(M) || throw(ArgumentError("M has to be a Coxeter matrix"))

        S = ["<"*string(i)*">" for i=1:size(M,1)]
        CG = new(M,S)
        return CG, [CoxeterElement(CG, [i], true) for i=1:length(S)]
    end

    function CoxeterGroup(M::MatElem, S::Array{String, 1})
        is_coxeter_matrix(M) || throw(ArgumentError("M has to be a Coxeter matrix"))
        length(S)==size(M,1) || throw(ArgumentError("For M of size "*string(size(M))*", S should have length "*string(size(M,1))))

        CG = new(M,S)
        return CG, [CoxeterElement(CG, [i], true) for i=1:length(S)]
    end
end


"""
    CoxeterGroup(M::Array{Int,2}, S::Array{String, 1})
    CoxeterGroup(M::Array{Int,2})

Returns a CoxeterGroup together with its generating set.
If S is not specified, the generating set will be ["<1>","<2>","<3>",…]
"""
CoxeterGroup(M::Array{Int,2}, S::Array{String, 1}) = CoxeterGroup(matrix(ZZ,M), S)
CoxeterGroup(M::Array{Int,2}) = CoxeterGroup(matrix(ZZ, M))


"""
    CoxeterGroup(rank::Int)
    CoxeterGroup(rank::Int, S::Array{String,1})

Returns a CoxeterGroup for the CoxeterMatrix of size rank, with 1's on the diagonal, and 2's everywhere else.
If S is not specified, the generating set will be ["<1>", "<2>", "<3>", …]
"""
function CoxeterGroup(rank::Int)
    rank>=0 || throw(ArgumentError("rank has to be positive"))
    M = matrix(ZZ, fill(ZZ(2), rank, rank))
    for i = 1:rank
        M[i,i] = 1
    end
    return CoxeterGroup(M)
end
function CoxeterGroup(rank::Int, S::Array{String,1})
    rank>=0 || throw(ArgumentError("rank has to be positive"))
    rank==length(S) || throw(ArgumentError("S has to have length rank"))
    M = matrix(ZZ, fill(ZZ(2), rank, rank))
    for i = 1:rank
        M[i,i] = 1
    end
    return CoxeterGroup(M, S)
end


"""
    CoxeterElement

A **Coxeter element** is an element of a Coxeter group. We use CoxeterElement in order to facilitate the calculation in a given Coxeter group.
A CoxeterElement is always represented by its InverseShortLex form.
"""
struct CoxeterElement <: CoxElt
    G::CoxeterGroup
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
    function CoxeterElement(G::CoxeterGroup, w::Array{Int,1}, isNormal=false::Bool)
        if isNormal
            new(G, w)
        else
            new(G, NF(w,G))
        end
    end
end

#This function is only for convenience. It can only be called to produce the identity element
function CoxeterElement(G::CoxeterGroup, w::Array{Any,1}, isNormal=false::Bool)
    isempty(w) || throw(ArgumentError("w has to be of type Array{Int,1}; I received an element of type: "*string(typeof(w))))
    return CoxeterElement(G, Int[], true)
end


function Base.show(io::IO, mime::MIME"text/plain", w::CoxeterElement)
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

function Base.show(io::IO, mime::MIME"text/plain", C::CoxeterGroup)
    println(io, "Coxeter Group with Coxeter Matrix:")
    show(io, mime, C.M)
    println(io, "\nand generating Set:")
    println(io, C.S)
end

"The identity element of the Coxeter group."
function Base.one(C::CoxeterGroup)
    return CoxeterElement(C, Int[], true)
end

"The Coxeter generators."
function generators(C::CoxeterGroup)
    return [CoxeterElement(C, [i], true) for i=1:rank(C)]
end

function Base.size(w::CoxeterElement)
    return size(w.w)
end

function Base.length(w::CoxeterElement)
    return length(w.w)
end

function Base.getindex(w::CoxeterElement, i::Int)
    return getindex(w.w, i)
end

function Base.insert!(w::CoxeterElement, i::Int, s::Int)
    insert!(w.w, i, s)
end

function Base.deleteat!(w::CoxeterElement, i::Int)
    deleteat!(w.w, i)
end

function Base.copy(w::CoxeterElement)
    return CoxeterElement(w.G, copy(w.w), true)
end

function Base.:(==)(w::CoxeterElement, x::CoxeterElement)
    return w.G == x.G && w.w == x.w
end

function Base.hash(w::CoxeterElement, h::UInt)
    return hash(w.w, h)
end

"""
    <(x::CoxeterElement, y::CoxeterElement)

Returns true if x<y according to the InverseShortLex ordering.
"""
function <(x::CoxeterElement, y::CoxeterElement)
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
function >(x::CoxeterElement, y::CoxeterElement)
    return <(y,x)
end

"""
    *(x::CoxeterElement, y::CoxeterElement)

The operation of the Coxeter group. Returns the CoxeterElement in InverseShortLex form of x⋅y.
"""
function *(x::CoxeterElement, y::CoxeterElement)
    x.G.M==y.G.M || throw(ArgumentError("x and y do not belong to the same Coxeter group"))
    return NF(x,y)
end

"""
    setRelation!(CG::CoxeterGroup, i::Int, j::Int, m::Int)

Sets the relations m_{i,j} and m_{j,i} to m.
Use 0 as placeholder if m_{s,t} = infinity .
"""
function setRelation!(CG::CoxeterGroup, i::Int, j::Int, m::Int)
    (i<=rank(CG) && j<=rank(CG)) || throw(BoundsError("attempt to acces "*string(rank(CG))*"x"*string(rank(CG))*" Matrix at index ["*string(i)*","*string(j)*"]"))
    m>=0 || throw(ArgumentError("m has to be positive"))
    m==1 || i!=j || throw(ArgumentError("m_{i,i} has to equal 1"))
    CG.M[i,j] = m
    CG.M[j,i] = m
end

"""
    getRelation(CG::CoxeterGroup, i::Int, j::Int, m::Int)

Returns the relation m_{i,j} such that (<i><j>)^m_{i,j} = 1
"""
function getRelation(CG::CoxeterGroup, i::Int, j::Int)
    (i>0 && i<=rank(CG) && j>0 && j<=rank(CG)) || throw(BoundsError("attempt to acces "*string(rank(CG))*"x"*string(rank(CG))*" Matrix at index ["*string(i)*","*string(j)*"]"))
    return CG.M[i,j]
end

"""
    getCoxeterMatrix(CG::CoxeterGroup)

Returns the associated Coxeter matrix for CG
"""
function getCoxeterMatrix(CG::CoxeterGroup)
    return deepcopy(CG.M)
end

"""
    getCoxeterGroup(w::CoxeterElement)

Returns the Coxeter group in which w currently resides.
"""
function getCoxeterGroup(w::CoxeterElement)
    return w.G
end

"""
    rank(CG::CoxeterGroup)

Returns the size of the generating set for CG.
"""
function rank(CG::CoxeterGroup)
  return length(CG.S)
end

"""
    is_coxeter_matrix(M::MatElem)

returns true if and only if M is a **Coxeter matrix**.
We use 0 as placeholder if m_{s,t} = infinity
"""
function is_coxeter_matrix(M::MatElem)
    #is square
    if size(M,1) != size(M,2)
        return false
    end
    #main diagonal = 1
    for i = 1:size(M,1)
        if M[i,i] != 1
            return false
        end
    end
    #M[i,j]=M[j,i]  and  M[i,j]≥0
    for i = 1:size(M,1)
        for j = 1:i-1
            if M[i,j]!=M[j,i] || (M[i,j]<0)
                return false
            end
        end
    end
    return true
end

#returns the normal form of S[x] in CG
function NF(x::Array{Int,1}, CG::CoxeterGroup)
    if length(x)>1
        return NF(x[1:end-1], CoxeterElement(CG, [x[end]], true))
    else
        return CoxeterElement(CG, x, true)
    end
end

#returns the normal form of x*y
function NF(x::CoxeterElement, y::CoxeterElement)
    result = copy(y)
    for i = length(x):-1:1
        result = NF!(x[i], result)
    end
    return result
end

#returns the normal form of x*y
function NF(x::Array{Int,1}, y::CoxeterElement)
    result = copy(y)
    for i = length(x):-1:1
        result = NF!(x[i], result)
    end
    return result
end


function NF(s::Int, w::CoxeterElement)
    return NF!(s, copy(w))
end


function NF!(s::Int, w::CoxeterElement)
    n = length(w)
    l = 1
    r = 1
    σ = s
    while r <= n
        c = w[r]
        r = r + 1
        if c >= σ
            σ = c
            t = Exchange(s, CoxeterElement(w.G, w.w[l:r-1], true))
            if t != nothing
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


function Exchange(s::Int, w::CoxeterElement)
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
        sy = CoxeterElement(CG, pushfirst!(w[1:n-1], s), true)
        z = NF(w[1], sy) #s_1sy
        if length(z) < n
            m = 1
            while m<length(z) && z[m+1] == sy[m+1]
                m = m + 1
            end
            ω = NF(vcat(w[1], s, w[1:m-1]), CoxeterElement(CG, w[m+1:n], true))

            for t in ω
                if t > ω[end]
                    return ω[end]
                end
            end
        end
        return nothing
    end
end

end # module
