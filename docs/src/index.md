# CoxeterGroups.jl Documentation

This package is for creating and working with Coxeter groups.


## Creating a group

There are various implementations of Coxeter groups, each represented by a concrete type which is a subtype of the abstract [`CoxGrp`](@ref).
Use one of the constructors below to create a group, depending on which implementation is desired.
For example:

```julia
# A Coxeter group of type A2, implemented using minimal roots.
W, (s, t) = coxter_group_min([2, -1; -1, 2])

# The symmetric group S3, implemented via permutations.
S, (u, v) = symmetric_group(3)
```

Here is a summary of the creation functions - there is more information on specific implementations in the [Implementations](@ref) section.
```@docs
coxeter_group_min
symmetric_group
coxeter_group_recursive
```


## Operations on the group

A Coxeter group is a subtype of the abstract [`CoxGrp`](@ref) type.
Coxeter groups compare equal by reference.
The following operations are supported on all group implementations:

```@docs
CoxeterGroups.rank(::CoxGrpMin)
CoxeterGroups.generators(::CoxGrpMin)
Base.one(::CoxGrpMin)
CoxeterGroups.coxeter_matrix(::CoxGrpMin)
CoxeterGroups.is_finite(::CoxGrpMin)
CoxeterGroups.longest_element(::CoxGrp)
```

## Operations on group elements

A Coxeter group element is a subtype of the abstract [`CoxElt`](@ref) type.
Coxeter group elements are immutable value types: they compare equal if they represent the same group element, they can be used as keys in hash tables, and so on.
The following operations are supported on all Coxeter group element implementations:

```@docs
Base.parent(::CoxEltMin)
Base.isone(::CoxEltMin)
Base.:(*)(::CoxEltMin, ::CoxEltMin)
Base.inv(::CoxEltMin)
Base.:(^)(::CoxElt, ::Integer)
CoxeterGroups.length(::CoxEltMin)
Base.sign(::CoxElt)
CoxeterGroups.is_right_descent(::CoxEltMin, ::Integer)
CoxeterGroups.is_left_descent(::Integer, ::CoxEltMin)
CoxeterGroups.right_multiply(::CoxEltMin, ::Integer)
CoxeterGroups.left_multiply(::Integer, ::CoxEltMin)
CoxeterGroups.short_lex(::CoxEltMin)
CoxeterGroups.inverse_short_lex(::CoxEltMin)
```

## Internals

A Coxeter group implementation consists of a new pair of concrete types subtyping [`CoxGrp`](@ref) and [`CoxElt`](@ref) respectively.

```@docs
CoxGrp
CoxElt
```

When adding a new concrete implementation, a small set of functions must be implemented explicitly, then generic implementations of other functions will begin working.
The concrete implementation should specialise these generic implementations where it is possible to take advantage of specific internal structure.

- The group type must implement a constructor such as those in [Creating a group](@ref), and the functions [`coxeter_matrix`](@ref) [`rank`](@ref), [`one`](@ref), [`generators`](@ref), and [`is_finite`](@ref).
- The element type must implement [`parent`](@ref), `==`, `hash`, [`is_left_descent`](@ref), [`is_right_descent`](@ref), [`left_multiply`](@ref), and [`right_multiply`](@ref).

The following functions will then be defined:

- For the group type: [`longest_element`](@ref).
- For the element type: [`isone`](@ref), [`length`](@ref), [`sign`](@ref), [`short_lex`](@ref), [`inverse_short_lex`](@ref), [`*`](@ref), [`inv`](@ref), [`^`](@ref).


## Code conventions

Code conventions are (currently) as follows:

- Library users should access data via methods, not via properties.
- A Coxeter group implementation called `foo` should define a concrete subtype `CoxGrpFoo` of `CoxGrp`, and a concrete subtype `CoxEltFoo` of `CoxElt`.

## Mathematical conventions

A *generalised Cartan matrix* or *GCM* is a square matrix ``A \in \operatorname{Mat}_I(\mathbb{Z})`` over the integers which satisfies the following properties:

- *(C1)* All diagonal elements are equal to ``2``,
- *(C2)* All off-diagonal entries are zero or negative, and
- *(C3)* ``A_{ij} = 0`` if and only if ``A_{ji} = 0``.

A *Coxeter-Cartan matrix* is a square matrix ``A`` over the reals which satisfies the conditions C1, C2, C3 of a GCM, and additionally

- *(C4)* The product ``A_{ij} A_{ji}`` is either equal to ``4 \cos^2(\pi/m)`` for some ``m \in \{2, 3, \ldots\}``, or ``A_{ij} A_{ji} \geq 4``.


A *based root datum* over the PID ``R`` (we will only be thinking of ``R = \mathbb{Z}`` or ``R = \mathbb{R}``) is the data of:

- Two free ``R``-modules ``X`` and ``X^\vee``,
- A perfect pairing ``\langle -, - \rangle \colon X^\vee \times X \to R``,
- A collection of *simple roots* ``\Delta = \{\alpha_i\}_{i \in I} \in X``,
- A collection of *simple coroots* ``\Delta^\vee = \{\alpha_i^\vee\}_{i \in I} \in X^\vee``,

such that the matrix ``A_{ij} := \langle \alpha_i^\vee, \alpha_j \rangle`` is a GCM.
We then say that the based root datum is a *realisation* of the GCM ``A``.
Define the reflections ``r_i \colon X \to X`` and ``r_i^\vee \colon X^\vee \to X^\vee`` by the formulas

```math
r_i(\lambda) = \lambda - \langle \alpha_i^\vee, \lambda \rangle \alpha_i, \quad
r_i^\vee(\nu) = \nu - \langle \nu, \alpha_i \rangle \alpha_i^\vee.
```

Then the reflections each define a faithful representation of a Coxeter group ``(W, I)`` on each of ``X`` and ``X^\vee``, which are dual in the sense that ``\langle w \nu, \lambda \rangle = \langle \nu, w \lambda \rangle``.

Every GCM ``A`` admits two canonical based root data.
The *simply-connected* realisation is the root datum where ``X^\vee`` is the free ``R``-module with basis ``(\alpha_i^\vee)_{i \in I}``, ``X = \operatorname{Hom}(X^\vee, R)`` is the dual, the pairing ``\langle -, - \rangle`` is the canonical pairing between a module and its dual, and the simple roots ``\alpha_i`` are uniquely determined by the condition that ``\langle \alpha_i^\vee, \alpha_j \rangle = A_{ij}``.
The *adjoint* realisation is the dual of the simply-connected realisation.

In the reflection representation implementation, we use the adjoint realisation of a GCM.
Concretely, the space ``X`` is the space of column vectors, with each simple root ``\alpha_i`` a coordinate vector, the space ``X^\vee`` is the space of row vectors, with each simple coroot ``\alpha_i^\vee`` equal to the ``i``th row of the Cartan matrix, and the pairing ``\langle -, - \rangle`` is the dot product.



## Implementations

### Minimal roots implementation

The _minimal roots_ implementation of a Coxeter group represents group elements by words in the Coxeter generators in ShortLex normal form.
Operations on these words are done with the assistance of a _minimal root reflection table_, which needs to be constructed when the group is created.

Currently this implementation can handle any Coxeter group defined by a generalised Cartan matrix, but this is only because generating the minimal root reflection table is more difficult for general Coxeter groups.
The same underlying algorithms will work once given a new table.

Performance characteristics for a Coxeter group of rank ``r``:

- The group itself takes ``O(r M)`` space, where ``M`` is the number of minimal roots.
  For finite types, the minimal roots are exactly the positive roots, and for affine types, there are twice as many minimal roots as there are positive roots for the corresponding finite system.
  The number of roots in a finite system tends to be ``O(r^2)`` in the rank ``r``.
- Length is an ``O(1)`` operation.
- ShortLex normal form should be ``O(1)``, but is ``O(l(w))`` since it returns a defensive copy of the internal word.
- Testing a left or right descent is ``O(l(w))``.
- Multiplication on the right ``xs`` by a simple generator ``s`` is ``O(l(x))``.
- Multiplication on the left ``sx`` by a simple generator ``s`` is ``O(l(x)^2)``.
- Multiplication ``xy`` is ``O(l(y)(l(x) + l(y)))``.
- Inversion ``x^{-1}`` is ``O(l(x)^2)``.

```@docs
CoxGrpMin
CoxeterGroups.create_reflection_table_gcm
```

### Symmetric groups

The symmetric group ``S_n`` on ``n`` letters is the group of permutations of the set ``[n] = \{1, \ldots, n\}``, whose standard Coxeter generators are ``s_1, \ldots, s_{n-1}`` where ``s_i = (i, i+1)`` is an adjacent transposition.
The symmetric group may be constructed using the [`symmetric_group`](@ref) function.


A permutation ``\pi \colon [n] \to [n]`` is represented as a pair of arrays recording both ``\pi`` and its inverse ``\pi^{-1}``:
```math
    \pi \mapsto (\mathtt{fwd} = [\pi(1), \ldots, \pi(n)], \mathtt{inv} = [\pi^{-1}(1), \ldots, \pi^{-1}(n)]).
```

Complexity of operations in terms of ``n``:

- Creation of the group is ``O(1)``.
- Memory usage for each element is ``O(n)``.
- Testing a left or right descent is ``O(1)``.
- Multiplication is ``O(n)``.
- Calculating the sign of a permutation is ``O(n)``
- Taking Coxeter length takes ``O(n^2)`` time (this could be improved to ``O(n \log n)`` by a standard divide-and-conquer trick to count inversions).
