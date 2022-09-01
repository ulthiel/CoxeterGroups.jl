# CoxeterGroups.jl Documentation

Coxeter groups and their elements are subtypes of the abstract types [`CoxGrp`](@ref) and [`CoxElt`](@ref).
Currently there are two (incomplete) implementations, which can be constructed using either [`CoxeterGroup`](@ref) or [`coxeter_group_min`](@ref).
The constructor takes a matrix, which can be either a generalised Cartan matrix, or a Coxeter matrix (currently `coxeter_group_min` cannot accept a Coxeter matrix), and returns a group object and list of generators.
For example:

```
W, (s, t) = coxeter_group_min([2 -1; -1 2])
s*t*s == t*s*t
```

The interface for Coxeter group objects is that they are subtypes of [`CoxGrp`](@ref), and implement the functions

- `CoxeterGroups.rank(W)` returns the rank, i.e. number of generators.
- `CoxeterGroups.generators(W)` returns a list of generators.
- `Base.one(W)` returns the group identity.
- Coxeter group objects compare equal `==` only if they are exactly the same object in memory. Isomorphic groups are not equal.

The interface for Coxeter group elements is they are subtypes of [`CoxElt`](@ref), and implement the functions

- `Base.parent(w)` returns the parent Coxeter group.
- `Base.isone(w)` to check if an element is the identity.
- `Base.length(w)` to return the Coxeter length.
- Coxeter group elements compare equal `==` if and only if they represent the same element of the same parent. Coxeter group elements are immutable and work like values: they can be hashed and used as keys in a dictionary.
- `*`, `inv`, and `w^k` should work on group elements.
- `CoxeterGroups.is_right_descent(w, s)` returns true if the integer `s`, which represents a simple generator, is in the right descent set of `w`.
- `CoxeterGroups.is_left_descent(s, w)` works similarly.

```@docs
CoxGrp
CoxElt
coxeter_group_min
CoxeterGroup
```


## Minimal roots implementation

The _minimal roots_ implementation of a Coxeter group represents group elements by words in the Coxeter generators in ShortLex normal form.
Operations on these words are done with the assistance of a _minimal root reflection table_, which needs to be constructed when the group is created.

Currently this implementation can handle any Coxeter group defined by a generalised Cartan matrix, but this is only because generating the minimal root reflection table is more difficult for general Coxeter groups.
The same underlying algorithms will work once given a new table.

```@docs
CoxeterGroups.CoxGrpMin
CoxeterGroups.create_reflection_table_gcm
```
