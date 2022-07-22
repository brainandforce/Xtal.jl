# Getting started

# Getting Xtal.jl

As of this release, Xtal.jl is not in the Julia package registry. You'll need to manually add this
repo to your Julia environment (which should be at least v1.6):

```julia-repl
(@v1.7 pkg)> add https://github.com/brainandforce/Xtal.jl
```

# Crystals

The `AbstractCrystal` supertype contains two concrete types, `Crystal{D}` and 
`CrystalWithDatasets{D,K,V}`. The `CrystalWithDatasets` type is a combination of a `Crystal{D}`
with a `Dict{K,V}`. You can access datasets with `getindex()`, so a `CrystalWithDatasets` behaves
similarly to a dictionary. It should be noted that autocompletion is not currently supported for
`CrystalWithDatasets`, so tabbing not be useful.

# Lattices

Xtal.jl provides types for lattices of all dimensionalities. All types listed here take a type
parameter `D` which represents the dimensionality. We aim for this to be useful in calculations
on highly complex structures with incommensurate modulation that may be best described in a higher-
dimensional superspace.

## Basis vectors

The `BasisVectors` type is used to store basis vectors.

If a system is not periodic, the zero element is used to represent this by convention.

### Why not use `SMatrix` for basis vectors?

While the `SMatrix` type seems to make sense for a collection of basis vectors, it poses one major 
problem: the calculation of type parameters.

The `SMatrix{D1,D2,T,L}` type requires four type parameters - `D1` and `D2` are each of the matrix
dimensions, and `T` is the element type of the matrix. However, one last parameter is needed: `L`, 
the length of the `NTuple` that backs the `SMatrix`.

Julia currently does not allow for the calculation of type parameters from other type parameters,
which poses a serious problem in the declaration of structs. Technically, a type like
`SMatrix{3,3,Float64}` is an abstract type, as the `L` parameter is undeclared, even though the 
value of `L` can be inferred from `D1` and `D2`. By declaring a struct to have this type, there
seems to be a significant performance drop.

The `BasisVectors{D}` type wraps an `SVector{D,SVector{D,Float64}}`. This only requires a single 
type parameter `D`, and allows for fully concrete struct declarations.

Other packages like [Crystalline.jl](https://github.com/thchr/Crystalline.jl) face this issue and
have devised their own solution. In the future, we may switch back to `SMatrix` as the default type
for basis vectors by using custom macros, or a package like
[ComputedFieldTypes.jl](https://github.com/vtjnash/ComputedFieldTypes.jl) which permits structs to
have computed type parameters.

### Real and reciprocal lattices

Xtal.jl provides real and reciprocal lattices as their own types: `RealLattice` and 
`ReciprocalLattice`. Both of these are instances of `AbstractLattice`, and they may be freely
converted between each other.

`AbstractLattice` types consist of pairs of `BasisVectors`, a primitive set of basis vectors which
may be accessed with `prim()` and a conventional set of basis vectors which may be accessed with
`conv()`.

# Datasets

Xtal.jl supports a good number of different data types, including:
  * Real space datagrids
  * Reciprocal space data by HKL index
  * Band structures and densities of states
  * k-point lists and grids
  * Data by atomic position
  * Spherical harmonic coefficients
