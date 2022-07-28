# REAL SPACE
#-------------------------------------------------------------------------------------------------#

"""
    RealSpaceDataGrid{D,T} <: AbstractRealSpaceData{D}

A data grid defined in real space, containing data of type T.
"""
struct RealSpaceDataGrid{D,T} <: AbstractRealSpaceData{D}
    # Basis vectors defining the lattice
    latt::BasisVectors{D}
    # Shift of the origin from the lattice
    orig::SVector{D,Float64}
    # The actual data grid
    grid::Array{T,D}
    # Inner constructor
    function RealSpaceDataGrid(
        latt::BasisVectors{D},
        orig::AbstractVector{<:Real},
        grid::Array{T,D}
    ) where {D,T}
        @assert length(orig) == D "Origin vector does not have the lattice dimensionality"
        # Make sure the shift values lie in (-0.5, 0.5]
        orig = orig - floor.(orig)
        return new{D,T}(latt, orig, grid)
    end
end

"""
    RealSpaceDataGrid(
        latt::AbstractLattice{D},
        orig::AbstractVector{<:Real},
        grid::AbstractArray{T,D};
        prim=true
    )

Creates a real space data grid using lattice information from an `AbstractLattice`. By default,
data is assumed to be given in terms of the primitive lattice (as is usually the case for
computational data).
"""
function RealSpaceDataGrid(
    latt::AbstractLattice{D},
    orig::AbstractVector{<:Real},
    grid::AbstractArray{T,D};
    prim=true
) where {D,T}
    # Conversion for safety
    l = RealLattice(latt)
    if prim
        return RealSpaceDataGrid(prim(l), orig, grid)
    else
        return RealSpaceDataGrid(conv(l), orig, grid)
    end
end

"""
    RealSpaceDataGrid(f, g::RealSpaceDataGrid)

Applies a function `f` elementwise to the grid elements of a `RealSpaceDataGrid` and returns a new
`RealSpaceDataGrid`.
"""
function RealSpaceDataGrid(f, g::RealSpaceDataGrid)
    return RealSpaceDataGrid(basis(g), shift(g), f.(grid(g)))
end

# getindex supports arbitrary integer indices for RealSpaceDataGrid
function Base.getindex(g::RealSpaceDataGrid, inds...)
    # Perform modulo math to get the indices
    # WARNING: Julia % is the remainder function, not modulo!
    imod = mod.(inds .- 1,  gridsize(g)) .+ 1
    return getindex(grid(g), imod...)
end

# Iterator definitions: pass through matrix iteration
Base.iterate(g::RealSpaceDataGrid) = iterate(grid(g))
Base.iterate(g::RealSpaceDataGrid, state) = iterate(grid(g), state)

"""
    basis(g::RealSpaceDataGrid{D,T}) -> BasisVectors{D}

Gets the basis vectors of a `RealSpaceDataGrid`.
"""
basis(g::RealSpaceDataGrid) = g.latt
shift(g::RealSpaceDataGrid) = g.orig

"""
    grid(g::RealSpaceDataGrid{D,T}) -> Array{T,D}

Gets the array that backs a `RealSpaceDataGrid{D,T}`, which is an `Array{T,D}`.
"""
grid(g::RealSpaceDataGrid) = g.grid
# Size of the data grid in entries per dimension
# TODO: should we overload Base.size() as well?

"""
    gridsize(g::RealSpaceDataGrid{D,T}) -> NTuple{D,Int}

Gets the dimensions of the backing array corresponding to a `RealSpaceDataGrid`.
"""
gridsize(g::RealSpaceDataGrid) = size(g.grid)
#Base.size(g::RealSpaceDataGrid) = gridsize(g)

"""
    volume(g::RealSpaceDataGrid) -> Float64

Gets the crystal volume associated with a `RealSpaceDataGrid`.

By default, units are assumed to be cubic angstroms.
"""
volume(g::RealSpaceDataGrid) = volume(basis(g))

"""
    voxelsize(g::RealSpaceDataGrid) -> Float64

Gets the size of a single voxel of a `RealSpaceDataGrid`.

By default, units are assumed to be cubic angstroms.
"""
voxelsize(g::RealSpaceDataGrid) = volume(g) / prod(gridsize(g))

"""
    Xtal.grid_check(g1::RealSpaceDataGrid, g2::RealSpaceDataGrid)

Performs a check on two `RealSpaceDataGrid`s to ensure that the basis, origin shift, and grid
dimensions are the same before performing mathematical operations.
"""
function grid_check(g1::RealSpaceDataGrid, g2::RealSpaceDataGrid)
    @assert basis(g1) === basis(g2) "Grid basis vectors for each grid are not identical."
    @assert shift(g1) === shift(g2) "Grid shifts from origin are not identical."
    @assert size(grid(g1)) === size(grid(g2)) "Grid sizes are different."
    return nothing
end

# TODO: ensure that type inference works
# We might be able to remove some type parameter tests as well
function Base.:+(g1::RealSpaceDataGrid{D,T1}, g2::RealSpaceDataGrid{D,T2}) where {D,T1,T2}
    # Check that the grids are identical
    grid_check(g1, g2)
    # Add the two datagrids elementwise
    newgrid = grid(g1) + grid(g2)
    return RealSpaceDataGrid(basis(g1), shift(g1), newgrid)
end

function Base.:*(g1::RealSpaceDataGrid{D,T1}, g2::RealSpaceDataGrid{D,T2}) where {D,T1,T2}
    # Check that the grids are identical
    grid_check(g1, g2)
    # Add the two datagrids elementwise
    newgrid = grid(g1) .* grid(g2)
    return RealSpaceDataGrid(basis(g1), shift(g1), newgrid)
end

function Base.:*(s::Number, g::RealSpaceDataGrid{D,T}) where {D,T}
    newgrid = s * grid(g)
    return RealSpaceDataGrid(basis(g), shift(g), newgrid)
end

Base.:*(g::RealSpaceDataGrid, s::Number) = s * g

# Multiply everything in the datagrid by -1
Base.:-(g::RealSpaceDataGrid) = RealSpaceDataGrid(basis(g), shift(g), -grid(g))
# Subtract two datagrids
Base.:-(g1::RealSpaceDataGrid, g2::RealSpaceDataGrid) = +(g1, -g2)

"""
    integrate(g::RealSpaceDataGrid{D,T<:Number}) -> <:Number

Performs an integration across all voxels, returning a scalar value.
"""
function integrate(g::RealSpaceDataGrid{D,T}) where {D,T<:Number}
    return sum(grid(g)) * voxelsize(g)
end

"""
    integrate(f, g::RealSpaceDataGrid{D,T<:Number}) -> <:Number

Applies the function `f` pointwise to the elements of a datagrid, then integrates the grid.
"""
function integrate(f, g::RealSpaceDataGrid{D,T}) where {D,T<:Number}
    return sum(f.(grid(g))) * voxelsize(g)
end

"""
    fft(g::RealSpaceDataGrid{D,<:Number}; maxhkl=zeros(Int,D)) -> HKLData{D,<:Complex}

Performs a fast Fourier transform on the data in a `RealSpaceDataGrid{D,<:Number}` and
generates an `HKLData{D,<:Complex}`.
"""
function FFTW.fft(
    g::RealSpaceDataGrid{D,<:Number};
    maxhkl::AbstractVector{<:Integer} = zeros(Int, D)
) where D 
    # Generate the bounds needed for the HKLdata
    bounds = [range(div(sz, 2, RoundUp) - sz, length = sz) for sz in gridsize(g)]
    @debug string("Bounds:\t", bounds)
    # Calculate the shift factors to put the values in the right places
    shifts = [last(b) for b in bounds]
    # Take the grid fft
    # Permute the elements to get the indexing working
    f = circshift(fft(grid(g)), shifts .+ 1)
    # With defined bounds, truncate the data
    if all(!iszero, maxhkl)
        ind = [(-x:x) .- first(b) .+ 1 for (x,b) in zip(abs.(maxhkl), bounds)]
        f = f[ind...]
        bounds = [-x:x for x in abs.(maxhkl)]
    end
    return HKLData(f, bounds)
end

function d_spacing(g::RealSpaceDataGrid, miller::AbstractVector{<:Integer})
    return d_spacing(basis(g), miller)
end

#=
"""
    interpolate(g::RealSpaceDataGrid{D,T}, inds...)

Evaluates a datagrid at a fractional index by linearly interpolating the values between nearest
set of points.
"""
function interpolate(g::RealSpaceDataGrid{D,T}, inds::Vararg{<:Real,D}) where {D,T}
    # Get the floor
    f = floor.(inds)
    c = ceil.(inds)
    # If they're the same, it's probably because the indices were integers
    # Just index the array
    if f == c
        return g[inds...]
    end
    # Get the set of points needed to perform the interpolation
end

"""
    interpolate(g::RealSpaceDataGrid{D,T}, r::AbstractVector{<:Real})

Evaluates a datagrid at the reduced coordinate `r` using linear interpolation.
"""
function interpolate(g::RealSpaceDataGrid, r::AbstractVector{<:Real})
    inds = reduced .* size
end
=#

# RECIPROCAL SPACE
#-------------------------------------------------------------------------------------------------#

"""
    KPointGrid{D} <: AbstractKPoints{D}

Contains a grid used to generate k-points during a calculation.

The grid itself is given as an `SMatrix{D,D,Int}`, and can be interpreted as a set of `D` vectors 
given in terms of the primitive basis. These vectors can alternatively be used to construct a 
supercell.

The shift of the k-point mesh off Γ is given as an `SVector{D,Float64}`.
"""
struct KPointGrid{D} <: AbstractKPoints{D}
    grid::SMatrix{D,D,Int}
    orig::SVector{D,Float64}
    function KPointGrid{D}(grid::AbstractMatrix{<:Integer}, orig::AbstractVector{<:Real}) where D
        # only allow positive values in the grid matrix
        @assert all(x -> x > 0, grid) "negative values are disallowed in the grid matrix"
        # Keep shift inside the Brillouin zone
        orig = orig -  round.(orig)
        return new(grid, orig)
    end
end

"""
    KPointList{D} <: AbstractKPoints{D}

Contains a list of k-points. This is useful for describing an ordered list of k-points that are 
not associated with a mesh - for instance, those used in band structure calculations. This can
also be used to store lists of k-points that are generated by a `KPointGrid{D}`.

Also included is a list of weights. When k-points are constructed from a grid, some k-points that
are placed on lattice symmetry operations are counted multiple times, and the weights compensate
for this. If no weights are provided, they are assumed to be 1.
"""
struct KPointList{D} <: AbstractKPoints{D}
    points::Vector{SVector{D,Float64}}
    weights::Vector{Float64}
    function KPointList(
        points::AbstractVector{<:SVector{D,<:Real}},
        weights::AbstractVector{<:Real}
    ) where D
        @assert length(points) == length(weights) "Number of k-points and weights do not match."
        return new{D}(points, weights)
    end
end

function KPointList{D}(
    points::AbstractVector{<:AbstractVector{<:Real}},
    weights::AbstractVector{<:Real}
) where D
    @assert length.(points) == D "k-points have the wrong dimensionality"
    svpoints = [SVector{D,Float64}(v) for v in points]
    return KPointList(svpoints, weights)
end

function KPointList{D}(points::AbstractVector{<:AbstractVector{<:Real}}) where D
    return KPointList{D}(points, ones(length(points)))
end

function KPointList(points::AbstractVector{<:SVector{D,<:Real}}) where D
    return KPointList(points, ones(length(points)))
end

# Index like the vector it is internally
Base.getindex(k::KPointList, i) = k.list[i]

function Base.setindex!(k::KPointList, v::AbstractVector, i)
    k.points[i] = v
end

"""
    nkpt(k::KPointList{D}) -> Int

Gets the number of k-points in a `KPointList`.
"""
nkpt(k::KPointList) = length(k.points)
Base.length(k::KPointList) = length(k.points)

#= TODO: figure out how to get a k-point list from a grid
#  This would require also getting the correct k-point weights
function KPointList{D}(k::KPointGrid)

end
=#

"""
    BandAtKPoint

Stores information about a band's energy and its occupancy at a specific k-point.
"""
struct BandAtKPoint
    # Energies
    e::Vector{Float64}
    # Occupancy
    occ::Vector{Float64}
    function BandAtKPoint(e::AbstractVector{<:Real}, occ::AbstractVector{<:Real})
        @assert length(e) == length(occ) "Size of energy and occupancy arrays do not match."
        return new(e, occ)
    end
end

"""
    BandAtKPoint(eocc::AbstractVector{NTuple{2,<:Real}})

Constructs a new `BandAtKPoint` from a vector containing tuples of energy and occupancy data
(in that order).
"""
function BandAtKPoint(eocc::AbstractVector{<:NTuple{2,<:Real}})
    return BandAtKPoint([x[1] for x in eocc], [x[2] for x in eocc])
end

# Access any pair of energy and occupancy with indexing
function Base.getindex(b::BandAtKPoint, inds...)
    return (b.e[inds...], b.occ[inds...])
end

"""
    nband(b::BandAtKPoint) -> Int

Returns the number of bands associated with a k-point.
"""
nband(b::BandAtKPoint) = length(b.e)

"""
    BandStructure{D}

Stores information about an electronic band structure, including the list of k-points used to
generate the data (as a `KPointList{D}`)and the band information at every k-point (as a 
`Vector{BandAtKPoint}`).
"""
struct BandStructure{D} <: AbstractReciprocalSpaceData{D}
    # k-points for which band data is defined
    kpts::KPointList{D}
    # Set of energy and occupancy data
    bands::Vector{BandAtKPoint}
    function BandStructure{D}(kpts::KPointList{D}, bands::Vector{BandAtKPoint}) where D
        @assert nkpt(kpts) == length(bands) "Incorrect number of k-points or band datasets."
        @assert _allsame(length(bands)) "Number of bands is inconsistent."
        return new(kpts, bands)
    end
end

"""
    BandStructure{D}(kpts::AbstractKPoints{D}, bands::AbstractVector{<:BandAtKPoint}) where D

Generates a new band structure from k-point information and a vector containing band information
at each k-point.
"""
function BandStructure{D}(kpts::AbstractKPoints{D}, bands::AbstractVector{<:BandAtKPoint}) where D
    return BandStructure{D}(kpts, bands)
end

# Get the pair of a k-point and associated band data
function Base.getindex(b::BandStructure{D}, inds...) where D
    return (b.kpts[inds...], b.bands[inds...])
end

nkpt(b::BandStructure{D}) where D = nkpt(b.kpts)
nband(b::BandStructure{D}) where D = nband(b.bands[1])

"""
    FatBands{D} <: AbstractReciprocalSpaceData{D}

Stores information relevant to plotting fatbands.

- FatBands.bands: matrix of energies at each [kpt, band].
- FatBands.projband: array of lm-decomposed band structure. [orbital, ion, band, kpt].
- FatBands.cband: array of complex-valued contributions to band structure.
"""
struct FatBands{D} <: AbstractReciprocalSpaceData{D}
    bands::Matrix{Float64}
    projband::Array{Float64,4}
    cband::Array{Complex{Float64},4}
end

"""
    HKLData{D,T} <: AbstractReciprocalSpaceData{D}

Stores information associated with specific sets of reciprocal lattice vectors. Data can be
accessed and modified using regular indexing, where indices may be negative.
"""
struct HKLData{D,T} <: AbstractHKL{D,T}
    # the actual data
    data::Array{T,D}
    # the bounds in each dimension
    # mutable since the dimensions of Array{D,T} can be changed, in principle
    bounds::MVector{D,UnitRange{Int}}
    function HKLData(
        data::AbstractArray{T,D},
        bounds::AbstractVector{<:AbstractRange{<:Integer}}
    ) where {D,T}
        # The size of the array should match the bounds given
        # For instance, a HKLData with bounds  [-10:10, -10:10, -10:10] should Be
        # [21, 21, 21]
        @assert [s for s in size(data)] == [length(r) for r in bounds] string(
            "Array size incompatible with bounds."
        )
        return new{D,T}(data, bounds)
    end
end

# Needed because HKLData will nearly always have unexpected indices
Base.has_offset_axes(g::HKLData) = true

"""
    Xtal.shiftbounds(g::HKLData{D,T}, inds) -> NTuple{D,<:Integer}

Checks that integer array indices used to access data in an `HKLData` are valid and shifts them to 
access the correct portions of the backing array.
"""
function shiftbounds(hkl::HKLData{D,T}, inds) where {D,T}
    # Check that all the indices are in bounds
    if !mapreduce((i,r) -> i in r, &, inds, bounds(hkl))
        # TODO: does this produce a reasonable error message with correct bounds?
        throw(BoundsError(hkl, inds))
    end
    # Adjust the indices to match the array
    # Subtract the minimum index then add 1
    # So if the range is -10:10, an index of 0 should be 0 - -10 + 1 = 11
    i = inds .- minimum.(bounds(hkl)) .+ 1
    return i
end

"""
    grid(hkl::HKLData{D,T}) -> Array{T,D}

Returns the array that contains the reciprocal space data. Note that this causes information
about the index offset to be lost!
"""
grid(hkl::HKLData) = hkl.data

"""
    bounds(hkl::HKLData{D,T}) -> MVector{D,UnitRange{Int64}}

Returns a range of minimum and maximum indices along each dimension of an `HKLData`
"""
bounds(hkl::HKLData) = hkl.bounds

# HKLData now supports indexing by Miller index
function Base.getindex(g::HKLData{D,T}, inds...) where {D,T}
    i = shiftbounds(g, inds)
    return g.data[i...]
end

function Base.setindex!(g::HKLData{D,T}, x::T, inds...) where {D,T}
    i = shiftbounds(g, inds)
    g.data[i...] = x
end

function HKLData(
    a::AbstractArray{T,D},
    bounds::Vararg{AbstractUnitRange{<:Integer},D}
) where {D,T}
    return HKLData(a, bounds)
end

function Base.zeros(
    ::Type{HKLData{D,T}},
    bounds::Vararg{AbstractUnitRange{<:Integer}, D}
) where {D,T}
    data = zeros(T, length.(bounds))
    return HKLData(data, MVector{D,UnitRange{Int}}(bounds))
end

Base.abs(hkl::HKLData) = HKLData(abs.(grid(hkl)), bounds(hkl))
Base.abs2(hkl::HKLData) = HKLData(abs2.(grid(hkl)), bounds(hkl))

"""
    HKLDict{D,T}

An alternative to `HKLData` uses a dictionary instead of an array as a backing field.

This is a more space-efficient alternative to `HKLData` in the case of reciprocal space data with
a large number of zero components. For wavefunction data, which is often specified to some energy
cutoff that corresponds to a distance in reciprocal space, there are many zero valued elements to
the array. Unspecified elements in an `HKLDict` are assumed to be zero.
"""
struct HKLDict{D,T} <: AbstractHKL{D,T}
    dict::Dict{SVector{D,Int},T}
end

Base.has_offset_axes(hkl::HKLDict) = true

function Base.getindex(hkl::HKLDict{D,T}, inds...) where {D,T}
    v = SVector{D,Int}(inds...)
    if haskey(hkl.dict, v)
        return hkl.dict[v]
    else
        # Return a zero element of some kind by default
        return zero(T)
    end
end

function Base.setindex!(hkl::HKLDict{D,T}, value::T, inds...) where {D,T}
    hkl.dict[SVector{D,Int}(inds...)] = value
end

Base.keys(hkl::HKLDict) = keys(hkl.dict)

Base.iterate(hkl::HKLDict) = iterate(hkl.dict)
Base.iterate(hkl::HKLDict, i) = iterate(hkl.dict, i)

"""
    vectors(hkl::HKLDict)

Returns the set of vectors in an `HKLDict` for which values have been defined.
"""
function vectors(hkl::HKLDict)
    return keys(hkl.dict)
end

function HKLDict(hkl::HKLData{D,T}) where {D,T<:Union{<:Number,<:AbstractArray{Number}}}
    dict = Dict{SVector{D,Int},T}()
    # Get the offset for the indices
    offset = minimum.(hkl.bounds) .- 1
    # Iterate through the matrix and get its coordinates
    for ind in CartesianIndices(hkl.data)
        # Only add nonzero elements
        if hkl.data[ind] != zero(T)
            dict[Tuple(ind) .+ offset] = hkl.data[ind]
        end
    end
    return HKLDict(dict)
end

function HKLData(hkl::HKLDict{D,T}) where {D,T<:Union{<:Number,<:AbstractArray{Number}}}
    # Find the bounds
    bounds = MVector(UnitRange(extrema(v[n] for v in keys(hkl.dict)...)) for n in 1:D)
    data = zeros(T, length.(bounds)...)
    # Loop through the dictionary
    for (k,v) in hkl.dict
        data[k...] = v
    end
    return HKLData(data, bounds)
end

"""
    ReciprocalWavefunction{D,T<:Real} <: AbstractReciprocalSpaceData{D}

Contains a wavefunction stored by k-points and bands in a planewave basis. Used to store data in
VASP WAVECAR files. Each k-point is expected to have the same number of bands.

Every band has associated data containing coefficients of the constituent planewaves stored in a 
`HKLData{D,Complex{T}}`. Unlike most data structures provided by this package, the type of
complex number used does not default to `Float64`: wavefunction data is often supplied as a 
`Complex{Float32}` to reduce the size of the data.
"""
struct ReciprocalWavefunction{D,T<:Real} <: AbstractReciprocalSpaceData{D}
    # Reciprocal lattice on which the k-points are defined
    rlatt::BasisVectors{D}
    # k-points used to construct the wavefunction
    kpts::KPointList{D}
    # Planewave coefficients: a Matrix (size nkpt*maxnband) of HKLData
    waves::Matrix{HKLData{D,Complex{T}}}
    function ReciprocalWavefunction(
        rlatt::BasisVectors{D},
        kpts::AbstractKPoints{D},
        waves::AbstractMatrix{HKLData{D,Complex{T}}}
    ) where {D,T<:Real}
        @assert length(kpts) == size(waves, 1) string(
            "k-point list length inconsistent with number of wavefunction entries"
        )
        return new{D,T}(rlatt, kpts, waves)
    end
end

function ReciprocalWavefunction(
    latt::AbstractLattice{D},
    kpts::AbstractKPoints{D},
    waves::AbstractMatrix{HKLData{D,Complex{T}}};
    primitive::Bool = true
) where {D,T<:Real}
    M = primitive ? prim(latt) : conv(latt)
    return ReciprocalWavefunction(M, kpts, waves)
end

# Getting indices should pull from the waves struct: wf[kpt, band]
Base.getindex(wf::ReciprocalWavefunction{D,T}, inds...) where {D,T} = wf.waves[inds...]

nkpt(wf::ReciprocalWavefunction{D,T}) where {D,T} = size(wf.waves, 1)
nband(wf::ReciprocalWavefunction{D,T}) where {D,T} = size(wf.waves, 2)

# DENSITY OF STATES
#-------------------------------------------------------------------------------------------------#

"""
    DensityOfStates

Contains the total density of states information.
"""
struct DensityOfStates <: AbstractDensityOfStates
    # Fermi energy
    fermi::Float64
    # Energy at each point
    energy::Vector{Float64}
    # Number of states at each energy
    dos::Vector{Float64}
    # Integrated density of states (electron count up to that energy)
    int::Vector{Float64}
    function DensityOfStates(
        fermi::Real,
        energy::AbstractVector{<:Real},
        dos::AbstractVector{<:Real},
        int::AbstractVector{<:Real}
    )
        @assert size(energy) == size(dos) == size(int) string(
            "Number of energy and DOS entries do not match."
        )
        return new(fermi, energy, dos, int)
    end
end

# Create a DOS curve without integrated DOS data
"""
    DensityOfStates(fermi::Real, energy::AbstractVector{<:Real}, dos::AbstractVector{<:Real})

Creates a `DensityOfStates` object without integrated DOS information by integrating the state 
density.
"""
function DensityOfStates(
    fermi::Real,
    energy::AbstractVector{<:Real},
    dos::AbstractVector{<:Real}
)
    # Generate an integrated DOS vector
    # TODO: perhaps we can improve this with a trapezoidal approximation?
    int = [sum(dos[1:n]) for n in 1:length(dos)]
    return DensityOfStates(fermi, energy, dos, int)
end

# Getting an index produces a tuple of the energy, 
Base.getindex(d::DensityOfStates, ind) = (d.energy[ind], d.dos[ind], d.int[ind])

# Setting up for gaussian smearing of DOS curve
function gaussian(
    dos::DensityOfStates,
    sigma::Real,
)
    # Center gaussian within energy range
    # Define gaussian
    # Permute gaussian to prevent shifting of smear output
    ctr = dos.energy[div(length(dos.energy), 2)]
    g = exp.(-((dos.energy .- ctr)/sigma).^2)
    return circshift(g, 1 - findmax(g)[2])
end

"""
    smear(dos::DensityOfStates, sigma::Real)

Smears the DOS function by convoluting it with a gaussian, sigma will determine how
much the DOS is smeared. 
"""
function smear(
    dos::DensityOfStates,
    sigma::Real
)
    # Convolute DOS with gaussian using Fourier transforms
    gsmear = real(ifft(fft(dos.dos) .* fft(gaussian(dos, sigma))))
    # Normalize gaussian smearing so that sum(gaussian smear) = sum(DOS)
    return DensityOfStates(fermi(dos), energies(dos), (gsmear)*(sum(dos.dos)/sum(gsmear)), dos.int)
end

"""
    ProjectedDensityOfStates

Contains projected density of states information.
"""
struct ProjectedDensityOfStates
    # Fermi energy
    fermi::Float64
    # Range of energies
    energy::Vector{Float64}
    # Matrix containing projected DOS data
    # Columns are the components, rows are the energies
    dos::Matrix{Float64}
    # Integrated DOS in the same format
#    int::Matrix{Float64}
    function ProjectedDensityOfStates(
        fermi::Real,
        energy::AbstractVector{<:Real},
        dos::AbstractMatrix{<:Real}
        #int::AbstractMatrix{<:Real}
    )
        @assert length(energy) == size(dos,2) "Number of energy and DOS entries do not match."
        return new(fermi, energy, dos) #, int)
    end
end

#="""
    ProjectedDensityOfStates(
        fermi::Real,
        energy::AbstractVector{<:Real},
        dos::AbstractVector{<:Real}
    )

Creates a `ProjectedDensityOfStates` object without integrated DOS information by integrating
the state density.
"""
function ProjectedDensityOfStates(
    fermi::Real,
    energy::AbstractVector{<:Real},
    dos::AbstractMatrix{<:Real}
)
    # Generate an integrated DOS matrix
    # TODO: perhaps we can improve this with a trapezoidal approximation?
    int = hcat(vec(sum(M[:,1:n], dims=2)) for n in 1:size(M,2))
    return ProjectedDensityOfStates(fermi, energy, dos, int)
end=#

"""
    fermi(d::AbstractDensityOfStates) -> Float64

Gets the Fermi energy from DOS data. There are no guarantees on the unit of energy used!
"""
fermi(d::AbstractDensityOfStates) = d.fermi

"""
    energies(d::AbstractDensityOfStates; usefermi=false) -> Vector{Float64}

Gets the range of energies in the dataset. If `usefermi` is set to true, the energies returned will
be adjusted such that the Fermi energy is set to zero.

There are no guarantees on the unit of energy used!
"""
energies(d::AbstractDensityOfStates; usefermi=false) = d.energy .- (usefermi * d.fermi)

"""
    nelectrons(d::DensityOfStates)

Gets the approximate number of electrons that are needed to reach the Fermi level.
"""
function nelectrons(d::DensityOfStates)
    # Get the nearest entries to the Fermi energy
    E = energies(d, usefermi=true)
    # Find the closest pair of energies to the Fermi energy
    inds = partialsortperm(abs.(E), 1:2)
    # Interpolate linearly between the nearest points
    # TODO: there's probably a nicer way to do this.
    (e1, e2) = [d[i][1] for i in inds]
    (i1, i2) = [d[i][3] for i in inds]
    m = (i2 - i1)/(e2 - e1)
    b = i1 - m*e1
    return m*fermi(d) + b
end

# ATOMIC DATA
#-------------------------------------------------------------------------------------------------#

"""
    AtomicData{D,T}

Data associated with individual atoms in a structure.

This is a type alias for `Dict{AtomPosition{D},T}`. Keys are `AtomPosition` entries, and the values
may be of any type.
"""
const AtomicData{D,T} = Dict{AtomPosition{D},T} where {D,T}

# Note: @computed structs cannot be documented normally
# Use an @doc after the struct, like such
@computed struct SphericalComponents{Lmax}
    v::NTuple{(Lmax+1)^2,Float64}
    # Default constructor without parameters takes numbers directly
    function SphericalComponents(x::Vararg{<:Real,N}) where N
        L = sqrt(length(x)) - 1
        if isinteger(L)
            Lmax = Int(L)
        else
            throw(ArgumentError(string("Cannot determine Lmax from number of arguments.")))
        end
        return new{Lmax}(x)
    end
    # Default constructor with parameters takes any iterator
    function SphericalComponents{Lmax}(v) where Lmax
        @assert length(v) == (Lmax+1)^2 "For Lmax == $Lmax, iterator have length $((Lmax+1)^2)"
        return new{Lmax}(Tuple(v))
    end
end

@doc """
    SphericalComponents{Lmax}

Real spherical harmonic components up to `Lmax`. This can be used to describe atomic orbitals or
projections of data onto atomic sites.
""" SphericalComponents

SphericalComponents(v::SVector{N,<:Real}) where N = SphericalComponents(v...)
SphericalComponents(t::NTuple{N,<:Real}) where N = SphericalComponents(t...)

"""
    Xtal.sc_ind(l::Integer, m::Integer) -> Int

Gets the associated linear index for a pair of (l,m) values used in `SphericalComponents`.
"""
sc_ind(l::Integer, m::Integer) = l^2 + l + 1 + m

# TODO: finish the inverse of the above function
# sc_ind(x) = 

function Base.getindex(s::SphericalComponents{Lmax}, l::Integer, m::Integer) where Lmax
    abs(m) <= l || error("|m| must be less than l")
    l <= Lmax || error("l exceeds lmax ($Lmax)")
    return s.v[sc_ind(l, m)]
end
