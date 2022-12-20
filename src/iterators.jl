<<<<<<< HEAD
# this file defines iterators used for looping over a grid
abstract type AbstractGridIterator end
Base.IteratorSize(::Type{T})   where {T<:AbstractGridIterator} = Base.HasLength() # this is default in Base
Base.IteratorEltype(::Type{T}) where {T<:AbstractGridIterator} = Base.HasEltype() # this is default in Base
Base.eltype(::Type{T})         where {T<:AbstractGridIterator} = T

=======
# This file defines iterators used for looping over a grid
>>>>>>> master

struct UpdateFlags
    nodes::Bool
    coords::Bool
    dofs::Bool
end

UpdateFlags(; nodes::Bool=true, coords::Bool=true, dofs::Bool=true) =
    UpdateFlags(nodes, coords, dofs)


###############
## CellCache ##
###############

"""
<<<<<<< HEAD
struct CellIterator{dim,C,T,DH<:Union{AbstractDofHandler,Nothing}} <: AbstractGridIterator
=======
    CellCache(grid::Grid)
    CellCache(dh::AbstractDofHandler)

Create a cache object with pre-allocated memory for the nodes, coordinates, and dofs of a
cell. The cache is updated for a new cell by calling `reinit!(cache, cellid)` where
`cellid::Int` is the cell id.

**Struct fields of `CellCache`**
 - `cc.nodes :: Vector{Int}`: global node ids
 - `cc.coords :: Vector{<:Vec}`: node coordinates
 - `cc.dofs :: Vector{Int}`: global dof ids (empty when constructing the cache from a grid)

**Methods with `CellCache`**
 - `reinit!(cc, i)`: reinitialize the cache for cell `i`
 - `cellid(cc)`: get the cell id of the currently cached cell
 - `getnodes(cc)`: get the global node ids of the cell
 - `getcoordinates(cc)`: get the coordinates of the cell
 - `celldofs(cc)`: get the global dof ids of the cell
 - `reinit!(fev, cc)`: reinitialize [`CellValues`](@ref) or [`FaceValues`](@ref)

See also [`CellIterator`](@ref).
"""
struct CellCache{X,G<:AbstractGrid,DH<:Union{AbstractDofHandler,Nothing}}
>>>>>>> master
    flags::UpdateFlags
    grid::G
    # Pretty useless to store this since you have it already for the reinit! call, but
    # needed for the CellIterator(...) workflow since the user doesn't necessarily control
    # the loop order in the cell subset.
    cellid::ScalarWrapper{Int}
    nodes::Vector{Int}
    coords::Vector{X}
    dh::DH
    dofs::Vector{Int}
end

function CellCache(grid::Grid{dim,C,T}, flags::UpdateFlags=UpdateFlags()) where {dim,C,T}
    N = nnodes_per_cell(grid)
    nodes = zeros(Int, N)
    coords = zeros(Vec{dim,T}, N)
    return CellCache(flags, grid, ScalarWrapper(-1), nodes, coords, nothing, Int[])
end

function CellCache(dh::Union{DofHandler{dim,T},MixedDofHandler{dim,T}}, flags::UpdateFlags=UpdateFlags()) where {dim,T}
    N = nnodes_per_cell(dh.grid)
    nodes = zeros(Int, N)
    coords = zeros(Vec{dim,T}, N)
    n = ndofs_per_cell(dh)
    celldofs = zeros(Int, n)
    return CellCache(flags, dh.grid, ScalarWrapper(-1), nodes, coords, dh, celldofs)
end

# TODO: Can always resize and combine the two reinit! methods maybe?
function reinit!(cc::CellCache, i::Int)
    cc.cellid[] = i
    if cc.flags.nodes
        cellnodes!(cc.nodes, cc.grid, i)
    end
    if cc.flags.coords
        cellcoords!(cc.coords, cc.grid, i)
    end
    if cc.dh !== nothing && cc.flags.dofs
        @assert cc.dh isa DofHandler
        celldofs!(cc.dofs, cc.dh, i)
    end
    return cc
end
function reinit!(cc::CellCache{<:Any,<:AbstractGrid,<:MixedDofHandler}, i::Int)
    @assert cc.dh isa MixedDofHandler
    cc.cellid[] = i
    if cc.flags.nodes
        resize!(cc.nodes, nnodes_per_cell(cc.dh, i))
        cellnodes!(cc.nodes, cc.dh, i)
    end
    if cc.flags.coords
        resize!(cc.coords, nnodes_per_cell(cc.dh, i))
        cellcoords!(cc.coords, cc.dh, i)
    end
    if cc.flags.dofs
        resize!(cc.dofs, ndofs_per_cell(cc.dh, i))
        celldofs!(cc.dofs, cc.dh, i)
    end
    return cc
end

# reinit! FEValues with CellCache
reinit!(cv::CellValues, cc::CellCache) = reinit!(cv, cc.coords)
reinit!(fv::FaceValues, cc::CellCache, f::Int) = reinit!(fv, cc.coords, f)

# Accessor functions (TODO: Deprecate? We are so inconsistent with `getxx` vs `xx`...)
getnodes(cc::CellCache) = cc.nodes
getcoordinates(cc::CellCache) = cc.coords
celldofs(cc::CellCache) = cc.dofs
cellid(cc::CellCache) = cc.cellid[]

# TODO: This can definitely be deprecated
celldofs!(v::Vector, cc::CellCache) = copyto!(v, cc.dofs) # celldofs!(v, cc.dh, cc.cellid[])

# TODO: These should really be replaced with something better...
nfaces(cc::CellCache) = nfaces(eltype(cc.grid.cells))
onboundary(cc::CellCache, face::Int) = cc.grid.boundary_matrix[face, cc.cellid[]]

##################
## CellIterator ##
##################

const IntegerCollection = Union{Set{<:Integer}, AbstractVector{<:Integer}}

"""
    CellIterator(grid::Grid, cellset=1:getncells(grid))
    CellIterator(dh::AbstractDofHandler, cellset=1:getncells(dh))

Create a `CellIterator` to conveniently iterate over all, or a subset, of the cells in a
grid. The elements of the iterator are [`CellCache`](@ref)s which are properly
`reinit!`ialized. See [`CellCache`](@ref) for more details.

Looping over a `CellIterator`, i.e.:
```julia
for cc in CellIterator(grid, cellset)
    # ...
end
```
is thus simply convenience for the following equivalent snippet:
```julia
cc = CellCache(grid)
for idx in cellset
    reinit!(cc, idx)
    # ...
end
```
"""
struct CellIterator{CC<:CellCache, IC<:IntegerCollection}
    cc::CC
    set::IC
end

function CellIterator(gridordh::Union{Grid,AbstractDofHandler},
                      set::Union{IntegerCollection,Nothing}=nothing,
                      flags::UpdateFlags=UpdateFlags())
    if set === nothing
        grid = gridordh isa AbstractDofHandler ? gridordh.grid : gridordh
        set = 1:getncells(grid)
    end
    if gridordh isa MixedDofHandler
        # TODO: Since the CellCache is resizeable this is not really necessary to check
        #       here, but might be useful to catch slow code paths?
        _check_same_celltype(gridordh.grid, set)
    end
    return CellIterator(CellCache(gridordh, flags), set)
end
function CellIterator(gridordh::Union{Grid,AbstractDofHandler}, flags::UpdateFlags)
    return CellIterator(gridordh, nothing, flags)
end

# Iterator interface
function Base.iterate(ci::CellIterator, state_in...)
    it = iterate(ci.set, state_in...)
    it === nothing && return nothing
    cellid, state_out = it
    reinit!(ci.cc, cellid)
    return (ci.cc, state_out)
end
Base.IteratorSize(::Type{<:CellIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:CellIterator}) = Base.HasEltype()
Base.eltype(::Type{<:CellIterator{CC}}) where CC = CC
Base.length(ci::CellIterator) = length(ci.set)


function _check_same_celltype(grid::AbstractGrid, cellset)
    celltype = typeof(grid.cells[first(cellset)])
    if !all(typeof(grid.cells[i]) == celltype for i in cellset)
        error("The cells in the cellset are not all of the same celltype.")
    end
end

struct BoundaryFaceIterator{CI<:CellIterator} <: AbstractGridIterator
    faces::Vector{Int}
    current_faceid::ScalarWrapper{Int}
    ci::CI
end

function BoundaryFaceIterator(dh::AbstractDofHandler, faces::Vector{Int}, cells::Vector{Int}, args...)
    if length(faces)!=length(cells)
        msg = "faces and cells have different lengths: $(length(faces)) vs $(length(cells))"
        throw(DimensionMismatch(msg))
    end
    return BoundaryFaceIterator(faces, ScalarWrapper(0), CellIterator(dh, cells, args...))
end

function BoundaryFaceIterator(dh::AbstractDofHandler, faceset, cellset=nothing, args...)
    cells, faces = _get_cells_and_faces(faceset, cellset)
    return BoundaryFaceIterator(faces, ScalarWrapper(0), CellIterator(dh, cells, args...))
end

function _get_cells_and_faces(faceset, ::Nothing)
    tuple((collect([faceindex[j] for faceindex in faceset]) for j in 1:2)...)
end

function _get_cells_and_faces(faceset, cellset)
    tuple((collect([faceindex[j] for faceindex in faceset if faceindex[1] in cellset]) for j in 1:2)...)
end

@inline Base.length(fi::BoundaryFaceIterator)  = length(fi.ci)
function Base.iterate(fi::BoundaryFaceIterator, state = 1)
    if state > length(fi)
        return nothing
    else
        return (reinit!(fi, state), state+1)
    end
end

# Use functions from CellIterator, except nfaces
# New functions: faceid, faceindex
function reinit!(fi::BoundaryFaceIterator, i::Int)
    reinit!(fi.ci, i)
    fi.current_faceid[] = fi.faces[i]
    return fi
end

for op = (:getnodes, :getcoordinates, :cellid, :celldofs)
    eval(quote
        function Ferrite.$op(fi::BoundaryFaceIterator, args...; kwargs...)
            return Ferrite.$op(fi.ci, args...; kwargs...)
        end
    end)
end
@inline faceid(fi::BoundaryFaceIterator) = fi.current_faceid[]
@inline celldofs!(v::Vector, fi::BoundaryFaceIterator) = celldofs!(v, fi.ci)
@inline onboundary(fi::BoundaryFaceIterator) = onboundary(fi.ci, faceid(fi))
@inline faceindex(fi::BoundaryFaceIterator) = FaceIndex(cellid(fi), faceid(fi))
@inline function reinit!(fv::FaceValues, fi::BoundaryFaceIterator)
    reinit!(fv, fi.ci, faceid(fi))
end