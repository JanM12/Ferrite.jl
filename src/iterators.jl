# This file defines iterators used for looping over a grid

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


"""
    FaceCache(grid::Grid)
    FaceCache(dh::AbstractDofHandler)

Create a cache object with pre-allocated memory for the nodes, coordinates, and dofs of a
cell suitable for looping over faces in a grid. This struct stores the current face number,
in addition to an underlying [`CellCache`](@ref). 
The cache is updated for a new face by calling `reinit!(cache, faceindex::FaceIndex)`

**Methods with `fc::FaceCache`**
 - `reinit!(fc, faceindex)`: reinitialize the cache for face `faceindex::FaceIndex`
 - `Ferrite.faceindex(fc)`: get the `FaceIndex` of the currently cached face
 - `Ferrite.faceid(fc)`: get the current faceid (`faceindex(fc)[2]`)
 - `cellid(fc)`: get the current cellid (`faceindex(fc)[1]`)
 - `getnodes(fc)`: get the global node ids of the cell
 - `getcoordinates(fc)`: get the coordinates of the cell
 - `celldofs(fc)`: get the global dof ids of the cell
 - `reinit!(fv, fc)`: reinitialize [`FaceValues`](@ref)
 
See also [`FaceIterator`](@ref).
"""
struct FaceCache{CC<:CellCache}
    cc::CC  # const for julia>1.8 
    current_facenr::ScalarWrapper{Int} 
end
FaceCache(args...) = FaceCache(CellCache(args...), ScalarWrapper(0))

function reinit!(fc::FaceCache, face::FaceIndex)
    cellnr, facenr = face
    reinit!(fc.cc, cellnr)
    fc.current_facenr[] = facenr
    return nothing
end

for op = (:getnodes, :getcoordinates, :cellid, :celldofs)
    eval(quote
        function Ferrite.$op(fc::FaceCache, args...; kwargs...)
            return Ferrite.$op(fc.cc, args...; kwargs...)
        end
    end)
end
@inline faceid(fc::FaceCache) = fc.current_faceid[]
@inline celldofs!(v::Vector, fc::FaceCache) = celldofs!(v, fc.cc)
@inline onboundary(fc::FaceCache) = onboundary(fc.cc, faceid(fc))
@inline faceindex(fc::FaceCache) = FaceIndex(cellid(fc), faceid(fc))
@inline function reinit!(fv::FaceValues, fc::FaceCache)
    reinit!(fv, fc.cc, faceid(fc))
end


###################
## GridIterators ##
###################
# CellIterator
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
@inline _getset(ci::CellIterator) = ci.set 
@inline _getcache(ci::CellIterator) = ci.cc

# FaceIterator
# Leaving flags undocumented as for CellIterator
"""
    FaceIterator(gridordh::Union{Grid,AbstractDofHandler}, faceset::Set{FaceIndex})

Iterate over the faces in `set`. 
Create a `FaceIterator` to conveniently iterate over the faces in `faceset`. 
The elements of the iterator are [`FaceCache`](@ref)s which are properly
`reinit!`ialized. See [`FaceCache`](@ref) for more details.

Looping over a `FaceIterator`, i.e.:
```julia
for fc in FaceIterator(grid, faceset)
    # ...
end
```
is thus simply convenience for the following equivalent snippet:
```julia
fc = FaceCache(grid)
for faceindex in faceset
    reinit!(fc, faceindex)
    # ...
end
"""
struct FaceIterator{FC<:FaceCache, IC<:Set{FaceIndex}}
    fc::FC
    set::F
end

function FaceIterator(gridordh::Union{Grid,AbstractDofHandler}, 
                      set, flags::UpdateFlags=UpdateFlags())
    if gridordh isa MixedDofHandler
        # TODO: Since the CellCache is resizeable this is not really necessary to check
        #       here, but might be useful to catch slow code paths?
        _check_same_celltype(gridordh.grid, set)
    end
    return FaceIterator(FaceCache(gridordh, flags), set)
end

@inline _getcache(fi::FaceIterator) = fi.fc
@inline _getset(fi::FaceIterator) = fi.fc

# Iterator interface
const GridIterators{C} = Union{CellIterator{C},FaceIterator{C}}

function Base.iterate(iterator::GridIterators, state_in...)
    it = iterate(_getset(iterator), state_in...)
    it === nothing && return nothing
    cellid, state_out = it
    cache = _getcache(iterator)
    reinit!(cache, cellid)
    return (cache, state_out)
end
Base.IteratorSize(::Type{<:GridIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:GridIterator}) = Base.HasEltype()
Base.eltype(::Type{<:GridIterator{C}}) where C = C
Base.length(iterator::GridIterator) = length(_getset(iterator))


function _check_same_celltype(grid::AbstractGrid, cellset::IntegerCollection)
    celltype = typeof(grid.cells[first(cellset)])
    if !all(typeof(grid.cells[i]) == celltype for i in cellset)
        error("The cells in the cellset are not all of the same celltype.")
    end
end

function _check_same_celltype(grid::AbstractGrid, faceset::Set{FaceIndex})
    celltype = typeof(grid.cells[first(first(cellset))])
    if !all(typeof(grid.cells[first(face)]) == celltype for face in faceset)
        error("The cells in the cellset are not all of the same celltype.")
    end
end