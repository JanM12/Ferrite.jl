"""
    Interpolation{ref_shape, order}()

Abstract type for interpolations defined on `ref_shape`
(see [`AbstractRefShape`](@ref)).
`order` corresponds to the order of the interpolation.
The interpolation is used to define shape functions to interpolate
a function between nodes.

The following interpolations are implemented:

* `Lagrange{RefLine,1}`
* `Lagrange{RefLine,2}`
* `Lagrange{RefQuadrilateral,1}`
* `Lagrange{RefQuadrilateral,2}`
* `Lagrange{RefTriangle,1}`
* `Lagrange{RefTriangle,2}`
* `Lagrange{RefTriangle,3}`
* `Lagrange{RefTriangle,4}`
* `Lagrange{RefTriangle,5}`
* `BubbleEnrichedLagrange{RefTriangle,1}`
* `CrouzeixRaviart{RefTriangle, 1}`
* `Lagrange{RefHexahedron,1}`
* `Lagrange{RefHexahedron,2}`
* `Lagrange{RefTetrahedron,1}`
* `Lagrange{RefTetrahedron,2}`
* `Lagrange{RefPrism,1}`
* `Lagrange{RefPrism,2}`
* `Serendipity{RefQuadrilateral,2}`
* `Serendipity{RefHexahedron,2}`

# Examples
```jldoctest
julia> ip = Lagrange{RefTriangle, 2}()
Lagrange{RefTriangle, 2}()

julia> getnbasefunctions(ip)
6
```
"""
abstract type Interpolation{shape #=<: AbstractRefShape=#, order, unused} end

const InterpolationByDim{dim} = Interpolation{<:AbstractRefShape{dim}}

abstract type ScalarInterpolation{      refshape, order} <: Interpolation{refshape, order, Nothing} end
abstract type VectorInterpolation{vdim, refshape, order} <: Interpolation{refshape, order, Nothing} end

# Number of components for the interpolation.
n_components(::ScalarInterpolation)                    = 1
n_components(::VectorInterpolation{vdim}) where {vdim} = vdim
# Number of components that are allowed to prescribe in e.g. Dirichlet BC
n_dbc_components(ip::Interpolation) = n_components(ip)
# n_dbc_components(::Union{RaviartThomas,Nedelec}) = 1

# TODO: Remove: this is a hotfix to apply constraints to embedded elements.
edges(ip::InterpolationByDim{2}) = faces(ip)
edgedof_indices(ip::InterpolationByDim{2}) = facedof_indices(ip)
edgedof_interior_indices(ip::InterpolationByDim{2}) = facedof_interior_indices(ip)
facedof_indices(ip::InterpolationByDim{1}) = vertexdof_indices(ip)

# TODO: Add a fallback that errors if there are multiple dofs per edge/face instead to force
#       interpolations to opt-out instead of silently do nothing.
"""
    adjust_dofs_during_distribution(::Interpolation)

This function must return `true` if the dofs should be adjusted (i.e. permuted) during dof
distribution. This is in contrast to i) adjusting the dofs during [`reinit!`](@ref) in the
assembly loop, or ii) not adjusting at all (which is not needed for low order
interpolations, generally).
"""
adjust_dofs_during_distribution(::Interpolation) = false

"""
    InterpolationInfo

Gathers all the information needed to distribute dofs for a given interpolation. Note that
this cache is of the same type no matter the interpolation: the purpose is to make
dof-distribution type-stable.
"""
struct InterpolationInfo
    nvertexdofs::Vector{Int}
    nedgedofs::Vector{Int}
    nfacedofs::Vector{Int}
    ncelldofs::Int
    reference_dim::Int
    adjust_during_distribution::Bool
    n_copies::Int
    function InterpolationInfo(interpolation::InterpolationByDim{3})
        n_copies = 1
        if interpolation isa VectorizedInterpolation
            n_copies = get_n_copies(interpolation)
            interpolation = interpolation.ip
        end
        new(
            [length(i) for i ∈ vertexdof_indices(interpolation)],
            [length(i) for i ∈ edgedof_interior_indices(interpolation)],
            [length(i) for i ∈ facedof_interior_indices(interpolation)],
            length(celldof_interior_indices(interpolation)),
            3,
            adjust_dofs_during_distribution(interpolation),
            n_copies,
        )
    end
    function InterpolationInfo(interpolation::InterpolationByDim{2})
        n_copies = 1
        if interpolation isa VectorizedInterpolation
            n_copies = get_n_copies(interpolation)
            interpolation = interpolation.ip
        end
        new(
            [length(i) for i ∈ vertexdof_indices(interpolation)],
            Int[],
            [length(i) for i ∈ facedof_interior_indices(interpolation)],
            length(celldof_interior_indices(interpolation)),
            2,
            adjust_dofs_during_distribution(interpolation),
            n_copies,
        )
    end
    function InterpolationInfo(interpolation::InterpolationByDim{1})
        n_copies = 1
        if interpolation isa VectorizedInterpolation
            n_copies = get_n_copies(interpolation)
            interpolation = interpolation.ip
        end
        new(
            [length(i) for i ∈ vertexdof_indices(interpolation)],
            Int[],
            Int[],
            length(celldof_interior_indices(interpolation)),
            1,
            adjust_dofs_during_distribution(interpolation),
            n_copies
        )
    end
end

# Some redundant information about the geometry of the reference cells.
nfaces(::Interpolation{RefHypercube{dim}}) where {dim} = 2*dim
nfaces(::Interpolation{RefTriangle}) = 3
nfaces(::Interpolation{RefTetrahedron}) = 4
nfaces(::Interpolation{RefPrism}) = 5

nedges(::Interpolation{RefLine}) = 0
nedges(::Interpolation{RefQuadrilateral}) = 0
nedges(::Interpolation{RefHexahedron}) = 12
nedges(::Interpolation{RefTriangle}) = 0
nedges(::Interpolation{RefTetrahedron}) = 6
nedges(::Interpolation{RefPrism}) = 9

nvertices(::Interpolation{RefHypercube{dim}}) where {dim} = 2^dim
nvertices(::Interpolation{RefTriangle}) = 3
nvertices(::Interpolation{RefTetrahedron}) = 4
nvertices(::Interpolation{RefPrism}) = 6

Base.copy(ip::Interpolation) = ip

"""
    Ferrite.getdim(::Interpolation)

Return the dimension of the reference element for a given interpolation.
"""
@inline getdim(::Interpolation{shape}) where {dim, shape <: AbstractRefShape{dim}} = dim

"""
    Ferrite.getrefshape(::Interpolation)::AbstractRefShape

Return the reference element shape of the interpolation.
"""
@inline getrefshape(::Interpolation{shape}) where {shape} = shape

"""
    Ferrite.getorder(::Interpolation)

Return order of the interpolation.
"""
@inline getorder(::Interpolation{shape,order}) where {shape,order} = order


#####################
# Utility functions #
#####################

"""
    Ferrite.getnbasefunctions(ip::Interpolation)

Return the number of base functions for the interpolation `ip`.
"""
getnbasefunctions(::Interpolation)

# The following functions are used to distribute the dofs. Definitions:
#   vertexdof: dof on a "corner" of the reference shape
#   facedof: dof in the dim-1 dimension (line in 2D, surface in 3D)
#   edgedof: dof on a line between 2 vertices (i.e. "corners") (3D only)
#   celldof: dof that is local to the element

"""
    shape_value(ip::Interpolation, ξ::Vec, i::Int)

Evaluate the value of the `i`th shape function of the interpolation `ip`
at a point `ξ` on the reference element. The index `i` must
match the index in [`vertices(::Interpolation)`](@ref), [`faces(::Interpolation)`](@ref) and
[`edges(::Interpolation)`](@ref).

For nodal interpolations the indices also must match the
indices of [`reference_coordinates(::Interpolation)`](@ref).
"""
shape_value(ip::Interpolation, ξ::Vec, i::Int)

"""
    shape_gradient(ip::Interpolation, ξ::Vec, i::Int)

Evaluate the gradient of the `i`th shape function of the interpolation `ip` in
reference coordinate `ξ`.
"""
function shape_gradient(ip::Interpolation, ξ::Vec, i::Int)
    return Tensors.gradient(x -> shape_value(ip, x, i), ξ)
end

function shape_gradient_and_value(ip::Interpolation, ξ::Vec, i::Int)
    return gradient(x -> shape_value(ip, x, i), ξ, :all)
end


"""
    reference_coordinates(ip::Interpolation)

Returns a vector of coordinates with length [`getnbasefunctions(::Interpolation)`](@ref)
and indices corresponding to the indices of a dof in [`vertices`](@ref), [`faces`](@ref) and
[`edges`](@ref).

    Only required for nodal interpolations.
    
    TODO: Separate nodal and non-nodal interpolations.
"""
reference_coordinates(::Interpolation)

"""
    vertexdof_indices(ip::Interpolation)

A tuple containing tuples of local dof indices for the respective vertex in local
enumeration on a cell defined by [`vertices(::Cell)`](@ref). The vertex enumeration must
match the vertex enumeration of the corresponding geometrical cell.

!!! note
    The dofs appearing in the tuple must be continuous and increasing! The first dof must be
    the 1, as vertex dofs are enumerated first.
"""
vertexdof_indices(ip::Interpolation) = ntuple(_ -> (), nvertices(ip))

"""
    edgedof_indices(ip::Interpolation)

A tuple containing tuples of local dof indices for the respective edge in local enumeration
on a cell defined by [`edges(::Cell)`](@ref). The edge enumeration must match the edge
enumeration of the corresponding geometrical cell.

The dofs are guaranteed to be aligned with the local ordering of the entities on the oriented edge.
Here the first entries are the vertex dofs, followed by the edge interior dofs.
"""
edgedof_indices(::Interpolation)

"""
    edgedof_interior_indices(ip::Interpolation)

A tuple containing tuples of the local dof indices on the interior of the respective edge in
local enumeration on a cell defined by [`edges(::Cell)`](@ref). The edge enumeration must
match the edge enumeration of the corresponding geometrical cell. Note that the vertex dofs
are included here.

!!! note
    The dofs appearing in the tuple must be continuous and increasing! The first dof must be
    computed via "last vertex dof index + 1", if edge dofs exist.
"""
edgedof_interior_indices(::Interpolation)

"""
    facedof_indices(ip::Interpolation)

A tuple containing tuples of all local dof indices for the respective face in local
enumeration on a cell defined by [`faces(::Cell)`](@ref). The face enumeration must match
the face enumeration of the corresponding geometrical cell.
"""
facedof_indices(::Interpolation)

"""
    facedof_interior_indices(ip::Interpolation)

A tuple containing tuples of the local dof indices on the interior of the respective face in
local enumeration on a cell defined by [`faces(::Cell)`](@ref). The face enumeration must
match the face enumeration of the corresponding geometrical cell. Note that the vertex and
edge dofs are included here.

!!! note
    The dofs appearing in the tuple must be continuous and increasing! The first dof must be
    the computed via "last edge interior dof index + 1", if face dofs exist.
"""
facedof_interior_indices(::Interpolation) 

"""
    celldof_interior_indices(ip::Interpolation)

Tuple containing the dof indices associated with the interior of the cell.

!!! note
    The dofs appearing in the tuple must be continuous and increasing! Celldofs are
    enumerated last.
"""
celldof_interior_indices(::Interpolation) = ()

# Some helpers to skip boilerplate
edgedof_indices(ip::InterpolationByDim{3}) = ntuple(_ -> (), nedges(ip))
edgedof_interior_indices(ip::InterpolationByDim{3}) = ntuple(_ -> (), nedges(ip))
facedof_indices(ip::Union{InterpolationByDim{2}, InterpolationByDim{3}}) =  ntuple(_ -> (), nfaces(ip))
facedof_interior_indices(ip::Union{InterpolationByDim{2}, InterpolationByDim{3}}) =  ntuple(_ -> (), nfaces(ip))

"""
    boundarydof_indices(::Type{<:BoundaryIndex})

Helper function to generically dispatch on the correct dof sets of a boundary entity.
"""
boundarydof_indices(::Type{<:BoundaryIndex})

boundarydof_indices(::Type{FaceIndex}) = Ferrite.facedof_indices
boundarydof_indices(::Type{EdgeIndex}) = Ferrite.edgedof_indices
boundarydof_indices(::Type{VertexIndex}) = Ferrite.vertexdof_indices

#########################
# DiscontinuousLagrange #
#########################
# TODO generalize to arbitrary basis positionings.
"""
Piecewise discontinuous Lagrange basis via Gauss-Lobatto points.
"""
struct DiscontinuousLagrange{shape, order, unused} <: ScalarInterpolation{shape, order}
    function DiscontinuousLagrange{shape, order}() where {shape <: AbstractRefShape, order}
        new{shape, order, Nothing}()
    end
end

getlowerorder(::DiscontinuousLagrange{shape,order}) where {shape,order} = DiscontinuousLagrange{shape,order-1}()

getnbasefunctions(::DiscontinuousLagrange{shape,order}) where {shape,order} = getnbasefunctions(Lagrange{shape,order}())
getnbasefunctions(::DiscontinuousLagrange{shape,0}) where {shape} = 1

# This just moves all dofs into the interior of the element.
celldof_interior_indices(ip::DiscontinuousLagrange) = ntuple(i->i, getnbasefunctions(ip))

# Mirror the Lagrange element for now.
function reference_coordinates(ip::DiscontinuousLagrange{shape, order}) where {shape, order}
    return reference_coordinates(Lagrange{shape,order}())
end
function shape_value(::DiscontinuousLagrange{shape, order}, ξ::Vec{dim}, i::Int) where {dim, shape <: AbstractRefShape{dim}, order}
    return shape_value(Lagrange{shape, order}(), ξ, i)
end

# Excepting the L0 element.
function reference_coordinates(ip::DiscontinuousLagrange{RefHypercube{dim},0}) where dim
    return [Vec{dim, Float64}(ntuple(x->0.0, dim))]
end

function reference_coordinates(ip::DiscontinuousLagrange{RefTriangle,0})
    return [Vec{2,Float64}((1/3,1/3))]
end

function reference_coordinates(ip::DiscontinuousLagrange{RefTetrahedron,0})
   return [Vec{3,Float64}((1/4,1/4,1/4))]
end

function shape_value(ip::DiscontinuousLagrange{shape, 0}, ::Vec{dim, T}, i::Int) where {dim, shape <: AbstractRefShape{dim}, T}
    i > 1 && throw(ArgumentError("no shape function $i for interpolation $ip"))
    return one(T)
end

############
# Lagrange #
############
struct Lagrange{shape, order, unused} <: ScalarInterpolation{shape, order}
    function Lagrange{shape, order}() where {shape <: AbstractRefShape, order}
        new{shape, order, Nothing}()
    end
end

# Vertices for all Lagrange interpolations are the same
vertexdof_indices(::Lagrange{RefLine}) = ((1,),(2,))
vertexdof_indices(::Lagrange{RefQuadrilateral}) = ((1,),(2,),(3,),(4,))
vertexdof_indices(::Lagrange{RefHexahedron}) = ((1,),(2,),(3,),(4,),(5,),(6,),(7,),(8,))
vertexdof_indices(::Lagrange{RefTriangle}) = ((1,),(2,),(3,))
vertexdof_indices(::Lagrange{RefTetrahedron}) = ((1,),(2,),(3,),(4,))
vertexdof_indices(::Lagrange{RefPrism}) = ((1,), (2,), (3,), (4,), (5,), (6,))

getlowerorder(::Lagrange{shape,order}) where {shape,order} = Lagrange{shape,order-1}()
getlowerorder(::Lagrange{shape,1}) where {shape} = DiscontinuousLagrange{shape,0}()

############################
# Lagrange RefLine order 1 #
############################
getnbasefunctions(::Lagrange{RefLine,1}) = 2

function reference_coordinates(::Lagrange{RefLine,1})
    return [Vec{1, Float64}((-1.0,)),
            Vec{1, Float64}(( 1.0,))]
end

function shape_value(ip::Lagrange{RefLine, 1}, ξ::Vec{1}, i::Int)
    ξ_x = ξ[1]
    i == 1 && return (1 - ξ_x) * 0.5
    i == 2 && return (1 + ξ_x) * 0.5
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

############################
# Lagrange RefLine order 2 #
############################
getnbasefunctions(::Lagrange{RefLine,2}) = 3

facedof_indices(::Lagrange{RefLine,2}) = ((1,), (2,))
celldof_interior_indices(::Lagrange{RefLine,2}) = (3,)

function reference_coordinates(::Lagrange{RefLine,2})
    return [Vec{1, Float64}((-1.0,)),
            Vec{1, Float64}(( 1.0,)),
            Vec{1, Float64}(( 0.0,))]
end

function shape_value(ip::Lagrange{RefLine, 2}, ξ::Vec{1}, i::Int)
    ξ_x = ξ[1]
    i == 1 && return ξ_x * (ξ_x - 1) * 0.5
    i == 2 && return ξ_x * (ξ_x + 1) * 0.5
    i == 3 && return 1 - ξ_x^2
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

#####################################
# Lagrange RefQuadrilateral order 1 #
#####################################
getnbasefunctions(::Lagrange{RefQuadrilateral,1}) = 4

facedof_indices(::Lagrange{RefQuadrilateral,1}) = ((1,2), (2,3), (3,4), (4,1))

function reference_coordinates(::Lagrange{RefQuadrilateral,1})
    return [Vec{2, Float64}((-1.0, -1.0)),
            Vec{2, Float64}(( 1.0, -1.0)),
            Vec{2, Float64}(( 1.0,  1.0,)),
            Vec{2, Float64}((-1.0,  1.0,))]
end

function shape_value(ip::Lagrange{RefQuadrilateral, 1}, ξ::Vec{2}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 1 && return (1 - ξ_x) * (1 - ξ_y) * 0.25
    i == 2 && return (1 + ξ_x) * (1 - ξ_y) * 0.25
    i == 3 && return (1 + ξ_x) * (1 + ξ_y) * 0.25
    i == 4 && return (1 - ξ_x) * (1 + ξ_y) * 0.25
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

#####################################
# Lagrange RefQuadrilateral order 2 #
#####################################
getnbasefunctions(::Lagrange{RefQuadrilateral,2}) = 9

facedof_indices(::Lagrange{RefQuadrilateral,2}) = ((1,2, 5), (2,3, 6), (3,4, 7), (4,1, 8))
facedof_interior_indices(::Lagrange{RefQuadrilateral,2}) = ((5,), (6,), (7,), (8,))
celldof_interior_indices(::Lagrange{RefQuadrilateral,2}) = (9,)

function reference_coordinates(::Lagrange{RefQuadrilateral,2})
    return [Vec{2, Float64}((-1.0, -1.0)),
            Vec{2, Float64}(( 1.0, -1.0)),
            Vec{2, Float64}(( 1.0,  1.0)),
            Vec{2, Float64}((-1.0,  1.0)),
            Vec{2, Float64}(( 0.0, -1.0)),
            Vec{2, Float64}(( 1.0,  0.0)),
            Vec{2, Float64}(( 0.0,  1.0)),
            Vec{2, Float64}((-1.0,  0.0)),
            Vec{2, Float64}(( 0.0,  0.0))]
end

function shape_value(ip::Lagrange{RefQuadrilateral, 2}, ξ::Vec{2}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 1 && return (ξ_x^2 - ξ_x) * (ξ_y^2 - ξ_y) * 0.25
    i == 2 && return (ξ_x^2 + ξ_x) * (ξ_y^2 - ξ_y) * 0.25
    i == 3 && return (ξ_x^2 + ξ_x) * (ξ_y^2 + ξ_y) * 0.25
    i == 4 && return (ξ_x^2 - ξ_x) * (ξ_y^2 + ξ_y) * 0.25
    i == 5 && return (1 - ξ_x^2) * (ξ_y^2 - ξ_y) * 0.5
    i == 6 && return (ξ_x^2 + ξ_x) * (1 - ξ_y^2) * 0.5
    i == 7 && return (1 - ξ_x^2) * (ξ_y^2 + ξ_y) * 0.5
    i == 8 && return (ξ_x^2 - ξ_x) * (1 - ξ_y^2) * 0.5
    i == 9 && return (1 - ξ_x^2) * (1 - ξ_y^2)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

################################
# Lagrange RefTriangle order 1 #
################################
getnbasefunctions(::Lagrange{RefTriangle,1}) = 3

facedof_indices(::Lagrange{RefTriangle,1}) = ((1,2), (2,3), (3,1))

function reference_coordinates(::Lagrange{RefTriangle,1})
    return [Vec{2, Float64}((1.0, 0.0)),
            Vec{2, Float64}((0.0, 1.0)),
            Vec{2, Float64}((0.0, 0.0))]
end

function shape_value(ip::Lagrange{RefTriangle, 1}, ξ::Vec{2}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 1 && return ξ_x
    i == 2 && return ξ_y
    i == 3 && return 1. - ξ_x - ξ_y
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

################################
# Lagrange RefTriangle order 2 #
################################
getnbasefunctions(::Lagrange{RefTriangle,2}) = 6

facedof_indices(::Lagrange{RefTriangle,2}) = ((1,2,4), (2,3,5), (3,1,6))
facedof_interior_indices(::Lagrange{RefTriangle,2}) = ((4,), (5,), (6,))

function reference_coordinates(::Lagrange{RefTriangle,2})
    return [Vec{2, Float64}((1.0, 0.0)),
            Vec{2, Float64}((0.0, 1.0)),
            Vec{2, Float64}((0.0, 0.0)),
            Vec{2, Float64}((0.5, 0.5)),
            Vec{2, Float64}((0.0, 0.5)),
            Vec{2, Float64}((0.5, 0.0))]
end

function shape_value(ip::Lagrange{RefTriangle, 2}, ξ::Vec{2}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    γ = 1. - ξ_x - ξ_y
    i == 1 && return ξ_x * (2ξ_x - 1)
    i == 2 && return ξ_y * (2ξ_y - 1)
    i == 3 && return γ * (2γ - 1)
    i == 4 && return 4ξ_x * ξ_y
    i == 5 && return 4ξ_y * γ
    i == 6 && return 4ξ_x * γ
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

######################################
# Lagrange RefTriangle order 3, 4, 5 #
######################################
# see https://getfem.readthedocs.io/en/latest/userdoc/appendixA.html

const Lagrange2Tri345 = Union{
    Lagrange{RefTriangle,3},
    Lagrange{RefTriangle,4},
    Lagrange{RefTriangle,5},
}

adjust_dofs_during_distribution(::Lagrange2Tri345) = true

function getnbasefunctions(ip::Lagrange2Tri345)
    order = getorder(ip)
    return (order + 1) * (order + 2) ÷ 2
end

# Permutation to switch numbering to Ferrite ordering
const permdof2DLagrange2Tri345 = Dict{Int,Vector{Int}}(
    1 => [1, 2, 3],
    2 => [3, 6, 1, 5, 4, 2],
    3 => [4, 10, 1, 7, 9, 8, 5, 2, 3, 6],
    4 => [5, 15, 1, 9, 12, 14, 13, 10, 6, 2, 3, 4, 7, 8, 11],
    5 => [6, 21, 1, 11, 15, 18, 20, 19, 16, 12, 7, 2, 3, 4, 5, 8, 9, 10, 13, 14, 17],
)

function facedof_indices(ip::Lagrange2Tri345)
    order = getorder(ip)
    order == 1 && return ((1,2), (2,3), (3,1))
    order == 2 && return ((1,2,4), (2,3,5), (3,1,6))
    order == 3 && return ((1,2,4,5), (2,3,6,7), (3,1,8,9))
    order == 4 && return ((1,2,4,5,6), (2,3,7,8,9), (3,1,10,11,12))
    order == 5 && return ((1,2,4,5,6,7), (2,3,8,9,10,11), (3,1,12,13,14,15))

    throw(ArgumentError("Unsupported order $order for Lagrange on triangles."))
end

function facedof_interior_indices(ip::Lagrange2Tri345)
    order = getorder(ip)
    order == 1 && return ((), (), ())
    order == 2 && return ((4,), (5,), (6,))
    order == 3 && return ((4,5), (6,7), (8,9))
    order == 4 && return ((4,5,6), (7,8,9), (10,11,12))
    order == 5 && return ((4,5,6,7), (8,9,10,11), (12,13,14,15))
    throw(ArgumentError("Unsupported order $order for Lagrange on triangles."))
end

function celldof_interior_indices(ip::Lagrange2Tri345)
    order = getorder(ip)
    ncellintdofs = (order + 1) * (order + 2) ÷ 2 - 3 * order
    totaldofs = getnbasefunctions(ip)
    return ntuple(i->totaldofs-ncellintdofs+i, ncellintdofs)
end

function reference_coordinates(ip::Lagrange2Tri345)
    order = getorder(ip)
    coordpts = Vector{Vec{2, Float64}}()
    for k = 0:order
        for l = 0:(order - k)
            push!(coordpts, Vec{2, Float64}((l / order, k / order)))
        end
    end
    return permute!(coordpts, permdof2DLagrange2Tri345[order])
end

function shape_value(ip::Lagrange2Tri345, ξ::Vec{2}, i::Int)
    if !(0 < i <= getnbasefunctions(ip))
        throw(ArgumentError("no shape function $i for interpolation $ip"))
    end
    order = getorder(ip)
    i = permdof2DLagrange2Tri345[order][i]
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i1, i2, i3 = _numlin_basis2D(i, order)
    val = one(ξ_y)
    i1 ≥ 1 && (val *= prod((order - order * (ξ_x + ξ_y ) - j) / (j + 1) for j = 0:(i1 - 1)))
    i2 ≥ 1 && (val *= prod((order * ξ_x - j) / (j + 1) for j = 0:(i2 - 1)))
    i3 ≥ 1 && (val *= prod((order * ξ_y - j) / (j + 1) for j = 0:(i3 - 1)))
    return val
end

function _numlin_basis2D(i, order)
    c, j1, j2, j3 = 0, 0, 0, 0
    for k = 0:order
        if i <= c + (order + 1 - k)
            j2 = i - c - 1
            break
        else
            j3 += 1
            c += order + 1 - k
        end
    end
    j1 = order - j2 -j3
    return j1, j2, j3
end

###################################
# Lagrange RefTetrahedron order 1 #
###################################
getnbasefunctions(::Lagrange{RefTetrahedron,1}) = 4

facedof_indices(::Lagrange{RefTetrahedron,1}) = ((1,3,2), (1,2,4), (2,3,4), (1,4,3))
edgedof_indices(::Lagrange{RefTetrahedron,1}) = ((1,2), (2,3), (3,1), (1,4), (2,4), (3,4))

function reference_coordinates(::Lagrange{RefTetrahedron,1})
    return [Vec{3, Float64}((0.0, 0.0, 0.0)),
            Vec{3, Float64}((1.0, 0.0, 0.0)),
            Vec{3, Float64}((0.0, 1.0, 0.0)),
            Vec{3, Float64}((0.0, 0.0, 1.0))]
end

function shape_value(ip::Lagrange{RefTetrahedron, 1}, ξ::Vec{3}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    ξ_z = ξ[3]
    i == 1 && return 1.0 - ξ_x - ξ_y - ξ_z
    i == 2 && return ξ_x
    i == 3 && return ξ_y
    i == 4 && return ξ_z
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

###################################
# Lagrange RefTetrahedron order 2 #
###################################
getnbasefunctions(::Lagrange{RefTetrahedron,2}) = 10

facedof_indices(::Lagrange{RefTetrahedron,2}) = ((1,3,2,7,6,5), (1,2,4,5,9,8), (2,3,4,6,10,9), (1,4,3,8,10,7))
edgedof_indices(::Lagrange{RefTetrahedron,2}) = ((1,2,5), (2,3,6), (3,1,7), (1,4,8), (2,4,9), (3,4,10))
edgedof_interior_indices(::Lagrange{RefTetrahedron,2}) = ((5,), (6,), (7,), (8,), (9,), (10,))

function reference_coordinates(::Lagrange{RefTetrahedron,2})
    return [Vec{3, Float64}((0.0, 0.0, 0.0)),
            Vec{3, Float64}((1.0, 0.0, 0.0)),
            Vec{3, Float64}((0.0, 1.0, 0.0)),
            Vec{3, Float64}((0.0, 0.0, 1.0)),
            Vec{3, Float64}((0.5, 0.0, 0.0)),
            Vec{3, Float64}((0.5, 0.5, 0.0)),
            Vec{3, Float64}((0.0, 0.5, 0.0)),
            Vec{3, Float64}((0.0, 0.0, 0.5)),
            Vec{3, Float64}((0.5, 0.0, 0.5)),
            Vec{3, Float64}((0.0, 0.5, 0.5))]
end

# http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
# http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch10.d/AFEM.Ch10.pdf
function shape_value(ip::Lagrange{RefTetrahedron, 2}, ξ::Vec{3}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    ξ_z = ξ[3]
    i == 1  && return (-2 * ξ_x - 2 * ξ_y - 2 * ξ_z + 1) * (-ξ_x - ξ_y - ξ_z + 1)
    i == 2  && return ξ_x * (2 * ξ_x - 1)
    i == 3  && return ξ_y * (2 * ξ_y - 1)
    i == 4  && return ξ_z * (2 * ξ_z - 1)
    i == 5  && return ξ_x * (-4 * ξ_x - 4 * ξ_y - 4 * ξ_z + 4)
    i == 6  && return 4 * ξ_x * ξ_y
    i == 7  && return 4 * ξ_y * (-ξ_x - ξ_y - ξ_z + 1)
    i == 8  && return ξ_z * (-4 * ξ_x - 4 * ξ_y - 4 * ξ_z + 4)
    i == 9  && return 4 * ξ_x * ξ_z
    i == 10 && return 4 * ξ_y * ξ_z
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

##################################
# Lagrange RefHexahedron order 1 #
##################################
getnbasefunctions(::Lagrange{RefHexahedron,1}) = 8

facedof_indices(::Lagrange{RefHexahedron,1}) = ((1,4,3,2), (1,2,6,5), (2,3,7,6), (3,4,8,7), (1,5,8,4), (5,6,7,8))
edgedof_indices(::Lagrange{RefHexahedron,1}) = ((1,2), (2,3), (3,4), (4,1), (5,6), (6,7), (7,8), (8,5), (1,5), (2,6), (3,7), (4,8))

function reference_coordinates(::Lagrange{RefHexahedron,1})
    return [Vec{3, Float64}((-1.0, -1.0, -1.0)),
            Vec{3, Float64}(( 1.0, -1.0, -1.0)),
            Vec{3, Float64}(( 1.0,  1.0, -1.0)),
            Vec{3, Float64}((-1.0,  1.0, -1.0)),
            Vec{3, Float64}((-1.0, -1.0,  1.0)),
            Vec{3, Float64}(( 1.0, -1.0,  1.0)),
            Vec{3, Float64}(( 1.0,  1.0,  1.0)),
            Vec{3, Float64}((-1.0,  1.0,  1.0))]
end

function shape_value(ip::Lagrange{RefHexahedron, 1}, ξ::Vec{3}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    ξ_z = ξ[3]
    i == 1 && return 0.125(1 - ξ_x) * (1 - ξ_y) * (1 - ξ_z)
    i == 2 && return 0.125(1 + ξ_x) * (1 - ξ_y) * (1 - ξ_z)
    i == 3 && return 0.125(1 + ξ_x) * (1 + ξ_y) * (1 - ξ_z)
    i == 4 && return 0.125(1 - ξ_x) * (1 + ξ_y) * (1 - ξ_z)
    i == 5 && return 0.125(1 - ξ_x) * (1 - ξ_y) * (1 + ξ_z)
    i == 6 && return 0.125(1 + ξ_x) * (1 - ξ_y) * (1 + ξ_z)
    i == 7 && return 0.125(1 + ξ_x) * (1 + ξ_y) * (1 + ξ_z)
    i == 8 && return 0.125(1 - ξ_x) * (1 + ξ_y) * (1 + ξ_z)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end


##################################
# Lagrange RefHexahedron order 2 #
##################################
# Based on vtkTriQuadraticHexahedron (see https://kitware.github.io/vtk-examples/site/Cxx/GeometricObjects/IsoparametricCellsDemo/)
getnbasefunctions(::Lagrange{RefHexahedron,2}) = 27

facedof_indices(::Lagrange{RefHexahedron,2}) = (
    (1,4,3,2, 12,11,10,9, 21),
    (1,2,6,5, 9,18,13,17, 22),
    (2,3,7,6, 10,19,14,18, 23),
    (3,4,8,7, 11,20,15,19, 24),
    (1,5,8,4, 17,16,20,12, 25),
    (5,6,7,8, 13,14,15,16, 26),
)
facedof_interior_indices(::Lagrange{RefHexahedron,2}) = (
    (21,), (22,), (23,), (24,), (25,), (26,),
)

edgedof_indices(::Lagrange{RefHexahedron,2}) = (
    (1,2, 9),
    (2,3, 10),
    (3,4, 11),
    (4,1, 12),
    (5,6, 13),
    (6,7, 14),
    (7,8, 15),
    (8,5, 16),
    (1,5, 17),
    (2,6, 18),
    (3,7, 19),
    (4,8, 20),
)
edgedof_interior_indices(::Lagrange{RefHexahedron,2}) = (
    (9,), (10,), (11,), (12,), (13,), (14,), (15,), (16,), (17), (18,), (19,), (20,)
)

celldof_interior_indices(::Lagrange{RefHexahedron,2}) = (27,)

function reference_coordinates(::Lagrange{RefHexahedron,2})
           # vertex
    return [Vec{3, Float64}((-1.0, -1.0, -1.0)), #  1
            Vec{3, Float64}(( 1.0, -1.0, -1.0)), #  2
            Vec{3, Float64}(( 1.0,  1.0, -1.0)), #  3
            Vec{3, Float64}((-1.0,  1.0, -1.0)), #  4
            Vec{3, Float64}((-1.0, -1.0,  1.0)), #  5
            Vec{3, Float64}(( 1.0, -1.0,  1.0)), #  6
            Vec{3, Float64}(( 1.0,  1.0,  1.0)), #  7
            Vec{3, Float64}((-1.0,  1.0,  1.0)), #  8
            # edge
            Vec{3, Float64}(( 0.0, -1.0, -1.0)), #  9
            Vec{3, Float64}(( 1.0,  0.0, -1.0)),
            Vec{3, Float64}(( 0.0,  1.0, -1.0)),
            Vec{3, Float64}((-1.0,  0.0, -1.0)),
            Vec{3, Float64}(( 0.0, -1.0,  1.0)),
            Vec{3, Float64}(( 1.0,  0.0,  1.0)),
            Vec{3, Float64}(( 0.0,  1.0,  1.0)),
            Vec{3, Float64}((-1.0,  0.0,  1.0)),
            Vec{3, Float64}((-1.0, -1.0,  0.0)),
            Vec{3, Float64}(( 1.0, -1.0,  0.0)),
            Vec{3, Float64}(( 1.0,  1.0,  0.0)),
            Vec{3, Float64}((-1.0,  1.0,  0.0)), # 20
            Vec{3, Float64}(( 0.0,  0.0, -1.0)),
            Vec{3, Float64}(( 0.0, -1.0,  0.0)),
            Vec{3, Float64}(( 1.0,  0.0,  0.0)),
            Vec{3, Float64}(( 0.0,  1.0,  0.0)),
            Vec{3, Float64}((-1.0,  0.0,  0.0)),
            Vec{3, Float64}(( 0.0,  0.0,  1.0)), # 26
            # interior
            Vec{3, Float64}((0.0, 0.0, 0.0)),    # 27
            ]
end

function shape_value(ip::Lagrange{RefHexahedron, 2}, ξ::Vec{3, T}, i::Int) where {T}
    # Some local helpers.
    @inline φ₁(x::T) = -0.5*x*(1-x)
    @inline φ₂(x::T) = (1+x)*(1-x)
    @inline φ₃(x::T) = 0.5*x*(1+x)
    (ξ_x, ξ_y, ξ_z) = ξ
    # vertices
    i == 1 && return φ₁(ξ_x) * φ₁(ξ_y) * φ₁(ξ_z)
    i == 2 && return φ₃(ξ_x) * φ₁(ξ_y) * φ₁(ξ_z)
    i == 3 && return φ₃(ξ_x) * φ₃(ξ_y) * φ₁(ξ_z)
    i == 4 && return φ₁(ξ_x) * φ₃(ξ_y) * φ₁(ξ_z)
    i == 5 && return φ₁(ξ_x) * φ₁(ξ_y) * φ₃(ξ_z)
    i == 6 && return φ₃(ξ_x) * φ₁(ξ_y) * φ₃(ξ_z)
    i == 7 && return φ₃(ξ_x) * φ₃(ξ_y) * φ₃(ξ_z)
    i == 8 && return φ₁(ξ_x) * φ₃(ξ_y) * φ₃(ξ_z)
    # edges
    i ==  9 && return φ₂(ξ_x) * φ₁(ξ_y) * φ₁(ξ_z)
    i == 10 && return φ₃(ξ_x) * φ₂(ξ_y) * φ₁(ξ_z)
    i == 11 && return φ₂(ξ_x) * φ₃(ξ_y) * φ₁(ξ_z)
    i == 12 && return φ₁(ξ_x) * φ₂(ξ_y) * φ₁(ξ_z)
    i == 13 && return φ₂(ξ_x) * φ₁(ξ_y) * φ₃(ξ_z)
    i == 14 && return φ₃(ξ_x) * φ₂(ξ_y) * φ₃(ξ_z)
    i == 15 && return φ₂(ξ_x) * φ₃(ξ_y) * φ₃(ξ_z)
    i == 16 && return φ₁(ξ_x) * φ₂(ξ_y) * φ₃(ξ_z)
    i == 17 && return φ₁(ξ_x) * φ₁(ξ_y) * φ₂(ξ_z)
    i == 18 && return φ₃(ξ_x) * φ₁(ξ_y) * φ₂(ξ_z)
    i == 19 && return φ₃(ξ_x) * φ₃(ξ_y) * φ₂(ξ_z)
    i == 20 && return φ₁(ξ_x) * φ₃(ξ_y) * φ₂(ξ_z)
    # faces
    i == 21 && return φ₂(ξ_x) * φ₂(ξ_y) * φ₁(ξ_z)
    i == 22 && return φ₂(ξ_x) * φ₁(ξ_y) * φ₂(ξ_z)
    i == 23 && return φ₃(ξ_x) * φ₂(ξ_y) * φ₂(ξ_z)
    i == 24 && return φ₂(ξ_x) * φ₃(ξ_y) * φ₂(ξ_z)
    i == 25 && return φ₁(ξ_x) * φ₂(ξ_y) * φ₂(ξ_z)
    i == 26 && return φ₂(ξ_x) * φ₂(ξ_y) * φ₃(ξ_z)
    # interior
    i == 27 && return φ₂(ξ_x) * φ₂(ξ_y) * φ₂(ξ_z)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end


#############################
# Lagrange RefPrism order 1 #
#############################
# Build on https://defelement.com/elements/examples/prism-Lagrange-1.html
getnbasefunctions(::Lagrange{RefPrism,1}) = 6

facedof_indices(::Lagrange{RefPrism,1}) = ((1,3,2), (1,2,5,4), (3,1,4,6), (2,3,6,5), (4,5,6))
edgedof_indices(::Lagrange{RefPrism,1}) = ((2,1), (1,3), (1,4), (3,2), (2,5), (3,6), (4,5), (4,6), (6,5))

function reference_coordinates(::Lagrange{RefPrism,1})
    return [Vec{3, Float64}((0.0, 0.0, 0.0)),
            Vec{3, Float64}((1.0, 0.0, 0.0)),
            Vec{3, Float64}((0.0, 1.0, 0.0)),
            Vec{3, Float64}((0.0, 0.0, 1.0)),
            Vec{3, Float64}((1.0, 0.0, 1.0)),
            Vec{3, Float64}((0.0, 1.0, 1.0))]
end

function shape_value(ip::Lagrange{RefPrism,1}, ξ::Vec{3}, i::Int)
    (x,y,z) = ξ
    i == 1 && return 1-x-y -z*(1-x-y)
    i == 2 && return x*(1-z)
    i == 3 && return y*(1-z)
    i == 4 && return z*(1-x-y)
    i == 5 && return x*z
    i == 6 && return y*z
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

#############################
# Lagrange RefPrism order 2 #
#############################
# Build on https://defelement.com/elements/examples/prism-Lagrange-2.html .
# This is simply the tensor-product of a quadratic triangle with a quadratic line.
getnbasefunctions(::Lagrange{RefPrism,2}) = 18

facedof_indices(::Lagrange{RefPrism,2}) = (
    #Vertices| Edges  | Face 
    (1,3,2  , 8,10,7         ),
    (1,2,5,4, 7,11,13,9,   16), 
    (3,1,4,6, 8,9,14,12,   17),
    (2,3,6,5, 10,12,15,11, 18),
    (4,5,6  , 13,15,14       ),
)
facedof_interior_indices(::Lagrange{RefPrism,2}) = (
    #Vertices| Edges  | Face 
    (), 
    (16,), 
    (17,), 
    (18,), 
    (),
)
edgedof_indices(::Lagrange{RefPrism,2}) = (
    #Vert|Edge
    (2,1, 7),
    (1,3, 8),
    (1,4, 9),
    (3,2, 10),
    (2,5, 11),
    (3,6, 12),
    (4,5, 13),
    (4,6, 14),
    (6,5, 15),
)
edgedof_interior_indices(::Lagrange{RefPrism,2}) = (
    #Vert|Edge
    (7,),
    (8,),
    (9,),
    (10,),
    (11,),
    (12,),
    (13,),
    (14,),
    (15,),
)

function reference_coordinates(::Lagrange{RefPrism,2})
    return [Vec{3, Float64}((0.0, 0.0, 0.0)),
            Vec{3, Float64}((1.0, 0.0, 0.0)),
            Vec{3, Float64}((0.0, 1.0, 0.0)),
            Vec{3, Float64}((0.0, 0.0, 1.0)),
            Vec{3, Float64}((1.0, 0.0, 1.0)),
            Vec{3, Float64}((0.0, 1.0, 1.0)),
            Vec{3, Float64}((1/2, 0.0, 0.0)),
            Vec{3, Float64}((0.0, 1/2, 0.0)),
            Vec{3, Float64}((0.0, 0.0, 1/2)),
            Vec{3, Float64}((1/2, 1/2, 0.0)),
            Vec{3, Float64}((1.0, 0.0, 1/2)),
            Vec{3, Float64}((0.0, 1.0, 1/2)),
            Vec{3, Float64}((1/2, 0.0, 1.0)),
            Vec{3, Float64}((0.0, 1/2, 1.0)),
            Vec{3, Float64}((1/2, 1/2, 1.0)),
            Vec{3, Float64}((1/2, 0.0, 1/2)),
            Vec{3, Float64}((0.0, 1/2, 1/2)),
            Vec{3, Float64}((1/2, 1/2, 1/2)),]
end

function shape_value(ip::Lagrange{RefPrism, 2}, ξ::Vec{3}, i::Int)
    (x,y,z) = ξ
    x² = x*x
    y² = y*y
    z² = z*z
    i == 1  && return 4*x²*z² - 6x²*z +2x² +8x*y*z² -12x*y*z +4x*y -6x*z² +9x*z -3x +4y²*z² -6y²*z + 2y² -6y*z² +9y*z -3*y +2z² -3z +1
    i == 2  && return x*(4x*z² -6x*z +2x -2z² +3z -1)
    i == 3  && return y*(4y*z² -6y*z +2y -2z² +3z -1)
    i == 4  && return z*(4x²*z -2x² + 8x*y*z -4x*y -6x*z +3x +4y²*z -2y² -6y*z +3y +2z -1)
    i == 5  && return x*z*(4x*z -2x -2z +1)
    i == 6  && return y*z*(4y*z -2y -2z +1)
    i == 7  && return 4x*(-2x*z² +3x*z -x -2*y*z² +3y*z -y +2z² -3z +1)
    i == 8  && return 4y*(-2x*z² +3x*z -x -2*y*z² +3y*z -y +2z² -3z +1)
    i == 9  && return 4z*(-2x²*z +2x² -4x*y*z +4x*y +3x*z -3x -2y²*z +2y² +3y*z -3y -z +1)
    i == 10 && return 4x*y*(2z² -3z +1)
    i == 11 && return 4x*z*(-2x*z +2x +z -1)
    i == 12 && return 4y*z*(-2y*z +2y +z -1)
    i == 13 && return 4x*z*(-2x*z +x -2y*z +y +2z -1)
    i == 14 && return 4y*z*(-2x*z +x -2y*z +y +2z -1)
    i == 15 && return 4x*y*z*(2z -1)
    i == 16 && return 16x*z*(x*z -x +y*z -y -z +1)
    i == 17 && return 16y*z*(x*z -x +y*z -y -z +1)
    i == 18 && return 16x*y*z*(1 -z)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

###################
# Bubble elements #
###################
"""
Lagrange element with bubble stabilization.
"""
struct BubbleEnrichedLagrange{shape, order, unused} <: ScalarInterpolation{shape, order}
    function BubbleEnrichedLagrange{shape, order}() where {shape <: AbstractRefShape, order}
        new{shape, order, Nothing}()
    end
end

#######################################
# Lagrange-Bubble RefTriangle order 1 #
#######################################
# Taken from https://defelement.com/elements/bubble-enriched-lagrange.html
getnbasefunctions(::BubbleEnrichedLagrange{RefTriangle,1}) = 4

vertexdof_indices(::BubbleEnrichedLagrange{RefTriangle,1}) = ((1,), (2,), (3,))
facedof_indices(::BubbleEnrichedLagrange{RefTriangle,1}) = ((1,2), (2,3), (3,1))
celldof_interior_indices(::BubbleEnrichedLagrange{RefTriangle,1}) = (4,)

function reference_coordinates(::BubbleEnrichedLagrange{RefTriangle,1})
    return [Vec{2, Float64}((1.0, 0.0)),
            Vec{2, Float64}((0.0, 1.0)),
            Vec{2, Float64}((0.0, 0.0)),
            Vec{2, Float64}((1/3, 1/3)),]
end

function shape_value(ip::BubbleEnrichedLagrange{RefTriangle, 1}, ξ::Vec{2}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 1 && return ξ_x*(9ξ_y^2 + 9ξ_x*ξ_y - 9ξ_y + 1)
    i == 2 && return ξ_y*(9ξ_x^2 + 9ξ_x*ξ_y - 9ξ_x + 1)
    i == 3 && return 9ξ_x^2*ξ_y + 9ξ_x*ξ_y^2 - 9ξ_x*ξ_y - ξ_x - ξ_y + 1
    i == 4 && return 27ξ_x*ξ_y*(1 - ξ_x - ξ_y)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

###############
# Serendipity #
###############
struct Serendipity{shape, order, unused} <: ScalarInterpolation{shape,order}
    function Serendipity{shape, order}() where {shape <: AbstractRefShape, order}
        new{shape, order, Nothing}()
    end
end

# Vertices for all Serendipity interpolations are the same
vertexdof_indices(::Serendipity{RefQuadrilateral}) = ((1,),(2,),(3,),(4,))
vertexdof_indices(::Serendipity{RefHexahedron}) = ((1,),(2,),(3,),(4,),(5,),(6,),(7,),(8,))

########################################
# Serendipity RefQuadrilateral order 2 #
########################################
getnbasefunctions(::Serendipity{RefQuadrilateral,2}) = 8
getlowerorder(::Serendipity{RefQuadrilateral,2}) = Lagrange{RefQuadrilateral,1}()

facedof_indices(::Serendipity{RefQuadrilateral,2}) = ((1,2,5), (2,3,6), (3,4,7), (4,1,8))
facedof_interior_indices(::Serendipity{RefQuadrilateral,2}) = ((5,), (6,), (7,), (8,))

function reference_coordinates(::Serendipity{RefQuadrilateral,2})
    return [Vec{2, Float64}((-1.0, -1.0)),
            Vec{2, Float64}(( 1.0, -1.0)),
            Vec{2, Float64}(( 1.0,  1.0)),
            Vec{2, Float64}((-1.0,  1.0)),
            Vec{2, Float64}(( 0.0, -1.0)),
            Vec{2, Float64}(( 1.0,  0.0)),
            Vec{2, Float64}(( 0.0,  1.0)),
            Vec{2, Float64}((-1.0,  0.0))]
end

function shape_value(ip::Serendipity{RefQuadrilateral,2}, ξ::Vec{2}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 1 && return (1 - ξ_x) * (1 - ξ_y) * 0.25(-ξ_x - ξ_y - 1)
    i == 2 && return (1 + ξ_x) * (1 - ξ_y) * 0.25( ξ_x - ξ_y - 1)
    i == 3 && return (1 + ξ_x) * (1 + ξ_y) * 0.25( ξ_x + ξ_y - 1)
    i == 4 && return (1 - ξ_x) * (1 + ξ_y) * 0.25(-ξ_x + ξ_y - 1)
    i == 5 && return 0.5(1 - ξ_x * ξ_x) * (1 - ξ_y)
    i == 6 && return 0.5(1 + ξ_x) * (1 - ξ_y * ξ_y)
    i == 7 && return 0.5(1 - ξ_x * ξ_x) * (1 + ξ_y)
    i == 8 && return 0.5(1 - ξ_x) * (1 - ξ_y * ξ_y)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

#####################################
# Serendipity RefHexahedron order 2 #
#####################################
# Note that second order serendipity hex has no interior face indices.
getnbasefunctions(::Serendipity{RefHexahedron,2}) = 20
getlowerorder(::Serendipity{RefHexahedron,2}) = Lagrange{RefHexahedron,1}()

facedof_indices(::Serendipity{RefHexahedron,2}) = (
    (1,4,3,2, 12,11,10,9),
    (1,2,6,5, 9,18,13,17),
    (2,3,7,6, 10,19,14,18),
    (3,4,8,7, 11,20,15,19),
    (1,5,8,4, 17,16,20,12),
    (5,6,7,8, 13,14,15,16)
)
edgedof_indices(::Serendipity{RefHexahedron,2}) = (
    (1,2, 9),
    (2,3, 10),
    (3,4, 11),
    (4,1, 12),
    (5,6, 13),
    (6,7, 14),
    (7,8, 15),
    (8,5, 16),
    (1,5, 17),
    (2,6, 18),
    (3,7, 19),
    (4,8, 20),
)

edgedof_interior_indices(::Serendipity{RefHexahedron,2}) = (
    (9,), (10,), (11,), (12,), (13,), (14,), (15,), (16,), (17), (18,), (19,), (20,)
)

function reference_coordinates(::Serendipity{RefHexahedron,2})
    return [Vec{3, Float64}((-1.0, -1.0, -1.0)),
            Vec{3, Float64}(( 1.0, -1.0, -1.0)),
            Vec{3, Float64}(( 1.0,  1.0, -1.0)),
            Vec{3, Float64}((-1.0,  1.0, -1.0)),
            Vec{3, Float64}((-1.0, -1.0,  1.0)),
            Vec{3, Float64}(( 1.0, -1.0,  1.0)),
            Vec{3, Float64}(( 1.0,  1.0,  1.0)),
            Vec{3, Float64}((-1.0,  1.0,  1.0)),
            Vec{3, Float64}((0.0, -1.0, -1.0)),
            Vec{3, Float64}((1.0, 0.0, -1.0)),
            Vec{3, Float64}((0.0, 1.0, -1.0)),
            Vec{3, Float64}((-1.0, 0.0, -1.0)),
            Vec{3, Float64}((0.0, -1.0, 1.0)),
            Vec{3, Float64}((1.0, 0.0, 1.0)),
            Vec{3, Float64}((0.0, 1.0, 1.0)),
            Vec{3, Float64}((-1.0, 0.0, 1.0)),
            Vec{3, Float64}((-1.0, -1.0, 0.0)),
            Vec{3, Float64}((1.0, -1.0, 0.0)),
            Vec{3, Float64}((1.0, 1.0, 0.0)),
            Vec{3, Float64}((-1.0, 1.0, 0.0)),]
end

function shape_value(ip::Serendipity{RefHexahedron, 2}, ξ::Vec{3}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    ξ_z = ξ[3]
    i == 1 && return 0.125(1 - ξ_x) * (1 - ξ_y) * (1 - ξ_z) - 0.5(shape_value(ip, ξ, 12) + shape_value(ip, ξ, 9) + shape_value(ip, ξ, 17))
    i == 2 && return 0.125(1 + ξ_x) * (1 - ξ_y) * (1 - ξ_z) - 0.5(shape_value(ip, ξ, 9) + shape_value(ip, ξ, 10) + shape_value(ip, ξ, 18))
    i == 3 && return 0.125(1 + ξ_x) * (1 + ξ_y) * (1 - ξ_z) - 0.5(shape_value(ip, ξ, 10) + shape_value(ip, ξ, 11) + shape_value(ip, ξ, 19))
    i == 4 && return 0.125(1 - ξ_x) * (1 + ξ_y) * (1 - ξ_z) - 0.5(shape_value(ip, ξ, 11) + shape_value(ip, ξ, 12) + shape_value(ip, ξ, 20))
    i == 5 && return 0.125(1 - ξ_x) * (1 - ξ_y) * (1 + ξ_z) - 0.5(shape_value(ip, ξ, 16) + shape_value(ip, ξ, 13) + shape_value(ip, ξ, 17))
    i == 6 && return 0.125(1 + ξ_x) * (1 - ξ_y) * (1 + ξ_z) - 0.5(shape_value(ip, ξ, 13) + shape_value(ip, ξ, 14) + shape_value(ip, ξ, 18))
    i == 7 && return 0.125(1 + ξ_x) * (1 + ξ_y) * (1 + ξ_z) - 0.5(shape_value(ip, ξ, 14) + shape_value(ip, ξ, 15) + shape_value(ip, ξ, 19))
    i == 8 && return 0.125(1 - ξ_x) * (1 + ξ_y) * (1 + ξ_z) - 0.5(shape_value(ip, ξ, 15) + shape_value(ip, ξ, 16) + shape_value(ip, ξ, 20))
    i == 9 && return 0.25(1 - ξ_x^2) * (1 - ξ_y) * (1 - ξ_z)
    i == 10 && return 0.25(1 + ξ_x) * (1 - ξ_y^2) * (1 - ξ_z)
    i == 11 && return 0.25(1 - ξ_x^2) * (1 + ξ_y) * (1 - ξ_z)
    i == 12 && return 0.25(1 - ξ_x) * (1 - ξ_y^2) * (1 - ξ_z)
    i == 13 && return 0.25(1 - ξ_x^2) * (1 - ξ_y) * (1 + ξ_z)
    i == 14 && return 0.25(1 + ξ_x) * (1 - ξ_y^2) * (1 + ξ_z)
    i == 15 && return 0.25(1 - ξ_x^2) * (1 + ξ_y) * (1 + ξ_z)
    i == 16 && return 0.25(1 - ξ_x) * (1 - ξ_y^2) * (1 + ξ_z)
    i == 17 && return 0.25(1 - ξ_x) * (1 - ξ_y) * (1 - ξ_z^2)
    i == 18 && return 0.25(1 + ξ_x) * (1 - ξ_y) * (1 - ξ_z^2)
    i == 19 && return 0.25(1 + ξ_x) * (1 + ξ_y) * (1 - ξ_z^2)
    i == 20 && return 0.25(1 - ξ_x) * (1 + ξ_y) * (1 - ξ_z^2)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end


#############################
# Crouzeix–Raviart Elements #
#############################
"""
Classical non-conforming Crouzeix–Raviart element.

For details we refer to the original paper:
M. Crouzeix and P. Raviart. "Conforming and nonconforming finite element 
methods for solving the stationary Stokes equations I." ESAIM: Mathematical Modelling 
and Numerical Analysis-Modélisation Mathématique et Analyse Numérique 7.R3 (1973): 33-75.
"""
struct CrouzeixRaviart{shape, order, unused} <: ScalarInterpolation{shape, order}
    CrouzeixRaviart{RefTriangle, 1}() = new{RefTriangle, 1, Nothing}()
end

getnbasefunctions(::CrouzeixRaviart) = 3

facedof_indices(::CrouzeixRaviart) = ((1,), (2,), (3,))
facedof_interior_indices(::CrouzeixRaviart) = ((1,), (2,), (3,))

function reference_coordinates(::CrouzeixRaviart)
    return [Vec{2, Float64}((0.5, 0.5)),
            Vec{2, Float64}((0.0, 0.5)),
            Vec{2, Float64}((0.5, 0.0))]
end

function shape_value(ip::CrouzeixRaviart, ξ::Vec{2}, i::Int)
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 1 && return 2*ξ_x + 2*ξ_y - 1.0
    i == 2 && return 1.0 - 2*ξ_x
    i == 3 && return 1.0 - 2*ξ_y
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

##################################################
# VectorizedInterpolation{<:ScalarInterpolation} #
##################################################

struct VectorizedInterpolation{vdim, refshape, order, SI <: ScalarInterpolation{refshape, order}} <: VectorInterpolation{vdim, refshape,order}
    ip::SI
    function VectorizedInterpolation{vdim}(ip::SI) where {vdim, refshape, order, SI <: ScalarInterpolation{refshape, order}}
        return new{vdim, refshape, order, SI}(ip)
    end
end

# Vectorize to reference dimension by default
function VectorizedInterpolation(ip::ScalarInterpolation{shape}) where {refdim, shape <: AbstractRefShape{refdim}}
    return VectorizedInterpolation{refdim}(ip)
end

Base.:(^)(ip::ScalarInterpolation, vdim::Int) = VectorizedInterpolation{vdim}(ip)
function Base.literal_pow(::typeof(^), ip::ScalarInterpolation, ::Val{vdim}) where vdim
    return VectorizedInterpolation{vdim}(ip)
end

# Helper to get number of copies for DoF distribution
get_n_copies(::VectorizedInterpolation{vdim}) where vdim = vdim

function getnbasefunctions(ipv::VectorizedInterpolation{vdim}) where vdim
    return vdim * getnbasefunctions(ipv.ip)
end
function shape_value(ipv::VectorizedInterpolation{vdim, shape}, ξ::Vec{refdim, T}, I::Int) where {vdim, refdim, shape <: AbstractRefShape{refdim}, T}
    i0, c0 = divrem(I - 1, vdim)
    i = i0 + 1
    c = c0 + 1
    v = shape_value(ipv.ip, ξ, i)
    return Vec{vdim, T}(j -> j == c ? v : zero(v))
end

# vdim == refdim
function shape_gradient_and_value(ipv::VectorizedInterpolation{dim, shape}, ξ::Vec{dim}, I::Int) where {dim, shape <: AbstractRefShape{dim}}
    return invoke(shape_gradient_and_value, Tuple{Interpolation, Vec, Int}, ipv, ξ, I)
end
# vdim != refdim
function shape_gradient_and_value(ipv::VectorizedInterpolation{vdim, shape}, ξ::V, I::Int) where {vdim, refdim, shape <: AbstractRefShape{refdim}, T, V <: Vec{refdim, T}}
    # Load with dual numbers and compute the value
    f = x -> shape_value(ipv, x, I)
    ξd = Tensors._load(ξ, Tensors.Tag(f, V))
    value_grad = f(ξd)
    # Extract the value and gradient
    val = Vec{vdim, T}(i -> Tensors.value(value_grad[i]))
    grad = zero(MMatrix{vdim, refdim, T})
    for (i, vi) in pairs(value_grad)
        p = Tensors.partials(vi)
        for (j, pj) in pairs(p)
            grad[i, j] = pj
        end
    end
    return SMatrix(grad), val
end

reference_coordinates(ip::VectorizedInterpolation) = reference_coordinates(ip.ip)
