# Defines CellScalarValues and CellVectorValues and common methods
"""
    CellScalarValues([::Type{T}], quad_rule::QuadratureRule, func_interpol::Interpolation, [geom_interpol::Interpolation])
    CellVectorValues([::Type{T}], quad_rule::QuadratureRule, func_interpol::Interpolation, [geom_interpol::Interpolation])

A `CellValues` object facilitates the process of evaluating values of shape functions, gradients of shape functions,
values of nodal functions, gradients and divergences of nodal functions etc. in the finite element cell. There are
two different types of `CellValues`: `CellScalarValues` and `CellVectorValues`. As the names suggest, `CellScalarValues`
utilizes scalar shape functions and `CellVectorValues` utilizes vectorial shape functions. For a scalar field, the
`CellScalarValues` type should be used. For vector field, both subtypes can be used.

**Arguments:**
* `T`: an optional argument (default to `Float64`) to determine the type the internal data is stored as.
* `quad_rule`: an instance of a [`QuadratureRule`](@ref)
* `func_interpol`: an instance of an [`Interpolation`](@ref) used to interpolate the approximated function
* `geom_interpol`: an optional instance of a [`Interpolation`](@ref) which is used to interpolate the geometry

**Common methods:**
* [`reinit!`](@ref)
* [`getnquadpoints`](@ref)
* [`getdetJdV`](@ref)

* [`shape_value`](@ref)
* [`shape_gradient`](@ref)
* [`shape_symmetric_gradient`](@ref)
* [`shape_divergence`](@ref)

* [`function_value`](@ref)
* [`function_gradient`](@ref)
* [`function_symmetric_gradient`](@ref)
* [`function_divergence`](@ref)
* [`spatial_coordinate`](@ref)
"""
CellValues

# Scalar fields:
#  - IP     <: ScalarInterpolation{refdim, refshape, order}
#  - N_t    <: Number
#  - dNdx_t <: Vec{spacedim}
#  - dNdξ_t <: Vec{refdim}
#  - T      <: Number
#  - M_t    <: Number
#  - dMdξ_t <: Vec{refdim}

# Vector fields:
#  - IP     <: VectorInterpolation{refdim, refshape, order}
#  - N_t    <: Vec{valdim}
#  - dNdx_t <: Tensor{valdim, spacedim}
#  - dNdξ_t <: Tensor{valdim, refdim}
#  - T      <: Number
#  - M_t    <: Number
#  - dMdξ_t <: Vec{refdim}

struct ReferenceMapping{M_t, dMdξ_t, W_t, QR_t, GIP <: ScalarInterpolation}
    M::Matrix{M_t}
    dMdξ::Matrix{dMdξ_t}
    detJdV::Vector{W_t}
    qr::QR_t
    interpolation::GIP
end

# Default the number type to the quad_rule type
ReferenceMapping(qr::QuadratureRule{<:Any,<:Any,T}, ip::ScalarInterpolation) where T = ReferenceMapping(T, qr, ip)

function ReferenceMapping(::Type{T}, qr::QuadratureRule{refdim, refshape},
        ip::ScalarInterpolation{refdim, refshape}) where {T <: Number, refdim, refshape}
    n_qpoints = length(getweights(qr))
    n_basefunc = getnbasefunctions(ip)
    M    = fill(zero(T)             * T(NaN), n_basefunc, n_qpoints)
    dMdξ = fill(zero(Vec{refdim,T}) * T(NaN), n_basefunc, n_qpoints)
    for (qp, ξ) in pairs(qr.points)
        for i in 1:n_basefunc
            dMdξ[i, qp], M[i, qp] = gradient(ξ -> value(ip, i, ξ), ξ, :all)
        end
    end
    detJdV = fill(T(NaN), n_qpoints)
    ReferenceMapping{eltype(M), eltype(dMdξ), eltype(detJdV), typeof(qr), typeof(ip)}(
        M, dMdξ, detJdV, qr, ip,
    )
end

struct CellValues{N_t, dNdx_t, dNdξ_t, RM, IP, GIP, QR_t} <: FEValues
    N::Matrix{N_t}
    dNdx::Matrix{dNdx_t}
    dNdξ::Matrix{dNdξ_t}
    func_interp::IP
    qr::QR_t
    reference_mapping::RM
end

const CellScalarValues = CellValues{<:Number}
const CellVectorValues = CellValues{<:AbstractVector}

# _valuetype(::CellValues{N_t}) where N_t = N_t
# _gradienttype(::CellValues{<:Any,dNdx_t}) where dNdx_t = dNdx_t

function CellValues(
        qr::QuadratureRule{refdim, refshape, T},
        ip::ScalarInterpolation{refdim, refshape, order},
        reference_mapping = ReferenceMapping(qr, Lagrange{refdim, refshape, 1}()),
    ) where {refdim, refshape, T, order}
    @assert qr === reference_mapping.qr
    # Allocate buffers
    n_qpoints = length(getweights(qr))
    n_func_basefuncs = getnbasefunctions(ip)
    N    = fill(zero(T)                          * T(NaN), n_func_basefuncs, n_qpoints)
    dNdx = fill(zero(Vec{#=spacedim=#refdim, T}) * T(NaN), n_func_basefuncs, n_qpoints)
    dNdξ = fill(zero(Vec{refdim, T}            ) * T(NaN), n_func_basefuncs, n_qpoints)
    # Compute value and gradient of shape functions
    for (qp, ξ) in pairs(qr.points)
        for i in 1:n_func_basefuncs
            dNdξ[i, qp], N[i, qp] = gradient(ξ -> value(ip, i, ξ), ξ, :all)
        end
    end
    CellValues{eltype(N), eltype(dNdx), eltype(dNdξ), typeof(ip), typeof(reference_mapping)}(N, dNdx, dNdξ, ip, qr, reference_mapping)
end

# CellScalarValues
# struct CellScalarValues{dim,T<:Real,refshape<:AbstractRefShape} <: CellValues{dim,T,refshape}
#     N::Matrix{T}
#     dNdx::Matrix{Vec{dim,T}}
#     dNdξ::Matrix{Vec{dim,T}}
#     detJdV::Vector{T}
#     M::Matrix{T}
#     dMdξ::Matrix{Vec{dim,T}}
#     qr::QuadratureRule{dim,refshape,T}
#     # The following fields are deliberately abstract -- they are never used in
#     # performance critical code, just stored here for convenience.
#     func_interp::Interpolation{dim,refshape}
#     geo_interp::Interpolation{dim,refshape}
# end

# function CellScalarValues(quad_rule::QuadratureRule, func_interpol::Interpolation,
#         geom_interpol::Interpolation=func_interpol)
#     CellScalarValues(Float64, quad_rule, func_interpol, geom_interpol)
# end

# function CellScalarValues(::Type{T}, quad_rule::QuadratureRule{dim,shape}, func_interpol::Interpolation,
#         geom_interpol::Interpolation=func_interpol) where {dim,T,shape<:AbstractRefShape}

#     @assert getdim(func_interpol) == getdim(geom_interpol)
#     @assert getrefshape(func_interpol) == getrefshape(geom_interpol) == shape
#     n_qpoints = length(getweights(quad_rule))

#     # Function interpolation
#     n_func_basefuncs = getnbasefunctions(func_interpol)
#     N    = fill(zero(T)          * T(NaN), n_func_basefuncs, n_qpoints)
#     dNdx = fill(zero(Vec{dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)
#     dNdξ = fill(zero(Vec{dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)

#     # Geometry interpolation
#     n_geom_basefuncs = getnbasefunctions(geom_interpol)
#     M    = fill(zero(T)          * T(NaN), n_geom_basefuncs, n_qpoints)
#     dMdξ = fill(zero(Vec{dim,T}) * T(NaN), n_geom_basefuncs, n_qpoints)

#     for (qp, ξ) in enumerate(quad_rule.points)
#         for i in 1:n_func_basefuncs
#             dNdξ[i, qp], N[i, qp] = gradient(ξ -> value(func_interpol, i, ξ), ξ, :all)
#         end
#         for i in 1:n_geom_basefuncs
#             dMdξ[i, qp], M[i, qp] = gradient(ξ -> value(geom_interpol, i, ξ), ξ, :all)
#         end
#     end

#     detJdV = fill(T(NaN), n_qpoints)

#     CellScalarValues{dim,T,shape}(N, dNdx, dNdξ, detJdV, M, dMdξ, quad_rule, func_interpol, geom_interpol)
# end

# # CellVectorValues
# struct CellVectorValues{#=vdim,=#dim,T<:Real,refshape<:AbstractRefShape,M} <: CellValues{dim,T,refshape}
#     N::Matrix{Vec{dim,T}}
#     dNdx::Matrix{Tensor{2,dim,T,M}}
#     dNdξ::Matrix{Tensor{2,dim,T,M}}
#     detJdV::Vector{T}
#     M::Matrix{T}
#     dMdξ::Matrix{Vec{dim,T}}
#     qr::QuadratureRule{dim,refshape,T}
#     # The following fields are deliberately abstract -- they are never used in
#     # performance critical code, just stored here for convenience.
#     func_interp::#=Vector=#Interpolation{#=vdim,=#dim,refshape}
#     geo_interp::Interpolation{dim,refshape}
# end

# # This helper is used to derive a corresponding geometric interpolation when not passed to
# # the constructor explicitly. For now this is only used to de-vectorize a
# # VectorizedInterpolation since we want to use the scalar version internally and let the
# # coordinates vectorize the interpolation instead.
# derive_geometric_interpolation(ip::VectorizedInterpolation) = ip.ip
# derive_geometric_interpolation(ip) = ip

# function CellVectorValues(quad_rule::QuadratureRule, func_interpol::Interpolation,
#         geom_interpol::Interpolation=derive_geometric_interpolation(func_interpol))
#     CellVectorValues(Float64, quad_rule, func_interpol, geom_interpol)
# end

# # TODO: Maybe deprecate this auto-vectorizing method?
# function CellVectorValues(::Type{T}, quad_rule::QuadratureRule, func_interpol::ScalarInterpolation,
#         geom_interpol::Interpolation=func_interpol) where {T}
#     return CellVectorValues(T, quad_rule, VectorizedInterpolation(func_interpol), geom_interpol)
# end

# function CellVectorValues(::Type{T}, quad_rule::QuadratureRule{dim,shape}, func_interpol::VectorInterpolation{vdim,dim,shape},
#         geom_interpol::Interpolation=derive_geometric_interpolation(func_interpol)) where {vdim,dim,T,shape<:AbstractRefShape}

#     @assert getdim(func_interpol) == getdim(geom_interpol)
#     @assert getrefshape(func_interpol) == getrefshape(geom_interpol) == shape
#     n_qpoints = length(getweights(quad_rule))

#     # Function interpolation
#     n_func_basefuncs = getnbasefunctions(func_interpol)
#     N    = fill(zero(Vec{dim,T})      * T(NaN), n_func_basefuncs, n_qpoints)
#     dNdx = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)
#     dNdξ = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)

#     # Geometry interpolation
#     n_geom_basefuncs = getnbasefunctions(geom_interpol)
#     M    = fill(zero(T)          * T(NaN), n_geom_basefuncs, n_qpoints)
#     dMdξ = fill(zero(Vec{dim,T}) * T(NaN), n_geom_basefuncs, n_qpoints)

#     for (qp, ξ) in enumerate(quad_rule.points)
#         for basefunc in 1:n_func_basefuncs
#             dNdξ[basefunc, qp], N[basefunc, qp] = gradient(ξ -> value(func_interpol, basefunc, ξ), ξ, :all)
#         end
#         for basefunc in 1:n_geom_basefuncs
#             dMdξ[basefunc, qp], M[basefunc, qp] = gradient(ξ -> value(geom_interpol, basefunc, ξ), ξ, :all)
#         end
#     end

#     detJdV = fill(T(NaN), n_qpoints)
#     MM = Tensors.n_components(Tensors.get_base(eltype(dNdx)))

#     CellVectorValues{dim,T,shape,MM}(N, dNdx, dNdξ, detJdV, M, dMdξ, quad_rule, func_interpol, geom_interpol)
# end

function reinit!(cv::CellValues{dim}, x::AbstractVector{Vec{dim,T}}) where {dim,T}
    n_geom_basefuncs = getngeobasefunctions(cv)
    n_func_basefuncs = getnbasefunctions(cv)
    length(x) == n_geom_basefuncs || throw_incompatible_coord_length(length(x), n_geom_basefuncs)

    @inbounds for i in 1:length(cv.qr.weights)
        w = cv.qr.weights[i]
        fecv_J = zero(Tensor{2,dim})
        for j in 1:n_geom_basefuncs
            fecv_J += x[j] ⊗ cv.dMdξ[j, i]
        end
        detJ = det(fecv_J)
        detJ > 0.0 || throw_detJ_not_pos(detJ)
        cv.detJdV[i] = detJ * w
        Jinv = inv(fecv_J)
        for j in 1:n_func_basefuncs
            cv.dNdx[j, i] = cv.dNdξ[j, i] ⋅ Jinv
        end
    end
end
