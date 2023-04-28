# Reference cells

The reference cells are used to i) define grid cells, ii) define shape functions, and iii)
define quadrature rules. The numbering of vertices, edges, faces are visualized below. See also
[`FerriteViz.elementinfo`](https://ferrite-fem.github.io/FerriteViz.jl/dev/api/#FerriteViz.elementinfo).

## Numbering and directions of entities

The local numbering of vertices, edges, and faces, on the reference cells, specifies, for
example, the cell node order, and the order in which interpolations distribute their DoFs.
It is important the the same convention is used everywhere to not cause mismatches. The
convention adopted by Ferrite.jl is documented below.

##### Vertices
Vertices (nodes) are numbered based on their coordinate in the reference system in *reverse*
lexicographic ordering, i.e. the inner-most numbering loop is for the ``ξ₁`` direction. See
e.g. [`RefQuadrilateral`](@ref Ferrite.RefQuadrilateral) for an example.

##### Edges
Edges are defined by the 2-tuple of vertices it connects, starting with the smallest vertex
number. To enumerate the edges of the reference cell lexicographic ordering is used on the
vertex tuples. See e.g. [`RefQuadrilateral`](@ref Ferrite.RefQuadrilateral) for an example.

##### Faces
Faces are defined by the n-tuple of vertices constructing the face in anti-clockwise order
(viewing the face from the outside of the reference cell), starting with the smallest vertex
number. To enumerate the faces of the reference cell, lexicographic ordering is used on the
n-tuples. See e.g. ([`RefHexahedron`](@ref Ferrite.RefHexahedron)) for an example.

### `AbstractRefShape` subtypes

```@docs
Ferrite.AbstractRefShape
Ferrite.RefLine
Ferrite.RefTriangle
Ferrite.RefQuadrilateral
Ferrite.RefTetrahedron
Ferrite.RefHexahedron
Ferrite.RefPrism
```
