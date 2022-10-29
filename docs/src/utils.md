# Utils
```@meta
CurrentModule = Utils
```

```@docs
a_matrix(nstates::Integer=2)
on(op::AbstractMatrix{<:Number}, q::Integer, n::Integer)
kron_concat(ops::AbstractVector{<:AbstractMatrix{<:Number}})
kron_concat(ops::AbstractMatrix{<:Number}, n::Integer)
algebra(
    n::Integer,
    m::Integer=2;
    basis::Union{AbstractMatrix{<:Number},Nothing}=nothing,
)
```
