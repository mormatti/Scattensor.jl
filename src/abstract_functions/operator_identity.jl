"""
    operator_identity(args...)

Abstract interface for constructing an identity operator.

Concrete implementations are provided by backends (e.g. dense/sparse matrices, MPOs). Typical method
signatures in this package include:
- `operator_identity(MatrixType, n::Integer)` for matrix identities
- `operator_identity(MPO, d::Integer, L::Integer)` for an identity MPO on a length-`L` chain
"""
function operator_identity end