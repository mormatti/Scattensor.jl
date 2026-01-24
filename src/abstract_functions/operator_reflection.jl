"""
    operator_reflection(args...)

Abstract interface for constructing a reflection (parity) operator.

The reflection operator `R` reverses the order of sites (e.g. `ABCDE -> EDCBA`). Concrete
implementations are provided by backends (e.g. matrices or MPOs).
"""
function operator_reflection end