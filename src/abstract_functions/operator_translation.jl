"""
    operator_translation(args...)

Abstract interface for constructing a translation operator.

The translation operator `T` performs a cyclic shift by one site (typically to the right, e.g.
`ABCDE -> EABCD`). Concrete implementations are provided by backends (e.g. matrices or MPOs).
"""
function operator_translation end