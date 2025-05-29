"""
        operator_reflection(args...)

    Return the translation operator from specified arguments.
    A reflection operator T shifts cyclically on the right (e.g. ABCDE -> EABCD) the system.
    It must have a concrete implementation in a specific type.
    """
function operator_reflection end