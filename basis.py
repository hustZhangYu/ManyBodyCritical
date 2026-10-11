"""Particle-number basis labels without constructing a dense number operator."""


def nlabel(L, n):
    """Return the ascending basis indices used by the original Kronecker basis.

    A zero bit denotes the occupied state because the local number operator is
    diag(1, 0). Thus particle number is L minus the number of set bits. The
    ordering and out-of-range-particle-number empty result remain unchanged.
    Only the returned list is stored; no 2**L by 2**L matrix is allocated.
    """
    return [index for index in range(2**L)
            if L - bin(index).count("1") == n]
