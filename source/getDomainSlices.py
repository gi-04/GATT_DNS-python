def get_domain_slices(n, p):
    """
    This function computes the size of each slice the same way the 2decomp library does.
    It is used to correctly distribute the info across the processes.
    The distribution is as even as possible. If an uneven distribution is needed, extra nodes are placed first in the last slices.
    For example, if 10 nodes are divided across 3 slices, the division would be 3 3 4.
    """
    import numpy as np

    if n == 1:
        return np.array([[1], [1]])

    n_points_base = n // p
    n_ceil = n - n_points_base * p

    n_points = np.ones(p, dtype=int) * n_points_base
    n_points[-n_ceil:] = n_points_base + 1

    slices = np.zeros((2, p), dtype=int)
    slices[1, :] = np.cumsum(n_points)
    slices[0, 0] = 1
    slices[0, 1:] = slices[1, :-1] + 1

    return slices



