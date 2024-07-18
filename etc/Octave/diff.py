import numpy as np
def diff(x, n=None, dim=None):
    if dim is None:
        dim = next((i for i, s in enumerate(x.shape) if s > 1), None)
        if dim is None:
            return np.array([])
    
    if n is None:
        n = 1
    
    if n >= x.shape[dim]:
        return np.array([])
    
    return np.diff(x, n=n, axis=dim)

