import numpy as np

def unique(A, op1=None, op2=None):
    if op1 is None:
        op1 = ''
    if op2 is None:
        op2 = ''

    if op1 == 'sorted':
        op1 = 'first'
    if op2 == 'sorted':
        op2 = 'first'

    if op1 != 'stable' and op2 != 'stable':
        if not op1:
            C, ia, ic = np.unique(A, return_index=True, return_inverse=True)
        elif not op2:
            C, ia, ic = np.unique(A, return_index=True, return_inverse=True, method=op1)
        else:
            C, ia, ic = np.unique(A, return_index=True, return_inverse=True, method=op1, axis=op2)
        return C, ia, ic

    if op1 == 'rows':
        _, ia, ic = np.unique(A, axis=0, return_index=True, return_inverse=True)
        ia.sort()
        icnew = np.zeros_like(ic)
        for i, idx in enumerate(ia):
            icnew[ic == idx] = i + 1
        ic = icnew
        C = A[ia, :]
    else:
        _, ia, ic = np.unique(A, return_index=True, return_inverse=True)
        ia.sort()
        icnew = np.zeros_like(ic)
        for i, idx in enumerate(ia):
            icnew[ic == idx] = i + 1
        ic = icnew
        C = A[ia]

    return C, ia, ic

def unique_builtin(x, *args):
    if not (isinstance(x, (np.ndarray, np.bool_, str, list))):
        raise ValueError("unique: X must be an array or cell array of strings")

    optrows = any(arg == 'rows' for arg in args)
    optfirst = any(arg == 'first' for arg in args)
    optlast = any(arg == 'last' for arg in args)
    if optfirst and optlast:
        raise ValueError('unique: cannot specify both "first" and "last"')
    elif optfirst + optlast + optrows != len(args):
        raise ValueError("unique: invalid option")

    if optrows and isinstance(x, list):
        print('unique: "rows" is ignored for cell arrays')
        optrows = False

    if isinstance(x, np.ndarray) and x.ndim == 2 and not optrows and len(args) <= 1:
        if np.count_nonzero(x) < x.size:
            y = np.unique(np.concatenate(([0], x.ravel())), *args)
        else:
            y = np.unique(x.flatten(), *args)
        return y, [], []

    if optrows:
        n = x.shape[0]
        dim = 1
    else:
        n = x.size
        dim = 2 if x.ndim == 1 else 1

    y = x.copy()
    if n == 0:
        if not optrows and x.size > 0:
            if isinstance(y, list):
                y = []
            else:
                y = np.zeros(0, x.dtype)
        return y, [], []
    elif n == 1:
        return y, [1], [1]

    if optrows:
        if len(args) > 1:
            y, i = np.unique(y, axis=0, return_index=True)
        else:
            y = np.unique(y, axis=0)
        match = np.all(y[:-1] == y[1:], axis=1)
        y = y[~match]
    else:
        y = y.flatten()
        if len(args) > 1:
            y, i = np.unique(y, return_index=True)
        else:
            y = np.unique(y)
        match = y[:-1] == y[1:]
        y = y[~match]

    if len(args) > 2:
        j = i.copy()
        if dim == 1:
            j[i] = np.cumsum(np.concatenate(([1], ~match)))
        else:
            j[i] = np.cumsum(np.concatenate(([1], ~match)))
    else:
        j = []

    if len(args) > 1:
        idx = np.where(match)[0]
        if optfirst:
            idx += 1
        i = np.delete(i, idx)

    return y, i, j


