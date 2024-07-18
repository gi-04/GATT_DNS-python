import scipy.io as sio
import sys

def save(*varargin):
    args = list(varargin)
    
    for i in range(1, len(args)):
        if not args[i].startswith('-'):
            # Get the variable from the caller's namespace
            caller_locals = sys._getframe(1).f_locals
            args[i] = caller_locals.get(args[i], args[i])
        
        if args[i] == '-v7.3':
            args[i] = '-v7'
    
    # Use scipy.io.savemat to save the data
    sio.savemat(args[0], {k: v for k, v in zip(args[1::2], args[2::2])}, do_compression=True)



