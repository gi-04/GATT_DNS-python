import numpy as np

def getMatrixTypeBlocks(typeMap, p_row, p_col):
    """
    This function transforms an array of derivative types into a list of coordinates for each type.
    For example, for derivatives in X, the input matrix will be ny x nz and each entry is the type of derivative to be applied at each row.
    The output is a list with 5 columns. The first is the derivative type, columns 2 to 5 are the limits for each region in Y and Z.
    One list is created for each MPI process.
    """
    
    maxBlockSize = 128  # Maximum number of rows per block, larger blocks will be divided
    
    blocks = [None] * (p_row * p_col)
    
    J, K = typeMap.shape
    
    # Get all blocks
    allBlocks = []
    for k in range(K):
        starts = np.concatenate(([0], np.where(np.diff(typeMap[:, k]) != 0)[0] + 1))
        ends = np.concatenate((starts[1:] - 1, [J - 1]))
        allBlocks.extend(np.column_stack((typeMap[starts, k], starts, ends, np.full(len(starts), k), np.full(len(starts), k))))
    
    allBlocks = np.array(allBlocks)
    
    # Merge blocks
    i = 0
    while i < allBlocks.shape[0]:
        j = i + 1
        while j < allBlocks.shape[0]:
            if np.array_equal(allBlocks[i, :3], allBlocks[j, :3]) and allBlocks[i, 4] + 1 == allBlocks[j, 3]:
                allBlocks[i, 4] = allBlocks[j, 4]
                allBlocks = np.delete(allBlocks, j, axis=0)
            else:
                j += 1
        i += 1
    
    # Divide blocks for processors
    domainSlicesY = getDomainSlices(J, p_row)
    domainSlicesZ = getDomainSlices(K, p_col)
    
    for j in range(p_row):
        for k in range(p_col):
            nProc = k + j * p_col
            
            bL = allBlocks.copy()
            
            Ji = domainSlicesY[0, j]
            Jf = domainSlicesY[1, j]
            Ki = domainSlicesZ[0, k]
            Kf = domainSlicesZ[1, k]
            
            bL = bL[(bL[:, 1] <= Jf) & (bL[:, 2] >= Ji) & (bL[:, 3] <= Kf) & (bL[:, 4] >= Ki)]
            
            bL[:, 1] = np.maximum(bL[:, 1], Ji)
            bL[:, 2] = np.minimum(bL[:, 2], Jf)
            bL[:, 3] = np.maximum(bL[:, 3], Ki)
            bL[:, 4] = np.minimum(bL[:, 4], Kf)
            
            blocks[nProc] = bL
    
    # Reduce block sizes if needed
    if not np.isinf(maxBlockSize):
        for nProc in range(p_col * p_row):
            bL = blocks[nProc]
            bLnew = np.zeros((0, 5), dtype=int)
            for j in range(bL.shape[0]):
                iSize = bL[j, 2] - bL[j, 1] + 1
                jSize = bL[j, 4] - bL[j, 3] + 1
                bSize = iSize * jSize
                nSlices = int(np.ceil(bSize / maxBlockSize))
                if nSlices == 1:
                    bLnew = np.vstack((bLnew, bL[j, :]))
                else:
                    iSlices = int(np.ceil(nSlices / jSize))
                    jSlices = int(np.floor(nSlices / iSlices))
                    
                    iSlicesInd = getDomainSlices(iSize, iSlices)
                    jSlicesInd = getDomainSlices(jSize, jSlices)
                    
                    iSlicesInd = bL[j, 1] + iSlicesInd.T - 1
                    jSlicesInd = bL[j, 3] + jSlicesInd.T - 1
                    
                    for ii in range(iSlices):
                        for jj in range(jSlices):
                            bLnew = np.vstack((bLnew, [bL[j, 0], iSlicesInd[ii, 0], iSlicesInd[ii, 1], jSlicesInd[jj, 0], jSlicesInd[jj, 1]]))
            blocks[nProc] = bLnew
    
    return blocks

def getDomainSlices(size, numSlices):
    """
    Helper function to divide a domain into slices.
    """
    sliceSize = size // numSlices
    remainder = size % numSlices
    slices = np.zeros((2, numSlices), dtype=int)
    start = 0
    for i in range(numSlices):
        end = start + sliceSize + (1 if i < remainder else 0)
        slices[:, i] = [start, end - 1]
        start = end
    return slices



