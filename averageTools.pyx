from __future__ import division

import numpy as np
cimport numpy as np

DTYPEi = np.int
DTYPEf = np.float64
ctypedef np.int_t DTYPE_i
ctypedef np.float64_t DTYPE_f

cimport cython
@cython.boundscheck(False)
def jackknife(np.ndarray[DTYPE_f,ndim=1] data ,
        np.ndarray[DTYPE_f,ndim=1] data2 = None,
        np.ndarray[DTYPE_f,ndim=1] data3 = None):
    ''' 
    Return jackknife average (accounting for bias) and error.
    '''
    
    cdef DTYPE_i numBins = int(len(data))
    cdef DTYPE_f dataAve
    
    cdef np.ndarray[DTYPE_f,ndim=1] jkTerms = np.zeros(numBins, dtype=float)
    #cdef np.ndarray[DTYPE_f,ndim=1] newsamp1
    #cdef np.ndarray[DTYPE_f,ndim=1] newsamp2
    #cdef np.ndarray[DTYPE_f,ndim=1] newsamp3

    if data2 == None:     # in case of single 
        for i in xrange(numBins):
            #newsamp1 = np.delete(data,i)
            jkTerms[i] = np.mean(np.delete(data,i))
        dataAve = np.mean(data)
    
    else:               # in case of specific heat
        for i in xrange(numBins):
            #newsamp1 = np.delete(data,i)
            #newsamp2 = np.delete(data2,i)
            #newsamp3 = np.delete(data3,i)
            jkTerms[i] = (np.mean(np.delete(data,i))
                    - np.mean(np.delete(data2,i))**2 
                    - np.mean(np.delete(data3,i)) )
            if i%1000==0:
                print i
        dataAve = ( np.mean(data) - np.mean(data2)**2 - np.mean(data3) )
    
    jkAve = np.mean(jkTerms)
    ActAve = (1.0*numBins*dataAve - 1.0*(numBins-1)*jkAve)
    jkVar = np.mean(jkTerms**2) - jkAve**2
    jkErr = np.sqrt((numBins-1)*jkVar)

    return ActAve, jkErr

# =============================================================================
