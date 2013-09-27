import numpy as np
import glob,argparse

def parseCMD():
    ''' Parse the command line. '''
    parser = argparse.ArgumentParser(description='pulls down lots of files.')
    parser.add_argument('fileNames', help='Data File Name.', nargs='+')
    parser.add_argument('-t', '--typeOfAverage', type=str,
            default='jackknife',
            help='NOT WORKING YET: Do you want jackknife or bootstrap?')
    parser.add_argument('-s', '--skip', type=int,
            default=1000,
            help='Number of bins to skip')

    return parser.parse_args()

def jackknife(data,data2=None,data3=None):
    ''' Return jackknife average (accounting for bias) and error.'''
    numBins = int(len(data))
    jkTerms = np.zeros(numBins)

    if data2==None:     # in case of single 
        for i in range(numBins):
            jkTerms[i] = np.mean(np.delete(data,i))
        dataAve = np.mean(data)
    else:               # in case of specific heat
        for i in range(numBins):
            jkTerms[i] = ( np.mean(np.delete(data,i))
                    - np.mean(np.delete(data2,i))**2 
                    - np.mean(np.delete(data3,i)) )
        dataAve = ( np.mean(data) - np.mean(data2)**2 
                - np.mean(data3) )
    
    jkAve = np.mean(jkTerms)
    ActAve = 1.0*numBins*dataAve - 1.0*(numBins-1)*jkAve
    jkVar = np.mean(jkTerms**2) - jkAve**2
    jkErr = np.sqrt((numBins-1)*jkVar)

    return ActAve, jkErr

# =============================================================================
