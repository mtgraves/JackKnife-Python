import pylab as pl
import glob,argparse,sys
import averagingTools as aTools

def main():

    args = aTools.parseCMD()
   
    # Check if our data file exists, if not: write one.
    # Otherwise, open the file and plot.
    check = glob.glob('*JackKnifeData_Cv.dat*')
    fileNames = args.fileNames
    skip = args.skip

    print fileNames
    
    if check == []:
        
        temps,Cvs,CvsErr = pl.array([]),pl.array([]),pl.array([])
        Es, EsErr = pl.array([]), pl.array([])
   
        # open new data file, write headers
        fout = open('JackKnifeData_Cv.dat', 'w')
        fout.write('#%15s\t%16s\t%16s\t%16s\t%16s\n'% (
            'T', 'E', 'Eerr', 'Cv', 'CvErr'))
        
        # perform jackknife analysis of data, writing to disk
        if args.Crunched:   # check if we have combined data
            tempList = aTools.getHeadersFromFile(fileNames[0])
            print fileNames[0]
            n = 0
            for temp in tempList:
                print n
                temps = pl.append(temps, float(temp))
                E, EEcv, Ecv, dEdB = pl.loadtxt(fileNames[0],\
                        unpack=True,usecols=(n,n+1,n+2,n+3),delimiter=',')
                EAve, Eerr = aTools.jackknife(E[skip:])
                jkAve, jkErr = aTools.jackknife(
                        EEcv[skip:],Ecv[skip:],dEdB[skip:])
                print 'T = ',float(temp),':'
                print '<E>  = ',EAve,' +/- ',Eerr
                print '<Cv> = ',jkAve,' +/- ',jkErr
                Es      = pl.append(Es, EAve)
                Cvs     = pl.append(Cvs, jkAve)
                EsErr   = pl.append(EsErr, Eerr)
                CvsErr  = pl.append(CvsErr, jkErr)
                fout.write('%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\n' %(
                    float(temp), EAve, Eerr, jkAve, jkErr)) 
                n += 4
        else:       # otherwise just read in individual (g)ce-estimator files
            for fileName in fileNames:
                temp = float(fileName[13:19])
                temps = pl.append(temps, temp)
                E, EEcv, Ecv, dEdB = pl.loadtxt(fileName, unpack=True, 
                        usecols=(4,11,12,13))
                jkAve, jkErr = aTools.jackknife(
                        EEcv[skip:],Ecv[skip:],dEdB[skip:])
                EAve, Eerr = aTools.jackknife(E[skip:])
                print 'T = ',temp
                print '<Cv> = ',jkAve,' +/- ',jkErr
                print '<E>  = ',EAve,' +/- ',Eerr
                Es      = pl.append(Es, EAve)
                Cvs     = pl.append(Cvs, jkAve)
                EsErr   = pl.append(EsErr, Eerr)
                CvsErr  = pl.append(CvsErr, jkErr)
                fout.write('%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\n' %(
                    float(temp), EAve, Eerr, jkAve, jkErr)) 
        
        fout.close()

    else:
        print 'Found existing data file in CWD.'
        temps, Es, EsErr, Cvs,CvsErr = pl.loadtxt('JackKnifeData_Cv.dat', 
                unpack=True)
   
    errCheck = True
    if errCheck:
        EsNorm, EsErrNorm = pl.array([]), pl.array([])
        for fileName in args.fileNames:
            #Ecv,Eth = pl.loadtxt(fileName, unpack=True, usecols=(4,-5))
            Ecv = pl.loadtxt(fileName, unpack=True, usecols=(4,))
            EsNorm = pl.append(EsNorm,pl.average(Ecv))
            #ET = pl.append(ET, pl.average(Eth))
            EsErrNorm = pl.append(EsErrNorm, pl.std(Ecv)/pl.sqrt(float(len(Ecv))))
            #ETerr = pl.append(ETerr, pl.std(Eth)/pl.sqrt(float(len(Eth))))

        pl.scatter(temps, EsErrNorm, label='Standard Error', color='Navy')
        pl.scatter(temps, EsErr, label='Jackknife Error', color='Orange')
        pl.grid()
        pl.legend()
        pl.show()

    QHO = True
    if QHO:
        # analytical solutions for 1D QHO with one particle
        tempRange = pl.arange(0.01,1.0,0.01)
        Eanalytic = 0.5/pl.tanh(1.0/(2.0*tempRange))
        CvAnalytic = 1.0/(4.0*(tempRange*pl.sinh(1.0/(2.0*tempRange)))**2)

    ShareAxis=True      # shared x-axis for Cv and Energy
    # plot the specific heat vs. temperature
    if ShareAxis:
        ax1 = pl.subplot(211)
    else:
        pl.figure(1)
    if QHO: # plot analytic result
        pl.plot(tempRange,CvAnalytic, label='Exact')
    pl.errorbar(temps,Cvs,CvsErr, label='PIMC',color='Violet',fmt='o')
    if not ShareAxis:
        pl.xlabel('Temperature [K]')
    pl.ylabel('Specific Heat', fontsize=20)
    pl.grid(True)
    pl.legend(loc=2)
    
    # plot the energy vs. temperature
    if ShareAxis:
        pl.setp(ax1.get_xticklabels(), visible=False)
        ax2 = pl.subplot(212, sharex=ax1)
    else:
        pl.figure(2)
    if QHO: # plot analytic result
        pl.plot(tempRange,Eanalytic, label='Exact')
    pl.errorbar(temps,Es,EsErr, label='PIMC virial',color='Lime',fmt='o')
    pl.xlabel('Temperature [K]', fontsize=20)
    pl.ylabel('Energy [K]', fontsize=20)
    pl.grid(True)
    pl.legend(loc=2)

    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
