"""
23/04/2018 - JFE
This file holds the description of a simple carbon model inspired by DALEC.
The simple model simulates GPP using ACM as a function of LAI, air temperature 
and radiation. A fixed fraction of GPP is allocated to autotrophic respiration 
and the remainder is allocated to a single plant pool. Plant turnover releases
carbon in a single dead organic matter pool in which microbial decomposition
leads to heterotrophic respiration
"""

import numpy as np

def ACM(lat,LAI,tmn,tmx,swrd,co2,doy,ceff):
    """
    The aggregated canopy model is an emulator of the multi-layer SPA model
    which calculates canopy-level daily GPP in g C m-2 d-1as a function of 
    LAI and weather variables under well-watered conditions. Parameters are:
    - lat               latitude                        [deg N]
    - LAI               Leaf Area Index                 [m2 m-2]
    - tmn               Daily minimum air temperature   [deg C]
    - tmx               Daily maximum air temperature   [deg C]
    - dswr              Downward short-wave radiation   [MJ m-2 d-1]
    - co2               Atmospheric CO2 concentration   [ppmv]
    - doy               Day of Year                     [-]
    - ceff              Carbon efficiency parameter     [-]
    """

    #define an array with parameters
    acmpars        =   np.zeros(11,dtype='float64')
    acmpars[0]     =   LAI
    acmpars[1]     =   tmx
    acmpars[2]     =   tmn
    acmpars[3]     =   1.
    acmpars[4]     =   co2
    acmpars[5]     =   doy
    acmpars[6]     =   lat
    acmpars[7]     =   dswr
    acmpars[8]     =   -2.
    acmpars[9]     =   1.
    acmpars[10]    =   np.pi
    
    #define constants
    consts=np.array((ceff,0.0156935,4.22273,208.868,0.0453194,0.37836,7.19298, 0.011136,2.1001,0.789798),dtype='float64')

    gc=np.power(np.abs(acmpars[8]),consts[9])/(consts[5] * acmpars[9] + 0.5 * ( acmpars[1] - acmpars[2]))
    pp=acmpars[0]*acmpars[3]/gc*consts[0]*np.exp(consts[7]*acmpars[1])
    qq=consts[2]-consts[3]
    ci=0.5*(acmpars[4]+qq-pp+np.power(np.power(acmpars[4]+qq-pp,2)-4*(acmpars[4]*qq-pp*consts[2]),0.5))
    e0=consts[6]*np.power(acmpars[0],2)/(np.power(acmpars[0],2)+consts[8])
    dec=-23.4*np.cos((360.*(acmpars[5]+10.)/365.)*acmpars[10]/180.)*acmpars[10]/180.
    mult=np.tan(acmpars[6]*acmpars[10]/180)*np.tan(dec)

    if mult>=1 :
        dayl=24. 
    elif(mult<=-1):
        dayl=0.
    else:
        dayl=24.*np.arccos(-mult) / acmpars[10]

    cps=e0*acmpars[7]*gc*(acmpars[4]-ci)/(e0*acmpars[7]+gc*(acmpars[4]-ci))
    GPP=cps*(consts[1]*dayl+consts[4])
    return GPP

def simple_carbon(lat,LAI,metdrivers,pars,cveg_0 = 1.,csom_0 = 1.):
    """
    This is the simple carbon model running at a location at latitude lat
    """

    #define the number of time steps
    nsteps = LAI.size
    # create arrays to store fluxes and pools at end of time step 
    gpp     = np.zeros(n,dtype='float64')
    ra      = np.zeros(n,dtype='float64')
    rh      = np.zeros(n,dtype='float64')

    cveg    = np.zeros(n+1, dtype = 'float64')     
    csom    = np.zeros(n+1, dtype = 'float64')     
    #store initial conditions
    cveg[0] = cveg_0
    csom[0] = csom_0

    #assign dummy variables to parameters
    ceff, frau, kveg, ksom, q10 = pars



    #loop over nsteps
    for ii in range(nsteps):
        #create dummy variables for readability        
        tmn = metdrivers[:,1]
        tmx = metdrivers[:,2]
        dswr= metdrivers[:,3]
        co2 = metdrivers[:,4]
        doy = metdrivers[:,5]

        gpp[ii]     = ACM(lat,LAI[ii],tmn,tmx,dswr,co2,doy,ceff)
        ra[ii]      = gpp[ii]*frau

        #calculate temperature response function assuming a ref temperature of 15C
        tavg    = (tmx+tmn)/2.
        ft      = q10**((tavg-15.)/10)

        # calculate turnovers
        cveg_turnover = cveg[ii]*kveg*ft
        rh[ii] = csom[ii]*ksom*ft

        # apply fluxes and update pools
        cveg[ii+1]    = cveg[ii]+gpp[ii]-ra[ii]-cveg_turnover
        csom[ii+1]    = csom[ii]+cveg_turnover-rh[ii]

    # returns time series of land atmosphere fluxes and pools
    return gpp,ra,rh,cveg,csom


        




