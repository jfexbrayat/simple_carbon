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
    pars        =   np.zeros(11,dtype='float64')
    pars[0]     =   LAI
    pars[1]     =   tmx
    pars[2]     =   tmn
    pars[3]     =   1.
    pars[4]     =   co2
    pars[5]     =   doy
    pars[6]     =   lat
    pars[7]     =   dswr
    pars[8]     =   -2.
    pars[9]     =   1.
    pars[10]    =   np.pi
    
    #define constants
    consts=np.array((ceff,0.0156935,4.22273,208.868,0.0453194,0.37836,7.19298, 0.011136,2.1001,0.789798),dtype='float64')

    gc=np.power(np.abs(pars[8]),consts[9])/(consts[5] * pars[9] + 0.5 * ( pars[1] - pars[2]))
    pp=pars[0]*pars[3]/gc*consts[0]*np.exp(consts[7]*pars[1])
    qq=consts[2]-consts[3]
    ci=0.5*(pars[4]+qq-pp+np.power(np.power(pars[4]+qq-pp,2)-4*(pars[4]*qq-pp*consts[2]),0.5))
    e0=consts[6]*np.power(pars[0],2)/(np.power(pars[0],2)+consts[8])
    dec=-23.4*np.cos((360.*(pars[5]+10.)/365.)*pars[10]/180.)*pars[10]/180.
    mult=np.tan(pars[6]*pars[10]/180)*np.tan(dec)

    if mult>=1 :
        dayl=24. 
    elif(mult<=-1):
        dayl=0.
    else:
        dayl=24.*np.arccos(-mult) / pars[10]

    cps=e0*pars[7]*gc*(pars[4]-ci)/(e0*pars[7]+gc*(pars[4]-ci))
    GPP=cps*(consts[1]*dayl+consts[4])
    return GPP


