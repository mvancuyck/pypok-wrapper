import astropy.units as u
import numpy as np

def cmn2cb( map_b, map_cmn, cb, sigma_b, sigma_b_1):

    w  = where( map_b != None)
    m1 = map_b[w]
    m2 = map_cmn[w]

    h         = histogram( m1, bin=1.0d0, reverse_ind=R)
    cb        = double(h*0.0d0)
    sigma_b   = double(h*0.0d0)
    sigma_b_1 = double(h*0.0d0)
    for i in range(0, n_elements(h)-1):

        if(R[i+1]-R[i] == 2):
            cb[i]        = avg(    m2[ R[ R[i] : R[i+1]-1]])
            sigma_b_1[i] = stddev( m2[ R[ R[i] : R[i+1]-1]])
            sigma_b[i]   = stddev( m2[ R[ R[i] : R[i+1]-1]])/sqrt((R[i+1]-R[i]))

        if(R[i+1]-R[i] == 1):
            cb[i]        = m2[ R[ R[i]]]
            sigma_b_1[i] = !values.f_nan
            sigma_b[i]   = !values.f_nan
            
    return
