import numpy as np 
import astropy.units as u

def give_mape_k(res, map_t):
 
    lmap_x = map_t.shape[0]*res
    lmap_y = map_t.shape[1]*res

    map_kx = np.zeros(map_t.shape)
    map_ky = np.zeros(map_t.shape)
    map_k  = np.zeros(map_t.shape)

    for m in range(0,nx):
        if(m <= nx/2): m1 = m
        else: m1 = m - nx
        for n in range(0,ny):
            if(n <= ny/2): n1 = n
            else: n1 = m - nx
            
            kx = m1/lmap_x
            ky = n1/lmap_y
            
            map_kx[m,n] = kx
            map_ky[m,n] = ky
            map_k[m,n] = sqrt( kx^2 + ky^2)

    return map_k, map_kx, map_ky

    
