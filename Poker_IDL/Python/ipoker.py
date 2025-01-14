import astropy.units as u 
import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.io import fits
from IPython import embed
from scipy import signal
import datetime
import matplotlib.pyplot as plt
import subprocess
subprocess.run([fct,list de arg])
def run_ipoker(parfile, data_map, mask, res, scale, nx_large, ny_large, beta, delta_l_over_l, log_binning, remove_1st_bin, keep_avg, apod_lenght = 0, dkmin = 0, include_pix_window_function = True,  bypass_mtt_bb = False, dirr = None , a_beam = None, b_beam = None,  map1 = None, radius_ratio = None,   ): 


    k_nyquist = np.pi / res.to(u.rad) #rad**-1

    #Instrumental beam to include in mtt_bb
    if(a_beam is not None):
        s = a_beam.shape
        if(s[0] != nx_large or s[1] != ny_large):
            print("Incompatible sizes : a_beam and mask_large")
    
    else: a_beam = np.ones((nx_large, ny_large))


    if(b_beam is not None):
        s = b_beam.shape
        if (s[0] != nx_large or s[1] != ny_large):
            print("Incompatible sizes : b_beam and mask_large")
            
    else: b_beam = a_beam

    if(mask is None): mask = np.ones(data_map.shape)

    

    mask_w8, patch, wp, k_map, k_mapx, k_mapy, bintab, bin_width, map_k_binning, map_b, x_mtt_bb = ipoker_arrays(parfile, data_map, mask, res, scale, beta, nx_large, ny_large, k_nyquist, a_beam, b_beam,  dkmin, apod_lenght, delta_l_over_l, log_binning = log_binning, bypass_mtt_bb=bypass_mtt_bb )


    #Embed data in the large patch
    if(keep_avg  == False): Map = data_map.copy() -  data_map[np.where(mask != 0)].mean()
    else: Map = data_map.copy()
    map_large     = np.ones(mask_w8.shape)
    map_large[wp] = Map ##??

    if(map1 is not None):
        map1_large = np.zeros(mask_w8.shape)
        if(keep_avg == False): data_map1 = map1 - map1.mean()
        else: data_map1 = map1
        map1_large[wp] = data_map1 ###????
    else:
        map1_large = map_large


    #Compute the data(+noise) pseudo-spectrum
    pb_mes =  my_p2(map_large,  mask_w8.shape[0] * mask_w8.shape[1], res, map1) 

    #Discard DC bin to improve Mtt_bb's conditioning
    if(remove_1st_bin  == True):
        pb_mes            = pb_mes[1:]
        noise_pseudo_spec = noise_pseudo_spec[1:]
        sigma_pb_mes      = sigma_pb_mes[1:] ##????


    #Binned k^beta*P(k) output power spectrum
    nk = k_out.shape[0]
    pb_out = np.reshape(x_mtt_bb_m1)##(pb_mes-noise_pseudo_spec), nk)  #????
    sigma_pk_out  = np.reshape(x_mtt_bb_m1)##sigma_pb_mes) #???

    #Effective power spectrum at bin average assuming the process is described by a
    #power law of index=beta
    pk_out       = pb_out * k_out**(-beta)
    sigma_pk_out = sigma_pk_out * k_out**(-beta)


    return  k_out, pk_out, sigma_pk_out, k_nyquist





def ipoker_arrays(parfile, data_map, mask, res, scale, beta, nx_large, ny_large, k_nyquist, a_beam, b_beam, dk_min = 0,  apod_lenght = 0, delta_l_over_l = 0,  kmin = -1, log_binning = False, bypass_mtt_bb = False,  bintab = None,  clean_poker = False, bin_width = None, ):



    
    mask_w8, patch, ix0, iy0 = new_make_mask( data_map.shape[0], data_map.shape[1], scale, nx_large, ny_large, 20, 0, 20 )
    mask_w8[ix0:ix0+mask.shape[0], iy0:iy0+mask.shape[1]] = mask *  mask_w8[ix0:ix0+mask.shape[0], iy0:iy0+mask.shape[1]]

    #Binning parameters
    k_map, kx_map, ky_map = give_mape_k(res, mask_w8) #rad-1
    k_map = k_map * 2 * np.pi 

    #Bins
    if(kmin == -1): kmin = np.min(k_map[np.nonzero(k_map)]) #rad-1

    if(dk_min <= 0):dk_min = 2 * 2 * np.pi / ( data_map.shape[0] * res.to(u.rad)) #rad-1


    if(bintab is None):
        k_range = [kmin, k_nyquist.value] #rad**-1
        bintab, bin_width = make_bintab(k_range, dk_min.value, delta_l_over_l, log_binning=log_binning)
        bintab[0] = bintab[0]*0.99
        bintab    = np.insert(bintab,0,0)



    #Init header

    hdr = {} 

    #CRPIX1  =      507.90474039695 / Pixel coordinate of reference point            
    #CRPIX2  =      503.02270076207 / Pixel coordinate of reference point               
    hdr["CDELT1"]  =   (res.to(u.deg)).value
    hdr["CDELT2"]  =   (res.to(u.deg)).value  
    hdr["CUNIT1"]  = 'deg'              
    hdr["CUNIT2"]  = 'deg'             
    hdr["CTYPE1"]  = 'RA---TAN'          
    hdr["CTYPE2"]  = 'DEC--TAN'          
    #CRVAL1  =     0.70395330810396 / [deg] Coordinate value at reference point      
    #CRVAL2  =     0.69721933449752 / [deg] Coordinate value at reference point           
    #LONPOLE =                180.0 / [deg] Native longitude of celestial pole       
    #LATPOLE =     0.69721933449752 / [deg] Native latitude of celestial pole        
    #DATEREF = '1858-11-17'         / ISO-8601 fiducial time                         
    #MJDREFI =                  0.0 / [d] MJD of fiducial time, integer part         
    #MJDREFF =                  0.0 / [d] MJD of fiducial time, fractional part      
    #RADESYS = 'ICRS'               / Equatorial coordinate system                  
    hdr["DATE"]    =  str(datetime.datetime.now())             
    hdr["COMMENT"] = "Ipoker"                 

    hdu = fits.ImageHDU(data_map,header = hdr)

    hdu = fits.HDUList([hdu])
    hdu.writeto("")
    
    #write patch
    poker_writefits("../patch.fits", patch, hdr)

    #write Beam transfer function (no longer complex)
    #poker_writefits("../beamA.fits",  a_beam, hdr)
    #poker_writefits("../beamB.fits",  b_beam, hdr)

    # Derive binning laws
    map_k_binning = k_map**beta 
    #poker_writefits("../map_k_binning.fits", map_k_binning, hdr)

    #Pass bintab to F90 (reduce units)

    """f = fits.PrimaryHDU(bintab)
    hdu = fits.HDUList([f])
    hdr = hdu[0].header
    poker_writefits("../input_bintab.fits", bintab/(2*np.pi/res), hdr)"""

    #Get mode mixing matrix
    #map_b_junk, xbintab = igive_map_b( nx_large, ny_large, bintab/(2*np.pi/res) ) 

    
    if(bypass_mtt_bb == True):
        x_mtt_bb = np.identity( nbins)  
        map_b    = map_b_junk 
    else:
        params2ascii( pars, parfile )

        ####################################### srun -c nbr coeur poker..
        '''
        if (!arch eq "franklin") or (!arch eq "parallel") then begin
        print, ""
        print, "=================="
        print, "ready to run poker_count_task and qsub run_poker.txt."
        print, "Then press .c to continue."
        print, "=================="
        stop
        endif else begin
        spawn, "poker_mbb "+parfile
        print, "" & print, "" & print, "MBB COMPUTED."
        endelse
        '''
        #######################################
        """
        map_b    = fits.getdata( pars.file_map_b)  #######################################
        x_mtt_bb = fits.getdata( pars.file_mtt_bb_x) #######################################
        """

    '''
    #Compute bin addresses
    cmn2cb, map_b, k_map, k_out
    
    poker_writefits, pars.k_out, k_out, /silent

    #Discard DC bin to improve Mtt_bb's conditioning
    if(remove_1st_bin):
        x_mtt_bb_m1 = invert( x_mtt_bb[1:*,1:*])
        k_out       = k_out[1:*]
    else: x_mtt_bb_m1 = invert( x_mtt_bb)

    out_arrays = {mask:mask, w8:w8, patch:patch, wp:wp, k_map:k_map, $
              bintab:bintab, a_beam:a_beam, b_beam:b_beam, $
              map_k_binning:map_k_binning, $
              x_mtt_bb:x_mtt_bb, $
              x_mtt_bb_m1:x_mtt_bb_m1, map_b:map_b, $
              k_out:k_out, pk_out:k_out*0.0d0, pseudo_pk:dblarr(nbins), noise_pseudo_spec:dblarr(nbins), $
              xbintab:xbintab}

    if(clean_poker): poker_clean, out_params'''

    return  mask_w8, patch, indices_k, k_map, k_mapx, k_mapy, bintab, bin_width, map_k_binning, map_b, x_mtt_bb



def new_make_mask(nx, ny, scale, nx_large, ny_large,  radius_ratio,  nholes = 0, apod_length = 0, clean_border = False ):

    mask = np.zeros((nx_large,ny_large))
    patch = mask.copy()
    holes = np.ones((nx_large,ny_large))

    if(scale >= 1): 
        ix0 = int(( (scale - 1)*nx/2))
        iy0 = int(( (scale - 1)*ny/2))
    else:
        ix0 = 0
        iy0 = 0

    mask[ ix0:ix0+nx-1, iy0:iy0+ny-1] = 1
    patch[ix0:ix0+nx-1, iy0:iy0+ny-1] = 1

    if(nholes >0):
        radius = nx / radius_ratio
        list_xc = np.random.randint(0, nx, nholes)
        list_yc = np.random.randint(0, ny, nholes)
        for x in range(0,nx):
            for y in range(0,ny):
                for h in range(0,nholes):
                    if( np.sqrt( (x-list_xc[h])**2 + (y-list_yc[h])**2 ) <= radius):
                        mask[ ix0+x,iy0+y] = 0
                        holes[ ix0+x,iy0+y] = 0

    if(clean_border == True):
        mask[ix0,         iy0:iy0+ny-1] = 1
        mask[ix0+nx-1,    iy0:iy0+ny-1] = 1
        mask[ix0:ix0+nx-1, iy0] = 1
        mask[ix0:ix0+nx, iy0+ny-1] = 1

    
    if(apod_length != 0):
        
        cos_axis = signal.cosine(2 * apod_length)
        cos_axis_half = int((cos_axis.shape[0] /  2) )
        x_axis = np.ones((nx))
        x_axis[:cos_axis_half] = cos_axis[:cos_axis_half]
        x_axis[int(-1* cos_axis_half):] = cos_axis[int(-1*cos_axis_half):]
        y_axis = np.ones((ny))
        y_axis[:cos_axis_half] = cos_axis[:cos_axis_half]
        y_axis[int(-1*cos_axis_half):] = cos_axis[int(-1*cos_axis_half):]
        taper = np.dot(x_axis[:,np.newaxis], y_axis[np.newaxis,:]) 
        mask[ix0:ix0+nx, iy0:iy0+ny] = mask[ix0:ix0+nx, iy0:iy0+ny]*taper

    #apodizer les trous ? 
       
    return mask, patch, ix0, iy0
    
        






def give_mape_k(res, map_t):

    res = res.to(u.rad)
 
    nx = map_t.shape[0]
    ny = map_t.shape[1]
    lmap_x = map_t.shape[0]*res #rad
    lmap_y = map_t.shape[1]*res #rad
    map_kx = np.zeros(map_t.shape)
    map_ky = np.zeros(map_t.shape)
    map_k  = np.zeros(map_t.shape)

    for m in range(0,nx):
        if(m <= nx/2): m1 = m
        else: m1 = m - nx
        for n in range(0,ny):
            if(n <= ny/2): n1 = n
            else: n1 = n - ny
            
            kx = ( m1/lmap_x ).value
            ky = ( n1/lmap_y ).value
            
            map_kx[m,n] = kx
            map_ky[m,n] = ky
            map_k[m,n] = np.sqrt( kx**2 + ky**2)

    return map_k, map_kx, map_ky #rad**-1








    
def make_bintab(l, delta_l_min, dll, delta_l_max=0, log_binning=False, ): 

    lmax      = l[1]
    lmin      = l[0]
    l1 = lmin
    delta_l = 0 

    bintab = []
    bintab.append(lmin)

    bin_width = []

    while(l1 + delta_l <= lmax):
        
   
        if(log_binning ): 
            delta_l = np.minimum( np.maximum(l1*dll, delta_l_min) , (lmax - l1) )
            print("here1")
        else: 
            delta_l = delta_l_min
            print("here2")

        if(delta_l_max != 0): 
            delta_l = np.minimum(delta_l, delta_l_max)
            print("here3")
        print("here4")
        l1 = l1 +  delta_l
        bintab.append(l1)
        bin_width.append(delta_l)

    bintab = np.asarray(bintab)

    if( bintab.max() <= lmax):
        n = bintab.shape[0]
        bintab[n-1]    = lmax
        bin_width[n-1] = bintab[n-1]-bintab[n-2]


    return bintab, bin_width












def my_p2(data_map, npix, res,  map1):
    
    if(map1 is not None):
        ft  = np.fft.fft2( data_map)
        ft1 = np.fft.fft2( map1)
        p2  = 0.5*( ft*np.conj(ft1) + np.conj(ft)*ft1) * npix * res.to(u.rad)**2
    else: p2 = abs( np.fft.fft2( data_map)**2 * npix * res.to(u.rad)**2 )
    return p2






def poker_writefits(filename, data, Hdr, heap=None, append = None, compress = None, CheckSum = None, nanValue= None, silent = None):

        

    f = fits.PrimaryHDU(data, header=Hdr)
    hdu = fits.HDUList([f])
    hdu.writeto(filename, overwrite=True)




def igive_map_b( naxis1, naxis2, bintab, map_b,):
    
    npix = naxis1 * naxis2
    nbins = bintab.shape[0] - 1
    pix = np.linspace(0,npix-1, npix)
    map_b = np.zeros((naxis1, naxis2))

    #pix to mn
    n = pix/naxis1
    m = pix - n*naxis1

    #mn to k 
    kx = np.zeros(m.shape)
    ky = np.zeros(m.shape)
    k  = np.zeros(m.shape)


    if(m.any() <= nx/2):
        w = np.where( m <= nx/2 ) 
        kx[w] = m[w]/nx

    if(m.any() >= nx/2):
        w = np.where( m >= nx/2 ) 
        kx[w] = nx-m[w]/nx

    if(n.any() <= ny/2):
        w = np.where( n <= nx/2 ) 
        ky[w] = n[w]/ny
                
    if(n.any() >= ny/2):
        w = np.where( n >= ny/2 ) 
        kx[w] = ny-m[w]/ny


    k = sqrt( kx^2 + ky^2) #reduced units

    #k to bin

    b = np.zeros(k.shape)
    xbintab = np.zeros(bintab.shape)

    for i in range(O, nbins):
        if(k.any() <= bintab[i] and k.any() >= bintab[i+1] ):
                          w = np.where( k <= bintab[i] )
                          w = np.where(k[w] >= bintab[i+1])
                          b[w]       = i
                          xbintab[i] = 1
            

    return b,xbintab"""






def params2ascii(params, parfile):
    
    File = open(parfile, w)

    File.write(params)
    
    File.close()






"""def cmn2cb( map_b, map_cmn, cb, sigma_b, sigma_b_1):

    w  = where( map_b != None)
    m1 = map_b[w]
    m2 = map_cmn[w]

    h         = histogram( m1, bin=1.0d0, reverse_ind=R)
    cb        = double(h*0.0d0)
    sigma_b   = double(h*0.0d0)
    sigma_b_1 = double(h*0.0d0)
    for i in range(0, n_elements(h)-1):

        if(R[i+1]-R[i] >= 2):
            cb[i]        = avg(    m2[ R[ R[i] : R[i+1]-1]])
            sigma_b_1[i] = stddev( m2[ R[ R[i] : R[i+1]-1]]) # faire un mask? 
            sigma_b[i]   = stddev( m2[ R[ R[i] : R[i+1]-1]])/sqrt((R[i+1]-R[i]))

        if(R[i+1]-R[i] == 1):
            cb[i]        = m2[ R[ R[i]]]
            sigma_b_1[i] = !values.f_nan #si 1 element  ==> sigma = 0
            sigma_b[i]   = !values.f_nan
            
    return"""
