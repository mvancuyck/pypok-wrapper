import astropy.units as u 
import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.io import fits
from IPython import embed
from scipy import signal
from scipy import interpolate
import datetime
import matplotlib.pyplot as plt
#plt.style.use('classic')
import subprocess
from scipy.ndimage import gaussian_filter

from scipy.ndimage import gaussian_filter1d
from scipy.ndimage import generic_filter
import astropy.convolution as conv

 
def ipoker(data_map_large, res, k_out, beta,  map_b, map_k_binning, x_mtt_bb1, noise_pseudo_spec = np.zeros((1)), map_2_large = np.zeros((1)), remove_1st_bin = True  ): 

    
    #Compute the data and noise pseudo spectrum
    if(map_2_large.all() != 0): c_mn_t =  my_p2( data_map_large, map_2_large, res)
    else: c_mn_t =  my_p2( data_map_large, res) 

    pseudo_pk_mes, sigma_pb_mes, sigma_b = cmn2cb( map_b, map_k_binning * c_mn_t, return_sigmab1 = True, return_sigmab = True ) 

    if( noise_pseudo_spec.any() == 0): noise_pseudo_spec = np.zeros(pseudo_pk_mes.shape)

    #Discard DC bin to improve Mtt_bb's conditioning
    if(remove_1st_bin):
        pseudo_pk_mes     = pseudo_pk_mes[1:]
        noise_pseudo_spec = noise_pseudo_spec[1:]
        sigma_pb_mes      = sigma_pb_mes[1:] 

    #Binned k^beta*P(k) output power spectrum
    
    pb_out = np.absolute(np.dot(x_mtt_bb1, pseudo_pk_mes - noise_pseudo_spec))
    sigma_pk_out  = np.dot(x_mtt_bb1,sigma_pb_mes) 

    #Effective power spectrum at bin average assuming the process is described by a
    #power law of index=beta
    pk_out       = pb_out * k_out**(-beta)
    sigma_pk_out = sigma_pk_out * k_out**(-beta)

    return  pseudo_pk_mes, pk_out, sigma_pk_out



def poker_initialization(parfile, data_map, mask, res, scale, beta, a_beam, b_beam, delta_l_over_l = 0,  apod_length = 0,  kernel = "sincos" , no_mask = False, mask2 = np.ones((1)), keep_avg = False, remove_1st_bin = True, bypass_mtt_bb = False):

    #If no mask given:
    if(mask.any() != 1 and not no_mask):
        mask = make_mask(data_map.shape[0], data_map.shape[1], 3, 5 )

    #In case of cross correlation: the mask is the combination of the masks of the maps to cross correlate
    if(mask2.any() !=1): mask = mask * mask2
        
    #compute the apodized large mask
    data_map_large, mask_large = enlarges_masks_apodizes(data_map, mask, scale, apod_length, keep_avg = keep_avg)

    #compute the maximum value of k
    k_nyquist = np.pi / res.to(u.rad).value #rad**-1
    #compute the map of the radial k value
    k_map  =  give_map_k(res, data_map_large.shape[0] , data_map_large.shape[1]) #rad-1
    #Measuring from the map of k the miminum value of k
    kmin = np.min(k_map[np.nonzero(k_map)]) #rad-1 #pourquoi pas formule ?
    #Computing the minimum size of a k bin
    dk_min = 2 * 2 * np.pi / ( np.minimum(data_map.shape[0], data_map.shape[1]) * res.to(u.rad).value) #rad-1
    k_range = [kmin, k_nyquist] #rad**-1
    #Setting the binning of the k map
    k_bin_tab, k_bin_width = make_bintab(k_range, dk_min, delta_l_over_l) #rad-1
    map_k_binning = k_map**beta #pourquoiiiiii

    #Get mode mixing matrix
    map_k_binned, xbintab = igive_map_b( data_map_large.shape[0], data_map_large.shape[1], k_bin_tab/(2*np.pi/res.to(u.rad).value) ) #pourquoi reduced units 

    #If no mask nor beam is provided, save computation time                                               
    if(mask.any() != 1 and a_beam.any() != 1 and b_beam.any() != 1): bypass_mtt_bb = True

    if(bypass_mtt_bb == True):
        x_mtt_bb = np.identity( nbins)  #if no mask, mixing matrix = identity matrix
        map_binned = map_k_binned
    else:

        f = fits.PrimaryHDU(mask_large)
        hdu = fits.HDUList([f])
        hdu.writeto("mask.fits",  overwrite = True)

        f = fits.PrimaryHDU(a_beam)
        hdu = fits.HDUList([f])
        hdu.writeto("beam.fits",  overwrite = True)

        f = fits.PrimaryHDU(b_beam)
        hdu = fits.HDUList([f])
        hdu.writeto("beam1.fits",  overwrite = True)

        f = fits.PrimaryHDU(map_k_binning)
        hdu = fits.HDUList([f])
        hdu.writeto("map_k_binning.fits",  overwrite = True)

        f = fits.PrimaryHDU(k_bin_tab/(2*np.pi/res.to(u.rad).value)) #reduce units ###reduce units kezaco
        hdu = fits.HDUList([f])
        hdu.writeto("bintab.fits",  overwrite = True)

        #params2ascii( pars, parfile )
        #le .par est deja ecrit, je ferait la fonction quand j'aurais implementer un dictionnaire contenant tous les arguments.
        ####################################### 
        subprocess.run(["poker_mbb",parfile])
        print( "MBB COMPUTED. ")

        #######################################
        map_binned  = fits.getdata( "map_b.fits")  
        x_mtt_bb = fits.getdata( "poker_mtt_bb_out_x.fits") 

    #Compute bin addresses
    k_out = cmn2cb( map_binned, k_map )

    f = fits.PrimaryHDU(k_out) #reduce units ###reduce units kezaco
    hdu = fits.HDUList([f])
    hdu.writeto("k_out.fits",  overwrite = True)

    #Discard DC bin to improve Mtt_bb's conditioning
    if(remove_1st_bin):
        x_mtt_bb_m1 = np.linalg.inv( x_mtt_bb[1:,1:] )
        k_out       = k_out[1:]
    else: x_mtt_bb_m1 = np.linalg.inv( x_mtt_bb)

    return  mask_large,  k_nyquist, kmin,dk_min, k_map, k_bin_tab, k_bin_width, map_k_binning, map_binned, xbintab, x_mtt_bb_m1, k_out





def make_mask(ny, nx, nholes, radius_ratio, clean_border = True):  
    
    mask = np.ones((ny, nx))

    if(nholes >0):
        radius = nx / radius_ratio
        list_yc = np.random.randint(0, ny, nholes)
        list_xc = np.random.randint(0, nx, nholes)
        for y in range(0,ny):
            for x in range(0,nx):
                for h in range(0,nholes):
                    if( np.sqrt( (x-list_xc[h])**2 + (y-list_yc[h])**2 ) <= radius and 0 <= x <= nx-1 and 0 <= y <= ny-1): mask[y,x] = 0
    if(clean_border):
        mask[0,:] = 1
        mask[ny-1,:] = 1
        mask[:,0] = 1
        mask[:,nx-1] = 1

    return mask #,list_yc, list_xc



    
def give_map_k(res, ny, nx):
    res = res.to(u.rad)
    lmap_y = ny*res #rad
    lmap_x = nx*res #rad
    map_ky = np.zeros((ny, nx))
    map_kx = np.zeros((ny, nx))
    map_k  = np.zeros((ny, nx))

    for m in range(0,nx):
        if(m <= nx/2): m1 = m
        else: m1 = m - nx
        for n in range(0,ny):
            if(n <= ny/2): n1 = n
            else: n1 = n - ny
            
            kx = ( m1/lmap_x ).value
            ky = ( n1/lmap_y ).value
            
            map_kx[n,m] = kx
            map_ky[n,m] = ky
            map_k[n,m] = np.sqrt( kx**2 + ky**2)

    return map_k * 2 * np.pi #rad**-1
    

    
def make_bintab(l, delta_l_min, dll = 0, delta_l_max=0, ): 

    lmax      = l[1]
    lmin      = l[0]
    l1 = lmin
    delta_l = 0 

    bintab = []
    bintab.append(lmin)

    bin_width = []

    while(l1 + delta_l <= lmax):
        
   
        if(dll != 0): 
            delta_l = np.minimum( np.maximum(l1*dll, delta_l_min) , (lmax - l1) )
        else: 
            delta_l = delta_l_min

        if(delta_l_max != 0): 
            delta_l = np.minimum(delta_l, delta_l_max)
        l1 = l1 +  delta_l
        bintab.append(l1)
        bin_width.append(delta_l)

    bintab = np.asarray(bintab) 

    if( bintab.max() <= lmax):
        bintab[-1]    = lmax
        bin_width[-1] = bintab[-1]-bintab[-2]

    bintab[0] = bintab[0]*0.99
    bintab    = np.insert(bintab,0,0)

    return bintab, bin_width #rad-1


def igive_map_b( naxisY, naxisX, bintab, ):
    
    npix = naxisY * naxisX
    pix = np.round(np.linspace(0,npix-1, npix), 0).astype(int)

    #pix to mn
    n = np.floor(pix/naxisY)
    m = pix - n*naxisY

    #mn to k 
    kx = np.zeros(m.shape)
    ky = np.zeros(m.shape)
    k  = np.zeros(m.shape)


    for M in range(m.shape[0]):
        if(m[M] <= naxisX/2): kx[M] = m[M] / naxisX
        if(m[M] >= naxisX/2): kx[M] = (naxisX - m[M]) / naxisX

    for N in range(n.shape[0]):
        if(n[N] <= naxisY/2): ky[N] = n[N] / naxisY
        if(n[N] >= naxisY/2): ky[N] = (naxisY - n[N]) / naxisY

    k = np.sqrt( kx**2 + ky**2) #reduced units

    #k to bin

    b = np.ones(k.shape) * (-1)
    xbintab = np.zeros(bintab.shape)

    for i in range(xbintab.shape[0]-1):
        for M in range(b.shape[0]):
            if(k[M] >= bintab[i] and k[M] <= bintab[i+1] ):
                b[M]       = i
                xbintab[i] = 1

    b = b[:,np.newaxis]
    map_b = b.reshape( (naxisY, naxisX) )  
    

    return map_b, xbintab



def cmn2cb( map_b, map_cmn, return_sigmab1 = False, return_sigmab = False):

    map1 = map_b.copy()
    map2 = map_cmn.copy()

    map1[np.where( map_b == -32768.0 )] = 0  # -1 instead of -32768 in my igive_map_b
    map2[np.where( map_b == -32768.0 )] = 0  #change qqchose de 0 au lieu de prendre juste array sans undef?
                                             #si oui voir si on peut ignorer un pixel avec NaN
    cb = []
    sigma_b_1 = []  
    sigma_b = []
    
    hist, edges = np.histogram(map1)

    for i in range(0, int(map1.max())+1):
        index=np.where(map1 == i)
        cb.append(np.mean(map2[index]))
        if(return_sigmab1):
            if(index[0].shape[0] >= 1): sigma_b_1.append(np.std(map2[index]))
            else: sigma_b_1.append(np.NaN)
        if(return_sigmab): 
            if(index[0].shape[0] >= 1): sigma_b.append(np.std(map2[np.where(map1 == i)] / np.sqrt(index[0].shape[0])))
            else: sigma_b.append(np.NaN)
        
    cb = np.asarray(cb)
    sigma_b_1 = np.asarray(sigma_b_1)
    sigma_b = np.asarray(sigma_b)                                        
    
    if(return_sigmab1 and return_sigmab): return cb, sigma_b_1, sigma_b
    if(return_sigmab1):return cb, sigma_b_1
    if(return_sigmab): return cb, sigma_b
    if(not return_sigmab1 and not return_sigmab): return cb



def cls2map(l, clt, ny, nx, res,  l_cutoff=0,   ampl = 1,   fwhm_arcmin = 0, cu_t = np.zeros((1)),  index=None, force = True, cu_t_reset = True, ):
    
    lmap_rad_y = ny*res.to(u.rad).value
    lmap_rad_x = nx*res.to(u.rad).value
    #Generate gaussian amplitudes
    norm = np.sqrt( (nx/lmap_rad_x)*(ny/lmap_rad_y))

    noise = np.random.normal(loc=0, scale=1, size=(nx,ny))
    dmn1     = np.fft.fft2( noise )

    #Interpolate input power spectrum
    if(cu_t.any() == 0):

        k_map = give_map_k(res, ny, nx) 

        if(l_cutoff != 0): lmax = l_cutoff
        else:
            if(index is not None): lmax = k_map.max()
            else: lmax = l.max()

        if(fwhm_arcmin):
            sigma = fwhm_arcmin.to(u.rad) * gaussian_fwhm_to_sigma
            bl    = np.exp(-l*(l+1)*sigma**2)
            bk    = np.exp(-k_map*(k_map+1)*sigma**2)
        else:
            bl = np.ones(l.shape)
            bk = np.ones(k_map.shape)

        cu_t = np.zeros(k_map.shape)
        w = np.where( k_map <= lmax)#pourquoi pas kmin
        if(not w[0].any()): print("wrong k range")
        else:
            #Power law spectrum
            if(index != None): cu_t[w] = ampl * k_map[w]**index * bk[w]
            else:
                print("interpolate")
                f = interpolate.interp1d( l, clt*bl,  kind='linear')
                cu_t[w] = ampl * f(k_map[w])
                w1 = np.where( cu_t <= 0)
                if(w1[0].shape[0] != 0 and force): cu_t[w1] = 0

    #Fill amn_t
    amn_t = dmn1 * norm * np.sqrt( cu_t)
    
    #Output mapw
    map_t = np.real(np.fft.ifft2( amn_t ))

    return  map_t, cu_t



def my_p2(data_map, res,  map1 = None):

    npix = data_map.shape[0] * data_map.shape[1]

    if(map1 is not None):
        ft  = np.fft.fft2( data_map)
        ft1 = np.fft.fft2( map1)
        p2  = 0.5*( ft*np.conj(ft1) + np.conj(ft)*ft1) *  (res.to(u.rad).value)**2 / npix
    else: 
        p2 = np.absolute( (np.fft.fft2( data_map))**2 * (res.to(u.rad).value)**2 ) / npix

    return p2


def enlarges_masks_apodizes(Map, mask, scale,  apod_length = 0,  window = np.zeros((20,20)),  keep_avg = False, kernel = "sincos",  save_patch = True):

    ny_large = int( scale*Map.shape[0])
    nx_large = int( scale*Map.shape[1])

    iy0 = int(( (scale - 1)*Map.shape[0]/2))
    ix0 = int(( (scale - 1)*Map.shape[1]/2))
    
    if(window.any() == 0):
        return_window = True
        window = np.zeros((ny_large, nx_large))
        window[iy0:iy0+Map.shape[1], ix0:ix0+Map.shape[0]] = mask
        if(apod_length >= 0 ):
            print("before apod")
            window = window_apodization(window, apod_length)
    else: return_window = False

    if(keep_avg == False): 
        Map = Map -  Map[np.where(mask != 0)].mean()

    Map_large = np.zeros((ny_large, nx_large))
    Map_large[iy0:iy0+Map.shape[1], ix0:ix0+Map.shape[0]] = Map
    Map_large = Map_large * window

    patch = np.zeros((ny_large, nx_large))
    patch[iy0:iy0+Map.shape[1], ix0:ix0+Map.shape[0]] = 1
        
    if(save_patch):
        f = fits.PrimaryHDU(patch)
        hdu = fits.HDUList([f])
        hdu.writeto("patch.fits", overwrite = True)

    if(return_window): return Map_large, window
    else: return Map_large
        


        
def window_apodization(window, apod_length, kernel_type = "sincos" ):   

    kernelSize = 2*apod_length+1
    kernel = np.ones((kernelSize, kernelSize))


    if(kernel_type == "sincos"):
        y = np.linspace(0, kernelSize-1, kernelSize)
        y_axis= np.ones(y.shape)
        w = (y <= apod_length)
        y_axis[w] = y[w]/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*y[w]/apod_length)
        w = ( (kernelSize-1-y) <= apod_length)
        y_axis[w] = (kernelSize-1-y[w])/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*(kernelSize-1-y[w])/apod_length)
        for y in range(0, kernelSize):
            kernel[y,:] = kernel[y,:] * y_axis
        for x in range(0, kernelSize):
            kernel[:,x] = kernel[:,x] * y_axis



    if(kernel_type == "sin"):
        for y in range(0, kernelSize):
            kernel[y,:] = kernel[y,:] * (np.cos(np.linspace(0, kernelSize-1, kernelSize)/(kernelSize-1)*2*np.pi-np.pi)+1)/2
        for x in range(0, kernelSize):
            kernel[:,x] = kernel[:,x] * (np.cos(np.linspace(0, kernelSize-1, kernelSize)/(kernelSize-1)*2*np.pi-np.pi)+1)/2

    if(kernel_type == "han"):
        for y in range(0, kernelSize):
            kernel[y,:] = kernel[y,:] * signal.hann(kernelSize)
        for x in range(0, kernelSize):
            kernel[:,x] = kernel[:,x] * signal.hann(kernelSize) 

    if(kernel_type == "cos"):
        for y in range(0, kernelSize):
            kernel[y,:] = kernel[y,:] * signal.cosine(kernelSize)
        for x in range(0, kernelSize):
            kernel[:,x] = kernel[:,x] * signal.cosine(kernelSize) 
    
    sWindow = conv.convolve(window, kernel)
    innerWindow = np.zeros(sWindow.shape)
    innerWindow[np.where(sWindow >= 1 )] = 1

    final_window = conv.convolve(innerWindow, kernel)
    
    embed()

    plt.imshow(final_window, origin = "lower")
    plt.colorbar()
    plt.show()


    

    return final_window


def mc_reduce(tab_res, quick = False):
    
    nmc = tab_res.shape[0]
    nbins = tab_res.shape[1]
    tab_avg = np.mean(tab_res, axis = 0)
    sigma_avg = np.std(tab_res, axis = 0) / np.sqrt(nmc)

    if(quick == False):
        cov_mat = np.zeros((nbins,nbins))
        xcor = np.zeros((nbins, nbins))
        #Compute covariance matrix
        for b in range(0, nbins):
            for b1 in range(0, nbins):
                cov_mat[b,b1] = np.mean((tab_res[b,:] - tab_avg[b]) * (tab_res[b1,:]-tab_avg[b1]))
        #Cross correlation
        for b in range(0, nbins):
            for b1 in range(0, nbins):
                xcorr[b,b1] = cov_mat[b,b1]/np.sqrt(cov_mat[b,b]*cov_mat[b1,b1])
        return tab_avg, sigma_avg, cov_mat, xcorr

    else: return tab_avg, sigma_avg

