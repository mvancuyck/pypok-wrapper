import numpy as np
import astropy.unist as u
import astropy.io as fits 


def ipoker_arrays( data_map, mask, scale, res_pix, radius_ratio, apod_lenght, nholes, res_arcmin, a_beam, b_beam, pars, k_nyquist, bypass_mtt_bb = True, bintab = None, kmin = None,  clean_poker):

    mask_w8, patch = poker_make_mask(data_map.shape[0], data_map[1], scale, radius_ratio,  nholes, apod_lenght )
    
    indices = np.where(patch == 1) #add apodization to mask if required
    mask_w8[indices] =  mask * mask_weight[indices] ###################################

    #Binning parameters
    k_map, k_mapx, k_mapy = give_mape_k(res_pix, mask_w8) #######################################
    k_map = k_map * 2 * np.pi 
    wk = np.where(k_map != 0 )

    #Bins
    if(kmin is None): kmin = k_map[wk].min()

    if(bintab is None):
        k_range = [kmin, k_nyquist]
        bintab = make_bintab()
        bintab[0] = bintab[0]*0.99
        bintab    = [0, bintab]


    ###############################################
    # Init header
    '''header    = strarr(9)
    header[0] = "SIMPLE  =                    T / Written by IDL/ipoker.pro                      "
    header[1] = "BITPIX  =                  -64 / IEEE double precision floating point           "
    header[2] = "NAXIS   =                    2 /       /                                        "
    header[3] = "NAXIS1  =                  1   /Number of positions along axis 1                "    ; place holder
    header[4] = "NAXIS2  =                  1   /Number of positions along axis 2                "    ; place holder
    header[5] = "BLOCKED =                    T         /                                        "
    header[6] = "CDELT1  =        "+string( res_arcmin, "(F13.11)")+" /Resolution arcmin                               "
    header[7] = "CDELT2  =        "+string( res_arcmin, "(F13.11)")+" /Resolution arcmin                               "
    header[8] = "END        '''
    ####################
    

    #write patch
    #poker_writefits, pars.patch, patch, header, /silent #######################################

    #write Beam transfer function (no longer complex)
    #poker_writefits, pars.beam,  a_beam, /silent #######################################
    #poker_writefits, pars.beam1, b_beam, /silent #######################################


    # Derive binning laws
    map_k_binning = np.ones(k_map.shape)
    map_k_binning[wk] = k_map[wk]**beta 

    #write
    poker_writefits, pars.map_k_binning, map_k_binning, /silentb#######################################

    #Pass bintab to F90 (reduce units)
    poker_writefits, pars.input_bintab, bintab/(2.0d0*!dpi/pars.res_pix), /silent #######################################
    nbins = bintab.shape-1

    #Get mode mixing matrix
    igive_map_b, pars.nx_large, pars.ny_large, bintab/(2.0d0*!dpi/pars.res_pix), map_b_junk, xbintab #######################################

    if(bypass_mtt_bb == True):
        x_mtt_bb = identity( nbins)  #######################################
        map_b    = map_b_junk  #######################################
    else:
        parfile = "poker.par"
        params2ascii, pars, parfile  #######################################

    #######################################
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

    map_b    = fits.getdata( pars.file_map_b,/silent)  #######################################
    x_mtt_bb = readfits( pars.file_mtt_bb_x,/silent) #######################################


    return mask, mask_w8, patch, wp, k_map, bintab, a_beam, b_beam, map_k_binning, x_mtt_bb, x_mtt_bb_m1, map_b, k_out, pk_out, pseudo_pk, noise_speudo_pk, xbintab
