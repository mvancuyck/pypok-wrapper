import numpy as np
from astropy.io import fits
import astropy.units as u
from IPython import embed 
from astropy import wcs
import scipy.constants as cst
import powspec 
from pathlib import Path
from matplotlib import gridspec
from ipoker import *
from progress.bar import Bar
from matplotlib.gridspec import GridSpec
import matplotlib

def poker_main():
    
    parser = argparse.ArgumentParser(description="poker main",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    #options

    parser.add_argument('params_poker', help=".par file with poker params", default = None)

    #parser.add_argument('nnodes', help = "number of nods")
    parser.add_argument('--load_directly', help = "load directly Mbb", action="store_true")

    parser.add_argument('--non_iteractive', help = "when compute Mbb, deactivate matplotlib", action="store_true")

    args = parser.parse_args()

    if(args.non_iteractive): matplotlib.use("Agg")

    #Load user parameters for Poker
    P = load_params(args.params_poker)
    P = set_pars(P, load_mbb_directly = args.load_directly)

    print('here')
    embed()

    #Charging the data map
    data_map = fits.getdata(f"{P['map_path']}/{P['map_name']}")
    
    #Charging the multipols l and the dust and noise angular power spectra
    l         = fits.getdata(f"{P['l_file']}") 
    cl_signal = fits.getdata(f"{P['cl_signal_file']}")
    cl_noise = fits.getdata(P['cl_noise_file'])

    data_map_large = apodizing_a_map(data_map, P['mask_array'], patch = P["patch"], iy0=P["iy0"], ix0=P["ix0"], keep_avg = P["keep_avg"])
    
    _, signal_pk_out = ipoker(data_map_large, P['res'], P['beta'], P['l'],P['l_bins'], P['map_l_power_beta'],P['l_map'],  P['m_bb_m1'])

    #------------ Noise ----------------

    # if the data are noisy, computes an estimation of the mean noise pseudo power spectrum <N^> (Eq. 17)
    
    if(P['include_noise']):
        
        #Charging the noise map
        measured_noise = fits.getdata(P['noise_map']) #noise of the data
                
        total_data_map = data_map + measured_noise
        
        pseudo_pk_noise_list = []

        #init cl_map_noise
        noise_map_large, cl_map_noise = gaussian_random_field(l, cl_noise, P["ny_large"], P["nx_large"], P["res"], l_cutoff= P["l_max"],  sigma = P['sigma_beam'])
        
        bar = Bar('Monte Carlo simulations for noise', max=P['nmc'])

        for imc in range(0, P['nmc']):
            
            bar.next()

            noise_map_large, cl_map_noise = gaussian_random_field(l, cl_noise, P["ny_large"], P["nx_large"], P["res"], l_cutoff= P["l_max"],  sigma = P['sigma_beam'], cl_map = cl_map_noise)
            
            noise_map_large = noise_map_large * P["patch"]
            
            pseudo_pk, _ = ipoker(noise_map_large, P['res'], P['beta'], P['l'],P['l_bins'], P['map_l_power_beta'], P['l_map'],  P['m_bb_m1'])
            
            pseudo_pk_noise_list.append( pseudo_pk )
        
        bar.finish
 
        pseudo_pk_noise_list = np.asarray(pseudo_pk_noise_list)
        noise_pseudo_pk, sigma_avg = mc_reduce(pseudo_pk_noise_list, quick = True)

    else:
    
        total_data_map = data_map.copy()
        cl_map_noise =np.zeros( data_map_large.shape )
        noise_pseudo_pk = np.zeros(len(signal_pk_out) )

    
    #---------------------------------- Data power spectrum ------------------------------------------
    
    total_data_map_large = apodizing_a_map(total_data_map,P['mask_array'], patch = P["patch"], iy0=P["iy0"], ix0=P["ix0"],  keep_avg = P["keep_avg"])
    data_pseudo_pk, data_pk_out = ipoker(total_data_map_large, P['res'], P['beta'], P['l'],  P['l_bins'], P['map_l_power_beta'],P['l_map'], P['m_bb_m1'])

    #-------- Error bars ----------------------
    
    # Init signal cl_map_signal
    Map, cl_map_signal = gaussian_random_field(l, cl_signal, P["ny_large"], P["nx_large"], P["res"], l_cutoff= P["l_max"],  sigma = P['sigma_beam'])
    
    pk_out_list = []

    bar = Bar('Monte Carlo simulation for the signal', max=P['nmc'])
    for imc in range(0, P['nmc']):
        bar.next()
        #Signal
        map_t, cl_map_signal = gaussian_random_field(l, cl_signal, P["ny_large"], P["nx_large"], P["res"], l_cutoff= P["l_max"], ny_exit = P["ny"], nx_exit = P["nx"], sigma = P['sigma_beam'], cl_map = cl_map_signal)
        if(P["include_noise"]):
            #Add white noise

            noise, cl_map_noise =  gaussian_random_field(l, cl_noise, P["ny_large"], P["nx_large"], P["res"], l_cutoff= P["l_max"], ny_exit = P["ny"], nx_exit = P["nx"], sigma = P['sigma_beam'], cl_map = cl_map_noise)
            map_t = map_t + noise

        map_t = apodizing_a_map(map_t, P['mask_array'], patch = P["patch"], iy0=P["iy0"], ix0=P["ix0"],  keep_avg = P["keep_avg"])

        #Power spectrum
        
        _, pk_out = ipoker(map_t, P['res'], P['beta'], P['l'], P['l_bins'], P['map_l_power_beta'],P['l_map'], P['m_bb_m1'], noise_pseudo_pk = noise_pseudo_pk)

        pk_out_list.append( pk_out )

    pk_out_list = np.asarray(pk_out_list) 
    
    pk_final, sigma_pk_final, cov_mat, corr_mat = mc_reduce(pk_out_list)

    bar.finish
    print('')
    
    #---Plot---

    print("end")
    
    if(P['include_noise']): pk_powspec, k_edges = powspec.power_spectral_density(data_map+measured_noise, res=P['res'], bins = P["l_bins"] / (2*np.pi))
    else: pk_powspec, k_edges = powspec.power_spectral_density(data_map, res=P['res'], bins = P["l_bins"] / (2*np.pi))
    pk_powspec = pk_powspec
    f = interpolate.interp1d( l, cl_signal,  kind='linear')
    cl_out = f(P['l'])
    
    L = P['l']
    L_signal = l
    
    # Configure Figure A
    size = 7; lw = 2; mk=4
    plt.rc('font', size=size)
    plt.rc('axes', titlesize=size)
    plt.rc('axes', labelsize=size)

    fig = plt.figure(figsize=(7, 3), dpi=200)
    gs = fig.add_gridspec(2, 2,hspace=0 )

    # Plot in the first column, first row
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlim(L.min() - 100, L.max() + 100)
    ax1.errorbar(L, pk_final, yerr=sigma_pk_final, fmt='o', color="green", ecolor="green", label='POKER', markersize=mk)
    ax1.plot(L_signal, cl_signal, label="Input", lw=lw)
    ax1.plot(L_signal, (cl_signal + cl_noise), label="Total", lw=lw)
    ax1.plot(L, pk_powspec, ls='--', color="grey", label="Naif", lw=lw)
    ax1.set_yscale("log")
    ax1.set_ylabel("P(k) [$ \\rm Jy^2/sr$]")
    ax1.set_ylim(1e-12, 1e-5)
    ax1.set_xscale("log")
    ax1.legend(loc='upper right')
    ax1.tick_params(axis="both", which='minor', tickdir="inout", top=True, length=4, color='k')

    # Plot in the first column, second row
    ax = fig.add_subplot(gs[1, 0], sharex=ax1)
    line_e = ax.errorbar(L,  pk_final/cl_out-1, yerr= sigma_pk_final/cl_out  , fmt='o', color = "green", ecolor = "green", markersize=mk)
    line_f = ax.plot((L.min() - 100,L.max()+100),(0,0), color = "gray")
    ax.set_xscale("log")
    ax.set_xlabel("k [$\\rm arcmin^{-1}$]")
    ax.set_ylabel("relative error")
    ax.set_xlim(L.min() - 100, L.max()+100)
    ax.xaxis.grid(True, which='minor')
    ax.tick_params(axis = "both", which='minor', tickdir = "inout", top = True,  length=4,color='k')

    # Plot in the second column (only one plot)
    ax2 = fig.add_subplot(gs[:, 1])
    title = "Bin Bin correlation matrix"
    ax2 = fig.add_subplot(122, projection='3d')
    _x = np.arange(corr_mat.shape[0])
    _y = np.arange(corr_mat.shape[1])
    _xx, _yy = np.meshgrid(_x, _y)
    x, y = _xx.ravel(), _yy.ravel()
    top = corr_mat.reshape(corr_mat.shape[0] * corr_mat.shape[1]) 
    bottom = np.zeros(top.shape)
    width = depth = 1
    ax2.bar3d(x, y, bottom, width, depth, top, shade=True, color='cyan')
    ax2.view_init(50, 20)
    ax2.xaxis.pane.set_edgecolor('w')
    ax2.yaxis.pane.set_edgecolor('w')
    ax2.zaxis.pane.set_edgecolor('w')
    
    plt.tight_layout()
    plt.savefig("poker.pdf")
    plt.show()
    
if __name__ == '__main__':

    poker_main()
