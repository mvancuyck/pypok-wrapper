import numpy as np
import astropy.unist as u
import math 

def poker_make_mask(nx, ny, scale, radius_ratio,  nholes = 0, apod_lenght = 0, clean_border = True ):
    
    mask = np.zeros(( int(math.ceil(nx*scale)), int(math.ceil(ny*scale))))
    patch = mask.copy()
    holes = np.ones(( int(math.ceil(nx*scale)), int(math.ceil(ny*scale))))

    if scale_x gt 1 then ix0 = long( (scale_x - 1)*nx/2.) else ix0 = 0
    if scale_y gt 1 then iy0 = long( (scale_y - 1)*ny/2.) else iy0 = 0

    x0 = (scale - 1)*nx/2.)
    y0 = (scale - 1)*ny/2.)


    mask[ ix0:ix0+nx-1, iy0:iy0+ny-1] = 1
    patch[ix0:ix0+nx-1, iy0:iy0+ny-1] = 1

    if(nholes >0):
        radius = nx / radius_ratio
        list_xc = np.random.randint(0, nx, nholes)
        list_yc = np.random_integers(0, ny, nholes)
        for x,y in zip(range(0,nx), range(0,ny)):
            for h in range(0,nholes):
                if( np.sqrt( (x-xc[h])**2 + (y-yc[h])**2 ) <= radius):
                    mask[  ix0+i,iy0+j] = 0
                    holes[ ix0+i,iy0+j] = 0

    if(clean_border == True):
        mask[ix0,      iy0:iy0+ny-1] = 1
        mask[ix0+nx-1, iy0:iy0+ny-1] = 1
        mask[ix0:ix0+nx-1,      iy0] = 1
        mask[ix0:ix0+nx-1, iy0+ny-1] = 1


    if(apod_lenght != 0):

        x_axis = np.linspace(0,nx, nx+1).astype(int)
        y_axis = np.ones(x_axis.shape)
        w = np.where( x_axis <= apod_length)
        y_axis[w] = x_axis[w]/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*x_axis[w]/apod_length)
        w = np.where( (nx-1-x) <= apod_length)
        y_axis[w] = (nx-1-x_axis[w])/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*(nx-1-x_axis[w])/apod_length)
   

        x_axis = np.linspace(0,nx, nx+1).astype(int)
        y_axis_1 = np.ones(x_axis.shape)
        w = np.where( x_axis <= apod_length)
        y_axis_1[w] = x_axis[w]/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*x_axis[w]/apod_length)
        w = np.where( (nx-1-x) <= apod_length)
        y_axis_1[w] = (nx-1-x_axis[w])/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*(nx-1-x_axis[w])/apod_length)

        taper = y#y1   ###################

        mask[ix0:ix0+nx-1, iy0:iy0+ny-1] = mask[ix0:ix0+nx-1, iy0:iy0+ny-1]*taper

    return mask, patch
