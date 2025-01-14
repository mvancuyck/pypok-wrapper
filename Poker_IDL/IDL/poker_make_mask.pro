;; Copyright 2010, Nicolas Ponthieu, Julien Grain, Guilaine Lagache
;; 
;; This file is part of Poker.
;; 
;; Poker is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;; 
;; Poker is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;; 
;; You should have received a copy of the GNU General Public License
;; along with Poker.  If not, see <http://www.gnu.org/licenses/>.
;; 
;; For more information about Poker, see http://www.ias.u-psud.fr/poker
;; ======================================================================



pro poker_make_mask, nx, ny, scale, nholes, radius_ratio, mask, $
                     apod_length=apod_length, patch=patch, $
                     holes=holes, clean_border=clean_border, $
                     n_large=n_large, nx_large=nx_large, ny_large=ny_large, ix0=ix0, iy0=iy0

if n_params() lt 1 then begin
   message, /info, "calling sequence: "
   print, "poker_make_mask, nx, ny, scale, nholes, radius_ratio, mask, $"
   print, "                 apod_length=apod_length, patch=patch, $"
   print, "                 holes=holes, clean_border=clean_border, $"
   print, "                 n_large=n_large, nx_large=nx_large, ny_large=ny_large"
   return
endif


if not keyword_set(apod_length) then apod_length = 0
if not keyword_set(nx_large)    then nx_large = round( scale*nx)
if not keyword_set(ny_large)    then ny_large = round( scale*ny)

if keyword_set(n_large) then begin
   nx_large = n_large
   ny_large = n_large
endif

scale_x = float( nx_large)/nx
scale_y = float( ny_large)/ny

mask  = dblarr( nx_large, ny_large)
patch = dblarr( nx_large, ny_large)
holes = dblarr( nx_large, ny_large) + 1.0d0

if scale_x gt 1 then ix0 = long( (scale_x - 1)*nx/2.) else ix0 = 0
if scale_y gt 1 then iy0 = long( (scale_y - 1)*ny/2.) else iy0 = 0

mask[ ix0:ix0+nx-1, iy0:iy0+ny-1] = 1.0d0
patch[ix0:ix0+nx-1, iy0:iy0+ny-1] = 1.0d0

if nholes gt 0 then begin
   xc = randomu( seed, nholes)*nx
   yc = randomu( seed, nholes)*ny
   radius = float(nx)/radius_ratio
   for i=0, nx-1 do begin
      for j=0, ny-1 do begin
         for ih=0, nholes-1 do begin
            if sqrt( (i-xc[ih])^2 + (j-yc[ih])^2) le radius then begin
               mask[  ix0+i,iy0+j] = 0
               holes[ ix0+i,iy0+j] = 0
            endif
         endfor
      endfor
   endfor
endif

if keyword_set(clean_border) then begin
   mask[ix0,      iy0:iy0+ny-1] = 1.0d0
   mask[ix0+nx-1, iy0:iy0+ny-1] = 1.0d0
   mask[ix0:ix0+nx-1,      iy0] = 1.0d0
   mask[ix0:ix0+nx-1, iy0+ny-1] = 1.0d0
endif

if apod_length ne 0 then begin

   x = dindgen(nx)
   y = x*0.0d0 + 1.0d0
   w = where( x le apod_length)
   y[w] = x[w]/apod_length - 1.0d0/(2.0d0*!dpi)*sin(2.0d0*!dpi*x[w]/apod_length)
   w = where( (nx-1-x) le apod_length)
   y[w] = (nx-1-x[w])/apod_length - 1.0d0/(2.0d0*!dpi)*sin(2.0d0*!dpi*(nx-1-x[w])/apod_length)
   print, "mask"

   x = dindgen(ny)
   y1 = x*0.0d0 + 1.0d0
   w = where( x le apod_length)
   y1[w] = x[w]/apod_length - 1.0d0/(2.0d0*!dpi)*sin(2.0d0*!dpi*x[w]/apod_length)
   w = where( (ny-1-x) le apod_length)
   y1[w] = (ny-1-x[w])/apod_length - 1.0d0/(2.0d0*!dpi)*sin(2.0d0*!dpi*(ny-1-x[w])/apod_length)
   taper = y#y1
   mask[ix0:ix0+nx-1, iy0:iy0+ny-1] = mask[ix0:ix0+nx-1, iy0:iy0+ny-1]*taper
endif


end
