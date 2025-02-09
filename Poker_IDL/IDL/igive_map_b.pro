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



pro pix2mn, ipix, nx, ny, m, n
  n = long(ipix)/long(nx)
  m = long(ipix) - long(n)*long(nx)
end

;; pro mn2k, m, n, nx, ny, kx, ky, k
;;   if (long(m) le long(nx)/2) then begin
;;      kx = double(m)/nx
;;   endif else begin
;;      kx = double(nx-m)/nx
;;   endelse
;;   if (long(n) le long(ny)/2) then begin
;;      ky = double(n)/ny
;;   endif else begin
;;      ky = double(ny-n)/ny
;;   endelse
;;   k = sqrt( kx^2 + ky^2) ;reduced units
;; end

pro mn2k, m, n, nx, ny, kx, ky, k

  kx = m*0.0d0
  ky = m*0.0d0
  k  = m*0.0d0

  w = where( long(m) le long(nx)/2, nw)
  if nw ne 0 then kx[w] = double(m[w])/nx

  w = where( long(m) gt long(nx)/2, nw)
  if nw ne 0 then kx[w] = double(nx-m[w])/nx

  w = where( long(n) le long(ny)/2, nw)
  if nw ne 0 then ky[w] = double(n[w])/ny

  w = where( long(n) gt long(ny)/2, nw)
  if nw ne 0 then ky[w] = double(ny-n[w])/ny

  k = sqrt( kx^2 + ky^2) ;reduced units
end


;; pro k2bin, k, nbins, bintab, b
;;   test = 0
;;   b    = 0
;;   while ( (test eq 0) and (b lt nbins)) do begin
;;      if ( (bintab[b] le k) and ( k lt bintab[b+1])) then test = 1
;;      b = b+1
;;   endwhile
;; 
;;   ;correct for extra addition and the end of "do while" loop
;;   b = b-1
;; 
;;   ;if k is outside bintab bounds, set b to undef value
;;   if (test eq 0) then b = -32768
;; 
;; end

pro k2bin, k, nbins, bintab, b, xbintab
  
  b = lonarr( n_elements(k)) + !undef

  xbintab = lonarr( n_elements(bintab))

  for i=0L, nbins-1 do begin
     w = where( k ge bintab[i] and k lt bintab[i+1], nw)
     if nw ne 0 then begin
        b[w]       = i
        xbintab[i] = 1
     endif
  endfor
end


;; pro igive_map_b, naxis1, naxis2, bintab, map_b
;; 
;;   npix  = long(naxis1)*long(naxis2)
;;   map_b = dblarr( naxis1, naxis2)
;;   nbins = n_elements(bintab)-1
;; 
;;   for i=0L, npix-1 do begin
;;      pix2mn, i, naxis1, naxis2, m, n
;;      mn2k, m, n, naxis1, naxis2, kx, ky, k
;;      k2bin, k, nbins, bintab, b
;;      
;;      map_b[m,n] = b
;;   endfor
;; 
;; end

pro igive_map_b, naxis1, naxis2, bintab, map_b, xbintab

  if n_params() lt 1 then begin
     message, /info, "Calling sequence"
     print, "igive_map_b, naxis1, naxis2, bintab [reduced units], map_b, xbintab"
     return
  endif

  npix = long(naxis1)*long(naxis2)
  nbins = n_elements(bintab)-1
  pix  = lindgen( npix)
  map_b = dblarr( naxis1, naxis2) + !undef

  pix2mn, pix, naxis1, naxis2, m, n
  mn2k, m, n, naxis1, naxis2, kx, ky, k
  k2bin, k, nbins, bintab, b, xbintab

  map_b[*] = b

end
