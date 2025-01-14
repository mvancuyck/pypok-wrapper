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





;;Purpose
;; Computes Fourier wavevectors for a given array

pro give_map_k, res, map_t, map_k, map_kx, map_ky


if n_params() lt 1 then begin
   message, /info, "calling sequence:"
   print, "give_map_k, res, map_t, map_k, map_kx, map_ky"
   return
endif

s = size(map_t)
case s[0] of
   1:begin
      nx = s[1]
      ny = 1
   end
   2:begin
      nx = s[1]
      ny = s[2]
   end
endcase

lmap_x = double(nx*res)
lmap_y = double(ny*res)

map_kx = double( map_t)*0.0d0
map_ky = double( map_t)*0.0d0
map_k  = double( map_t)*0.0d0

for im=0L, nx-1 do begin
   if im le nx/2 then m1 = double(im) else m1 = double(im-nx)
   for in=0L, ny-1 do begin
      if in le ny/2 then n1 = double(in) else n1 = double(in-ny)
      kx = m1/lmap_x
      ky = n1/lmap_y
      map_kx[im,in] = kx
      map_ky[im,in] = ky
      map_k[ im,in] = sqrt( kx^2 + ky^2)
   endfor
endfor

end
