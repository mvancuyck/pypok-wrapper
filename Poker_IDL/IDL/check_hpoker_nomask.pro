
;; Demonstration script to use POKER in hierarchical mode
;;-------------------------------------------------------

;; Number of Monte-Carlo simulations
nmc = 30 ; 500

;; Size of the simulated patch (in pixels)
;; For this test, I take a convenient size close to Guilaine's true parameters
nx = 832
ny = nx
nx_large = 1024 ; for the input simulation
ny_large = 1024
res_arcmin = 16.d0/60.d0 ; guessed Herschel resolution ?

index_tt             = -1       ; to start easy, then try with -3
bypass               = 0        ; pour ne pas passer par le F90
remove_1st_bin       = 1        ; 1st bin (good practice)
apodization_fraction = 0.2

;; Input power spectrum
kmax = 2*!dpi*sqrt(2)/2.0d0/(res_arcmin*!arcmin2rad)
k = dindgen( 10000)/9999.*( kmax-1) + 1
pk_in = k^index_tt

;; Input map
cls2map, [0,1], [0,1], nx_large, ny_large, res_arcmin, map, cu_t, index=index_tt

;; ;; Generate mask
;; mask = dblarr(nx, ny) + 1.d0
;; nholes = 30
;; radius_ratio = 40
;; poker_make_mask, nx, ny, 1, nholes, radius_ratio, mask

;; Full map, degraded by a factor 8
nx1   = nx/8
ny1   = nx1
nx1_large = 128
ny1_large = 128
map1  = dblarr(nx1,ny1)
;; mask1 = dblarr(nx1,ny1)
for i1=0, nx1-1 do begin
   for j1=0, ny1-1 do begin
      map1[ i1,j1] = avg( map[ i1*8:(i1+1)*8-1, j1*8:(j1+1)*8-1])
      ;; mask1[i1,j1] = avg( mask[ i1*8:(i1+1)*8-1, j1*8:(j1+1)*8-1])
   endfor
endfor

;; Generate apodized mask to compute the power spectrum
apod_length = round(apodization_fraction*min([nx1,ny1]))
poker_make_mask, nx1, ny1, 1, 0, 1, mask_apod, apod_length=apod_length

;; Compute power spectrum and initialize arrays for the large map
ipoker, map1, res_arcmin/8, k1, pk1, mask=mask_apod, $
        beta=beta, remove=remove_1st_bin, out_params=params1, $
        out_arrays=arrays1, nx_large=nx1_large, ny_large=ny1_large

;; Full resolution, take the same nx and ny to save time on mbb computation
nx2 = 104
ny2 = 104
nx2_large = 128
ny2_large = 128
ipoker, map[0:nx2-1,0:ny2-1], res_arcmin, k2, pk2, mask=mask_apod, $
        beta=beta, remove=remove_1st_bin, out_params=params2, $
        out_arrays=arrays2, nx_large=nx2_large, ny_large=ny2_large

print, "minmax(k1): ", minmax(k1)
print, "minmax(k2): ", minmax(k2)
;; Not much overlap, but I start like this for a try

pk2_all = dblarr(8,8,n_elements(k2))
for i=0,7 do begin
   for j=0,7 do begin
      print, "i, j: ", i, j
      ipoker, map[i*nx2:(i+1)*nx2-1,j*ny2:(j+1)*ny2-1], res_arcmin, k2, pk2, $
              in_arrays=arrays2, in_params=params2
      pk2_al[i,j,*] = pk2
   endfor
endfor

print, "stopped here"
stop








;; Make simulations to compute the transfer function, the covariances etc...
map_temp = dblarr(nx1,ny1)
pk1_res = dblarr(nmc,n_elements(k1))

for imc=0, nmc-1 do begin
   message, /info, "fix me: I average map assuming i know the signal behind the mask"
   message, /info, "it's wrong, i should not assume this"
   print, imc
   cls2map, [0,1], [0,1], nx_large, ny_large, res_arcmin, map_junk, cu_t, /no_cu_t_reset
   for i1=0, nx1-1 do begin
      for j1=0, ny1-1 do begin
         map_temp[ i1,j1] = avg( map[ i1*8:(i1+1)*8-1, j1*8:(j1+1)*8-1])
      endfor
   endfor
   ipoker, map_temp, res_arcmin/8, k1, pk1, mask=mask_apod1*mask1, $
           beta=beta,  remove=remove_1st_bin, in_params=params1, $
           in_arrays=arrays1
   pk1_res[imc,*] = pk1












   
endfor
mc_reduce, pk1_res, pk1_out_avg, sigma_pk1_out, cov_mat_pk1_out, xcorr_pk1_out
new_pk_in = interpol(pk_in,k,k1)

wind, 1, 1, /free, /large
my_multiplot, 2, 2, pp, pp1, /rev
imview, map, title='Input map', position=pp1[0,*]
imview, mask, title='Input mask', position=pp1[1,*], /noerase
imview, map1, title='Full map degraded /8', position=pp1[2,*], /noerase
plot, k, pk_in, position=pp1[3,*], /noerase, /xlog, /ylog, xra=xra, yra=yra
oploterror, k1, pk1_out_avg, sigma_pk1_out, col=70
legendastro, ['Input pk', 'pk1_out_avg'], col=[!p.color,70], box=0, line=0
my_multiplot, /reset

;; Deal with submaps

stop

; pour les petites échelles et les 1/4 de cartes

ipoker, dblarr(nx/2, ny/2), res_arcmin, k2a, pk2a, mask= mask_apod*mask_in[0:nx/2-1,0:ny/2-1],  bypass=bypass, beta=beta,  $
        out_params=params_2a, out_arrays=arrays_2a, remove=remove ,nx_large = nx_large , ny_large=ny_large
ipoker, dblarr(nx/2, ny/2), res_arcmin, k2b, pk2b, mask= mask_apod*mask_in[nx/2:nx-1,0:ny/2-1],  bypass=bypass, beta=beta, $
        out_params=params_2b, out_arrays=arrays_2b, remove=remove ,nx_large = nx_large , ny_large=ny_large
ipoker, dblarr(nx/2, ny/2), res_arcmin, k2c, pk2c, mask= mask_apod*mask_in[0:nx/2-1,ny/2:ny-1],  bypass=bypass, beta=beta, $
        out_params=params_2c, out_arrays=arrays_2c, remove=remove ,nx_large = nx_large , ny_large=ny_large
ipoker, dblarr(nx/2, ny/2), res_arcmin, k2d, pk2d, mask= mask_apod*mask_in[nx/2:nx-1,ny/2:ny-1],  bypass=bypass, beta=beta, $
        out_params=params_2d, out_arrays=arrays_2d, remove=remove ,nx_large = nx_large , ny_large=ny_large

; pour les plus grandes echelles et la carte complète

poker_rebin, mask_in, 4, mask_in_deg
w=where(mask_in_deg ne 1 ,nw)
if nw ne 0 then mask_in_deg[w]=0d0
ipoker, dblarr(nx/2, ny/2), res_arcmin*2, k1, pk1, mask= mask_apod*mask_in_deg, $
        out_params=params_1, out_arrays=arrays_1, remove=remove, bypass=bypass, beta=beta,nx_large = nx_large , ny_large=ny_large

;; Pour la fonction de transfert
;; pas necessaire d'appliquer le masque, mais on laisse pour
;; eviter de reinitialiser une fois ipoker sans masque
;cls2map, [0,1], [0,1], nx, ny, res_arcmin, map,cu_t ,  index=index_tt
cls2map, [0,1], [0,1], long(scale_factor*nx), long(scale_factor*ny), res_arcmin, map,cu_t ,  index=index_tt
nmc1 = 100
for imc=0, nmc1-1 do begin 
   if (imc mod 10) eq 0 then print, imc

   ;; on cree la carte de depart
    ;cls2map, [0,1], [0,1], nx, ny, res_arcmin, map, cu_t , index=index_tt, /no_cu_t_reset
   cls2map, [0,1], [0,1], long(scale_factor*nx), long(scale_factor*ny), res_arcmin, map, cu_t, /no_cu_t_reset
   map = map[0:nx-1,0:ny-1]
   
   ;; Carte complète dégradée une fois pour grandes échelles ( 128*128 )
   poker_rebin, map, 4, map_deg4 ; 4 est le nombre de pixels concatenés dans la map source
   ipoker, map_deg4, res_arcmin*2, k_out, pk_out, in_params=params_1, in_arrays=arrays_1
   
   ;; init simulations de pk_out
   if imc eq 0 then pk_out_res=dblarr(nmc1, n_elements(k_out))
   
   ;; stockage pk_out
   pk_out_res[imc,*]=pk_out
   
endfor

mc_reduce , pk_out_res , pk_out_moy , sigma_pk_out , cov_mat_pk_out , xcorr_pk_out

;; Fonction de transfert du rebin
new_pk_in=interpol(pk_in , k , k_out)
F=pk_out_moy/new_pk_in

wind, 1, 1 , /free
outplot, file='transfer_func', /png
plot, k, pk_in , /xlog , /ylog, xra=[100, 2e4], /xs, /ys
oplot, k_out, pk_out_moy, psym=-8, col=250
oplot, k_out, pk_out_moy/F, psym=-8, col=55
legendastro, ["Spectre d'entree theorique",$
              "Sous-spectre de la carte degradee",$
              "Sous-spectre de la carte degradee avec correction par F"], col=[0,250,55], line=0 , /right
outplot, /close
stop

pk_out_res1a=dblarr(nmc, n_elements(k2a))
pk_out_res1b=dblarr(nmc, n_elements(k2a))
pk_out_res1c=dblarr(nmc, n_elements(k2a))
pk_out_res1d=dblarr(nmc, n_elements(k2a))

;cls2map, [0,1], [0,1], long(scale_factor*nx), long(scale_factor*ny), res_arcmin, map,cu_t ,  index=index_tt
for imc=0, nmc-1 do begin 
   
    ;; on cree la carte de depart
   cls2map, [0,1], [0,1], long(scale_factor*nx), long(scale_factor*ny), res_arcmin, map,cu_t,  index=index_tt, /no_cu_t_reset
   ;cls2map, [0,1], [0,1], nx, ny, res_arcmin, map, cu_t , index=index_tt, /no_cu_t_reset
   map=map[0:nx-1,0:ny-1]

    ;; Carte complète dégradée une fois pour grandes échelles ( 128*128 )
   poker_rebin, map, 4, map_deg4 ; 4 est le nombre de pixels concatenés dans la map source
   ipoker, map_deg4, res_arcmin*2, k_out, pk_out, in_params=params_1, in_arrays=arrays_1

    ;; correction du pk_out par la fonction de transfert
   pk_out_corrige=pk_out/F

    ;; 4 sous cartes non dégradée
   ipoker, map[0:nx/2-1,0:ny/2-1],   res_arcmin, k_out1a, pk_out1a, in_params=params_2a, in_arrays=arrays_2a
   ipoker, map[nx/2:nx-1,0:ny/2-1],  res_arcmin, k_out1b, pk_out1b, in_params=params_2b, in_arrays=arrays_2b
   ipoker, map[0:nx/2-1,ny/2:ny-1],  res_arcmin, k_out1c, pk_out1c, in_params=params_2c, in_arrays=arrays_2c
   ipoker, map[nx/2:nx-1,ny/2:ny-1], res_arcmin, k_out1d, pk_out1d, in_params=params_2d, in_arrays=arrays_2d
   
   pk_out_res1a[imc,*]=pk_out1a
   pk_out_res1b[imc,*]=pk_out1b
   pk_out_res1c[imc,*]=pk_out1c
   pk_out_res1d[imc,*]=pk_out1d
   
   ;; interpolation : pour les pk des grandes echelles avec les k des petites
   w2 = where( k_out1a ne 0 and k_out1a le max(k_out),nw2, compl=w22, ncompl=nw22)
   new_k_full = k_out1a[w2]
   new_pk_out = interpol( pk_out_corrige, k_out, new_k_full)
   
   ;; init vecteur des pi
   if imc eq 0 then pk_total= dblarr(nmc,nw2+1+4*(nw2+nw22))
   
   pk_total[imc, 0:nw2]                         = [pk_out_corrige[0], new_pk_out]
   pk_total[imc, nw2+1          : 2*nw2+nw22]   = pk_out1a
   pk_total[imc, 2*nw2+nw22+1   : 3*nw2+2*nw22] = pk_out1b
   pk_total[imc, 3*nw2+2*nw22+1 : 4*nw2+3*nw22] = pk_out1c
   pk_total[imc, 4*nw2+3*nw22+1 : 5*nw2+4*nw22] = pk_out1d

   if (imc mod 10) eq 0 then print, imc
endfor

mc_reduce, pk_total, pk_total_avg, sigma_pk_total, cov_mat_total, xcorrtot
mc_reduce, pk_out_res1a, pk_out1a_avg, sigma_pk_out1a, cov_mat_1a, xcorrtot1a
mc_reduce, pk_out_res1b, pk_out1b_avg, sigma_pk_out1b, cov_mat_1b, xcorrtot1b
mc_reduce, pk_out_res1c, pk_out1c_avg, sigma_pk_out1c, cov_mat_1c, xcorrtot1c
mc_reduce, pk_out_res1d, pk_out1d_avg, sigma_pk_out1d, cov_mat_1d, xcorrtot1d

my_multiplot, 2, 1, pp, pp1
xra=[100, 2e4]
wind, 1, 1, /free, xs=1200
outplot, file='sous_spectres', /png
plot, k, pk_in , /xlog, /ylog, /xs, /ys, xra=xra, position=pp1[0,*]
oplot, [k_out[0], new_k_full], pk_total_avg, psym=-8, col=200
oplot, k_out1a, pk_out1a_avg, psym=-8, col=10
oplot, k_out1b, pk_out1b_avg, psym=-8, col=50
oplot, k_out1c, pk_out1c_avg, psym=-8, col=150
oplot, k_out1d, pk_out1d_avg, psym=-8, col=250
legendastro, ["Spectre d'entree theorique","Sous-spectre de la carte degradee par Monte-Carlo",$
              "Sous-spectre de la premiere sous-carte par Monte-Carlo",$
              "Sous-spectre de la deuxieme sous-carte par Monte-Carlo",$
              "Sous-spectre de la troisieme sous-carte par Monte-Carlo",$
              "Sous-spectre de la quatrieme sous-carte par Monte-Carlo"], col=[0,200,10,50,150,250], line=0 , /right
outplot, /close
imview, alog( abs(cov_mat_total)), position=pp1[1,*], /noerase

imview, alog( abs(cov_mat_total)), png="cov_mat_total.png"
imview, xcorrtot, png="xcorrtot.png"

inv_C=invert(cov_mat_total)

wind, 1, 1, /free
imview, inv_c##cov_mat_total, title='inv_c##cov_mat_total'
stop


;; ;; si on neglige les cross-correlations entre sous-cartes
;; cov_approx = cov_mat_total
;; cov_approx[nw2+1:*,nw2+1:*] = 0.d0
;; for i=0, n_elements(cov_approx[*,0])-1 do cov_approx[i,i] = cov_mat_total[i,i]
;; inv_c = invert(cov_mat_total)

;;;;--------------------------
;;;; SVD ?
;;la_svd, cov_mat_total, s, U, V, /double, status=status
;;plot, s/s[0], /ylog
;;
;;w = where( s/s[0] lt 1e-5)
;;s[w] = 0.d0
;;w = where( s lt 1e-10, compl=w1)
;;s[w1] = 1.d0/s[w1] 
;;
;;inv_c = V##diag_matrix(s)##transpose(u)
;;imview, inv_c##cov_mat_total
;;;;--------------------------


bin_vec=[k_out[0],k_out1a]
nbins = n_elements(bin_vec)
nl = 1+n_elements(new_k_full)+4*n_elements(k_out1a)
A=dblarr(nbins,nl)

k_vec_tot=[k_out[0],new_k_full,k_out1a,k_out1b,k_out1c,k_out1d]
for i=0, nl-1  do begin
   for j=0, nbins-1  do begin
      if k_vec_tot[i] eq bin_vec[j]  then begin
         A[j,i]=1
      endif
   endfor
endfor

t_A=transpose(A)
bloc_mat1= t_A##inv_C##A
inv_bloc_mat1=invert(bloc_mat1)
vect2=t_A##inv_C##pk_total_avg
vec_SP=inv_bloc_mat1##vect2

;; moyenne des réalisations
vec_SPres=dblarr(nmc,n_elements(bin_vec))
for imc=0, nmc-1 do vec_SPres[imc,*]=inv_bloc_mat1##t_A##inv_C##pk_total[imc,*]

my_multiplot, 2, 1, pp, pp1
mc_reduce, vec_SPres, vec_SPres_moy , sigma_Pk_moy , cov_mat_pk , corr_mat_pk
wind, 1, 1, /free, xs=1200
;outplot, file='final_spec', /png
ploterror, bin_vec, vec_SPres_moy , sigma_Pk_moy , psym=-8 , /xlog , /ylog, $
           xtitle="k", ytitle="P(k)", title="Combined spectrum : indice = -3", position=pp1[0,*]
oplot, k , pk_in , col=250
oploterror, bin_vec, vec_SPres_moy , sigma_Pk_moy , psym=8
imview, abs(corr_mat_pk), position=pp1[1,*], /noerase, title='Abs( corr. matrix)', imrange=[0,1]
;outplot, /close
stop
stop

imview, abs(corr_mat_pk), png='corr_mat_pk.png'


wind, 1, 1, /free
ploterror, bin_vec, vec_spres_moy-interpol(pk_in,k,bin_vec), sigma_Pk_moy, /nodata
oploterror, bin_vec, vec_spres_moy-interpol(pk_in,k,bin_vec), sigma_Pk_moy, psym=-8, col=250, errcol=250
oplot, [0,2e4], [0,0]








wind, 1, 1, /free
imview, cov_mat_pk


wind, 1, 1, /free
imview, alog(abs(cov_mat_pk))

stop

mc_reduce, pk_res1full, pk_avg, sigma_pk_avg, cov_mat1, xcorr1
mc_reduce, pk_res2a, pk_1_avga, sigma_pk_1_avga, cov_mat2a, xcorr2
mc_reduce, pk_res2b, pk_1_avgb, sigma_pk_1_avgb, cov_mat2b
mc_reduce, pk_res2c, pk_1_avgc, sigma_pk_1_avgc, cov_mat2c
mc_reduce, pk_res2d, pk_1_avgd, sigma_pk_1_avgd, cov_mat2d

;init matrices
c1=dblarr(nw2,nw2)
c2=dblarr(nw2,nw2)
c3=dblarr(nw2,nw2)
c4=dblarr(nw2,nw2)
c5=dblarr(nw2,nw2)
c6=dblarr(nw2,nw2)
c7=dblarr(nw2,nw2)
c8=dblarr(nw2,nw2)

d1=dblarr(nw2,nw2)
d2=dblarr(nw2,nw2)
d3=dblarr(nw2,nw2)
d4=dblarr(nw2,nw2)
d5=dblarr(nw2,nw2)
d6=dblarr(nw2,nw2)
d7=dblarr(nw2,nw2)
d8=dblarr(nw2,nw2)


for imc=0, nmc-1 do begin
   
   for i=0, n_elements(new_k_full)-1 do begin

      for j=0, n_elements(new_k_full)-1 do begin

;matrice de covariance de la corrélation des petites cartes entre elles
        
        ; c1[i,j]=(pk_res2a[imc,i]-pk_1_avga[i])*(pk_res2b[imc,j]-pk_1_avgb[j])/sqrt(cov_mat2a[i,i]*cov_mat2b[j,j])
        ; c2[i,j]=(pk_res2b[imc,i]-pk_1_avgb[i])*(pk_res2c[imc,j]-pk_1_avgc[j])/sqrt(cov_mat2b[i,i]*cov_mat2c[j,j])
       ;  c3[i,j]=(pk_res2c[imc,i]-pk_1_avgc[i])*(pk_res2d[imc,j]-pk_1_avgd[j])/sqrt(cov_mat2c[i,i]*cov_mat2d[j,j])
       ;  c4[i,j]=(pk_res2d[imc,i]-pk_1_avgd[i])*(pk_res2a[imc,j]-pk_1_avga[j])/sqrt(cov_mat2d[i,i]*cov_mat2a[j,j])

         
         c1[i,j]=(pk_res2a[imc,i]-pk_1_avga[i])*(pk_res2b[imc,j]-pk_1_avgb[j])
         c2[i,j]=(pk_res2b[imc,i]-pk_1_avgb[i])*(pk_res2c[imc,j]-pk_1_avgc[j])
         c3[i,j]=(pk_res2c[imc,i]-pk_1_avgc[i])*(pk_res2d[imc,j]-pk_1_avgd[j])
         c4[i,j]=(pk_res2d[imc,i]-pk_1_avgd[i])*(pk_res2a[imc,j]-pk_1_avga[j])
         
;matrice de covariance entre les grandes echelles de cartes et
;les petites echelles de cartes

        ; c5[i,j]=(pk_res1full[imc,i]-pk_avg[i])*(pk_res2a[imc,j]-pk_1_avga[j])/sqrt(cov_mat1[i,i]*cov_mat2a[j,j])
       ;  c6[i,j]=(pk_res1full[imc,i]-pk_avg[i])*(pk_res2b[imc,j]-pk_1_avgb[j])/sqrt(cov_mat1[i,i]*cov_mat2b[j,j])
        ; c7[i,j]=(pk_res1full[imc,i]-pk_avg[i])*(pk_res2c[imc,j]-pk_1_avgc[j])/sqrt(cov_mat1[i,i]*cov_mat2c[j,j])
        ; c8[i,j]=(pk_res1full[imc,i]-pk_avg[i])*(pk_res2d[imc,j]-pk_1_avgd[j])/sqrt(cov_mat1[i,i]*cov_mat2d[j,j])

         c5[i,j]=(pk_res1full[imc,i]-pk_avg[i])*(pk_res2a[imc,j]-pk_1_avga[j])
         c6[i,j]=(pk_res1full[imc,i]-pk_avg[i])*(pk_res2b[imc,j]-pk_1_avgb[j])
         c7[i,j]=(pk_res1full[imc,i]-pk_avg[i])*(pk_res2c[imc,j]-pk_1_avgc[j])
         c8[i,j]=(pk_res1full[imc,i]-pk_avg[i])*(pk_res2d[imc,j]-pk_1_avgd[j])
      

; accumulateur des petites cartes entre elles
         d1[i,j]=d1[i,j]+c1[i,j]
         d2[i,j]=d2[i,j]+c2[i,j]
         d3[i,j]=d3[i,j]+c3[i,j]
         d4[i,j]=d4[i,j]+c4[i,j]

; accumulateur des grandes echelles de cartes correlees avec
;les petites echelles de cartes
         d5[i,j]=d5[i,j]+c5[i,j]
         d6[i,j]=d6[i,j]+c6[i,j]
         d7[i,j]=d7[i,j]+c7[i,j]
         d8[i,j]=d8[i,j]+c8[i,j]
         
      endfor

   endfor
   
endfor

; moyenne des matrices de covariance
d1=d1/nmc
d2=d2/nmc
d3=d3/nmc
d4=d4/nmc
d5=d5/nmc
d6=d6/nmc
d7=d7/nmc
d8=d8/nmc





print,'C12=', d1
print,'C23=', d2
print,'C34=', d3
print,'C41=', d4
print,'Ca=',  d5
print,'Cb=',  d6
print,'Cc=',  d7
print,'Cd=',  d8


nk = n_elements( new_k_full)
xcorr5 = dblarr(nk,nk)
for i=0, nk-1 do begin
   for j=0, nk-1 do begin
      ;; covariance de la grande carte avec la petite sous carte
      xcorr5[i,j] = d5[i,j]/sqrt( d5[i,i]*d5[j,j])
   endfor
endfor


my_multiplot, 4, 4, /rev, pp, pp1
wind, 1, 1, /free, /large
imview, cov_mat2a, title='Covar Petite carte', position=pp[0,0,*], /noerase
imview, cov_mat1, title='Covar Grande carte', position=pp[0,1,*], /noerase
imview, d1, title='Covar Petite_a x Petite_b', position=pp[1,0,*], /noerase
imview, d2, title='Covar Petite_b x Petite_c', position=pp[2,0,*], /noerase
imview, d3, title='Covar Petite_c x Petite_d', position=pp[3,0,*], /noerase
imview, d4, title='Covar Petite_d x Petite_a', position=pp[0,2,*], /noerase
imview, d5, title='Covar Petite_a x Grande', position=pp[1,1,*], /noerase
imview, d6, title='Covar Petite_b x Grande', position=pp[2,1,*], /noerase
imview, d7, title='Covar Petite_c x Grande', position=pp[1,2,*], /noerase
imview, d8, title='Covar Petite_d x Grande', position=pp[2,2,*], /noerase
imview, xcorr2, title='XCorr Petite x Petite', position=pp[3,2,*], /noerase
imview, xcorr1, title='XCorr Grande x Grande', position=pp[2,3,*], /noerase
imview, xcorr5, title='XCorr Petite x Grande', position=pp[3,3,*], /noerase

stop

xrange = minmax( k_out[where(k_out ne 0)])
;P(k) grandes échelles
plot, k_out, pk_out, /xlog, /ylog, $
      xrange=xrange, /xs, /nodata
;P(k) grandes échelles avec mc_reduce
oploterror, k_out, pk_avg, sigma_pk_avg, psym=4, col=250, errcol=250

;P(k) sans ipoker
oplot, k,  pk_in

;P(k) petites échelles avec mc_reduce
oploterror, k_out1a, pk_1_avga, sigma_pk_1_avga, col=5, errcol=5
oploterror, k_out1b, pk_1_avgb, sigma_pk_1_avgb, col=55, errcol=55
oploterror, k_out1c, pk_1_avgc, sigma_pk_1_avgc, col=100, errcol=100
oploterror, k_out1d, pk_1_avgd, sigma_pk_1_avgd, col=150, errcol=150

end
