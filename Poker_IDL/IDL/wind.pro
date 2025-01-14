
;; PRO wind, nwin, pos, free=free, xsize=xsize, ysize=ysize, title=title, large=large
;; 
;; if not keyword_set( xsize) then xsize = 600 ; 700 ; 425
;; if not keyword_set( ysize) then ysize = 500 ; 600 ; 350
;; 
;; if keyword_set(large) then begin
;;    xsize=1000
;;    ysize=900
;; endif
;; 
;; xpos = 0
;; ypos = 500
;; 
;; if defined(pos) then xpos = xsize*(pos-1)
;; 
;; if !my_window eq -1 then begin
;;    window, nwin, free=free, xpos = xpos, ypos = ypos, xsize = xsize, ysize = ysize, title=title
;; endif else begin
;;    wset, !my_window
;; endelse
;; 
;; END
; FXD: added /iconic keyword
PRO wind, nwin, pos, free=free, xsize=xsize, ysize=ysize, $
          title=title, large=large, xpos=xpos, ypos=ypos, $
          xlarge=xlarge, ylarge=ylarge, iconic = iconic

  if !d.name EQ 'NULL' then begin
     nwin=10
     pos=[0,0,0,0]
     return
  endif
  
my_screen_size = get_screen_size()

if not keyword_set( xsize) then xsize = 0.36*my_screen_size[0]
if not keyword_set( ysize) then ysize = 0.48*my_screen_size[1]

if keyword_set(xlarge) then begin
   if strupcase(!host) eq "NIKA2B" then begin
      xsize = 1000
   endif else begin
      xsize = 0.8*my_screen_size[0]
   endelse
endif
if keyword_set(ylarge) then begin
   if strupcase(!host) eq "NIKA2B" then begin
      ysize = 700
   endif else begin
      ysize = 0.85*my_screen_size[1]
   endelse
endif
if keyword_set(large) then begin
   if strupcase(!host) eq "NIKA2B" then begin
      xsize = 1000
      ysize = 700
   endif else begin
      xsize = 0.8*my_screen_size[0]
      ysize = 0.85*my_screen_size[1]
   endelse
endif

;; Leave some margin
xsize = xsize < (0.95*my_screen_size[0])
ysize = ysize < (0.95*my_screen_size[1])

if not keyword_set(xpos) then xpos = 0
if not keyword_set(ypos) then ypos = my_screen_size[1]-ysize
if defined(pos) then xpos = float(pos-1)*0.2*my_screen_size[0]

if !my_window eq -1 then begin
   window, nwin, free=free, xpos = xpos, ypos = ypos, xsize = xsize, ysize = ysize,title=title
;; wanted to keep track for the !d.window for blink purpose but it does not work yet
;;   if keyword_set(title) then begin
;;      ;; recreate the window now that I have its number
;;      title=strtrim(!d.window,2)+"/ "+title
;;      window, !d.window, xpos = xpos, ypos = ypos, xsize = xsize, ysize = ysize, title=title
;;      print, !d.window
;;      stop
;;   endif

endif else begin
   print, "wseting !my_window = "+strtrim(!my_window,2)
   wset, !my_window
endelse
wshow, !d.window, iconic = iconic

;; to be compatible with cgplot and restore ordinary white/black background
device, decompose=0;, static_gray=1

END
