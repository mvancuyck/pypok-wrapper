
pro poker_init

defsysv, "!fwhm2sigma", 1.00d0/sqrt(8.000d0*alog(2.000d0))
defsysv, "!arcmin2rad", 1.000d0/60.000d0*!dtor
defsysv, "!h2d", 360.0d0/24.0d0
defsysv, "!undef", -32768
defsysv, '!scratch', getenv('SCRATCH')
defsysv, '!arch',    getenv("MY_ARCH")
defsysv, "!arcsec2rad", 1.0d0/3600.0d0*!dtor

print, ""
print, "====> POKER now available"
end
