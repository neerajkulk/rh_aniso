@initrh
.r readall

nx = 4
ny = 4
nz = 82
nl = 568

arr = dblarr(nx,ny,nz,nl)
offset = 0.0
openr, lun, /GET_LUN, 'J.dat'
point_lun, lun,offset
readu, lun, arr
writefits, 'J00.fits', arr
free_lun,lun


arr = dblarr(nx,ny,nz,nl)
offset = 0.0
openr, lun, /GET_LUN, 'J20.out'
point_lun, lun,offset
readu, lun, arr
writefits, 'J20.fits', arr
free_lun,lun

arr = dblarr(nx,ny,nz,nl)
offset = 0.0
openr, lun, /GET_LUN, 'reJ21.out'
point_lun, lun,offset
readu, lun, arr
writefits, 'reJ21.fits', arr
free_lun,lun


arr = dblarr(nx,ny,nz,nl)
offset = 0.0
openr, lun, /GET_LUN, 'reJ22.out'
point_lun, lun,offset
readu, lun, arr
writefits, 'reJ22.fits', arr
free_lun,lun


arr = dblarr(nx,ny,nz,nl)
offset = 0.0
openr, lun, /GET_LUN, 'imJ21.out'
point_lun, lun,offset
readu, lun, arr
writefits, 'imJ21.fits', arr
free_lun,lun


arr = dblarr(nx,ny,nz,nl)
offset = 0.0
openr, lun, /GET_LUN, 'imJ22.out'
point_lun, lun,offset
readu, lun, arr
writefits, 'imJ22.fits', arr
free_lun,lun

