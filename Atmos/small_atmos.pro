; IDL script to make an atmosphere smaller


a = read3datmos('snapshot_B_0d0_252x126x252_le_00100.atmos')
nx = 50
ny = 50

openw, unit, /GET_LUN, 'small.atmos', /XDR
writeu, unit, long([nx,ny,a.Nz,a.Nhydr]) 
writeu, unit, long([1, 2])
writeu, unit, double(a.dx), double(a.dy)
writeu, unit, double(a.z)               
writeu, unit, double(a.t[0:nx-1,0:ny-1,*])
writeu, unit, double(a.n_elec[0:nx-1,0:ny-1,*])
writeu, unit, double(a.vturb[0:nx-1,0:ny-1,*])
writeu, unit, double(a.vx[0:nx-1,0:ny-1,*])
writeu, unit, double(a.vy[0:nx-1,0:ny-1,*])
writeu, unit, double(a.vz[0:nx-1,0:ny-1,*])
writeu, unit, double(a.nh[0:nx-1,0:ny-1,*])

free_lun, unit
