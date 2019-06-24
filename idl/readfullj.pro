FUNCTION readfullJ, L0=l0, LN=lN, J20, reJ21, imJ21, reJ22, imJ22

  DEFAULT_J_FILE   = "J.dat"
  DEFAULT_J20_FILE = "J20.out"
  DEFAULT_reJ21_FILE = "reJ21.out"
  DEFAULT_imJ21_FILE = "imJ21.out"
  DEFAULT_reJ22_FILE = "reJ22.out"
  DEFAULT_imJ22_FILE = "imJ22.out"

@geometry.common
@spectrum.common

  IF (NOT keyword_set(L0)) THEN l0 = 0
  IF (NOT keyword_set(LN)) THEN lN = spectrum.Nspect-1
  Nspect = lN - l0 + 1

  CASE geometryType OF
    "ONE_D_PLANE": BEGIN
      Jnu = dblarr(geometry.Ndep, Nspect)
      offset = 8L*l0 * geometry.Ndep
    END
    "TWO_D_PLANE": BEGIN
      Jnu = dblarr(geometry.Nx, geometry.Nz, Nspect)
      offset = 8L*l0 * (geometry.Nx*geometry.Nz)
    END
    "THREE_D_PLANE": BEGIN
      Jnu = dblarr(geometry.Nx, geometry.Ny, geometry.Nz, Nspect)
      offset = 8L*l0 * (geometry.Nx*geometry.Ny*geometry.Nz)
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      Jnu = dblarr(geometry.Nradius, Nspect)
      offset = 8L*l0 * geometry.Nradius
    END
  ENDCASE

  openr, lun, /GET_LUN, DEFAULT_J_FILE
  point_lun, lun, offset
  readu, lun, Jnu
  free_lun, lun
  
  J20 = Jnu
  openr, lun, /GET_LUN, DEFAULT_J20_FILE
  point_lun, lun, offset
  readu, lun, J20
  free_lun, lun

  reJ21 = Jnu
  openr, lun, /GET_LUN, DEFAULT_reJ21_FILE
  point_lun, lun, offset
  readu, lun, reJ21
  free_lun, lun


  imJ21 = Jnu
  openr, lun, /GET_LUN, DEFAULT_imJ21_FILE
  point_lun, lun, offset
  readu, lun, imJ21
  free_lun, lun


  reJ22 = Jnu
  openr, lun, /GET_LUN, DEFAULT_reJ22_FILE
  point_lun, lun, offset
  readu, lun, reJ22
  free_lun, lun

  imJ22 = Jnu
  openr, lun, /GET_LUN, DEFAULT_imJ22_FILE
  point_lun, lun, offset
  readu, lun, imJ22
  free_lun, lun
  
  
  return, Jnu
END

