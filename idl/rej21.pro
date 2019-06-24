FUNCTION readfullJ, L0=l0, LN=lN, reJ21

  DEFAULT_J_FILE   = "J.dat"
  DEFAULT_reJ21_FILE = "reJ21.out"

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
      Jnu = dblarr(geometry.Nx, geometry.Ny, geometry.Nz, Nspect) + 1.0
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

  IF (n_params() EQ 1) THEN BEGIN
    reJ21 = Jnu
;    openr, lun, /GET_LUN, DEFAULT_reJ21_FILE
;    point_lun, lun, offset
    readu, lun, reJ21
    free_lun, lun
  ENDIF

  return, Jnu
END

