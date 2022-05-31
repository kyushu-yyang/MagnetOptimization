FUNCTION getIronPermeability( Model, n, H ) RESULT(mu)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: B, mu, H, Ba, pa, Bb, pb, La, Lb
  
  ! variables needed inside function
  Logical :: GotIt
  TYPE(ValueList_t), POINTER :: material

  ! get pointer on list for material
  material => GetMaterial()
  IF (.NOT. ASSOCIATED(material)) THEN
      CALL Fatal('getIronPermeability', 'No material found')
  END IF

  ! fit constants
  !Ba = 0.46628
  !Bb = 1.7218
  !pa = 9337.7
  !pb = 89.434

  ! data at 4K provided by TOSHIBA
  Ba = 1.330823
  Bb = 1.649819
  pa = 81610.10
  pb = 72.71

  ! permeability of iron at 4K
  La = dcosh(H/pa) / dsinh(H/pa) - (pa/H)
  Lb = dcosh(H/pb) / dsinh(H/pb) - (pb/H)

  B  = Ba*La + Bb * dsinh(dabs(H)/pb)/dcosh(dabs(H)/pb) * Lb
  mu = B / H
  
END FUNCTION getIronPermeability
