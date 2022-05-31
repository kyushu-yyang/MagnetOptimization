!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Module for solving magnetic vector potential in Cartesian and
! *  cylindrically symmetric 2D cases. In both cases the vector potential
! *  is reduced to a single component. 
! *
! *  Authors: Juha Ruokolainen, Mika Malinen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 30.11.2012
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! * Module derived from MagnetoDynamincs2D to calculate magnetic field
! * generated from source field without conductors.
! *
! * Authors: Ye Yang
! * Email:   yang.ye@qst.go.jp
! * Address: QST, Chiba, Japan
! *
! * Original Date: 10.11.2021
! *
!/******************************************************************************

!> \ingroup Solvers
!> \{

!------------------------------------------------------------------------------
SUBROUTINE StaticMagneticField2D_Init( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: HandleAsm, Found
  CHARACTER(*), PARAMETER :: Caller = 'StaticMagneticField2D_Init'
 
  Params => GetSolverParams()
  CALL ListAddInteger( Params, 'Variable Dofs',1 )
  CALL ListAddNewString( Params,'Variable','Potential')
  CALL ListAddNewLogical( Params,'Apply Mortar BCs',.TRUE.)
  CALL ListAddNewLogical( Params,'Use Global Mass Matrix',.TRUE.)

  HandleAsm = ListGetLogical( Params,'Handle Assembly',Found )
  
  IF( HandleAsm ) THEN
    IF( CurrentCoordinateSystem() == AxisSymmetric .OR. &
        CurrentCoordinateSystem() == CylindricSymmetric ) THEN 
      CALL Warn(Caller,'Handle assembly not yet available in axisymmetric case!')
      HandleAsm = .FALSE.
    END IF
    IF( ListGetLogicalAnyMaterial(Model, 'Zirka material') ) THEN
      CALL Warn(Caller,'Handle assembly not yet available for Zirka material!')
      HandleAsm = .FALSE.
    END IF
    IF( ListCheckPresentAnyBodyForce(Model, 'Lorentz velocity') ) THEN
      CALL Warn(Caller,'Handle assembly not yet available for "Lorentz velocity"!')
      HandleAsm = .FALSE.
    END IF
    IF( ListCheckPresentAnyBodyForce(Model, 'Angular Velocity') ) THEN
      CALL Warn(Caller,'Handle assembly not yet available for "Angular Velocity"!')
      HandleAsm = .FALSE.
    END IF
    IF(.NOT. HandleAsm ) THEN
      CALL Info(Caller,'Reverting to old bulk assembly routine!')
      CALL ListAddLogical( Params, 'Handle Assembly',.FALSE. )
    END IF
  END IF
  
!------------------------------------------------------------------------------
END SUBROUTINE StaticMagneticField2D_Init ! }}}
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve the magnetic vector potential expressed in terms of a single component.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE StaticMagneticField2D( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  USE CircuitUtils
  USE ZirkaUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t), POINTER :: Element

  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,n, nb, nd, t, istat, Active, NonlinIter, iter, tind

  TYPE(ValueList_t), POINTER :: BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
               NEWX(:), NEWY(:), POT(:)

  TYPE(Mesh_t),   POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  
  LOGICAL :: NewtonRaphson = .FALSE., CSymmetry, SkipDegenerate, &
      HandleAsm, MassAsm, ConstantMassInUse = .FALSE.
  LOGICAL :: InitHandles
  INTEGER :: CoupledIter
  TYPE(Variable_t), POINTER :: IterV, CoordVar

  TYPE(Matrix_t),POINTER::CM
  TYPE(GlobalHysteresisModel_t), POINTER :: ZirkaModel
  CHARACTER(*), PARAMETER :: Caller = 'StaticMagneticField2D'

  REAL(KIND=dp), ALLOCATABLE, SAVE :: MassValues(:)

  TYPE(TabulatedBasisAtIp_t), POINTER, SAVE :: BasisFunctionsAtIp(:)=>NULL()
  LOGICAL, SAVE :: BasisFunctionsInUse = .FALSE.

  ! new
  REAL(KIND=dp) :: Hsx, Hsy, x, y, curr
  REAL(KIND=dp), PARAMETER :: phi1 = -PI/8, phi2 = PI/8
  INTEGER :: nmax
  TYPE(Nodes_t) :: Nodes
  REAL(KIND=dp), ALLOCATABLE :: Hs(:,:)
  
!------------------------------------------------------------------------------

  CALL Info( Caller,'------------------------------------------------', Level=4 )
  CALL Info( Caller, 'Solving equation for magnetic vector potential', Level=4 )
  CALL Info( Caller,'------------------------------------------------', Level=4 )

  CALL DefaultStart()
  
  ! Allocate some permanent storage, this is done first time only:
  ! --------------------------------------------------------------
  NULLIFY(BC)
  Mesh => GetMesh()
  SolverParams => GetSolverParams()


  IF( ListGetLogical( SolverParams,'Store Basis Functions',Found ) ) THEN
    CALL TabulateBasisFunctions()
  END IF
  
  CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
  
  HandleAsm = ListGetLogical( SolverParams,'Handle Assembly',Found )
  IF( HandleAsm ) THEN
    CALL Info(Caller,'Performing handle version of bulk element assembly',Level=7)
  ELSE
    CALL Info(Caller,'Performing legacy version of bulk element assembly',Level=7)      
  END IF

  MassAsm = Transient
  IF( ConstantMassInUse ) MassAsm = .FALSE.
  
  NewtonRaphson = GetLogical(SolverParams, 'Newton-Raphson Iteration', Found)
  IF(GetCoupledIter()>1) NewtonRaphson = .TRUE.

  NonlinIter = GetInteger(SolverParams,'Nonlinear System Max Iterations',Found)
  IF(.NOT.Found) NonlinIter = 1

  SkipDegenerate = GetLogical(SolverParams, 'Skip Degenerate Elements',Found ) 
  
  CALL Info(Caller,'Initializing Zirka hysteresis models', Level=10)
  CALL InitHysteresis(Model, Solver)

  ! initializing the magnetic field derived from line current
  CALL Info( Caller, 'Initializing source field, Hs', Level=10 )
  curr = GetConstReal( Model % Constants, 'Operation Current', Found )
  IF (.NOT.FOUND) CALL Fatal('StaticMagneticField', 'The operation current is not defined.')

  nmax = GetInteger(Model % Constants, 'Field Order', Found )
  IF (.NOT.FOUND) nmax = 5

  PRINT *, '- operation current:', curr, 'A'
  PRINT *, '- order of calculation:', nmax

  Active = GetNOFActive()
  ALLOCATE(Hs(2,Active))

  DO t = 1, Active
    Element => GetActiveElement(t)
    n = GetElementNOFNodes(Element) ! n = 3 for triangluar element
    CALL GetElementNodes( Nodes, Element )

    x = SUM(Nodes % x(1:n))/n
    y = SUM(Nodes % y(1:n))/n

!    CALL GetArcSourceField( x, y, phi1, phi2, curr, Hsx, Hsy )
    CALL GetSourceField( x, y, curr, nmax, Hsx, Hsy )
    Hs(1,t) = Hsx
    Hs(2,t) = Hsy

    !IF ( MOD(t,50).EQ.0 ) THEN
      !PRINT *, t, x, y, Hsx*4*PI*1e-7, Hsy*4*PI*1e-7
    !END IF
  END DO

  ! debug
  !x = 0.0
  !y = 0.0
  !curr = 265.0
  !nmax = 6
  !CALL GetSourceField( x, y, curr, nmax, Hsx, Hsy )
  !PRINT *, '- Center field, Bx: ', Hsx*PI*4d-7, ' [Tesla] at 265A'
  !PRINT *, '- Center field, By: ', Hsy*PI*4d-7, ' [Tesla] at 265A'
  PRINT *, '- Length of Hs array: ', SIZE(Hs,1), 'x', SIZE(Hs,2)

  ! Nonlinear iteration loop
  ! -------------------------
  DO iter = 1,NonlinIter
    IF(Iter > 1) NewtonRaphson=.TRUE.
    ! System assembly:
    ! ----------------

    Active = GetNOFActive()
    CALL DefaultInitialize()
    IF( ConstantMassInUse ) THEN
      Solver % Matrix % MassValues = MassValues
    END IF

    InitHandles = .TRUE.
    tind = 0
    
    ! Domain element loop
    ! -------------------
!$omp parallel do private(Element,n,nd)   
    DO t=1,active
       Element => GetActiveElement(t)
       n  = GetElementNOFNodes(Element)
       nd = GetElementNOFDOFs(Element)
       nb = GetElementNOFBDOFs(Element)
       IF( SkipDegenerate .AND. DegenerateElement( Element ) ) THEN
         CALL Info(Caller,'Skipping degenerate element:'//TRIM(I2S(t)),Level=12)
         CYCLE
       END IF
       IF( HandleAsm ) THEN
         CALL LocalMatrixHandles(  Element, n, nd+nb, nb, InitHandles )
       ELSE
         CALL LocalMatrix(Element, n, nd, Hs(1,t), Hs(2,t))
       END IF
    END DO
!$omp end parallel do
    
    CALL DefaultFinishBulkAssembly()

    ! Boundary element loop
    ! ---------------------
    Active = GetNOFBoundaryElements()
!$omp parallel do private(Element, n, nd, BC,Found)
    DO t=1,active
      Element => GetBoundaryElement(t)
      BC=>GetBC( Element )
      IF(.NOT.ASSOCIATED(BC)) CYCLE

      IF(GetLogical(BC,'Infinity BC',Found)) THEN
         n  = GetElementNOFNodes( Element )
         nd = GetElementNOFDOFs( Element )
         CALL LocalMatrixInfinityBC(Element, n, nd)
      ELSE IF(GetLogical(BC,'Air Gap',Found)) THEN
         n  = GetElementNOFNodes( Element )
         nd = GetElementNOFDOFs( Element )
         CALL LocalMatrixAirGapBC(Element, BC, n, nd)
      END IF
    END DO
!$omp end parallel do

    CALL DefaultFinishAssembly()
    
    CALL SetMagneticFluxDensityBC()
    CALL DefaultDirichletBCs()

    IF( ListGetLogical( SolverParams,'Constant Mass Matrix',Found ) ) THEN
      IF( .NOT. ConstantMassInUse ) THEN
        ALLOCATE( MassValues( SIZE( Solver % Matrix % MassValues ) ) )
        MassValues = Solver % Matrix % MassValues 
        ConstantMassInUse = .TRUE.
        MassAsm = .FALSE.
      END IF
    END IF
    
    Norm = DefaultSolve()
 
    IF( DefaultConverged() ) EXIT
  END DO
  
  ! For cylindrical symmetry the model lumping has not been implemented
  IF( .NOT. CSymmetry ) THEN
    CALL CalculateLumpedTransient()
  END IF

  CALL DriveHysteresis(model, solver)

  ! This updates coordinates when using ElmerPost for visualization
  CoordVar => VariableGet(Mesh % Variables,'Coordinates')
  IF(ASSOCIATED(CoordVar)) THEN
    DO i=1,Mesh % NumberOfNodes
      j = 3*(CoordVar % Perm(i)-1)
      CoordVar % Values(j+1) = Mesh % Nodes % x(i)
      CoordVar % Values(j+2) = Mesh % Nodes % y(i)
      CoordVar % Values(j+3) = Mesh % Nodes % z(i)
    END DO
  END IF

  CALL DefaultFinish()

  CALL Info(Caller,'All done',Level=8)
  
CONTAINS


  !> Tabulate basis functions and their weights so that we do not need to compute them in the assembly
  !> process. This assumes that the geometry is not changing. This could be later moved to library
  !> but for now we use a local implementation. 
  !---------------------------------------------------------------------------------------------------
  SUBROUTINE TabulateBasisFunctions()

    INTEGER :: i, t, n, tind
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    LOGICAL :: stat
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp) :: detJ, Weight
    INTEGER :: Phase

    IF( BasisFunctionsInUse ) RETURN
    
    n = Mesh % MaxElementDofs
    ALLOCATE(Basis(n), dBasisdx(n,3))
    
    DO Phase = 0,1
      
      tind = 0
      
      DO i=1,GetNOFActive()
        Element => GetActiveElement(i)
        
        IP = GaussPointsAdapt( Element )      
        CALL GetElementNodes( Nodes, UElement=Element )
        n  = GetElementNOFNodes(Element)     
        
        DO t=1,IP % n
          tind = tind + 1
          
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )
          Weight = IP % s(t) * DetJ

          ! Only populate the table the 2nd round
          IF( Phase == 1 ) THEN
            ALLOCATE(BasisFunctionsAtIp(tind) % Basis(n))
            BasisFunctionsAtIp(tind) % Basis(1:n) = Basis(1:n)
            ALLOCATE(BasisFunctionsAtIp(tind) % dBasisdx(n,3))
            BasisFunctionsAtIp(tind) % dBasisdx(1:n,1:3) = dBasisdx(1:n,1:3)
            BasisFunctionsAtIp(tind) % Weight = Weight          
          END IF
        END DO
      END DO

      ! Allocate for the basis functions when the size has been computed
      IF(Phase == 0) THEN
        ALLOCATE( BasisFunctionsAtIp(tind) )
      END IF
    END DO
      
    DEALLOCATE(Basis, dBasisdx)    

    BasisFunctionsInUse = .TRUE.
        
    CALL Info(Caller,'Number of tabulated basis functions:'//TRIM(I2S(tind)),Level=5)
    
  END SUBROUTINE TabulateBasisFunctions
  

!------------------------------------------------------------------------------
! This is monolithic lumping routine that has been optimized for speed.
! This way we need to evaluate the Basis functions only once for each element.
! It is assumed that inertial moment is computed the 1st time only and it
! stays constant.
!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumpedTransient()
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: torq,TorqArea,IMoment,IA, &
       rinner,router,rmean,rdiff,ctorq,detJ,Weight,&
       Bp,Br,Bx,By,x,y,r,rho
   REAL(KIND=dp), ALLOCATABLE :: a(:),u(:),POT(:),dPOT(:), &
       pPot(:),Density(:)
   REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
   LOGICAL, ALLOCATABLE :: TorqueElem(:)
   INTEGER :: i,bfid,n,nd,nbf,NoSlices,PrevComm
   LOGICAL :: Found, Stat
   TYPE(ValueList_t),POINTER::Params
   TYPE(GaussIntegrationPoints_t) :: IP
   TYPE(Nodes_t) :: Nodes
   LOGICAL :: CalcTorque, CalcPot, CalcInert
   LOGICAL :: ThisTorque, ThisPot, ThisInert, Parallel, HaveRange, SliceAverage
   LOGICAL :: Visited = .FALSE.
   
   SAVE Visited, Nodes, Basis, dBasisdx, a, u, POT, dPOT, pPot, Density, Ctorq, TorqueElem
   
!------------------------------------------------------------------------------

   CALL Info(Caller,'Calculating lumped parameters',Level=8)

   ! Define whether we have something to compute
   rinner = ListGetCRealAnyBody( Model,'r inner',CalcTorque )
   IF( CalcTorque ) THEN
     router = ListGetCRealAnyBody( Model,'r outer')
     rmean = (rinner+router)/2
     rdiff = (router-rinner)
     HaveRange = .TRUE.     
   ELSE
     rmean = ListGetConstReal( CurrentModel % Simulation,'Rotor Radius',CalcTorque)
     rmean = ParallelReduction( rmean, 2 )
     IF(.NOT. CalcTorque .AND. rmean > EPSILON(rmean) ) THEN
       CALL ListAddConstReal( CurrentModel % Simulation,'Rotor Radius',rmean)
       CalcTorque = .TRUE.
     END IF
     rdiff = ListGetConstReal( CurrentModel % Simulation,'Rotor Air Gap Width',Found)
     IF(.NOT. Found ) rdiff = 1.0e-3 * rmean
     HaveRange = .FALSE.
   END IF
   
   CalcPot = ListGetLogicalAnyBodyForce( Model,'Calculate Potential' )
   CalcInert = CalcTorque .AND. .NOT. Visited 
   
   IF(.NOT. (CalcTorque .OR. CalcPot .OR. CalcInert) ) RETURN

   Parallel = ( ParEnv % PEs > 1 ) .AND. (.NOT. Mesh % SingleMesh ) 

   NoSlices = ListGetInteger( Model % Simulation,'Number Of Slices', SliceAverage ) 
   IF( NoSlices > 1 .AND. NoSlices < ParEnv % PEs ) THEN
     CALL Info(Caller,'Changing communicator for slice operation!',Level=5)
     PrevComm = ParEnv % ActiveComm
     ParEnv % ActiveComm = ParallelSlicesComm() 
   END IF
   
   
   nbf = Model % NumberOfBodyForces
   IF(.NOT. Visited ) THEN
     n = Model % Mesh % MaxElementDofs
     ALLOCATE( a(nbf), u(nbf), POT(n), dPOT(n), pPot(n), Density(n), &
         Basis(n), dBasisdx(n,3) )
   END IF

   IF( CalcTorque ) THEN
     torq = 0._dp
     TorqArea = 0._dp
   END IF
   IF( CalcInert ) THEN
     IMoment = 0._dp
     IA = 0.0_dp
   END IF
   IF( CalcPot ) THEN   
     U=0._dp
     a=0._dp
   END IF

   tind = 0

   IF(.NOT. Visited .AND. CalcTorque ) THEN
     ALLOCATE( TorqueElem( GetNOFActive() ) )
     TorqueElem = .FALSE.
     
     DO i=1,GetNOFActive()
       Element => GetActiveElement(i)
     
       ThisTorque = .FALSE.
       
       n  = GetElementNOFNodes(Element)     
       CALL GetElementNodes( Nodes, Element )

       IF( HaveRange ) THEN       
         ! We are given range in classical Arkkio style. 
         ! Check how the center lies with respect to the range.
         x = SUM(Nodes % x(1:n))/n
         y = SUM(Nodes % y(1:n))/n
         r = SQRT(x**2+y**2)
         IF (r >= rinner .AND. r <= router) THEN
           TorqueElem(i) = .TRUE.
         END IF
       ELSE
         ! We are not given a range. Just take any element
         ! which has even one node at the given radius. 
         DO j=1,n
           x = Nodes % x(j)
           y = Nodes % y(j)
           r = SQRT(x**2+y**2)
           IF( ABS(r-rmean) < rdiff / 2 ) THEN
             TorqueElem(i) = .TRUE.
             EXIT
           END IF
         END DO
       END IF
     END DO
              
     i = COUNT( TorqueElem )
     CALL Info(Caller,'Number of elements to compute torque: '//TRIM(I2S(i)))
   END IF

   
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)
     
     ThisTorque = .FALSE.
     ThisPot = .FALSE.
     ThisInert = .FALSE.
     
     IF( CalcPot ) THEN
       Params => GetBodyForce(Element)
       IF(ASSOCIATED(Params)) THEN         
         ThisPot = GetLogical(Params,'Calculate Potential',Found)
         IF( ThisPot ) THEN
           bfid = GetBodyForceId(Element)           
           CALL GetLocalSolution(POT, UElement=Element)
           CALL GetLocalSolution(pPOT,tstep=-1,UElement=Element)
           IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
             dPot = (POT - pPOT)/dt
           ELSE
             dPot = 1.5_dp*POT - 2*pPOT
             CALL GetLocalSolution(pPOT,tstep=-2,UElement=Element)
             dPot = (dPOT + 0.5_dp*pPOT)/dt
           END IF
         END IF
       END IF
     END IF
     
     IF( CalcInert ) THEN
       Params=>GetBodyParams(Element)
       IF(ASSOCIATED(Params)) THEN
         ThisInert = GetLogical(Params,'Calculate Inertial Moment',Found)
       END IF
       Density(1:n) = GetReal(GetMaterial(),'Density',Found,Element)
     END IF
     
     IF( CalcTorque ) THEN
       ThisTorque = TorqueElem(i)
       IF(ThisTorque .AND. .NOT. ThisPot ) THEN
         CALL GetLocalSolution(POT, UElement=Element)
       END IF
     END IF

     ! Only treat the element if we have something to compute
     IF( .NOT. (ThisPot .OR. ThisInert .OR. ThisTorque ) ) THEN
       ! If nothing to compute still update the counter for IP points
       IF( BasisFunctionsInUse ) THEN
         IP = GaussPoints(Element)
         tind = tind + IP % n
       END IF
       CYCLE
     END IF
       
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)     
     CALL GetElementNodes( Nodes, Element )
     
     ! Numerical integration:
     !-----------------------
     IP = GaussPoints(Element)
     
     DO t=1,IP % n
       
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       IF( BasisFunctionsInUse ) THEN      
         tind = tind + 1
         Basis => BasisFunctionsAtIp(tind) % Basis
         dBasisdx => BasisFunctionsAtIp(tind) % dBasisdx        
         Weight = BasisFunctionsAtIp(tind) % Weight 
       ELSE IF( ThisTorque ) THEN
         ! Only torque needs the derivatives of basis function
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )
         weight = IP % s(t) * detJ
       ELSE
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis )
         weight = IP % s(t) * detJ
       END IF

       ! Coordinates of the intergration point
       x = SUM(Nodes % x(1:n)*Basis(1:n))
       y = SUM(Nodes % y(1:n)*Basis(1:n))
       r = SQRT(x**2+y**2)
       
       IF(ThisPot ) THEN
         A(bfid) = A(bfid) + Weight
         U(bfid) = U(bfid) + Weight * SUM(dPot(1:nd)*Basis(1:nd))
       END IF
       
       IF( ThisTorque ) THEN                      
         Bx =  SUM(POT(1:nd)*dBasisdx(1:nd,2))
         By = -SUM(POT(1:nd)*dBasisdx(1:nd,1))
         Br =  x/r*Bx + y/r*By
         Bp = -y/r*Bx + x/r*By
         Torq = Torq + Weight * r*Br*Bp / (PI*4.0d-7*rdiff)
         TorqArea = TorqArea + Weight
       END IF
       
       IF( ThisInert ) THEN
         IF( r < rmean ) THEN
           rho = SUM( density(1:n) * Basis(1:n) ) 
           IF( rho > EPSILON( rho ) ) THEN
             IA = IA + Weight
             U = U + Weight * r * rho
           END IF
         END IF
       END IF
     END DO
   END DO
     

   ! Finally perform parallel reduction if needed, and
   ! store the results for saving by SaveScalars.
   !-------------------------------------------------------------------------   
   IF( CalcPot ) THEN
     IF( Parallel ) THEN
       DO i=1,nbf
         a(i) = ParallelReduction(a(i))
         u(i) = ParallelReduction(u(i))
       END DO
     END IF     
     DO i=1,nbf
       IF(a(i)>0) THEN
         CALL ListAddConstReal(Model % Simulation,'res: Potential / bodyforce ' &
             //TRIM(i2s(i)),u(i)/a(i))
         CALL ListAddConstReal(Model % Simulation,'res: area / bodyforce ' &
             //TRIM(i2s(i)),a(i))
       END IF
     END DO
   END IF
   
   IF( CalcTorque ) THEN   
     IF( Parallel ) THEN
       Torq = ParallelReduction(Torq)
       TorqArea = ParallelReduction(TorqArea)
     END IF

     ! Arkkios formula assumes that rinner and router are nicely aligned with elements.
     ! This may not the case, so the 1st time we make a geometric correction. 
     IF(.NOT. Visited ) THEN
       WRITE(Message,'(A,ES15.4)') 'Air gap initial torque:', Torq
       CALL Info(Caller,Message,Level=6)

       IF (TorqArea /= 0) THEN
         !Ctorq = PI*(router**2-rinner**2) / TorqArea
         Ctorq = 2 * PI * rmean * rdiff / TorqArea

         WRITE(Message,'(A,F8.4)') 'Air gap correction initial:', cTorq
         CALL Info(Caller,Message,Level=4)
         
         ! The correction factor also corrects for the number of periods.
         ! We don't want that - so let us take back that and the torque
         ! can be compared to inertial moment of the sector still. 
         i = ListGetInteger( CurrentModel % Simulation,'Rotor Periods',Found )
         i = NINT( ParallelReduction( 1.0_dp * i, 2 ) )
         IF( i > 1 ) THEN
           WRITE(Message,'(A,I0)') 'Air gap correction rotor periods: ',i
           CALL Info(Caller,Message,Level=4)
           Ctorq = Ctorq / i 
         END IF         
       ELSE
         Ctorq = 1.0_dp
       END IF
       
       WRITE(Message,'(A,F8.4)') 'Air gap correction:', cTorq
       CALL Info(Caller,Message,Level=4)
       !CALL ListAddConstReal(Model % Simulation,'res: air gap correction', cTorq)
     END IF
       
     Torq = Ctorq * Torq

     IF( SliceAverage ) THEN
       ! Save slice torque even for one slice since then the output for scalars is the same
       ! for any number of slices.       
       WRITE(Message,'(A,ES15.4)') 'Air gap torque for slice'//TRIM(I2S(ParEnv % MyPe))//':', Torq
       CALL Info(Caller,Message,Level=5)
       CALL ListAddConstReal(Model % Simulation,'res: air gap torque for slice', Torq)

       ! But the averaging makes sense only for more than one slice
       IF( NoSlices > 1 ) THEN
         Torq = ParallelReduction(Torq) / NoSlices
       END IF
     END IF
       
     WRITE(Message,'(A,ES15.4)') 'Air gap torque:', Torq
     CALL Info(Caller,Message,Level=5)
     CALL ListAddConstReal(Model % Simulation,'res: air gap torque', Torq)
   END IF
   
   IF( CalcInert ) THEN
     IF( Parallel ) THEN
       IMoment = ParallelReduction(IMoment)
       IA = ParallelReduction(IA)
     END IF

     WRITE(Message,'(A,ES15.4)') 'Inertial volume:', IA
     CALL Info(Caller,Message,Level=7)

     WRITE(Message,'(A,ES15.4)') 'Inertial moment:', Imoment
     CALL Info(Caller,Message,Level=7)
       
     CALL ListAddConstReal(Model % Simulation,'res: inertial volume', IA)
     CALL ListAddConstReal(Model % Simulation,'res: inertial moment', IMoment)
   END IF

   IF( SliceAverage ) THEN
     IF( NoSlices > 1 .AND. NoSlices < ParEnv % PEs ) THEN
       CALL Info(Caller,'Reverting communicator from slice operation!',Level=10)
       ParEnv % ActiveComm = PrevComm
     END IF
   END IF
     
   Visited = .TRUE.
   
!------------------------------------------------------------------------------
 END SUBROUTINE CalculateLumpedTransient
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Old style local matrix. 
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Element, n, nd, Hsx, Hsy)
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: Hsx, Hsy
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueListEntry_t), POINTER :: Lst
    TYPE(ValueList_t), POINTER :: Material, BodyForce
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueList_t), POINTER :: CompParams

    REAL(KIND=dp), POINTER :: Bval(:), Hval(:), Cval(:)
    REAL(KIND=dp), POINTER :: CubicCoeff(:) => NULL(), HB(:,:) => NULL()

    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), &
        LOAD(nd),R(2,2,n),C(n), mu,muder,Babs,POT(nd), &
        JAC(nd,nd),Agrad(3),C_ip,M(2,n),M_ip(2),x, y,&
        Lorentz_velo(3,nd), Velo(3), omega_velo
    REAL(KIND=dp) :: Bt(nd,2), Ht(nd,2)
    REAL(KIND=dp) :: nu_tensor(2,2)
    REAL(KIND=dp) :: B_ip(2), Alocal, H_ip(2)

    INTEGER :: i,p,q,t,siz

    LOGICAL :: Cubic, HBcurve, WithVelocity, WithAngularVelocity, Found, Stat
    LOGICAL :: CoilBody, StrandedCoil    

!$omp threadprivate(Nodes, CubicCoeff, HB)
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType

    ! Zirka related
    LOGICAL :: Zirka
    TYPE(Variable_t), POINTER :: hystvar
    TYPE(GlobalHysteresisModel_t), pointer :: zirkamodel

    ! new parameters for reduced vector potential formulation
    REAL(KIND=dp) :: Hs(2,n), Hs_ip(2)

!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,Element )
    STIFF = 0._dp
    JAC  = 0._dp
    FORCE = 0._dp
    IF(Transient) MASS = 0._dp

    Material => GetMaterial(Element)

    CALL GetConstRealArray( Material, HB, 'H-B curve', HBCurve )
    Zirka = ListGetLogical(Material, 'Zirka material', Zirka)

    siz = 0
    Cval => NULL()
    IF ( HBCurve ) THEN
      siz = SIZE(HB,1)
      IF(siz>1) THEN
        Bval=>HB(:,1)
        Hval=>HB(:,2)
        Cubic = GetLogical( Material, 'Cubic spline for H-B curve',Found)
        IF (Cubic.AND..NOT.ASSOCIATED(CubicCoeff)) THEN
          ALLOCATE(CubicCoeff(siz))
          CALL CubicSpline(siz,Bval,Hval,CubicCoeff)
        END IF
        Cval=>CubicCoeff
        HBCurve = .TRUE.
      END IF
    END IF

    IF(siz<=1) THEN
      Lst => ListFind(Material,'H-B Curve',HBcurve)
      IF(HBcurve) THEN
        Cval => Lst % CubicCoeff
        Bval => Lst % TValues
        Hval => Lst % FValues(1,1,:)
      END IF
    END IF

    if (zirka) then
      CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
      zirkamodel => GetZirkaPointer(Material)
      hystvar => GetZirkaVariable(Material)
    end if

    IF(HBcurve) THEN
      CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
      IF (.NOT. ASSOCIATED(Bval) ) CALL Fatal (Caller,'Bval not associated')
      IF (.NOT. ASSOCIATED(Hval) ) CALL Fatal (Caller,'Hval not associated')
    ELSE
      CALL GetReluctivity(Material,R,n,Element)
    END IF

    C = GetReal( Material, 'Electric Conductivity', Found, Element)

    M(1,:) = GetReal( Material, 'Magnetization 1', Found, Element)
    M(2,:) = GetReal( Material, 'Magnetization 2', Found, Element)

    Load = 0.0d0

    WithVelocity = .FALSE.
    WithAngularVelocity = .FALSE.
    
    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) THEN
      Load(1:n) = GetReal(BodyForce, 'Current Density', Found, Element)
      CALL GetRealVector(BodyForce, Lorentz_velo, 'Lorentz velocity', WithVelocity)
      omega_velo = ListGetCReal(BodyForce, 'Angular velocity', WithAngularVelocity) 
    END IF

    ! setup magnetic field derived from the current source
    ! =====================================================
    Hs(1,:) = Hsx
    Hs(2,:) = Hsy

    CoilBody = .FALSE.
    StrandedCoil = .FALSE.
    CompParams => GetComponentParams( Element )
    IF (ASSOCIATED(CompParams)) THEN
      CoilType = GetString(CompParams, 'Coil Type', Found)
      IF (Found) THEN
        SELECT CASE (CoilType)
        CASE ('stranded')
          CoilBody = .TRUE.
          StrandedCoil = .TRUE.
        CASE ('massive')
          CoilBody = .TRUE.
        CASE ('foil winding')
          CoilBody = .TRUE.
!          CALL GetElementRotM(Element, RotM, n)
        CASE DEFAULT
          CALL Fatal ('StaticMagneticField2D', 'Non existent Coil Type Chosen!')
        END SELECT
      END IF
    END IF

    
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      IF ( CSymmetry ) THEN
        x = SUM( Basis(1:n) * Nodes % x(1:n) )
        detJ = detJ * x
      END IF

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      nu_tensor = 0.0_dp

      IF(Zirka .OR. HBCUrve) THEN
        Agrad = 0.0_dp
        Agrad = MATMUL( POT,dBasisdx )
        Alocal = SUM( POT(1:nd) * Basis(1:nd) )
        ! Sign? This convention: \vec A = A u_z
        ! -----
        B_ip(1) = Agrad(2) 
        B_ip(2) = -Agrad(1)
        IF( CSymmetry ) THEN
          B_ip = -B_ip
          B_ip(2) = B_ip(2) + Alocal/x
        END IF
      END IF

      IF (HBcurve ) THEN
        ! -----
        Babs = MAX( SQRT(SUM(B_ip**2)), 1.d-8 )
        mu = InterpolateCurve(Bval,Hval,Babs,CubicCoeff=Cval)/Babs
        muder = (DerivateCurve(Bval,Hval,Babs,CubicCoeff=Cval)-mu)/Babs
        nu_tensor(1,1) = mu ! Mu is really nu!!! too lazy to correct now...
        nu_tensor(2,2) = mu
      ELSEIF(Zirka) THEN
        CALL GetZirkaHBAtIP(t, solver, element, hystvar, zirkamodel, B_ip, H_ip, nu_tensor)
      ELSE
        muder=0._dp
        DO p=1,2
          DO q=1,2
            nu_tensor(p,q) = SUM(Basis(1:n) * R(p,q,1:n))
          END DO
        END DO
      END IF

      C_ip = SUM( Basis(1:n) * C(1:n) )
      M_ip = MATMUL( M,Basis(1:n) )

      ! setup Hs
      !Hs_ip = MATMUL( MATMUL(1.0-4*PI*1e-7*nu_tensor,Hs), Basis(1:n) )
      !Hs_ip = MATMUL( (1.0-mu*PI*4d-7)*Hs, Basis(1:n) )
      Hs_ip = MATMUL( Hs, Basis(1:n) )

      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      IF (Transient .AND. C_ip/=0._dp .AND. .NOT. StrandedCoil ) THEN
        DO p=1,nd
          DO q=1,nd
            MASS(p,q) = MASS(p,q) + IP % s(t) * detJ * C_ip * Basis(q)*Basis(p)
          END DO
        END DO
      END IF

      ! Is the sign correct?
      !---------------------
      Bt(1:nd,1) =  dbasisdx(1:nd,2)  ! (curlw)_x
      Bt(1:nd,2) = -dbasisdx(1:nd,1)  ! (curlw)_y
      IF ( CSymmetry ) THEN
        Bt(1:nd,1:2) = -Bt(1:nd,1:2)
        Bt(1:nd,2) = Bt(1:nd,2) + Basis(1:nd)/x
      END IF

      ! nu * (curlw)
      DO p = 1,nd
        Ht(p,:) = MATMUL(nu_tensor, Bt(p,:))
      END DO

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             MATMUL(Ht(1:nd,:), TRANSPOSE(Bt(1:nd,:)))

      ! Csymmetry is not yet considered in the Newton linearization
      IF (HBcurve .AND. NewtonRaphson) THEN
!        DO p=1,nd
!          DO q=1,nd
!            JAC(p,q) = JAC(p,q) + IP % s(t) * DetJ * &
!              muder/babs*SUM(Agrad*dBasisdx(q,:))*SUM(Agrad*dBasisdx(p,:))
!          END DO
!        END DO
        DO p=1,nd
          DO q=1,nd
            JAC(p,q) = JAC(p,q) + IP % s(t) * DetJ * &
              muder/babs*SUM(B_ip*Bt(q,:))*SUM(B_ip*Bt(p,:))
          END DO
        END DO
      END IF

      IF (WithVelocity .OR. WithAngularVelocity ) THEN
        !
        ! Create an additional Lorentz effect so that the electric field
        ! has an added term v x curl A:
        !
        IF( WithVelocity ) THEN        
          Velo(1:2) = [ SUM(Basis(1:n)*Lorentz_velo(1,1:n)), &
              SUM(Basis(1:n)*Lorentz_velo(2,1:n)) ]
        ELSE
          x = SUM( Basis(1:n) * Nodes % x(1:n) )
          y = SUM( Basis(1:n) * Nodes % y(1:n) ) 
          
          ! Simplified omega \times r in 2D
          Velo(1) = -omega_velo * y
          Velo(2) = omega_velo * x
        END IF
          
        IF (CSymmetry) THEN
          DO p=1,nd
            STIFF(p,1:nd) = STIFF(p,1:nd) + IP % s(t) * DetJ * C_ip * Basis(p) * ( & 
                -Velo(2) * Bt(1:nd,1) + Velo(1) * Bt(1:nd,2) )
          END DO
        ELSE
          DO p=1,nd
            STIFF(p,1:nd) = STIFF(p,1:nd) + IP % s(t) * DetJ * C_ip * Basis(p) * ( & 
                Velo(2) * Bt(1:nd,1) - Velo(1) * Bt(1:nd,2) )
          END DO
        END IF
      END IF

      IF ( CSymmetry ) THEN
        FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (LoadAtip * Basis(1:nd) - &
            M_ip(1)*dBasisdx(1:nd,2)+M_ip(2)*(dBasisdx(1:nd,1) + Basis(1:nd)/x))
      ELSE
        !FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (LoadAtip * Basis(1:nd) + &
            !M_ip(1)*dBasisdx(1:nd,2)-M_ip(2)*dBasisdx(1:nd,1))

        FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (LoadAtip * Basis(1:nd) + &
            (M_ip(1)+Hs_ip(1))*dBasisdx(1:nd,2)-(M_ip(2)+Hs_ip(2))*dBasisdx(1:nd,1))
      END IF
      IF(zirka) then
        FORCE(1:nd) = FORCE(1:nd) - (H_ip(1)*Bt(1:nd,1) + H_ip(2)*Bt(1:nd,2)) * IP % s(t) * detJ
      END IF
    END DO

    IF (HBcurve .AND. NewtonRaphson) THEN
      STIFF = STIFF + JAC
      FORCE = FORCE + MATMUL(JAC,POT)
    END IF

    IF(Zirka) THEN
      FORCE = FORCE + MATMUL(STIFF, POT)
    END IF

    IF(Transient) THEN
      CALL Default1stOrderTime( MASS, STIFF, FORCE,UElement=Element, USolver=Solver )
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE,UElement=Element, USolver=Solver)

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Assembly using handles. A little faster even for linear triangles, maybe 20%. 
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixHandles( Element, n, nd, nb, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER, SAVE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:), POT(:)    
    REAL(KIND=dp), POINTER, SAVE :: CVal(:),BVal(:),HVal(:)
    REAL(KIND=dp) :: Nu0, Nu, weight, SourceAtIp, CondAtIp, DetJ, Mu, MuDer, Babs
    LOGICAL :: Stat,Found, HBCurve
    INTEGER :: i,j,t,p,q,dim,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueListEntry_t), POINTER :: Lst
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueList_t), POINTER :: Material, PrevMaterial => NULL()
    REAL(KIND=dp) :: B_ip(2), Ht(nd,2), Bt(nd,2), Agrad(2), JAC(nd,nd), Alocal
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: StrandedCoil
    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, CondCoeff_h, PermCoeff_h, &
        RelPermCoeff_h, RelucCoeff_h, Mag1Coeff_h, Mag2Coeff_h, &
        CoilType_h
    
    SAVE HBCurve, Nu0, PrevMaterial
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Current Density')
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
      CALL ListInitElementKeyword( PermCoeff_h,'Material','Permeability')
      CALL ListInitElementKeyword( RelPermCoeff_h,'Material','Relative Permeability')
      CALL ListInitElementKeyword( RelucCoeff_h,'Material','Reluctivity')
      CALL ListInitElementKeyword( Mag1Coeff_h,'Material','Magnetization 1')
      CALL ListInitElementKeyword( Mag2Coeff_h,'Material','Magnetization 2')
      Found = .FALSE.
      IF( ASSOCIATED( Model % Constants ) ) THEN
        Nu0 = ListGetCReal( Model % Constants,'Permeability of Vacuum',Found)
      END IF
      IF( .NOT. Found ) Nu0 = PI * 4.0d-7
      InitHandles = .FALSE.
      CALL ListInitElementKeyword( CoilType_h,'Component','Coil Type')
    END IF

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(MASS)) THEN
      m = Mesh % MaxElementDofs
      ALLOCATE(MASS(m,m), STIFF(m,m),FORCE(m), POT(m), STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
      IF(.NOT. BasisFunctionsInUse ) THEN
        ALLOCATE(Basis(m), dBasisdx(m,3))
      END IF
    END IF
      
    Material => GetMaterial(Element)
    IF( .NOT. ASSOCIATED( Material, PrevMaterial ) ) THEN
      PrevMaterial => Material           
      Lst => ListFind(Material,'H-B Curve',HBcurve)
      IF(HBcurve) THEN
        Cval => Lst % CubicCoeff
        Bval => Lst % TValues
        Hval => Lst % FValues(1,1,:)
        IF (.NOT. ASSOCIATED(Bval) ) CALL Fatal (Caller,'Bval not associated')
        IF (.NOT. ASSOCIATED(Hval) ) CALL Fatal (Caller,'Hval not associated')
      END IF      
    END IF


    StrandedCoil = .FALSE.
    CoilType = ListGetElementString(CoilType_h, Element, Found ) 
    IF( Found ) THEN
      IF( CoilType == 'stranded' ) THEN
        StrandedCoil = .TRUE.
      ELSE
        CALL Fatal(Caller,'Implemented only for stranded coils for now!')
      END IF
    END IF
        
    ! Initialize
    MASS  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp
    
    ! Integration rule
    IP = GaussPointsAdapt( Element )
      
    CALL GetElementNodes( Nodes, UElement=Element )

    IF(HBcurve) THEN
      CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
      JAC = 0.0_dp
    END IF
    
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      IF( BasisFunctionsInUse ) THEN      
        tind = tind + 1
        Basis => BasisFunctionsAtIp(tind) % Basis
        dBasisdx => BasisFunctionsAtIp(tind) % dBasisdx        
        Weight = BasisFunctionsAtIp(tind) % Weight 
      ELSE
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )
        Weight = IP % s(t) * DetJ
      END IF
        
      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      IF( HBCurve ) THEN
        Agrad(1:2) = MATMUL( POT(1:nd),dBasisdx(1:nd,1:2) )
        Alocal = SUM( POT(1:nd) * Basis(1:nd) )

        B_ip(1) = Agrad(2) 
        B_ip(2) = -Agrad(1)         
        Babs = MAX( SQRT(SUM(B_ip**2)), 1.d-8 )

        mu = InterpolateCurve(Bval,Hval,Babs,CubicCoeff=Cval)/Babs
        muder = (DerivateCurve(Bval,Hval,Babs,CubicCoeff=Cval)-mu)/Babs
      ELSE
        Nu = ListGetElementReal( RelPermCoeff_h, Basis, Element, Found, GaussPoint = t )
        IF( Found ) THEN
          Mu = 1.0_dp / (Nu0 * Nu)
        ELSE
          Nu = ListGetElementReal( PermCoeff_h, Basis, Element, Found, GaussPoint = t )
          IF( Found ) THEN
            Mu = 1.0_dp / Nu
          ELSE
            Mu = ListGetElementReal( RelucCoeff_h, Basis, Element, Found, GaussPoint = t )
          END IF
          IF(.NOT. Found ) THEN
            CALL Fatal(Caller,'Could not define reluctivity in any way!')
          END IF
        END IF
      END IF


      Bt(1:nd,1) =  dbasisdx(1:nd,2)
      Bt(1:nd,2) = -dbasisdx(1:nd,1)

      ! Here istrophy is assumed!
      Ht(1:nd,:) = mu * Bt(1:nd,:)
           
      IF ( HBCurve .AND. NewtonRaphson) THEN
        DO p=1,nd
          DO q=1,nd
            JAC(p,q) = JAC(p,q) + Weight * &
                muder/babs * SUM(B_ip(:) * Bt(q,:)) * SUM(B_ip(:)*Bt(p,:))
          END DO
        END DO
      END IF
            
      ! diffusive term: STIFF=STIFF+(a*grad(u),grad(v))   
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
          MATMUL(Ht(1:nd,:), TRANSPOSE(Bt(1:nd,:)))

      IF( MassAsm .AND. .NOT. StrandedCoil ) THEN
        CondAtIp = ListGetElementReal( CondCoeff_h, Basis, Element, Found )
        IF( Found ) THEN
          DO p=1,nd
            MASS(p,1:nd) = MASS(p,1:nd) + Weight * CondAtIp * Basis(1:nd) * Basis(p)
          END DO
        END IF
      END IF

      ! Current density source 
      SourceAtIP = ListGetElementReal( SourceCoeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * Basis(1:nd)
      END IF

      ! Magnetization source, weak form
      SourceAtIP = ListGetElementReal( Mag1Coeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * dBasisdx(1:nd,2)        
      END IF
      SourceAtIP = ListGetElementReal( Mag2Coeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) - Weight * SourceAtIP * dBasisdx(1:nd,1)        
      END IF     
    END DO

    IF (HBcurve .AND. NewtonRaphson) THEN
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + JAC(1:nd,1:nd)
      FORCE(1:nd) = FORCE(1:nd) + MATMUL(JAC(1:nd,1:nd),POT(1:nd))
    END IF
    
    IF( MassAsm ) THEN
      CALL DefaultUpdateMass(MASS,UElement=Element)
    END IF
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element) !, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixHandles
!------------------------------------------------------------------------------

  
!-------------------------------------------------------------------------------
! Calculates H and dHdB in 2D given B. This should be always inlined in LocalMatrix.
!-------------------------------------------------------------------------------
SUBROUTINE GetZirkaHBAtIP(i_IP, Solver, Element, HystVar, ZirkaModel, B_ip, H_ip, dHdB) ! {{{
!-------------------------------------------------------------------------------
  INTEGER, intent(in) :: i_IP
  TYPE(Solver_t) :: Solver
  TYPE(Element_t) :: Element
  TYPE(Variable_t), POINTER :: HystVar
  TYPE(GlobalHysteresisModel_t), POINTER :: ZirkaModel
  REAL(KIND=dp), INTENT(IN) :: B_ip(2)
  REAL(KIND=dp), INTENT(OUT) :: H_ip(2)
  REAL(KIND=dp), intent(INOUT) :: dHdB(2,2)
!-------------------------------------------------------------------------------
  INTEGER :: ipindex, n_dir, k,l
  REAL(KIND=dp) :: dH, B0(3)
!-------------------------------------------------------------------------------
  ipindex = getipindex(i_IP, usolver=solver, element=element, ipvar=hystvar)
  IF (ipindex == 0 ) RETURN

  H_ip = 0.0_dp
  DO n_dir = 1, UBOUND(zirkamodel % curves, 1)
    B0 = zirkamodel % curves(n_dir, ipindex) % B0
    ASSOCIATE(Bdir => SUM(B_ip*B0(1:2)))
      ! H_ip(1:2) = H_ip(1:2) + zirkamodel % curves(n_dir,ipindex) % &
      !     eval(sum(B_ip*B0(1:2)), cached = .true., dhdb=dH) * &
      !     B0(1:2)
      H_ip(1:2) = H_ip(1:2) + zirkamodel % curves(n_dir,ipindex) % &
          eval(Bdir, cached = .TRUE., dhdb=dH) * &
          B0(1:2)
    END ASSOCIATE
    DO k = 1,2
      DO l = 1,2
        dHdB(k,l) = dHdB(k,l) + dH*B0(k)*B0(l)
      END DO
    END DO
  END DO
  
END SUBROUTINE ! }}}
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixInfinityBC(Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), R(2,2,n), R_ip, &
            Inf_ip,Coord(3),Normal(3),mu,u,v

    TYPE(ValueList_t), POINTER :: Material

    TYPE(Element_t), POINTER :: Parent
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )
    STIFF = 0._dp
    FORCE = 0._dp

    Parent=>Element % BoundaryInfo % Left
    IF(.NOT.ASSOCIATED(Parent)) THEN
      Parent=>Element % BoundaryInfo % Right
    END IF
    IF(.NOT.ASSOCIATED(Parent)) RETURN

    Material => GetMaterial(Parent)
    CALL GetReluctivity(Material,R,n,Parent)

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                 IP % W(t), detJ, Basis )

      mu = SUM(Basis(1:n)*R(1,1,1:n)) ! We assume isotropic reluctivity here.

      Normal = NormalVector( Element, Nodes, u, v, .TRUE. )
      Coord(1) = SUM(Basis(1:n) * Nodes % x(1:n))
      Coord(2) = SUM(Basis(1:n) * Nodes % y(1:n))
      Coord(3) = SUM(Basis(1:n) * Nodes % z(1:n))

      IF( CSymmetry ) THEN
        detJ = detJ * Coord(1)
      END IF

      Inf_ip = mu * SUM(Coord*Normal)/SUM(Coord*Coord)

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + IP % s(t)*detJ*Inf_ip*Basis(q)*Basis(p)
        END DO
      END DO
    END DO
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixInfinityBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixAirGapBC(Element, BC, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat, Found
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), R(n), R_ip, &
            Inf_ip,Coord(3),Normal(3),mu,u,v, AirGapLength(nd), &
            AirGapMu(nd), AirGapL

    TYPE(ValueList_t), POINTER :: Material
    TYPE(ValueList_t), POINTER :: BC

    TYPE(Element_t), POINTER :: Parent
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )
    STIFF = 0._dp
    FORCE = 0._dp

    AirGapLength = ListGetConstReal( BC, 'Air Gap Length', UnfoundFatal = .TRUE.)

    AirGapMu = ListGetConstReal( BC, 'Air Gap Relative Permeability', Found)
    IF (.NOT. Found) AirGapMu=1.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      mu = 4*pi*1d-7*SUM(Basis(1:n)*AirGapMu(1:n))
      AirGapL = SUM(Basis(1:n)*AirGapLength(1:n))

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
          AirGapL/mu*MATMUL(dBasisdx, TRANSPOSE(dBasisdx))

    END DO
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixAirGapBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SetMagneticFluxDensityBC()
!------------------------------------------------------------------------------
! P. Lombard, G. Meunier, "A general purpose method for electric and magnetic 
! combined problems for 2D, axisymmetric and transient systems", IEEE Trans.
! magn. 29(2), p. 1737 - 1740, Mar 1993
! -ettaka- 
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: A
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), POINTER :: b(:)
    INTEGER :: i, n, j, k
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found
    REAL(KIND=dp) :: Bx(Solver % Mesh % MaxElementDofs), &
        By(Solver % Mesh % MaxElementDofs)
    REAL(KIND=dp) :: x, y
    INTEGER, POINTER :: Perm(:)

    Perm => Solver % Variable % Perm
    A => Solver % Matrix
    b => A % RHS
    
    DO i=1,GetNofBoundaryElements()
      Element => GetBoundaryElement(i)
      n = GetELementNofNodes()
      BC => GetBC()
      IF ( ASSOCIATED(BC)) THEN
        IF ( ListCheckPrefix( BC, 'Magnetic Flux Density') ) THEN
          Bx = 0._dp
          By = 0._dp

          Bx(1:n) = GetReal(BC, 'Magnetic Flux Density 1', Found)
          By(1:n) = GetReal(BC, 'Magnetic Flux Density 2', Found)

          DO j = 1,n
            k = Element % NodeIndexes(j)
            x = Mesh % Nodes % x(k)
            y = Mesh % Nodes % y(k)
            k = Perm(k)

            CALL UpdateDirichletDof( A, k, y * Bx(j) - x * By(j) )
          END DO 
        END IF  
      END IF  
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SetMagneticFluxDensityBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetSourceField( x, y, curr, nmax, Hsx, Hsy )
!------------------------------------------------------------------------------
! Calculate magetic field derived from a bundle of wires
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE

  REAL(KIND=dp) :: x, y, Hsx, Hsy, curr
  INTEGER :: nmax
!------------------------------------------------------------------------------
  INTEGER :: i
  INTEGER, PARAMETER :: nwire = 2140 
  !REAL(KIND=dp) :: curr = 23
  REAL(KIND=dp) :: xs(nwire), ys(nwire)

  LOGICAL :: FirstVisit = .TRUE.
  SAVE FirstVisit, xs, ys
!------------------------------------------------------------------------------

  ! do something only the first time the function is visited
  IF (FirstVisit) THEN
    OPEN(1,file='CoilWirePosition.txt',status='old',action='read')
    !WRITE(Message,'(2A15)') 'X[m]', 'Y[m]'
    !PRINT *, 'X[m]', 'Y[m]'

    DO i = 1, nwire
      READ(1,*) xs(i), ys(i)
      !WRITE(Message,'(2E15.6)') xs(i), ys(i)
      !PRINT *, xs(i), ys(i)
    END DO

    FirstVisit = .FALSE.
  END IF

  ! calculate magnetic field
  Hsx = 0.0d0
  Hsy = 0.0d0

  DO i = 1, nwire
    CALL GetFieldFromLineCurrent( xs(i), ys(i), x, y, curr, nmax, Hsx, Hsy )  
  END DO

!------------------------------------------------------------------------------
 END SUBROUTINE GetSourceField
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetArcSourceField( x, y, phi1, phi2, curr, Hsx, Hsy )
!------------------------------------------------------------------------------
! Calculate magetic field derived from a bundle of wires
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  REAL(KIND=dp) :: x, y, Hsx, Hsy, curr, phi1, phi2
!------------------------------------------------------------------------------
  INTEGER :: i
  INTEGER, PARAMETER :: NWIRE = 2140 
  INTEGER, PARAMETER :: NGAUSS= 38
  REAL(KIND=dp), PARAMETER :: RC = 1.890
  
  REAL(KIND=dp) :: xs(NWIRE), ys(NWIRE), xi(NGAUSS), wi(NGAUSS)
  REAL(KIND=dp) :: ps, dhr1, dhz1, dhr2, dhz2, Iop, hr, hz

  LOGICAL :: FirstVisit = .TRUE.
  SAVE FirstVisit, xs, ys, xi, wi
!------------------------------------------------------------------------------

  ! do something only the first time the function is visited
  IF (FirstVisit) THEN
    ! load conductor data
    OPEN(1,file='CoilWirePosition.txt',status='old',action='read')
    PRINT *, 'LOADING WIRE DATA'
    DO i = 1, NWIRE
      READ(1,*) xs(i), ys(i)
    END DO
    PRINT *, 'FINISHED LOADING WIRE DATA'

    ! calculate gauss points and weight
    PRINT *, 'CALCULATE GAUSS POINTS AND WEIGHTS'
    CALL  gauleg(NGAUSS, xi, wi)
    PRINT *, 'NUMBER OF GAUSS POINTS:', NGAUSS
    DO i = 1, NGAUSS
      PRINT *, xi(i), wi(i)
    END DO

    FirstVisit = .FALSE.
  END IF

  ! calculate magnetic field
  Hsx = 0.0d0
  Hsy = 0.0d0

  DO i = 1, NWIRE
    IF (xs(i).GT.0) THEN
        Iop =  1.0*curr
    ELSE
        Iop = -1.0*curr
    END IF

    hr = 0.0d0
    hz = 0.0d0

    ! gauss integral
    DO j=1, NGAUSS/2
      ps = (phi2-phi1)*0.5*xi(j) 
      
      CALL GetFieldFromArcCurrent( xs(i)+RC, ps,  ys(i), x+RC, y, dhr1, dhz1 ) 
      ! fliped wires to mirror the field 
      CALL GetFieldFromArcCurrent( xs(i)+RC, ps, -ys(i), x+RC, y, dhr2, dhz2 ) 

      hr = hr + wi(j)*(dhr1+dhr2) 
      hz = hz + wi(j)*(dhz1+dhz2)
    END DO
    
    Hsx = Hsx + hr*Iop
    Hsy = Hsy + hz*Iop
  END DO

  ! double and integral half of the field to accelerate calculation
  Hsx = Hsx * (phi2-phi1) / (4.0*PI)
  Hsy = Hsy * (phi2-phi1) / (4.0*PI)

!------------------------------------------------------------------------------
 END SUBROUTINE GetArcSourceField
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetFieldFromArcCurrent(rs, ps, zs, r, z, Hr, Hz )
!------------------------------------------------------------------------------
! Calculate magnetic field strength, Hs, from a line current source
! Reference: R.Gupta, Field Calculations and Computations
!
! rs, ps, zs: source point of line current in cylinderical coordinate
! r , p , z : observe point
! Hr, Hp, Hz: parameters of magnetic flux density
!
! Author: Ye Yang
! Date  : 11.24.2021
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: rs, ps, zs, r, z, Hr, Hz
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: gam, fi, parD
!------------------------------------------------------------------------------
    ! parameters
    gam = zs - z
    fi  = ps 
    parD= gam**2 + rs**2 + r**2 - 2.0*rs*r*COS(fi)

    ! calculate each component of field
    Hr  = -rs * parD**(-1.5) * gam * COS(fi)
    !Hp  = -rs * parD**(-1.5) * gam * SIN(fi)
    Hz  =  rs * parD**(-1.5) * (rs - r*COS(fi))

!------------------------------------------------------------------------------
 END SUBROUTINE GetFieldFromArcCurrent
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetFieldFromLineCurrent(xs, ys, x, y, curr, nmax, Hx, Hy )
!------------------------------------------------------------------------------
! Calculate magnetic field strength, Hs, from a line current source
! Reference: R.Gupta, Field Calculations and Computations
!
! xs, ys: source point of line current
! x , y : observe point
! curr  : operation current
! nmax  : maximum order of field calculation related to the precision
! Hx, Hy: magnetic field strength
!
! Author: Ye Yang
! Date  : 11.11.2021
!------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=dp) :: x, y, Hx, Hy, Hr, Ht, a, r, theta, phi, curr, Iop
    REAL(KIND=dp) :: xs, ys
    INTEGER :: nmax, n
!------------------------------------------------------------------------------
    Hr = 0.0d0
    Ht = 0.0d0

    ! coordinate transformation
    r    = SQRT(  x**2 +  y**2 )
    a    = SQRT( xs**2 + ys**2 )
    theta= DATAN2( y, x)
    !phi  = DATAN2(ys,xs)

    ! change the direction of current when x > 0
    IF (xs.GT.0) THEN
      Iop = -1 * curr
    ELSE
      Iop =  1 * curr
    END IF

    ! calculate field
    IF (r.LT.a) THEN
      phi = DATAN2(ys,xs)
      DO n = 1, nmax
        Hr  = Hr + (r/a)**(n-1) * SIN(n*(phi-theta))
        Ht  = Ht + (r/a)**(n-1) * COS(n*(phi-theta))
      END DO

      ! mirror the line current
      phi = DATAN2(-ys,xs)
      DO n = 1, nmax
        Hr  = Hr + (r/a)**(n-1) * SIN(n*(phi-theta))
        Ht  = Ht + (r/a)**(n-1) * COS(n*(phi-theta))
      END DO

      Hr =  Iop / (2.0*PI*a) * Hr
      Ht = -Iop / (2.0*PI*a) * Ht
    END IF

    IF (r.GT.a) THEN
      phi = DATAN2(ys,xs)
      DO n = 0, nmax
        Hr = Hr + (a/r)**(n+1) * SIN(n*(phi-theta))
        Ht = Ht + (a/r)**(n+1) * COS(n*(phi-theta))
      END DO

      ! mirror the line current
      phi = DATAN2(-ys,xs)
      DO n = 0, nmax
        Hr = Hr + (a/r)**(n+1) * SIN(n*(phi-theta))
        Ht = Ht + (a/r)**(n+1) * COS(n*(phi-theta))
      END DO

      Hr = Iop / (2.0*PI*a) * Hr
      Ht = Iop / (2.0*PI*a) * Ht
    END IF

    ! coordinate transformation for magnetic field
    IF (r.NE.a) THEN
      Hx = Hx + Hr*COS(theta) - Ht*SIN(theta)
      Hy = Hy + Hr*SIN(theta) + Ht*COS(theta) 
    END IF

!------------------------------------------------------------------------------
 END SUBROUTINE GetFieldFromLineCurrent
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivity(Material,Acoef,n,Element)
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(2,2,n)
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp), SAVE :: Avacuum
    LOGICAL :: Found
    LOGICAL, SAVE :: FirstTime = .TRUE.
    !$OMP THREADPRIVATE(Avacuum, FirstTime)
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF(.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF

    Acoef = GetTensor(Element, n, 2, 'Relative Permeability', 're', Found)

    IF ( Found ) THEN
      Acoef = Avacuum * Acoef
    ELSE
      Acoef = GetTensor(Element, n, 2, 'Permeability', 're', Found)
    END IF
    IF ( Found ) THEN
      Acoef = Get2x2TensorInverse(Acoef, n)
    ELSE
      Acoef = GetTensor(Element, n, 2, 'Reluctivity', 're', Found)
    END IF

    IF( .NOT. Found ) THEN
      CALL Warn('GetReluctivity',&
          'Could not get either > Reluctivity > or > Relative Permeability < !')
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivity
!------------------------------------------------------------------------------
  
!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula.  For detailed explanations finding weights & abscissas, see
!* "Numerical Recipes in Fortran */
!********************************************************************************
  SUBROUTINE  gauleg(ngp, xabsc, weig)
    !implicit none
    INTEGER  i, j, m
    REAL(KIND=dp) :: p1, p2, p3, pp, z, z1
    INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
    REAL(KIND=dp), INTENT(OUT) :: xabsc(ngp), weig(ngp)
    REAL(KIND=dp), PARAMETER :: EPS = 1.0d-15, M_PI = 3.141592654d0

    m = (ngp + 1) / 2

!* Roots are symmetric in the interval - so only need to find half of them  */
    do i = 1, m              ! Loop over the desired roots */

      z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */
100   p1 = 1.0d0
      p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */
        do j = 1, ngp
          p3 = p2
          p2 = p1
          p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
      pp = ngp*(z*p1-p2)/(z*z-1.0d0)
      z1 = z
      z  = z1 - p1/pp             ! Newton's Method  */

      if (abs(z-z1) .gt. EPS) GOTO  100

      xabsc(i) =  - z                     ! Roots will be bewteen -1.0 & 1.0 */
      xabsc(ngp+1-i) =  + z               ! and symmetric about the origin  */
      weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
      weig(ngp+1-i) = weig(i)             ! symmetric counterpart         */

    end do     ! i loop

  END SUBROUTINE gauleg

END SUBROUTINE StaticMagneticField2D
!------------------------------------------------------------------------------

