MODULE DriftDiffusion

    IMPLICIT NONE
    
    ! Some constants
    REAl(KIND=16), PARAMETER :: q     = 1.60217646e-19    !Electron charge (C)
    REAl(KIND=16), PARAMETER :: kb     = 1.3806503e-23        !Boltzmann constant (m2 Kg s-2 K-1)
    REAl(KIND=16), PARAMETER :: m0     = 9.10938188e-31    !Electron rest mass (kg)
    REAl(KIND=16), PARAMETER :: Epsi0 = 8.854187817e-12    !Vacuum permitivity (F/m)
    REAl(KIND=16), PARAMETER :: hp     = 6.626068e-34        !Planck's constant (m^2 Kg/s)
    REAl(KIND=16), PARAMETER :: Pi     = 3.14159265359
    !
    ! ---------------------------------------------------------------------------
    ! Everything included here are global variables available in the whole module
    ! ---------------------------------------------------------------------------
    ! Some useful variables
    REAl(KIND=16) :: T = 300                             !Temperature. Default room temperature.
    ! 
    ! Variable inputs and/or outputs. They will have M+1 elements.
    REAl(KIND=16), DIMENSION(0:6000) :: X                !Node possition (m)
    REAl(KIND=16), DIMENSION(0:6000) :: dX            !Node spacing (m)
    REAl(KIND=16), DIMENSION(0:6000) :: n, p            !Electron and hole densities (m-3)
    REAl(KIND=16), DIMENSION(0:6000) :: Rho            !Total charge density Rho = Nd+p-Nd-n (m-3)
    REAl(KIND=16), DIMENSION(0:6000) :: ni               !Carrier intrinsic densities (m-3)
    REAl(KIND=16), DIMENSION(0:6000) :: Nc, Nv        !Total effective density of states of electrons and holes (m-3)
    REAl(KIND=16), DIMENSION(0:6000) :: Nd, Na        !Density of ionised donors and acceptors (m-3).
    REAl(KIND=16), DIMENSION(0:6000) :: Fn, Fp        !Quasi-Fermi potential for electrons and holes (V)
    REAl(KIND=16), DIMENSION(0:6000) :: Psi            !Electrostatic potential (V)
    REAl(KIND=16), DIMENSION(0:6000) :: Eg            !Energy gap (eV)
    REAl(KIND=16), DIMENSION(0:6000) :: Xi            !Electron afinity (eV)
    REAl(KIND=16), DIMENSION(0:6000) :: Mun, Mup        !Mobilities of electrons and holes (m^2/Vs)
    REAl(KIND=16), DIMENSION(0:6000) :: Epsi            !Relative permitivity (-)
    REAl(KIND=16), DIMENSION(0:6000) :: Ncc, Nvhh, Nvlh        !Effective density of states of electrons and holes (m-3)
    REAl(KIND=16), DIMENSION(0:6000) :: tn, tp          !Lifetime of minority carriers in the SRH model
    REAl(KIND=16), DIMENSION(0:6000) :: Brad            !Radiative recombination coeficient
    REAl(KIND=16), DIMENSION(0:6000) :: CCn, CCp        !Auger recombination coeficients
    REAl(KIND=16), DIMENSION(0:6000) :: alfa            !Absorption coefficient. 
    REAl(KIND=16), DIMENSION(0:6000, 0:3000) :: AbsProfile    !Absorption coefficient as a function of wavelength.
    REAl(KIND=16), DIMENSION(0:6000) :: IQE, IQEsrh, IQErad, IQEaug, IQEsurb, IQEsurf ! Internal quantum efficiency of the device as a function of wavelength 
    !
    ! Some derived potentials useful for the calculation
    REAl(KIND=16), DIMENSION(0:6000) :: Vn, Vp        !Band edge potentials with respect certain reference
    REAl(KIND=16), DIMENSION(0:6000) :: Cn, Cp        !Modified electric potentials
    !
    !
    ! Bulk generation and recombination, including all processes
    REAl(KIND=16), DIMENSION(0:6000) :: GR            ! Generation-Recombination = Rsrh + Rrad + Raug - G
    REAl(KIND=16), DIMENSION(0:6000) :: Rrad            ! Radiative recombination    
    REAl(KIND=16), DIMENSION(0:6000) :: Rsrh            ! SRH recombinaiton
    REAl(KIND=16), DIMENSION(0:6000) :: Raug            ! Auger recombinaiton
    REAl(KIND=16), DIMENSION(0:6000) :: G                ! Generation
    REAl(KIND=16), DIMENSION(0:6000) :: vpoint            ! Voltage in an IV curve
    REAl(KIND=16), DIMENSION(0:6000) :: jpoint            ! Total current in an IV curve
    REAl(KIND=16), DIMENSION(0:6000) :: jsrhpoint        ! SRH current in an IV curve
    REAl(KIND=16), DIMENSION(0:6000) :: jradpoint        ! Radiative current in an IV curve
    REAl(KIND=16), DIMENSION(0:6000) :: jaugpoint        ! Auger current in an IV curve
    REAl(KIND=16), DIMENSION(0:6000) :: jsurpoint        ! Surface recombination current in an IV curve
    REAl(KIND=16), DIMENSION(0:6000) :: residual        ! residual in an IV curve
    INTEGER :: nvolt = 0

    REAl(KIND=16) :: PhotonFlux                        ! Photon flux
    REAl(KIND=16), DIMENSION(0:3000) :: PFspectrum = 0.0    ! Photon flux as a function of wavelength (same wl that the abs. coef.)
    INTEGER :: SRH, RAD, AUG, GEN                ! 1 = Included, 0 = Not included
    REAl(KIND=16) :: Jtot2, CurrentsBias(6), Currents(6)
    LOGICAL :: SingleWL = .FALSE.                ! Controls the generation in IQE running mode
    LOGICAL :: Dynamic = .FALSE.                ! Controls if there is dynamic meshing or not
    !
    ! Extra values for the boundary conditions    
!     REAl(KIND=16) :: Sn, Sp                            !Surface recombination velocity of minority carriers
    REAl(KIND=16) :: Snfront, Spfront, Snback, Spback    !Surface recombination velocity of minority carriers
    REAl(KIND=16) :: fneq = 1
    REAl(KIND=16) :: bneq = 1
    REAl(KIND=16) :: fpeq = 1
    REAl(KIND=16) :: bpeq = 1
    INTEGER :: FTYPE, BTYPE, FSUR, BSUR
    REAl(KIND=16) :: Vbarf, Vbarb
    INTEGER :: EQ, SC, OC, OCn, OCp
    !
    ! Reference values, tipically those at x = 0, on the left end of the device
    REAl(KIND=16) :: nir
    REAl(KIND=16) :: Munr, Mupr
    REAl(KIND=16) :: Xir
    REAl(KIND=16) :: Ncr, Nvr
    REAl(KIND=16) :: Egr
    REAl(KIND=16) :: Epsir
    !
    ! Doping in the device. We start asuming a simple pin junction
    REAl(KIND=16) :: Aceptors, Intrinsic, Donors        ! Doping of the p, i and n regions.  
    REAl(KIND=16) :: XD = 0.0
    !
    ! Other variables
    INTEGER :: M                                 ! The number of nodes -1
    REAl(KIND=16) :: MasterNodes(1000) = 0                ! Array containing the position of the Masternodes
    REAl(KIND=16) :: DML(1:1000, 20)                     ! DeviceMaterialsLibrary, array containing all the properties of the materials 
                                                ! used in the device.
    REAl(KIND=16) :: DoppingLibrary(200, 4)            ! An array containing all the constant doping profiles used in the device.
    REAl(KIND=16) :: AbsLibrary(-1:1000, 0:3000) = 0.0    ! An array containing all the absorption profiles used in the device.
    INTEGER :: MGrid = 1                         ! The number of grid lines (Max MGrid=200)
    INTEGER :: MReg  = 0                         ! The number of different material regions (Max MReg=200)

    ! Mesh variables
    INTEGER :: NumWL  = 2                         ! The number of wavelengths in the photon flux and the absroption coefficient
    REAl(KIND=16) :: Coarse, Fine, Ultrafine            ! The different mesh sizes
    REAl(KIND=16) :: Growth                            ! Growth parameter for the dynamic meshing
    ! The clamp for the variables Fn, Psi and Fp.
    REAl(KIND=16) :: clamp = 20
    REAl(KIND=16) :: ATol = 3.1622776601683796e-17    ! SQRT of machine epsilon at quadruple precission
    REAl(KIND=16) :: RTol = 1e-6
    INTEGER :: nitermax = 40
    !
    ! Voltage, current an series resistance information
    REAl(KIND=16) :: Vbi, Vi, Vap                            ! Built-in voltage, Vi = q*Vbi/kbT, applied voltage (used only in dep.aprox.)
    REAl(KIND=16) :: Voc, Isc, Vmax, Imax, Pmax, FF
    REAl(KIND=16) :: Rs = 0
    !
    ! Set of equations to be solved simultaneously [f]=0. For each internal node k = 1, M-2:
    !     - f(3k-1) corresponds to the continuty of Jp, associated to Fp
    !     - f(3k) corresponds to the Poisson equation, associated to Psi
    !     - f(3k+1) corresponds to the continuity of Jn, associated to Fn
    ! The total vector of equations with 3M-1 elements and auxiliary vector
    REAl(KIND=16), DIMENSION(18003) :: f, dsol
    ! The Jacobian matrix in compact form. It only contains the non-zero elements
    REAl(KIND=16), DIMENSION(18003,11) :: Jac
    !
    !
    ! Scaling factors
!     REAl(KIND=16) :: x0                                ! Max length scale x0 = total device thickness
    REAl(KIND=16) :: b                                ! Inverse of thermal voltage b = q/(kb*T)
    REAl(KIND=16) :: C0                                ! Maximum intrinsic concentration
    REAl(KIND=16) :: Mu0                                ! Maximum mobility
    REAl(KIND=16) :: D0                                ! D0 = Mu0/b
    REAl(KIND=16) :: G0                                ! Recombination-Generation G0 = D0*C0/x0**2
    REAl(KIND=16) :: t0                                ! t0 = X0/D0
    REAl(KIND=16) :: J0                                ! Current density J0 = q*D0*C0/X0
    
    ! Output file name
    CHARACTER(200) :: output
    INTEGER :: ou = 6
    LOGICAL :: make_log = .FALSE.

    ! ---------------------------------------------------------------------------
    ! End of the definition of global variables
    ! ---------------------------------------------------------------------------
    
CONTAINS
!-------------------------------------------------
    SUBROUTINE log_file(my_log_file)
        CHARACTER(200) :: my_log_file
        
        output = my_log_file
        make_log = .TRUE.
        
    END SUBROUTINE log_file
!-------------------------------------------------
    SUBROUTINE cancel_log()
        
        make_log = .FALSE.
        
    END SUBROUTINE cancel_log
!-------------------------------------------------
    SUBROUTINE open_log()
        LOGICAL :: exist_file, opened_unit

        IF (make_log) THEN
            ou = 2
            INQUIRE(file=output, exist=exist_file)
            INQUIRE(unit=ou, opened=opened_unit)
            
            IF (.NOT.opened_unit) THEN            
                IF (exist_file) THEN
                  OPEN(ou, file=output, status="old", position="append", action="write")
                ELSE
                  OPEN(ou, file=output, status="new", action="write")
                END IF
            END IF

        END IF
        
    END SUBROUTINE open_log
!-------------------------------------------------
    SUBROUTINE close_log()

        IF (make_log) THEN
            close(unit=ou)
        END IF
        ou = 6
        
    END SUBROUTINE close_log
!-------------------------------------------------
    SUBROUTINE version()                                    
        CHARACTER(10) :: ver
        
        ver = '0.5.0'
        
        CALL open_log()
        WRITE(ou,*) 'Fotran Poisson - DriftDiffusion version: ', ver
        CALL close_log()
        
    END SUBROUTINE version
!-------------------------------------------------     
    SUBROUTINE InitDevice(MM)
        INTEGER :: i, j, k
        INTEGER :: MM
        REAl(KIND=16) :: Nqw, Nbulk, Vconf

        CALL open_log()
        
        ! Scaling factors        
        b   = q/(kb*T)
        C0  = MAXVAL(DML(:, 16))
        Mu0 = MAX( MAXVAL(DML(:, 5)), MAXVAL(DML(:, 6)) )
        D0  = Mu0/b
        G0  = D0*C0/XD**2
        t0  = XD**2/D0
        J0  = q*D0*C0/XD

        !We create the mesh
        CALL CreateMesh(MM)
        
        ! Refine the mesh if appropiate and show the information
        WRITE(ou,*) 'CREATE MESH...'
        IF (MM <= 0) THEN
            WRITE(ou,*) 'Masternodes at (nm):'
            WRITE(ou,'(1f10.1)')( MasterNodes(i)/1e-9, i = 1, MGrid)
        END IF
        
        ! We fill the arrays with the material properties and doping as a function of position
        DO i = 1, MReg        ! Loop over the layers    
            DO j=0,M         ! Loop over the nodes
                IF ( (X(j)>=DML(i, 1)).AND.(X(j)<=DML(i, 2)) )  THEN
                    Eg(j)    = DML(i, 3)
                    Xi(j)    = DML(i, 4)
                    Mun(j)    = DML(i, 5)
                    Mup(j)    = DML(i, 6)
                    Nc(j)    = DML(i, 7)
                    Nv(j)    = DML(i, 8)
                    tn(j)    = DML(i, 10)
                    tp(j)    = DML(i, 11)        
                    Epsi(j)    = DML(i, 12)
                    Brad(j)    = DML(i, 13)
                    CCn(j)    = DML(i, 17)
                    CCp(j)    = DML(i, 18)
                    
                    ni(j)    = SQRT(Nc(j)*Nv(j)*EXP( -b*Eg(j) ) )                    

                    AbsProfile(j, 0:NumWL) = AbsLibrary(i, 0:NumWL)
                    
                    Na(j) = DoppingLibrary(i, 1)
                    Nd(j) = DoppingLibrary(i, 2)
                    
                END IF
            END DO
        END DO
        
        ! Set the reference values to the properties of the last point, X(M)
        Egr = Eg(M)
        Xir = Xi(M)
        Munr = Mun(M)
        Mupr = Mup(M)
        Ncr = Nc(M)
        Nvr = Nv(M)
        nir = ni(M)

        DO i = 0, M            
            CALL Bandedge(i)
        END DO

        ! Apply neutrality condition to find the initial values for the potential
        WHERE (Nd(0:M)>Na(0:M))
            n = 0.5*(Nd-Na) + 0.5*SQRT((Nd-Na)**2 + 4*ni**2)
            p = ni**2/n
        ELSEWHERE
            p = -0.5*(Nd-Na)+ 0.5*SQRT((Nd-Na)**2 + 4*ni**2)
            n = ni**2/p
        ENDWHERE
        
        Fn(0:M) = - LOG(n(M)/nir) + Vn(M)
        Fp(0:M) = Fn(0:M)
        
        Psi(0:M) = LOG(n(0:M)/nir) - Vn(0:M) + Fn(0:M)

        ! We smooth the potential to facilitate the initial convergence. We smooth 10 times
        Do j = 1, 10
            DO i = 1, M-1
                Psi(i) = (Psi(i-1) + Psi(i) + Psi(i+1)) / 3.0    
            END DO
        END DO

        DO j = 0, M            
            CALL ModPotential(j)
            CALL Carriers(j)
            Rho(j)=q*( p(j)-n(j)+Nd(j)-Na(j) )
        END DO
        
        fneq = n(0)
        bneq = n(M)

        WRITE(ou,*) ' '
        IF (Dynamic) THEN
            WRITE(ou,*) 'Initial number of nodes (M+1): ', M+1
            WRITE(ou,*) 'Refining mesh... '
        
            CALL DynamicMesh(1)
        
            WRITE(ou,*) '... Finished!'
        END IF
        WRITE(ou,*) 'Mesh with ', M+1, ' nodes.'
        WRITE(ou,*) '----------------------------------'
        
        CALL close_log()
        
    END SUBROUTINE InitDevice
!-------------------------------------------------    
    SUBROUTINE AddLayer(args, dum2)                        
        !External variables
        REAl(KIND=8) :: args(0:dum2)
        INTEGER :: dum2
        
        !Internal variables
        REAl(KIND=16) :: xini, xfin
        
        MReg = MReg + 1
        xini = XD
        xfin = XD + REAL(args(0),16)
        XD   = xfin
        
        DML(MReg, 1)  = xini
        DML(MReg, 2)  = xfin
        DML(MReg, 3)  = REAL(args(1),16)    ! Eg
        DML(MReg, 4)  = REAL(args(2),16)    ! Xi
        DML(MReg, 5)  = REAL(args(3),16)    ! Mun
        DML(MReg, 6)  = REAL(args(4),16)    ! Mup
        DML(MReg, 7)  = REAL(args(5),16)    ! Nc
        DML(MReg, 8)  = REAL(args(6),16)    ! Nv
        DML(MReg, 10) = REAL(args(7),16)    ! tn
        DML(MReg, 11) = REAL(args(8),16)    ! tp
        DML(MReg, 12) = REAL(args(9),16)    ! Epsi
        DML(MReg, 13) = REAL(args(10),16)    ! Brad        
        DML(MReg, 17) = REAL(args(11),16)     ! CCn
        DML(MReg, 18) = REAL(args(12),16)    ! CCp
        
        DoppingLibrary(MReg, 1) = REAL(args(13),16)        ! Acceptors
        DoppingLibrary(MReg, 2) = REAL(args(14), 16)    ! Donors
        
        CALL AddMasterNode(REAL(xini,16))        
        CALL AddMasterNode(REAL(xfin,16)-1e-10)    
        CALL AddMasterNode(REAL(xfin,16))

    END SUBROUTINE AddLayer
!-------------------------------------------------
    SUBROUTINE AddAbsorption(Ab, WL, dum)                
        ! Add the absorption coefficients to the structure. They MUST be added in the same order than the layers before initialise the structure.
        REAl(KIND=8) :: Ab(0:dum), WL(0:dum)
        INTEGER :: dum
        
        IF (NumWL==2) THEN        ! If the number of wavelengths is equal to 2, then this is the first call to this function. 
            NumWL = SIZE(WL)-1    
            AbsLibrary(-1, 0:NumWL) = REAL(WL(0:NumWL),16)    ! We asign the wavelength values 
        END IF
        AbsLibrary(MReg, 0:NumWL) = REAL(Ab(0:NumWL), 16)        
        
    END SUBROUTINE AddAbsorption
!-------------------------------------------------    
    SUBROUTINE set_generation(gen_profile, dum_m, dum_wl)
        REAl(KIND=8) :: gen_profile(-1:dum_m, 0:dum_wl)
        INTEGER :: dum_m, dum_wl
        NumWL = dum_wl
        AbsProfile(0:M, 0:NumWL) = REAL(gen_profile(0:dum_m, 0:dum_wl), 16)
        AbsLibrary(-1, 0:NumWL) = REAL(gen_profile(-1, 0:dum_wl), 16)
                    
    END SUBROUTINE set_generation
!-------------------------------------------------    
    SUBROUTINE AddMasterNode(newpoint)                    
        REAl(KIND=16) :: newpoint
        INTEGER :: i, j
        
        DO i = 1, MGrid
            IF (MasterNodes(i) > newpoint + 0.5e-10) THEN
                MasterNodes(i+1:MGrid+1) = MasterNodes(i:MGrid)
                MasterNodes(i) = newpoint 
                MGrid = MGrid + 1
                RETURN
            ELSE IF ( ABS(MasterNodes(i)-newpoint) < 0.5e-10 ) THEN
                RETURN
            END IF            
        END DO    

        IF (MasterNodes(MGrid) < newpoint - 0.5e-10) THEN
            MGrid = MGrid + 1
            MasterNodes(MGrid) = newpoint    
        END IF
        
    END SUBROUTINE AddMasterNode
!-------------------------------------------------
    SUBROUTINE CreateMesh(MM)                              
        INTEGER :: i, j, k, lc, lf, lu, extra
        REAl(KIND=16) :: delta, deltaF, deltaUF 
        REAl(KIND=16) :: TempMasterNodes(1000)
        INTEGER :: MM
        
        IF (MM>0) THEN
            M = MM
            DO i = 0, M
                X(i)  = i*XD/M
                dX(i) = XD/M
            END DO
            RETURN
        END IF        
        
        MasterNodes(MGrid-1) = MasterNodes(MGrid)
        MGrid = MGrid-1
        
        IF (MM < 0) Dynamic = .TRUE.
        
        j = 0
        ! Loop for the coarse mesh
        DO i = 1, MGrid-1
            
            X(j) = MasterNodes(i)
            lc = CEILING((MasterNodes(i+1)-MasterNodes(i))/Coarse) - 1
            
            delta = (MasterNodes(i+1)-MasterNodes(i))/(lc+1)
            
            IF (lc==0) THEN 
                extra=0
            ELSE
                extra=1
            END IF
            
            ! The fine mesh after a GridLine
            IF (Fine<Coarse) THEN
                lf = CEILING(delta/Fine) - 1
                deltaF = delta/(lf+1)
                
                ! The ultrafine mesh after a GridLine
                IF (Ultrafine<Fine) THEN
                    lu = CEILING(deltaF/Ultrafine) - 1
                    deltaUF = deltaF/(lu+1)
                    
                    ! Loop for the points of the ultrafine mesh
                    DO k = 1, lu+1
                        j = j + 1                
                        X(j) = X(j-1) + deltaUF
                    END DO
                    
                ELSE
                    ! The first point of the fine mesh
                    j = j + 1                
                    X(j) = X(j-1) + deltaF
                END IF
                
                ! Intermediate points of the fine mesh
                DO k = 1, lf-1+extra
                    j = j + 1                
                    X(j) = X(j-1) + deltaF
                END DO
            END IF

            ! Loop for the intermediate points of the coarse mesh
            DO k = 1, lc-1!+extra
                j = j + 1                
                X(j) = X(j-1) + delta
            END DO
            
            ! Loop for the fine mesh before a GridLine
            IF (Fine<Coarse) THEN
                IF (lc/=0) THEN                    
                    DO k = 1, lf
                        j = j + 1                
                        X(j) = X(j-1) + deltaF
                    END DO
                
                    ! Loop for the ultrafine mesh before a GridLine
                    IF (Ultrafine<Fine) THEN
                        DO k = 1, lu
                            j = j + 1                
                            X(j) = X(j-1) + deltaUF
                        END DO
                        j = j + 1
                    ELSE
                        j = j + 1
                    END IF
                ELSE 
                    IF (lf/=0) THEN
                        ! Loop for the ultrafine mesh before a GridLine
                        IF (Ultrafine<Fine) THEN
                            DO k = 1, lu
                                j = j + 1                
                                X(j) = X(j-1) + deltaUF
                            END DO
                            j = j + 1
                        ELSE
                            j = j + 1
                        END IF
                    END IF
                END IF
            ELSE
                j = j + 1
            END IF

        END DO
        
        X(j) = MasterNodes(MGrid)
        M = j    
        
        DO i = 0, M-1
            dX(i) = X(i+1)-X(i)
        END DO
        
    END SUBROUTINE CreateMesh
!-------------------------------------------------
    SUBROUTINE DynamicMesh(Initial)                        
        ! External variables
        INTEGER :: Initial
        ! Internal variables
        INTEGER :: i, j, k, l, frac, NodesPerRegion(1:MGrid), ilast
        REAl(KIND=16) :: Xtemp(0:6000), Dpot, Dmaj, Dpot2, Dmaj2, Dvar, Dvar2, RegionSize(1:MGrid), DGR, DGR2, minSize
        REAl(KIND=16) :: psiint, fpint, fnint, nint, pint, Gint, step
        LOGICAL :: Join, NotMasterNode(0:6000)
        
        ! Temporal material arrays
        REAl(KIND=16), DIMENSION(0:6000) :: TempFn, TempFp        !Quasi-Fermi potential for electrons and holes (V)
        REAl(KIND=16), DIMENSION(0:6000) :: TempPsi                !Electrostatic potential (V)
        REAl(KIND=16), DIMENSION(0:6000, 0:3000) :: TempAbsProfile!Absorption coefficient as a function of wavelength.
        
        l = 1
        RegionSize = 0.0
        NotMasterNode(0) = .FALSE.
        DO i = 1, M
            IF (X(i)>=MasterNodes(l+1)) THEN
                NotMasterNode(i) = .FALSE.
                RegionSize(l) = MasterNodes(l+1)-MasterNodes(l)
                l = l+1
            ELSE
                NotMasterNode(i) = .TRUE.
            END IF    
        END DO
        
        minSize = 2e-10
        l = 1
        j = 0
        ilast = 0
        Xtemp(0) = X(0)
        DO i = 1, M-1
            
            !print*, i, ilast
            psiint = Interpol(Xtemp(j), X(ilast), Psi(ilast), X(i), Psi(i))
            fpint = Interpol(Xtemp(j), X(ilast), Fp(ilast), X(i), Fp(i))
            fnint = Interpol(Xtemp(j), X(ilast), Fn(ilast), X(i), Fn(i))
            nint = EXP( Interpol(Xtemp(j), X(ilast), LOG(n(ilast) ), X(i), LOG( n(i)) ) )
            pint = EXP( Interpol(Xtemp(j), X(ilast), LOG(p(ilast) ), X(i), LOG( p(i)) ) )
!             Gint = EXP( Interpol(Xtemp(j), X(ilast), LOG(G(ilast) ), X(i), LOG( G(i)) ) )
                        
            Dpot = MAX( MAX( ABS(psiint-Psi(i)), ABS(fnint-Fn(i))  ), ABS(fpint-Fp(i))  ) 
            Dpot2 = MAX( MAX( ABS(psiint-Psi(i+1)), ABS(fnint-Fn(i+1))  ), ABS(fpint-Fp(i+1))  ) 
            
            Dmaj = MAX( ABS( LOG(nint/n(i)) ), ABS( LOG(pint/p(i)) ) )
            Dmaj2 = MAX( ABS( LOG(nint/n(i+1)) ), ABS( LOG(pint/p(i+1)) ) )
            
!             IF (Gint > 1e-7) THEN
!                 DGR = ABS( LOG(Gint/G(i)) )
!                 DGR2 = ABS( LOG(Gint/G(i+1)) )
!                 print*, DGR2, DGR, Dmaj, Growth
!             ELSE
!                 DGR = -1
!                 DGR2 = -1
!             END IF
            
!             Dvar = MAX(Dpot, MAX(Dmaj, DGR)) /2
!             Dvar2 = MAX(Dpot2, MAX(Dmaj2, DGR2)) /2
            
            Dvar = MAX(Dpot, Dmaj) / 2
            Dvar2 = MAX(Dpot2, Dmaj2) / 2
            
            ! Conditions for joining nodes.
            Join = Dvar2 < Growth                ! Not too much variation of the potentials or the carrier densities
            Join = Join.AND.NotMasterNode(i)    ! Not a master node
            Join = Join.AND.( (X(i+1)-Xtemp(j)) < Growth*(X(i+1)-MasterNodes(l))  )    ! Not too close to the left masternode and
            Join = Join.AND.( (X(i+1)-Xtemp(j)) < Growth*(MasterNodes(l+1)-X(i))  )    ! Nor to the rigth one
            Join = Join.AND.( (X(i+1)-Xtemp(j)) < Growth*RegionSize(l)/10.0)   ! Not too big for the region
            
            ! Divide if the variation of the potentials or of the carrier density is too large
            IF ( Dvar > Growth ) THEN                
                frac = CEILING(Dvar/Growth)
                IF (CEILING(dX(i-1)/minSize) < frac) frac = CEILING(dX(i-1)/minSize)    
                DO k = 1, frac-1
                    j = j+1
                    Xtemp(j) = X(i-1) + k*dX(i-1)/frac
                END DO
                j = j+1
                Xtemp(j) = X(i)
                ilast = i
            ! Divide if we are too close to the master node on the left 
            ELSE IF ( dX(i-1) >= Growth*(X(i)-MasterNodes(l)) ) THEN        
                frac = CEILING(dX(i-1)/Growth/(X(i)-MasterNodes(l)))
                IF (CEILING(dX(i-1)/minSize) < frac) frac = CEILING(dX(i-1)/minSize)
                DO k = 1, frac-1
                    j = j+1
                    Xtemp(j) = X(i-1) + k*dX(i-1)/frac
                END DO
                j = j+1
                Xtemp(j) = X(i)
                ilast = i
            ! Divide if we are too close to the master node on the rigth
            ELSE IF ( dX(i-1) >= Growth*(MasterNodes(l+1)-X(i)) ) THEN        
                frac = CEILING(dX(i-1)/Growth/(MasterNodes(l+1)-X(i)) )        
                IF (CEILING(dX(i-1)/minSize) < frac) frac = CEILING(dX(i-1)/minSize)    
                DO k = 1, frac-1
                    j = j+1
                    Xtemp(j) = X(i-1) + k*dX(i-1)/frac
                END DO
                j = j+1
                Xtemp(j) = X(i)
                ilast = i
            ! Divide if the element is too big for the region
            ELSE IF ( dX(i-1) >= Growth*RegionSize(l)/10.0 ) THEN        
                frac = CEILING(dX(i-1) / (Growth*RegionSize(l)/10.0) )
                IF (CEILING(dX(i-1)/minSize) < frac) frac = CEILING(dX(i-1)/minSize)        
                DO k = 1, frac-1
                    j = j+1
                    Xtemp(j) = X(i-1) + k*dX(i-1)/frac
                END DO
                j = j+1
                Xtemp(j) = X(i)
                ilast = i
            ! If we haven't deleted de node nor divided the element, we keep it as it is 
            ELSE IF (.NOT.Join) THEN
                j = j+1
                Xtemp(j) = X(i)
                ilast = i        
            END IF
            
            ! Check if we are about to change the region
            IF (.NOT.NotMasterNode(i)) l = l+1
    
        END DO
        
        ! The last interval
        frac = CEILING(dX(M-1)/Growth/(X(M)-X(M-1)) )            
        DO k = 1, frac-1
            j = j+1
            Xtemp(j) = X(M-1) + k*dX(M-1)/frac
        END DO
        
        j = j+1
        Xtemp(j) = X(M)
        
        ! Now it's time to interpolate all the other variables. 
        ! We just need to be carful with the masternodes to avoid mixing the properties of the regions. 
        
        ! First, we update the position of the masternodes
        dX(0:M) = 0.0
        l = 2
        NotMasterNode(0) = .FALSE.
        DO i = 1, j
            IF (Xtemp(i)>=MasterNodes(l)) THEN
                NotMasterNode(i) = .FALSE.
                l = l+1
            ELSE
                NotMasterNode(i) = .TRUE.
            END IF    
            dX(i-1) = Xtemp(i)-Xtemp(i-1)
        END DO
        
        ! We smooth the position of the nodes to avoid neibourgh elements too different. We smooth 10 times
        Do k = 1, 10
            DO i = 1, j
                IF (NotMasterNode(i)) THEN
                    Xtemp(i) = (Xtemp(i-1) + Xtemp(i) + Xtemp(i+1)) / 3.0    
                END IF    
                dX(i-1) = Xtemp(i)-Xtemp(i-1)
            END DO
        END DO
            
        k = -1
        DO i = 0, j
            IF (NotMasterNode(i)) THEN
                ! Between master nodes, we interpolate
                DO WHILE ( Xtemp(i)>X(k+1))
                    k = k+1        
                END DO

                TempFn(i)     = Interpol(Xtemp(i), X(k), Fn(k), X(k+1), Fn(k+1) )
                TempFp(i)     = Interpol(Xtemp(i), X(k), Fp(k), X(k+1), Fp(k+1) )
                TempPsi(i)     = Interpol(Xtemp(i), X(k), Psi(k), X(k+1), Psi(k+1) )
                
                DO l = 0, NumWL
                    TempAbsProfile(i, l) = Interpol(Xtemp(i), X(k), AbsProfile(k, l), X(k+1),AbsProfile(k+1,l) )
                END DO
                
            ELSE
                k = k + 1 
                DO WHILE ( Xtemp(i)>X(k))
                    k = k+1        
                END DO
                ! At the masternodes, we keep the same values
                TempFn(i)     = Fn(k)
                TempFp(i)     = Fp(k)
                TempPsi(i)     = Psi(k)
                TempAbsProfile(i, 0:NumWL) = AbsProfile(k, 0:NumWL)
            END IF
        END DO
        
        ! And update them with the new values
        M = j                
        X(0:M)      = Xtemp(0:M)        
        Fn(0:M)     = TempFn(0:M)    
        Fp(0:M)     = TempFp(0:M)    
        Psi(0:M)    = TempPsi(0:M)    
        AbsProfile(0:M, 0:NumWL) = TempAbsProfile(0:M, 0:NumWL)

        ! We fill the arrays with the material properties and doping as a function of position
        DO i = 1, MReg        ! Loop over the layers    
            DO j=0,M         ! Loop over the nodes
                IF ( (X(j)>=DML(i, 1)).AND.(X(j)<=DML(i, 2)) )  THEN
                    Eg(j)    = DML(i, 3)
                    Xi(j)    = DML(i, 4)
                    Mun(j)    = DML(i, 5)
                    Mup(j)    = DML(i, 6)
                    Nc(j)    = DML(i, 7)
                    Nv(j)    = DML(i, 8)
                    tn(j)    = DML(i, 10)
                    tp(j)    = DML(i, 11)        
                    Epsi(j)    = DML(i, 12)
                    Brad(j)    = DML(i, 13)
                    CCn(j)    = DML(i, 17)
                    CCp(j)    = DML(i, 18)
            
                    ni(j)    = SQRT(Nc(j)*Nv(j)*EXP( -b*Eg(j) ) )                    
            
                    Na(j) = DoppingLibrary(i, 1)
                    Nd(j) = DoppingLibrary(i, 2)
            
                END IF
            END DO
        END DO
        
        DO i = 0, M
            CALL Bandedge(i)
            CALL ModPotential(i)
            CALL Carriers(i)
            Rho(i)=q*( p(i)-n(i)+Nd(i)-Na(i))
        END DO
        CALL GR_sub        
        
    END SUBROUTINE DynamicMesh
!-------------------------------------------------    
    FUNCTION Interpol(xx, x1, y1, x2, y2)                
        REAl(KIND=16) :: xx, x1, y1, x2, y2, Interpol
        REAl(KIND=16) :: a, b
        
        a = (y2-y1)/(x2-x1)
        b = (y1*x2-y2*x1)/(x2-x1)
        Interpol = a*xx + b
        
    END FUNCTION Interpol
!-------------------------------------------------    
    SUBROUTINE Reset()                                    
        MGrid = 1
        Mreg = 0
        NumWL  = 2
        M = 0
        XD = 0
        fneq = 1
        bneq = 1
        
        X(:)     = 0.0
        Nd(:)     = 0.0
        Na(:)     = 0.0
        Nc(:)     = 0.0
        Nv(:)     = 0.0
        ni(:)     = 0.0
        Eg(:)     = 0.0
        Xi(:)     = 0.0
        Epsi(:) = 0.0
        Mun(:)     = 0.0
        Mup(:)    = 0.0
        Ncc(:)     = 0.0
        Nvhh(:) = 0.0
        Nvlh(:) = 0.0
        tn(:)     = 0.0
        tp(:)     = 0.0
        Brad(:) = 0.0
        alfa(:) = 0.0
        CCn(:)     = 0.0
        CCp(:)     = 0.0
        Fn(:)     = 0.0
        Fp(:)     = 0.0
        Psi(:)     = 0.0
        AbsProfile(:, :) = 0.0
        f(:) = 0
        dsol(:) = 0
        Jac(:,:) = 0
        MasterNodes(:) = 0
        AbsLibrary(:, :) = 0
        PFspectrum = 0.0
        G(:) = 0
        
        clamp = 20
        ATol = 3.1622776601683796e-17
        RTol = 1e-6
        nitermax = 40
        Rs = 0

        SingleWL = .FALSE.
    
    END SUBROUTINE Reset
!-------------------------------------------------    
    FUNCTION Get(VarName)                                
        REAl(KIND=8) :: Get (0:6000)
        INTEGER:: k
        CHARACTER(len=30) :: VarName
        
        SELECT CASE ( VarName )
        
        ! Material information (inputs)
            CASE ( 'x' )
                Get(0:M) = REAL(X(0:M), 8)
            CASE ( 'dx' )
                Get(0:M-1) = REAL(dX(0:M-1), 8)
            CASE ( 'eg' )
                Get(0:M) = REAL(Eg(0:M), 8)
            CASE ( 'xi' )
                Get(0:M) = REAL(Xi(0:M), 8)
            CASE ( 'na' )
                Get(0:M) = REAL(Na(0:M), 8)
            CASE ( 'nd' )
                Get(0:M) = REAL(Nd(0:M), 8)
            CASE ( 'mun' )
                Get(0:M) = REAL(Mun(0:M), 8)
            CASE ( 'mup' )
                Get(0:M) = REAL(Mup(0:M), 8)
            CASE ( 'epsi' )
                Get(0:M) = REAL(Epsi(0:M), 8)
            CASE ( 'tn' )
                Get(0:M) = REAL(tn(0:M), 8)
            CASE ( 'tp' )
                Get(0:M) = REAL(tp(0:M), 8)
            CASE ( 'brad' )
                Get(0:M) = REAL(Brad(0:M), 8)
            CASE ( 'ccn' )
                Get(0:M) = REAL(CCn(0:M), 8)
            CASE ( 'ccp' )
                Get(0:M) = REAL(CCp(0:M), 8)
            CASE ( 'ni' )
                Get(0:M) = REAL(ni(0:M), 8)
            CASE ( 'nc' )
                Get(0:M) = REAL(Nc(0:M), 8)
            CASE ( 'nv' )
                Get(0:M) = REAL(Nv(0:M), 8)
            CASE ( 'ncc' )
                Get(0:M) = REAL(ncc(0:M), 8)
            CASE ( 'nvhh' )
                Get(0:M) = REAL(nvhh(0:M), 8)
            CASE ( 'nvlh' )
                Get(0:M) = REAL(nvlh(0:M), 8)
                
        ! Carrier and charge densites densities
            CASE ( 'n' )
                Get(0:M) = REAL(n(0:M), 8)
            CASE ( 'p' )
                Get(0:M) = REAL(p(0:M), 8)
            CASE ( 'rho' )
                Get(0:M) = REAL(Rho(0:M), 8)
                
        ! Generation/Recombination
            CASE ( 'gr' )
                Get(0:M-1) = REAL(GR(0:M-1)/q/dX(0:M-1), 8)
            CASE ( 'rsrh' )
                Get(0:M-1) = REAL(Rsrh(0:M-1)/q/dX(0:M-1), 8)
            CASE ( 'rrad' )
                Get(0:M-1) = REAL(Rrad(0:M-1)/q/dX(0:M-1), 8)
            CASE ( 'raug' )
                Get(0:M-1) = REAL(Raug(0:M-1)/q/dX(0:M-1), 8)
            CASE ( 'g' )
                Get(0:M-1) = REAL(G(0:M-1)/q/dX(0:M-1), 8)
        
        ! Bandstructure
            CASE ( 'fp' )
                Get(0:M) = REAL(Fp(0:M)/b, 8)
            CASE ( 'fn' )
                Get(0:M) = REAL(Fn(0:M)/b, 8)
            CASE ( 'vn' )
                Get(0:M) = REAL(Vn(0:M)/b, 8)
            CASE ( 'vp' )
                Get(0:M) = REAL(Vp(0:M)/b, 8)
            CASE ( 'cn' )
                Get(0:M) = REAL(Cn(0:M)/b, 8)
            CASE ( 'cp' )
                Get(0:M) = REAL(Cp(0:M)/b, 8)
            CASE ( 'efh' )
                Get(0:M) = REAL((Fp(0)-Fp(0:M))/b, 8)
            CASE ( 'efe' )
                Get(0:M) = REAL((Fp(0)-Fn(0:M))/b, 8)
            CASE ( 'psi' )
                Get(0:M) = REAL(Psi(0:M)/b, 8)
            CASE ( 'ev' )
                Get(0:M) = REAL( (Fp(0)-Fp(0:M)+LOG( p(0:M)/Nv(0:M) )) /b, 8)
            CASE ( 'ec' )
                Get(0:M) = REAL( (-LOG( n(0:M)/Nc(0:M) ) +Fp(0)-Fn(0:M)  ) /b, 8)
            
        ! Voltage and current information
            CASE ( 'voc' )
                Get(0) = REAL(Voc, 8)
            CASE ( 'isc' )
                Get(0) = REAL(Isc, 8)    
            CASE ( 'vmax' )
                Get(0) = REAL(Vmax, 8)    
            CASE ( 'imax' )
                Get(0) = REAL(Imax, 8)
            CASE ( 'ff' )
                Get(0) = REAL(FF, 8)            
            CASE ( 'volt' )
                Get(1:nvolt) = REAL(vpoint(0:nvolt-1), 8)
            CASE ( 'jtot' )
                Get(1:nvolt) = REAL(jpoint(0:nvolt-1), 8)
            CASE ( 'jsrh' )
                Get(1:nvolt) = REAL(jsrhpoint(0:nvolt-1), 8)
            CASE ( 'jrad' )
                Get(1:nvolt) = REAL(jradpoint(0:nvolt-1), 8)
            CASE ( 'jaug' )
                Get(1:nvolt) = REAL(jaugpoint(0:nvolt-1), 8)
            CASE ( 'jsur' )
                Get(1:nvolt) = REAL(jsurpoint(0:nvolt-1), 8)
            CASE ( 'residual' )
                Get(1:nvolt) = REAL(residual(0:nvolt-1), 8)
        
        ! Internal quantum efficiency
            CASE ( 'iqe' )
                Get(0:NumWL) = REAL(iqe(0:NumWL), 8)
            CASE ( 'iqesrh' )
                Get(0:NumWL) = REAL(iqesrh(0:NumWL), 8)
            CASE ( 'iqerad' )
                Get(0:NumWL) = REAL(iqerad(0:NumWL), 8)
            CASE ( 'iqeaug' )
                Get(0:NumWL) = REAL(iqeaug(0:NumWL), 8)
            CASE ( 'iqesurf' )
                Get(0:NumWL) = REAL(iqesurf(0:NumWL), 8)
            CASE ( 'iqesurb' )
                Get(0:NumWL) = REAL(iqesurb(0:NumWL), 8)

        END SELECT
        
    END FUNCTION Get
!-------------------------------------------------
    SUBROUTINE Set(VarName, VarVal, index, index2)        
        CHARACTER(len=30) :: VarName
        REAl(KIND=8) :: VarVal
        INTEGER, OPTIONAL :: index, index2
        
        SELECT CASE ( VarName )
            CASE ( 't' )
                T = REAL(VarVal,16)
        
        ! Material information (inputs)
!             CASE ( 'sn' )
!                 Sn = REAL(VarVal,16)
!             CASE ( 'sp' )
!                 Sp = REAL(VarVal,16)
            CASE ( 'eg' )
                Eg(index) = REAL(VarVal,16)
            CASE ( 'xi' )
                Epsi(index) = REAL(VarVal,16)
            CASE ( 'epsi' )
                Mun(index) = REAL(VarVal,16)
            CASE ( 'mun' )
                Mun(index) = REAL(VarVal,16)
            CASE ( 'mup' )
                Mup(index) = REAL(VarVal,16)
            CASE ( 'ncc' )
                ncc(index) = REAL(VarVal,16)
            CASE ( 'nvhh' )
                nvhh(index) = REAL(VarVal,16)
            CASE ( 'nvlh' )
                nvlh(index) = REAL(VarVal,16)
            CASE ( 'tn' )
                tn(index) = REAL(VarVal,16)
            CASE ( 'tp' )
                tp(index) = REAL(VarVal,16)
            CASE ( 'Brad' )
                Brad(index) = REAL(VarVal,16)
            CASE ( 'ccn' )
                CCn(index) = REAL(VarVal,16)
            CASE ( 'ccp' )
                CCp(index) = REAL(VarVal,16)
            CASE ( 'absprofile' )
                AbsProfile(index, index2) = REAL(VarVal,16)    
                        
        ! Meshing and convergence        
            CASE ( 'coarse' )
                Coarse = REAL(VarVal,16)
            CASE ( 'fine' )
                Fine = REAL(VarVal,16)
            CASE ( 'ultrafine' )
                Ultrafine = REAL(VarVal,16)
            CASE ( 'clamp' )
                Clamp = REAL(VarVal,16)
            CASE ( 'atol' )
                ATol = REAL(VarVal,16)
            CASE ( 'rtol' )
                RTol = REAL(VarVal,16)
            CASE ( 'growth' )
                Growth = REAL(VarVal,16)
                
        ! Device inputs
            CASE ( 'rs' )
                Rs = REAL(VarVal,16)
        
        ! Others        
            CASE ( 'xd' )
                XD = REAL(VarVal,16)
                
        END SELECT

    END SUBROUTINE Set
!-------------------------------------------------
    SUBROUTINE Illumination(spectrum, dum)                
        REAl(KIND=8) :: spectrum(0:dum)
        INTEGER:: dum, i

        PFspectrum(0:NumWL) = REAL(spectrum(0:NumWL),16)
        
        PhotonFlux = 0
        DO i = 0, NumWL-1    
            PhotonFlux = PhotonFlux + (PFspectrum(i)+PFspectrum(i+1))/2 * ( AbsLibrary(-1,i+1)-AbsLibrary(-1, i) ) 
        END DO
        
    END SUBROUTINE Illumination
!-------------------------------------------------
    SUBROUTINE FrontBoundary(type, surfe, surfh, barrier)        
        !External variables
        CHARACTER(len=30) :: type
        REAl(KIND=8) :: surfe, surfh, barrier
        
        Snfront = REAL(surfe,16)
        Spfront = REAL(surfh,16)
        FSUR  = 1        
        FTYPE = 0
        
        SELECT CASE ( type )
            CASE ( 'ohmic' )
                FTYPE = 1
            CASE ( 'schottky' )
                print*, 'Not working yet. Chose Ohmic or Insulating.' 
                RETURN
                FTYPE = 2
                Vbarf = REAL(barrier,16)
            CASE ( 'insulator' )
                FTYPE = 3
                OC = 1
            CASE ( 'charged' )
                print*, 'Not working yet. Chose Ohmic or Insulating.' 
                RETURN
                FTYPE = 4
                OC = 1
                Vbarf = REAL(barrier,16)
        END SELECT
    
    END SUBROUTINE FrontBoundary
!-------------------------------------------------
    SUBROUTINE BackBoundary(type, surfe, surfh, barrier)        
        !External variables
        CHARACTER(len=30) :: type
        REAl(KIND=8) :: surfe, surfh, barrier
        
        Snback = REAL(surfe,16)
        Spback = REAL(surfh,16)
        BSUR  = 1        
        BTYPE = 0
        
        SELECT CASE ( type )
            CASE ( 'ohmic' )
                BTYPE = 1
            CASE ( 'schottky' )
                print*, 'Not working yet. Chose Ohmic or Insulating.' 
                RETURN
                BTYPE = 2
                Vbarb = REAL(barrier,16)
            CASE ( 'insulator' )
                BTYPE = 3
                OC = 1
            CASE ( 'charged' )
                print*, 'Not working yet. Chose Ohmic or Insulating.' 
                RETURN
                BTYPE = 4
                OC = 1
                Vbarb = REAL(barrier,16)
        END SELECT
    
    END SUBROUTINE BackBoundary    
!-------------------------------------------------    
! Functions to populate and update the potentials and carrier densities
!-------------------------------------------------    
    SUBROUTINE Carriers(i)                                
        ! Calculation of the carrier densities        
        INTEGER :: i
                
        n(i) = nir*EXP( Psi(i) - Fn(i) + Vn(i))
        p(i) = nir*EXP(-Psi(i) + Fp(i) + Vp(i))
        
    END SUBROUTINE Carriers
!-------------------------------------------------            
    SUBROUTINE Bandedge(i)                                 
        !Calculation of the bandegde potentials. Bandgap narrowing due to heavy dopping ignored
        INTEGER :: i
        
        Vn(i)=LOG(Nc(i)/Ncr) + b*(Xi(i)-Xir)         ! + b*DEgc(i)
        Vp(i)=LOG(Nv(i)/Nvr) - b*(Eg(i)-Egr) - b*(Xi(i)-Xir)         ! + b*DEgv(i)  
        
    END SUBROUTINE Bandedge
!-------------------------------------------------    
    SUBROUTINE ModPotential(i)                            
        ! Calculation of the modified electrostatic potential    
        INTEGER :: i
        
        Cn(i) = Psi(i) + Vn(i) + LOG(Mun(i)/Munr) 
        Cp(i) = Psi(i) - Vp(i) - LOG(Mup(i)/Mupr) 
                
    END SUBROUTINE ModPotential
!-------------------------------------------------    
! Auxiliary functions and their derivatives
!-------------------------------------------------
    FUNCTION Z(u)
        REAl(KIND=16) :: u
        REAl(KIND=16) :: Z
        REAl(KIND=16) :: tiny = 1.0e-6

        IF (abs(u) < tiny) THEN
            Z = 1.0 - 0.5*u
        ELSE
            Z = u/(EXP(u)-1.0)
        END IF

    END FUNCTION Z
!-------------------------------------------------
    FUNCTION Y(u)
        REAl(KIND=16) :: u
        REAl(KIND=16) :: Y
        REAl(KIND=16) :: tiny = 1.0e-6

        IF (abs(u) < tiny) THEN
            Y = 0.5 - 1.0/12.0*u
        ELSE
            Y = (1.0-Z(u))/u
        END IF

    END FUNCTION Y
!-------------------------------------------------
    FUNCTION dZ(u)
        REAl(KIND=16) :: u
        REAl(KIND=16) :: dZ
        REAl(KIND=16) :: tiny = 1.0e-6

        IF (abs(u) < tiny) THEN
            dZ = - 0.5 + 1.0/6.0*u
        ELSE
            dZ = 1.0/(EXP(u)-1.0) - EXP(u)*u/(EXP(u)-1.0)**2
        END IF

    END FUNCTION dZ
!-------------------------------------------------
    FUNCTION dY(u)
        REAl(KIND=16) :: u
        REAl(KIND=16) :: dY
        REAl(KIND=16) :: tiny = 1.0e-6

        IF (abs(u) < tiny) THEN
            dY = - 1.0/12.0 + 1.0/240.0*u**2
        ELSE
            dY = (Z(u)-1.0)/u**2 - dZ(u)/u

        END IF

    END FUNCTION dY
!-------------------------------------------------    
! The calculation of X such that F(X)=0, with F the poison and continuity equations and X the potentials.
!     - bandec:    calculates the LU decomposition of A
!     - bandbks:    calculates the X vector
!   - SolveLin: combines the above to obtain the solution of AX=B, with A a band matrix in compact form.
!   - backtracking: Apply the clamp and backtracks along the solution direction to find a solution that reduces norm(f). 
!   - SolveNonLin: Solves the no-linear sets of equations using the previous functions. 
!-------------------------------------------------
    SUBROUTINE bandec(a, nrow, m1, m2, np, mp, al, mp1, indx, d)    
        ! Subroutine reproduced from Numerical Recipies in Fortran. 2nd Ed.
        ! It calculate the LU decomposition of a band matrix stored in the compact form. The upper triangula rmatrix is stored in 
        ! a. 
        !     a(np, mp):    On input, is a band matrix in its compact form. On output is the upper diagonal matrix of the LU 
        !                decomposition.
        !    np, mp     : The physical dimensions of a. The "logical" or "useful" dimensions are nrow and m1+m2+1. 
        !     nrow    : The number of "useful" rows of a.
        !    m1, m2    : the number of diagonals below and above the main diagonal. 
        !    al(np, mp1) : Lower triangular matrix. 
        !    mp1     : Numer of columns of the lower triangular matrix. It mus be >= m1
        !    indx       : Vectros that stores the pivoting sequence
        !      d        : stores +1 or -1 depending on whether the number of row interchanges is even or odd, respectively. 
        !
        INTEGER    :: m1, m2, mp, mp1, nrow, np, indx(nrow)
        REAl(KIND=16)    :: d, a(np, mp), al(np, mp1)
        REAl(KIND=16), PARAMETER :: TINY = 1.0e-20    
        
        INTEGER :: i, j, k, l, mm
        REAl(KIND=16)  :: dum
        
        mm = m1+m2+1
        IF (mm>mp .or. m1>mp1 .or. nrow>np) STOP "Bad arguments at bandec."
        
        l = m1
        DO i = 1, m1        ! Rearrange the storage (I do not know why)
            DO j = m1+2-i, mm
                a(i, j-l) = a(i, j)
            END DO
            l = l-1
            DO j = mm-l, mm
                a(i, j) = 0.0
            END DO
        END DO

        d = 1.0
        l = m1
        DO k = 1, nrow                            ! For each row...
            dum = a(k, 1)
            i = k
            IF (l<nrow) l = l+1
            DO j = k+1, l                 
                IF (ABS(a(j, 1))>ABS(dum)) THEN ! ... find the pivot.
                    dum = a(j, 1)
                    i     = j                
                END IF
            END DO
            indx(k) = i
            
            ! The matrix is algorithmically singular, but proceed anyway with TINY pivot.
            IF (dum==0) a(k, 1) = TINY
            
            IF (i/=k) THEN            ! Interchange rows. Probably, this can be done in a more compact way. 
                d = -d
                DO j = 1, mm
                    dum     = a(k, j)
                    a(k,j)     = a(i, j)
                    a(i, j) = dum
                END DO
            END IF
            
            DO i = k+1, l             ! Do the elimination, filling the upper diagonal matrix in a and the lower diagonal in al. 
                dum = a(i, 1) / a(k, 1)
                
                al(k, i-k) = dum
                DO j = 2, mm
                    a(i, j-1) = a(i, j) - dum*a(k, j)
                END DO
                a(i, mm) = 0.0
            END DO
        END DO
        
    END SUBROUTINE bandec 
!-------------------------------------------------
    SUBROUTINE bandbks(a, nrow, m1, m2, np, mp, al, mp1, indx, b)    
        ! Subroutine reproduced from Numerical Recipies in Fortran. 2nd Ed.
        ! It solves the system of equations AX=B when A is a banded matrix, calling the subroutine bandec
        !     ain(np, mp):    The band matrix A in its compact form. It is not modified. 
        !    np, mp     : The physical dimensions of a. The "logical" or "useful" dimensions are nrow and m1+m2+1. 
        !     nrow    : The number of "useful" rows of a.
        !    m1, m2    : the number of diagonals below and above the main diagonal. 
        !      b        : On input, it is the B vector of the equation. On output, it is the solution vector X.  
        ! 
        INTEGER    :: nrow, m1, m2, np, mp, mp1, indx(nrow)  
        REAl(KIND=16)    :: a(np, mp), b(np), al(np, mp1)
        
        INTEGER :: i, k, l, mm 
        REAl(KIND=16)  :: dum 
        
        mm = m1+m2+1
        IF (mm>mp .or. m1>mp1 .or. nrow>np) STOP "Bad arguments at bandbks."

        l = m1
        DO k = 1, nrow                    ! Forward substitution, unescrambling the permutted rows as we go.
            i = indx(k)
            IF (i/=k) THEN
                dum  = b(k)
                b(k) = b(i)
                b(i) = dum
            END IF
            IF (l<nrow) l = l+1
            DO i = k+1, l
                b(i) = b(i)-al(k, i-k)*b(k)
            END DO
        END DO
        l = 1
        DO i = nrow, 1, -1                ! Backsubstitution
            dum = b(i)
            DO k = 2, l
                !print *, i, l, k
                dum = dum-a(i,k)*b(k+i-1)
            END DO
            b(i) = dum/a(i,1)
            IF (l<mm) l = l+1
        END DO
                                        
    END SUBROUTINE bandbks
!-------------------------------------------------
    SUBROUTINE SolveLin(ain, nrow, m1, m2, np, mp, bX)                
        ! It solves the system of equations AX=B when A is a banded matrix, calling the subroutine bandec and bandbks.
        ! 
        !     ain(np, mp):    The band matrix A in its compact form. It is not modified. 
        !    np, mp     : The physical dimensions of a. The "logical" or "useful" dimensions are "nrow" and "m1+m2+1". 
        !     nrow    : The number of "useful" rows of a.
        !    m1, m2    : the number of diagonals below and above the main diagonal. 
        !      bX        : On input, it is the B vector of the equation. On output, it is the solution vector X.  
        ! 
        INTEGER    :: nrow, m1, m2, np, mp 
        REAl(KIND=16)    :: ain(np, mp), bX(np)
        
        INTEGER :: i, k, l, mm, indx(nrow), info 
        REAl(KIND=16)  :: dum, a(np, mp), al(np, m1), xtemp(np), xtemp2(np), d
        
        info = 1
        mm = m1+m2+1
        IF (mm>mp .or. nrow>np) STOP "Bad arguments at SolveLin."
        
        a = ain
        
        ! We perform the LU decomposition
        CALL bandec(a, nrow, m1, m2, np, mp, al, m1, indx, d)
        
        !And solve the system.
        CALL bandbks(a, nrow, m1, m2, np, mp, al, m1, indx, bX)

    END SUBROUTINE SolveLin
!-------------------------------------------------
    SUBROUTINE backtracking(sum0, outBacktrack, MaxCorr)            
        REAl(KIND=16) :: CopyFp(0:M), CopyPsi(0:M), CopyFn(0:M), ModFp(0:M), ModPsi(0:M), ModFn(0:M)
        REAl(KIND=16) :: lambda, alfa, sum0, sum1, sum2, sum3, dum, MaxCorr, corr
        LOGICAL :: outBacktrack
        
        INTEGER :: i, l, k, j, maxi

        MaxCorr = MAXVAL(ABS(f(1:3*M-1)))
        !print*, MAXLOC(f(1:3*M-1)), REAL()
        
        ModFn = 0.0
        ModFp = 0.0
        ModPsi = 0.0
        CopyFn = 0.0
        CopyFp = 0.0
        CopyPsi = 0.0
        
        ! First, we copy the old potentials and the corrections
        DO l = 0, M
            CopyFp(l)  = Fp(l)
            CopyPsi(l) = Psi(l)
            CopyFn(l)  = Fn(l)
            
            ! Maximum change 
            ModFp(l)  = f(3*l+1)
            ModPsi(l) = f(3*l+2)  
            ModFn(l)  = f(3*l+3)
        END DO

        ! And now, we perform a line-search to have a solution that reduces the residual. 
        maxi = 10            ! If maxi=0, there is no line-search, only clamps. 
        alfa = 0.0001
        
        IF (MAXCORR<=clamp) THEN
            ! In this case, we don't apply any clamp and accept the solution as it is.
            DO l = 0,M
                Fp(l)  = CopyFp(l)  + ModFp(l)
                Psi(l) = CopyPsi(l) + ModPsi(l)
                Fn(l)  = CopyFn(l)  + ModFn(l)
            END DO
            CALL FillF
            RETURN
        ELSE
            ! In this case, we re-scale the solution vector using linesearch
            DO j = 0, maxi
                lambda = 0.5**j*clamp/MAXCORR
            
                DO l = 0,M
                    Fp(l)  = CopyFp(l)  + lambda*ModFp(l)
                    Psi(l) = CopyPsi(l) + lambda*ModPsi(l)
                    Fn(l)  = CopyFn(l)  + lambda*ModFn(l)
                END DO
        
                CALL FillF
                sum1 = 0
                sum2 = 0
                sum3 = 0
                DO k=0, M
                    sum1 = sum1 + ABS(f(3*k+1))**2
                    sum2 = sum2 + ABS(f(3*k+2)*XD/t0)**2    
                    sum3 = sum3 + ABS(f(3*k+3))**2
                END DO
                sum1  = SQRT(sum1+sum2+sum3)

                IF (outBacktrack) WRITE(ou, '(1I12, 3g14.4)') j, sum0, sum1, lambda
                IF (sum1 < sum0*(1-alfa)) THEN
                    sum0 = sum1
                    RETURN
                END IF 
        
            END DO
            
        END IF

    END SUBROUTINE backtracking
!-------------------------------------------------
    SUBROUTINE SolveNonLin(sum, sum1, sum2, sum3, niter, info, OutputLevel)        
        
        REAl(KIND=16) :: Jtot, Jsrh, Jrad, Jaug, Jsur, Jn, Jp
        REAl(KIND=16) :: sum1, sum2, sum3, sum, sumOld, VeryOldSum, MaxCorrection
        INTEGER :: niter, info, HIconv
        INTEGER :: k
        INTEGER :: OutputLevel
        
        ! We calculate the residual of the starting condition.
!         IF (AccountShifts) CALL BGN
        CALL FillF
        sum1 = 0
        sum2 = 0
        sum3 = 0
        DO k=0, M
            sum1 = sum1 + ABS(f(3*k+1))**2
            sum2 = sum2 + ABS(f(3*k+2)*XD/t0)**2    
            sum3 = sum3 + ABS(f(3*k+3))**2
        END DO
        sum  = SQRT(sum1+sum2+sum3)
        sum1 = SQRT(sum1)
        sum2 = SQRT(sum2)
        sum3 = SQRT(sum3)
        
           niter = 0
        info = 0
        Jtot = 0
        MaxCorrection = 0

IF(OutputLevel>=2)WRITE(ou,*)' '
IF(OutputLevel>=2)WRITE(ou,*)'     nit  Jtot (A/m^2)          sum          sum1         sum2         sum3        MaxCorr     Info'
IF(OutputLevel>=2)WRITE(ou,'(1I10,6g14.4,1I10)') niter, Jtot, sum, sum1, sum2, sum3, MaxCorrection, info
        DO WHILE (info==0)   ! The loop to solve the DD equations
            
            sumOld = sum
            niter = niter + 1
        
            ! The core of the solver
            CALL FillF
            CALL FillJacobian
            CALL SolveLin(Jac, 3*M+3, 5, 5, 18003, 11, f)
            CALL backtracking(sum, .FALSE., MaxCorrection)    

            ! We also calculate the partial residuals for the poisson equation and the two continuity equations. 
            sum1 = 0
            sum2 = 0
            sum3 = 0
            DO k=0, M
                sum1 = sum1 + ABS(f(3*k+1))**2
                sum2 = sum2 + ABS(f(3*k+2)*XD/t0)**2    
                sum3 = sum3 + ABS(f(3*k+3))**2
            END DO
            sum  = SQRT(sum1+sum2+sum3)
            sum1 = SQRT(sum1)
            sum2 = SQRT(sum2)
            sum3 = SQRT(sum3)

            ! Evaluate the conditions to finalize the loop.
            IF (niter >= nitermax) THEN
                info = 1
            ELSE IF (sum < Atol) THEN
                info = 2
            ELSE IF ((ABS((sum-sumOld)/sumOld) < Rtol).AND.(niter>1)) THEN
                info = 3
            END IF    
            
            ! And we update the current
            Currents(2) = 0
            Currents(3) = 0
            Currents(4) = 0
            
            ! Depending on having p-on-n or n-on-p, the surface recombination is a bit different
            IF (Nd(M)>Na(M)) THEN    ! p-on-n        
                Currents(5) = q*Snfront*( n(0)-fneq )
                Currents(6) = +q*Spback*( p(M)-ni(M)**2/bneq )
            ELSE IF (Nd(M)<Na(M)) THEN    ! n-on-p
                Currents(5) = +q*Spfront*( p(0)-ni(0)**2/fneq ) 
                Currents(6) = q*Snback*( n(M)-bneq )
            ELSE
                Currents(5) = 0
                Currents(6) = 0
            END IF
            
            Currents(1) = Currents(5) + Currents(6)
            DO k = 0, M-1
                Currents(2) = Rsrh(k) + Currents(2)
                Currents(3) = Rrad(k) + Currents(3)
                Currents(4) = Raug(k) + Currents(4)
                Currents(1) = Currents(1) + GR(k)
            END DO

            IF (OutputLevel>=2) WRITE(ou, '(1I10, 6g14.4, 1I10)') niter, Currents(1), sum, sum1, sum2, sum3, MaxCorrection, info

          END DO
        
    END SUBROUTINE SolveNonLin
!-------------------------------------------------
! Calculation of the functions (f), the jacobian (Jac), and the generation/recombination current
!-------------------------------------------------
    SUBROUTINE FillF                                
    ! This subroutine calculates the continuity equations (f) at each node.
    ! They should be zero in the steady state.
    !
    INTEGER :: k
    REAl(KIND=16) :: xx, T1, T2, T3
    REAl(KIND=16) :: funcp(0:M), funcn(0:M)

    !Update the potentials, carrier densities and the recombination rate
    DO k = 0, M            
        CALL ModPotential(k)
        CALL Carriers(k)
        Rho(k)=q*( p(k)-n(k)+Nd(k)-Na(k) )
    END DO
    CALL GR_sub    

    DO k = 0, M-1
        funcn(k) = -q*Mun(k)/b/dX(k)*n(k+1)*Z(Cn(k+1)-Cn(k))*(EXP(Fn(k+1)-Fn(k)) - 1) - GR(k)
        funcp(k) = -q*Mup(k)/b/dX(k)*p(k)  *Z(Cp(k+1)-Cp(k))*(EXP(Fp(k+1)-Fp(k)) - 1) + GR(k)
    END DO

    ! At each internal node, currents of adjacet elements must balance
    DO k = 1, M-1
        ! Hole current equation
        f(3*k+1) = funcp(k-1) - funcp(k) - GR(k-1)

        ! Poisson equation
        T1 = (Epsi(k+1)+Epsi(k))/dX(k)   * (Psi(k+1)-Psi(k))
        T2 = (Epsi(k)+Epsi(k-1))/dX(k-1) * (Psi(k)-Psi(k-1))
        T3 = b*Rho(k)*(dX(k)+dX(k-1))

         f(3*k+2) = T1 - T2 + T3

        ! Electron current equation
        f(3*k+3) = funcn(k-1) - funcn(k) + GR(k-1)

!         f(3*k+1) = -f(3*k+1)    ! In "SolveLin" subroutine, f enters negative; we change it now.
!         f(3*k+2) = -f(3*k+2)
!         f(3*k+3) = -f(3*k+3)
    END DO

    ! And finally, the boundary conditions
    ! At front surface
    f(1) = -funcp(0) - (1-EQ)*(1-OCn)*q*Spfront*( p(0)-fpeq ) 
    f(2) = (Epsi(1)+Epsi(0))/dX(0) * ((1-SC)*Psi(1) + SC*Vap - Psi(0))
    f(3) = -funcn(0) + (1-EQ)*(1-OCp)*q*Snfront*( n(0)-fneq ) 
    
!     f(1) = -f(1)
!     f(2) = -f(2)
!     f(3) = -f(3)
    
    ! At back surface
    f(3*M+1) = funcp(M-1) - GR(k-1) - (1-EQ)*(1-OCp)*q*Spback*( p(M)-bpeq )    
    f(3*M+2) = -(Epsi(M)+Epsi(M-1))/dX(M-1) * Psi(M)
    f(3*M+3) = funcn(M-1) + GR(k-1) + (1-EQ)*(1-OCn)*q*Snback*( n(M)-bneq )
    
!     f(3*M+1) = -f(3*M+1)
!     f(3*M+2) = -f(3*M+2)
!     f(3*M+3) = -f(3*M+3)
    
    f(1:3*M+3) = -f(1:3*M+3)
    
    END SUBROUTINE FillF
!-------------------------------------------------
    SUBROUTINE FillJacobian                           
    ! This subroutine calculates the Jacobian of "f" at each node.
    !
        INTEGER :: i, j, k
        REAl(KIND=16) :: dfuncp(6,0:M), dfuncn(6,0:M), dgrec(6,0:M)
        REAl(KIND=16)    :: urec, urad, term
                                
        DO k = 0, M-1

!             First we calculate the derivatives of the GR term
!             Calculate urec at k
            IF (k == 0) THEN
                term = tn(k)*(p(k)+ni(k)) + tp(k)*(n(k)+ni(k))
                urec = (n(k)*p(k)-ni(k)**2)/term
            END IF

            ! Derivatives w.r.t. (k)
            ! dg/dfp(k)
            dgrec(1,k) =               SRH*q*dX(k)/2.0 * (  n(k)*p(k)-tn(k)*p(k)*urec  )/term
            dgrec(1,k) = dgrec(1,k) + RAD*q*dX(k)/2.0 * Brad(k)*n(k)*p(k)
            dgrec(1,k) = dgrec(1,k) + AUG*q*dX(k)/2.0 * ( CCp(k)*p(k)*(n(k)*p(k)-ni(k)**2) + (CCn(k)*n(k)+CCp(k)*p(k))*n(k)*p(k) )
            ! dg/dPsi(k)
            dgrec(2,k) =               SRH*q*dX(k)/2.0 * (  tn(k)*p(k)-tp(k)*n(k) )*urec/term
            dgrec(2,k) = dgrec(2,k) + AUG*q*dX(k)/2.0 *  (CCn(k)*n(k)-CCp(k)*p(k))*(n(k)*p(k)-ni(k)**2)
            ! dg/dfn(k)
            dgrec(3,k) =               SRH*q*dX(k)/2.0 * (  tp(k)*n(k)*urec-n(k)*p(k)  )/term
            dgrec(3,k) = dgrec(3,k) - RAD*q*dX(k)/2.0 * Brad(k)*n(k)*p(k)
            dgrec(3,k) = dgrec(3,k) - AUG*q*dX(k)/2.0 * ( CCn(k)*n(k)*(n(k)*p(k)-ni(k)**2) + (CCn(k)*n(k)+CCp(k)*p(k))*n(k)*p(k) )

!             Now calculate urec at k+1
            term = tn(k+1)*(p(k+1)+ni(k+1)) + tp(k+1)*(n(k+1)+ni(k+1))
            urec = (n(k+1)*p(k+1)-ni(k+1)**2)/term

            ! Derivatives w.r.t. (k+1)
            ! dg/dfp(k+1)
            dgrec(4,k) =               SRH*q*dX(k)/2.0 * (  n(k+1)*p(k+1)-tn(k+1)*p(k+1)*urec  )/term
            dgrec(4,k) = dgrec(4,k) + RAD*q*dX(k)/2.0 * Brad(k+1)*n(k+1)*p(k+1)
            dgrec(4,k) = dgrec(4,k) + AUG*q*dX(k)/2.0 * (CCp(k+1)*p(k+1)*(n(k+1)*p(k+1)-ni(k+1)**2) &
                                    + (CCn(k+1)*n(k+1)+CCp(k+1)*p(k+1))*n(k+1)*p(k+1) )
            ! dg/dPsi(k+1)
            dgrec(5,k) =               SRH*q*dX(k)/2.0 * (  tn(k+1)*p(k+1)-tp(k+1)*n(k+1) )*urec/term
            dgrec(5,k) = dgrec(5,k) + AUG*q*dX(k)/2.0 *  (CCn(k+1)*n(k+1)-CCp(k+1)*p(k+1))*(n(k+1)*p(k+1)-ni(k+1)**2)
            ! dg/dfn(k+1)
            dgrec(6,k) = SRH*q*dX(k)/2.0 * (  tp(k+1)*n(k+1)*urec-n(k+1)*p(k+1)  )/term
            dgrec(6,k) = dgrec(6,k) - RAD*q*dX(k)/2.0 * Brad(k+1)*n(k+1)*p(k+1)
            dgrec(6,k) = dgrec(6,k) - AUG*q*dX(k)/2.0 * (CCn(k+1)*n(k+1)*(n(k+1)*p(k+1)-ni(k+1)**2) &
                                    + (CCn(k+1)*n(k+1)+CCp(k+1)*p(k+1))*n(k+1)*p(k+1) )

            ! Calculate derivatives of "funcn" and "funcp" w.r.t. the potentials
            ! funcn(k) = -q*Mun(k)/b/dX(k)*n(k+1)*Z(Cn(k+1)-Cn(k))*(EXP(Fn(k+1)-Fn(k)) - 1) - GR(k)*Y(Cn(k+1)-Cn(k))
            ! Derivatives w.r.t. (k)
            ! dfuncn/dfp(k)
            dfuncn(1, k) = - dgrec(1,k)*Y(Cn(k+1)-Cn(k))
            ! dfuncn/dPsi(k)
            dfuncn(2, k) = q*Mun(k)/b/dX(k)*n(k+1)*dZ(Cn(k+1)-Cn(k))*(EXP(Fn(k+1)-Fn(k))-1)
            dfuncn(2, k) = dfuncn(2, k) - dgrec(2,k)*Y(Cn(k+1)-Cn(k)) + GR(k)*dY(Cn(k+1)-Cn(k))
            ! dfuncn/dfn(k)
            dfuncn(3, k) = q*Mun(k)/b/dX(k)*n(k+1)*Z(Cn(k+1)-Cn(k))*EXP(Fn(k+1)-Fn(k)) - dgrec(3,k)*Y(Cn(k+1)-Cn(k))
            ! Derivatives w.r.t. (k+1)
            ! dfuncn/dfp(k+1)
            dfuncn(4, k) = - dgrec(4,k)*Y(Cn(k+1)-Cn(k))
            ! dfuncn/dPsi(k+1)
            dfuncn(5, k) = - q*Mun(k)/b/dX(k)*n(k+1)*( Z(Cn(k+1)-Cn(k))+dZ(Cn(k+1)-Cn(k)) )*(EXP(Fn(k+1)-Fn(k))-1)
            dfuncn(5, k) = dfuncn(5, k) - dgrec(5,k)*Y(Cn(k+1)-Cn(k)) - GR(k)*dY(Cn(k+1)-Cn(k))
            ! dfuncn/dfn(k+1)
            dfuncn(6, k) = - q*Mun(k)/b/dX(k)*n(k+1)*  Z(Cn(k+1)-Cn(k)) - dgrec(6,k)*Y(Cn(k+1)-Cn(k))

            ! funcp(k) = -q*Mup(k)/b/dX(k)*p(k)  *Z(Cp(k+1)-Cp(k))*(EXP(Fp(k+1)-Fp(k)) - 1) + GR(k)*Y(Cp(k)-Cp(k+1))
            ! Derivatives w.r.t. (k)
            ! dfuncp/dfp(k)
            dfuncp(1, k) = q*Mup(k)/b/dX(k)*p(k)*  Z(Cp(k+1)-Cp(k))    + dgrec(1,k)*Y(Cp(k)-Cp(k+1))
            ! dfuncp/dPsi(k)
            dfuncp(2, k) = q*Mup(k)/b/dX(k)*p(k)*( Z(Cp(k+1)-Cp(k))+dZ(Cp(k+1)-Cp(k)) )*(EXP(Fp(k+1)-Fp(k))-1)
            dfuncp(2, k) = dfuncp(2, k) + dgrec(2,k)*Y(Cp(k)-Cp(k+1)) + GR(k)*dY(Cp(k)-Cp(k+1))
            ! dfuncp/dfn(k)
            dfuncp(3, k) = dgrec(3,k)*Y(Cp(k)-Cp(k+1))
            ! Derivatives w.r.t. (k+1)
            ! dfuncp/dfp(k+1)
            dfuncp(4, k) = - q*Mup(k)/b/dX(k)*p(k)* Z(Cp(k+1)-Cp(k))*EXP(Fp(k+1)-Fp(k)) + dgrec(4,k)*Y(Cp(k)-Cp(k+1))
            ! dfuncp/dPsi(k+1)
            dfuncp(5, k) = - q*Mup(k)/b/dX(k)*p(k)*dZ(Cp(k+1)-Cp(k))*(EXP(Fp(k+1)-Fp(k))-1)
            dfuncp(5, k) = dfuncp(5, k) + dgrec(5,k)*Y(Cp(k)-Cp(k+1)) - GR(k)*dY(Cp(k)-Cp(k+1))
            ! dfuncp/dfn(k+1)
            dfuncp(6, k) = dgrec(6,k)*Y(Cp(k)-Cp(k+1))

        END DO

        ! For each internal node, calculate derivatives of hole current function, 
        ! Poisson's equation & electron current function w.r.t. each of nine potentials
 
        DO k = 1, M-1
            i = 3*k+1

            ! Derivatives of hole current continuity equation
            !   f(3*k+1) = funcp(k-1) - funcp(k) - GR(k-1)
            ! w.r.t. (k-1)
            Jac(i, 3) = dfuncp(1, k-1) - dgrec(1, k-1)
            Jac(i, 4) = dfuncp(2, k-1) - dgrec(2, k-1)
            Jac(i, 5) = dfuncp(3, k-1) - dgrec(3, k-1)
            ! w.r.t. (k)
            Jac(i, 6) = dfuncp(4, k-1) - dfuncp(1, k) - dgrec(4, k-1)
            Jac(i, 7) = dfuncp(5, k-1) - dfuncp(2, k) - dgrec(5, k-1)
            Jac(i, 8) = dfuncp(6, k-1) - dfuncp(3, k) - dgrec(6, k-1)
            ! w.r.t. (k+1)
            Jac(i, 9)  =                - dfuncp(4, k) 
            Jac(i, 10) =               - dfuncp(5, k) 
            Jac(i, 11) =               - dfuncp(6, k) 
            
            i = i + 1
            ! Derivatives of Poisson's equation
            !
!             T1 = (Epsi(k+1)+Epsi(k))/dX(k)   * (Psi(k+1)-Psi(k))
!             T2 = (Epsi(k)+Epsi(k-1))/dX(k-1) * (Psi(k)-Psi(k-1))
!             T3 = b*Rho(k)*(dX(k)+dX(k-1))
!
!              f(3*k+2) = T1 - T2 + T3            

!             ! w.r.t. (k-1)
            Jac(i, 2) = 0.0
            Jac(i, 3) = (Epsi(k)+Epsi(k-1))/dX(k-1)
            Jac(i, 4) = 0.0
            ! w.r.t. (k)
            Jac(i, 5) = b*q*p(k)*(dX(k)+dX(k-1))
            Jac(i, 6) = -(Epsi(k+1)+Epsi(k))/dX(k)-(Epsi(k)+Epsi(k-1))/dX(k-1)-b*q*(p(k)+n(k))*(dX(k)+dX(k-1))
            Jac(i, 7) = b*q*n(k)*(dX(k)+dX(k-1))
            ! w.r.t. (k+1)
            Jac(i, 8)  = 0.0
            Jac(i, 9)  = (Epsi(k+1)+Epsi(k))/dX(k)
            Jac(i, 10) = 0.0
            
        
            i = i + 1
            ! Derivatives of electron continuity equation
            !    f(3*k+3) = funcn(k-1) - funcn(k) + GR(k-1)
            ! w.r.t. (k-1)
            Jac(i, 1) = dfuncn(1, k-1) + dgrec(1, k-1)
            Jac(i, 2) = dfuncn(2, k-1) + dgrec(2, k-1)
            Jac(i, 3) = dfuncn(3, k-1) + dgrec(3, k-1)
            ! w.r.t. (k)
            Jac(i, 4) = dfuncn(4, k-1) - dfuncn(1, k) + dgrec(4, k-1)
            Jac(i, 5) = dfuncn(5, k-1) - dfuncn(2, k) + dgrec(5, k-1)
            Jac(i, 6) = dfuncn(6, k-1) - dfuncn(3, k) + dgrec(6, k-1)
            ! w.r.t. (k+1)
            Jac(i, 7) =                - dfuncn(4, k) 
            Jac(i, 8) =                - dfuncn(5, k) 
            Jac(i, 9) =                - dfuncn(6, k) 
            
        END DO
        
        
        ! Finally, we calculate the entries of the boundary conditions
        ! At front surface
        !     f(1) = -funcp(0) - (1-EQ)*(1-OC)*q*Spfront*( p(0)-ni(0)**2/fneq )
        i = 1
        ! w.r.t. (k)
        Jac(i, 6) = - dfuncp(1, 0) - (1-EQ)*(1-OCn)*q*Spfront*p(0) 
        Jac(i, 7) = - dfuncp(2, 0) + (1-EQ)*(1-OCn)*q*Spfront*p(0)
        Jac(i, 8) = - dfuncp(3, 0)
        ! w.r.t. (k+1)
        Jac(i, 9)  = - dfuncp(4, 0)
        Jac(i, 10) = - dfuncp(5, 0)
        Jac(i, 11) = - dfuncp(6, 0)
        
        !    f(2) = (Epsi(1)+Epsi(0))/dX(0) * ((1-SC)*Psi(1) + SC*Vap - Psi(0))
        i = 2
        ! w.r.t. (k)
        Jac(i, 5) = 0.0
        Jac(i, 6) = -(Epsi(1)+Epsi(0))/dX(0)
        Jac(i, 7) = 0.0
        ! w.r.t. (k+1)
        Jac(i, 8)  = 0.0
        Jac(i, 9)  = (1-SC)*(Epsi(1)+Epsi(0))/dX(0)
        Jac(i, 10) = 0.0
        
        !    f(3) = -funcn(0) + (1-EQ)*(1-OC)*q*Snfront*( n(0)-fneq )
        i = 3
        ! w.r.t. (k)
        Jac(i, 4) = - dfuncn(1, 0)
        Jac(i, 5) = - dfuncn(2, 0) + (1-EQ)*(1-OCp)*q*Snfront*n(0) 
        Jac(i, 6) = - dfuncn(3, 0) - (1-EQ)*(1-OCp)*q*Snfront*n(0)
        ! w.r.t. (k+1)
        Jac(i, 7) = - dfuncn(4, 0)
        Jac(i, 8) = - dfuncn(5, 0)
        Jac(i, 9) = - dfuncn(6, 0)


    
        ! At back surface
        !     f(3*M+1) = funcp(M-1) - GR(M-1) - (1-EQ)*(1-OC)*q*Spback*( p(M)-ni(M)**2/bneq )    
        i = 3*M+1
        ! w.r.t. (k-1)
        Jac(i, 3) = dfuncp(1, M-1) - dgrec(1, M-1)
        Jac(i, 4) = dfuncp(2, M-1) - dgrec(2, M-1)
        Jac(i, 5) = dfuncp(3, M-1) - dgrec(3, M-1) 
        ! w.r.t. (k)
        Jac(i, 6) = dfuncp(4, M-1) - dgrec(4, M-1) - (1-EQ)*(1-OCp)*q*Spback*p(M) 
        Jac(i, 7) = dfuncp(5, M-1) - dgrec(5, M-1) + (1-EQ)*(1-OCp)*q*Spback*p(M) 
        Jac(i, 8) = dfuncp(6, M-1) - dgrec(6, M-1)

        !    f(3*M+2) = -(Epsi(M)+Epsi(M-1))/dX(M-1) * (Psi(M)-OC*Psi(M-1)) 
        i = 3*M+2
        ! w.r.t. (k-1)
        Jac(i, 2) = 0.0
        Jac(i, 3) = 0.0
        Jac(i, 4) = 0.0
        ! w.r.t. (k)
        Jac(i, 5) = 0.0
        Jac(i, 6) = -(Epsi(M)+Epsi(M-1))/dX(M-1)
        Jac(i, 7) = 0.0

        !    f(3*M+3) = funcn(M-1) + GR(M-1) + (1-EQ)*(1-OC)*q*Snback*( n(M)-bneq )
        i = 3*M+3
        ! w.r.t. (k-1)
        Jac(i, 1) = dfuncn(1, M-1) + dgrec(1, M-1)
        Jac(i, 2) = dfuncn(2, M-1) + dgrec(2, M-1)
        Jac(i, 3) = dfuncn(3, M-1) + dgrec(3, M-1)
        ! w.r.t. (k)
        Jac(i, 4) = dfuncn(4, M-1) + dgrec(4, M-1) 
        Jac(i, 5) = dfuncn(5, M-1) + dgrec(5, M-1) + (1-EQ)*(1-OCn)*q*Snback*n(M)
        Jac(i, 6) = dfuncn(6, M-1) + dgrec(6, M-1) - (1-EQ)*(1-OCn)*q*Snback*n(M)


    END SUBROUTINE FillJacobian
!-------------------------------------------------
    SUBROUTINE GR_sub                                
    
        INTEGER :: k, j
        REAl(KIND=16) :: urec, urad, uaug
        
        DO k = 0, M-1
            
            IF (k == 0) THEN    
                urec = (n(k)*p(k)-ni(k)**2) / (tn(k)*(p(k)+ni(k)) + tp(k)*(n(k)+ni(k)))
                urad = n(k)*p(k)-ni(k)**2
                uaug = (CCn(k)*n(k)+CCp(k)*p(k)) * (n(k)*p(k)-ni(k)**2)
            END IF
            Rsrh(k) = urec 
            Rrad(k) = urad
            Raug(k) = uaug

            ! Now, calculate urec and urad at k+1
            urec = (n(k+1)*p(k+1)-ni(k+1)**2) / (tn(k+1)*(p(k+1)+ni(k+1)) + tp(k+1)*(n(k+1)+ni(k+1)))
            urad = n(k+1)*p(k+1)-ni(k+1)**2 
            uaug = (CCn(k+1)*n(k+1)+CCp(k+1)*p(k+1)) * (n(k+1)*p(k+1)-ni(k+1)**2)
        
            ! The average sheet R rate = average of volume R rate x element width
            Rsrh(k) = SRH* q*(Rsrh(k) +  urec)*dX(k)/2.0            ! Average sheet SRH recombination
            Rrad(k) = RAD* q*Brad(k)*(Rrad(k) + urad)*dX(k)/2.0        ! Average sheet radiative recombination
            Raug(k) = AUG* q*(Raug(k) +  uaug)*dX(k)/2.0            ! Average sheet Auger recombination        
            
            GR(k) =  Rsrh(k) +  Rrad(k) + Raug(k) - G(k)
            
        END DO
                    
    END SUBROUTINE GR_sub
!-------------------------------------------------
    SUBROUTINE Generation(wl)                        
        REAl(KIND=16) :: photonfluxini, PFspectrumini(0:NumWL), PF, PFWL(0:NumWL), TempG, sigma, GG(0:M)
        INTEGER :: i, j, k
        INTEGER, OPTIONAL :: wl
        
        ! No absorption
        IF (GEN == 0) THEN
            G(0:M) = 0.0
        
        ! Broadband absorption 
        ELSE IF (.NOT.SingleWL) THEN
            
            PFWL = q*PFspectrum(0:NumWL)
            G(0:M) = 0.0
            
            DO j = 0, NumWl-1
                GG(0:M) = (AbsProfile(0:M,j) * PFWL(j) + AbsProfile(0:M,j+1) * PFWL(j+1)) * (AbsLibrary(-1,j+1)-AbsLibrary(-1,j))/2
                G(0:M) = G(0:M) + GG(0:M)          
            END DO 
            
            G(0:M-1) = (G(0:M-1) + G(1:M)) * dX(0:M-1) / 2                 
        
        ! Single wavelength
        ELSE
            PF = q*PhotonFlux
        
            G(0:M) = 0.0    
            G(0:M-1) = PF*(AbsProfile(0:M-1,wl) + AbsProfile(1:M,wl)) * dX(0:M-1) / 2      
        END IF
    
    END SUBROUTINE Generation
!-------------------------------------------------    
! Running modes
!-------------------------------------------------
    FUNCTION OutputCode(info)                        
        INTEGER :: info
        CHARACTER(50) :: OutputCode
        
        IF (info==1) THEN
            OutputCode = 'Not converging: Reached Maximum iterations.'
            ERROR STOP OutputCode
        ELSE IF (info==2) THEN
            OutputCode = 'Reached Absolute Tolerance.'
        ELSE IF (info==3) THEN
            OutputCode = 'Reached Relative Tolerance.'
        END IF

    END FUNCTION OutputCode
!-------------------------------------------------
    SUBROUTINE Equilibrium(OutputLevel)                
        ! Solve the DD equations at 0V in the Dark         
    
        REAl(KIND=16) :: start_time, end_time
        REAl(KIND=16) :: sum, Jtot
        INTEGER :: info, GENtemp
        REAl(KIND=16) :: dum1, dum2, dum3, dum4, dum5, dum7    ! Dummy variables that must be different. 
        INTEGER :: dum6
        INTEGER :: OutputLevel

        CALL open_log()
                
        IF (OutputLevel>=1) THEN
                CALL CPU_TIME (start_time)
            WRITE(ou,*) ' '
            WRITE(ou,*) 'Starting EQUILIBRIUM... '
        END IF

        GENtemp = GEN
        GEN = 0
        EQ = 1
        SC = 0
        OC = 0
        OCp = 0
        OCn = 0
        
        CALL SolveNonLin(sum, dum3, dum4, dum5, dum6, info, OutputLevel)
        
        IF (Dynamic) THEN 
            CALL DynamicMesh(0)
            WRITE(ou,*) 'Remeshing...  M+1 = ',  M+1, ' nodes.'
            CALL SolveNonLin(sum, dum3, dum4, dum5, dum6, info, OutputLevel)
        END IF
        
        GEN = GENtemp
        
        IF (OutputLevel>=1) THEN
            CALL CPU_TIME (end_time)
            WRITE(ou, * ) 'EQUILIBRIUM Output Code: ', OutputCode(info)
            WRITE(ou, * ) '    Res: ', REAL(sum,4)
            WRITE(ou, * ) 'Elapsed time = ', REAL((end_time - start_time), 4) , 's'    
            WRITE(ou,*) ' '
        END IF
        
        fneq = n(0)
        bneq = n(M)
        fpeq = p(0)
        bpeq = p(M)
        Vbi  = Psi(0)
        EQ = 0
        
        CALL close_log()
    
    END SUBROUTINE Equilibrium
!-------------------------------------------------
    SUBROUTINE LightSC(OutputLevel, qmode)            
        ! Solve the DD equations at 0V and with illumination.         
    
        REAl(KIND=16) :: start_time, end_time
        REAl(KIND=16) :: sum, Jtot, Jsrh, Jrad, Jaug, Jsur
        REAl(KIND=16) :: photonfluxini, PFspectrumini(0:NumWL), PF, PFWL(0:NumWL)!, TempG 
        REAl(KIND=16), DIMENSION(0:6000) :: TempG    
        INTEGER :: info, i, j, k, maxsteps, qmode
        REAl(KIND=16) :: dum1, dum2, dum3, dum4, dum5, dum7     ! Dummy variables that must be different. 
        INTEGER :: dum6
        INTEGER :: OutputLevel, SWL
        
        CALL open_log()
                
        IF (OutputLevel>=1) THEN
                CALL CPU_TIME (start_time)
            WRITE(ou,*) ' '
            WRITE(ou,*) 'Starting LIGHTSC... '
        END IF

        GEN = 1
        SC = 0
        OC = 0
        OCp = 0
        OCn = 0
        IF (qmode==1)  THEN
            SC = 1
            Vap = Vbi
        ELSE IF (qmode==-1) THEN
            IF (Nd(M)>Na(M)) THEN  ! p-on-n
                OCn=1
                OCp=0
            ELSE            ! n-on-p
                OCn=0
                OCp=1
            END IF
            OC = 1
        END IF

        CALL Generation
        TempG(0:M) = G(0:M)
        
        maxsteps = CEILING(LOG10(PhotonFlux))

        IF (OutputLevel>=1) WRITE(ou,*) '     step  Jtot (A/m^2)          Res          Res-h    Res Poisson     Res-e        Info'
        DO i = 1, maxsteps 
            G = GEN* TempG*10**(i-REAL(maxsteps,4))
            CALL SolveNonLin(sum, dum3, dum4, dum5, dum6, info, OutputLevel)
            IF (OutputLevel>=1) WRITE(ou, '(1I10, 5g14.4, 1I10)') i, Currents(1), sum, dum3, dum4, dum5, info
            
        END DO
        
        Jtot = Currents(1)

        IF (Nd(M)>Na(M)) THEN 
            Voc = (Fp(0)-Fn(M))/b
        ELSE
            Voc = (Fn(0)-Fp(M))/b
        END IF
        Isc = Jtot

        IF (OutputLevel>=1) THEN
            CALL CPU_TIME (end_time)
            WRITE(ou, * ) 'LIGHTSC Output Code: ', OutputCode(info)
            WRITE(ou, * ) '    Res: ', REAL(sum,4)
            WRITE(ou, * ) '    J: ', REAL(Isc,4), ' A/m2'
            WRITE(ou, * ) '    V: ', REAL(Voc,4), ' V'
!             WRITE(ou, * ) '    IQE: ', REAL(-Jtot/q/PhotonFlux*100,4), ' %'
            WRITE(ou, * ) 'Elapsed time = ', REAL((end_time - start_time), 4) , 's'    
            WRITE(ou,*) ' '
        END IF
        
        CALL close_log()
    
    END SUBROUTINE LightSC
!-------------------------------------------------
    SUBROUTINE RunIV(Vfin, Vstep, OutputLevel, escape)            
        ! Solve the DD equations as a function of voltage in the dark.         
    
        REAl(KIND=16) :: start_time, end_time
        REAl(KIND=8) :: Vini, Vfin, Vstep
        REAl(KIND=16) :: Vapp, Vreal, step, Vend, factor
        REAl(KIND=16) :: sum, sum1, sum2, sum3, Jtot, Jsrh, Jrad, Jaug, Jsur
        REAl(KIND=16) :: intertot(3000)
        INTEGER :: info, niter
        INTEGER :: OutputLevel, LS, ESC
        INTEGER, OPTIONAL :: escape 
        INTEGER :: GENtemp
        LOGICAL :: continue_loop = .TRUE.
        REAl(KIND=16) :: xx, yy, sol
        
        INTEGER :: i
        
        CALL open_log()
                
        continue_loop = .TRUE.

        IF (PRESENT(escape)) THEN
            ESC = escape
        ELSE
            ESC = 0
        END IF
        
        IF (OutputLevel>=1) THEN
                CALL CPU_TIME (start_time)
            WRITE(ou,*) ' '
            WRITE(ou,*) 'Starting RUNIV... '
        END IF
        
        Vapp = 0.0
        Vreal = Vapp
        Vap = Vbi
        step = REAL(Vstep,16)
        Vend = REAL(Vfin,16) + 0.1*step
        nvolt = 0
        factor = 1
        SC = 1
        continue_loop = .TRUE.
        
        IF (OutputLevel>=1) WRITE(ou,*) 'nvolt   nit   info    Vapp       Jtot (A/m^2)       Res'
        
        DO WHILE ( continue_loop )   

            ! We create a new suitable initial condition based on the one obtained for the previous voltage.
            IF (nVolt/=0) THEN                            
                IF (Na(M)<Nd(M)) THEN
                    Fp(0:M) = Fp(0:M) + b*step*factor
                ELSE
                    Fn(0:M) = Fn(0:M) + b*step*factor
                END IF
                Psi(0:M) = Psi(0:M) * Vap/Psi(0)    
            END IF    

!            IF ((Dynamic).AND.(FLOOR(2*nvolt*ABS(step))/=FLOOR(2*(nvolt+1)*ABS(step)))) THEN
!                CALL DynamicMesh(0)
!                WRITE(ou,*) 'Remeshing...  M+1 = ',  M+1, ' nodes.'
!                IF (GEN/=0) CALL Generation
!            END IF
            
            CALL SolveNonLin(sum, sum1, sum2, sum3, niter, info, OutputLevel)

            WHERE (ABS(Currents) < Atol) Currents = 0
                
            ! We fill the arrays. 
            vpoint(nVolt)        = Vapp !+ Currents(1)*Rs        ! The external voltage depends on the series resistance
            jpoint(nVolt)        = Currents(1)
            jsrhpoint(nVolt)     = Currents(2)
            jradpoint(nVolt)    = Currents(3)
            jaugpoint(nVolt)    = Currents(4)
            jsurpoint(nVolt)    = Currents(5) + Currents(6)
            residual(nVolt)        = sum
            
            IF (OutputLevel>=1) WRITE(ou, '(3I6, 6g14.6)') nvolt, niter, info, Vapp + Currents(1)*Rs, jpoint(nVolt), sum  
            
            IF (((GEN==1).AND.(ESC/=0)).AND.( ABS(jpoint(nVolt)*vpoint(nVolt)) >= Pmax) ) THEN
                Imax = -jpoint(nVolt)
                Vmax = vpoint(nVolt)
                Pmax = Imax*Vmax
            END IF
            
            IF (((GEN==1).AND.(ESC/=0)).AND.(Currents(1)>0)) THEN
                
                Voc = vpoint(nVolt-1) - jpoint(nVolt-1) * (vpoint(nVolt-1) - vpoint(nVolt)) / (jpoint(nVolt-1) - jpoint(nVolt))
                FF = Pmax/(-Isc*Voc)*100
                
                WRITE(ou, * ) '     '
                EXIT
            END IF
            
!            ! To take into account corrections due to series resistence
!            IF (GEN/=0) THEN
!                factor = 1/MIN(20.0, 1.0 + 20*ABS(LOG( ABS( jpoint(nVolt)/Isc ) ) ) )
!            ELSE
!                IF (nvolt>0) factor = (step/( vpoint(nVolt) - vpoint(nVolt-1) ) )**2
!            END IF
            
            nvolt = nvolt+1            
            Vapp = Vapp + step*factor
            Vap = Vbi + Vapp*b
                    
            IF (Vstep>0) THEN
                Vreal = MAX(Vapp, vpoint(nVolt-1) )  ! In case series resistance is too high
                continue_loop = ( Vreal <= Vend )
            ELSE 
                Vreal = MIN(Vapp, vpoint(nVolt-1) )  ! In case series resistance is too high
                continue_loop = ( Vreal >= Vend )
            END IF
            
        END DO   
        
        IF (OutputLevel>=1) THEN
            CALL CPU_TIME (end_time)
            WRITE(ou, * ) 'Elapsed time = ', REAL((end_time - start_time), 4) , 's'    
            WRITE(ou,*) ' '
        END IF
    
        CALL close_log()
            
    END SUBROUTINE RunIV
!-------------------------------------------------
    SUBROUTINE RunIQE(OutputLevel)        

        ! External variables
        INTEGER :: OutputLevel
    !         REAl(KIND=8), OPTIONAL :: Vfin, Vstep
    
        ! Internal variables
        REAl(KIND=16) :: start_time, end_time
        REAl(KIND=16) :: sumsum, Jtot, PF, Jsrh, Jrad, Jaug, Jsur
        REAl(KIND=16) :: dum1, dum2, dum3, dum4, dum5, dum7     ! Dummy variables that must be different. 
        INTEGER :: dum6
        REAl(KIND=16) :: Jbias
        INTEGER :: info, niter
        INTEGER :: i, k, maxsteps
        REAl(KIND=16), DIMENSION(0:6000) :: TempG    

        CALL open_log()
            
        SC = 1

        TempG(0:M) = G(0:M)
        CurrentsBias = Currents
        PhotonFlux = PhotonFlux*1e-8
    
        SingleWL = .TRUE.    
    
        IF (OutputLevel>=1) THEN
               CALL CPU_TIME (start_time)
            WRITE(ou,*) ' '
            WRITE(ou,*) 'Starting RunIQE... '
        END IF
    
        IF (OutputLevel>=1) WRITE(ou,*) 'WLindex  nit   info    Wavelength (nm)    IQE (%)       Res'
        DO i = 0, NumWL-1
            CALL Generation(i)
            G(0:M) = G(0:M)+TempG(0:M)
        
            CALL SolveNonLin(sumsum, dum4, dum5, dum7, niter, info, OutputLevel)

            IQE(i) = -(Currents(1)-CurrentsBias(1))/q/PhotonFlux
            IQEsrh(i) = (Currents(2)-CurrentsBias(2))/q/PhotonFlux
            IQErad(i) = (Currents(3)-CurrentsBias(3))/q/PhotonFlux
            IQEaug(i) = (Currents(4)-CurrentsBias(4))/q/PhotonFlux
            IQEsurf(i) = (Currents(5)-CurrentsBias(5))/q/PhotonFlux
            IQEsurb(i) = (Currents(6)-CurrentsBias(6))/q/PhotonFlux
            IF (OutputLevel>=1) WRITE(ou,'(3I6,6g14.6)') i, niter, info, AbsLibrary(-1,i)/1e-9,IQE(i)*100,sumsum
        END DO    
    
        IF (OutputLevel>=1) THEN
            CALL CPU_TIME (end_time)
            WRITE(ou, * ) 'Elapsed time = ', REAL((end_time - start_time), 4) , 's'    
            WRITE(ou,*) ' '
        END IF
    
        CALL close_log()
            
    END SUBROUTINE RunIQE
!-------------------------------------------------    
END MODULE


