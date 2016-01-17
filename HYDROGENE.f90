!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAME PRINCIPAL
program HYDROGENE
  implicit none
  integer,parameter           :: R = 8
  real(kind = 8), parameter   :: hbar = 1.054571726E-34, m = 9.109E-31, e = 1.602E-19, a0 = 0.52917721092, EPS0 = 8.85418782E-12, &
                                 pi = 3.1415926535
  integer                     :: NUMS, INDX, NCRL, CORL, NMBR, LIMT
  real(kind = R)              :: AOLD, ASOL = .75*a0, DIST, ALPH = 2*a0, BETA, PASM, VRND(6), V, T, SINI, SFIN, STEP, FR12, PHI,  &
                                 PHI0, R12, R1L, R1R, R2L, R2R, THT1, THT2, PHI1, PHI2, LPH1, LPH2, LPCf, LPHT, GRD1, GRD2, &
                                 EVAR, U, R1(3), R2(3), RL(3) = 0., R1NEW(3), R2NEW(3)
  real(kind = R), allocatable :: EPSN(:)

!! OUVERTURE DU FICHIER D'ENTREE ET LECTURE DES DONNEES
    open(100,file = 'INPUT.dat')
    read(100,*) SINI, SFIN, STEP, BETA, CORL, NMBR, PASM
    close(100)

!! ALLOCATION DU VECTEUR EPSN CONTENANT LES ENERGIES DE CHAQUE CONFIGURATION
    allocate(EPSN(NMBR))

!! OUVERTURE DES FICHIERS DE SORTIE (SOLUTION DE L'EQUATION INTRINSEQUE, POSITION DES ELECTRONS ET COURBE DE MORSE)
    open(101,file = 'SOLUTION.dat')
    open(102,file = 'COORDONNEES.dat')
    open(103,file = 'MORSE.dat')

!! INITIALISATION DU GENERATEUR ALEATOIRE
    call init_random_seed()

!! BOUCLE SUR LA DISTANCE INTER-ATOMIQUE
    do DIST= SINI, SFIN, STEP
!!!! CALCUL DE LA SOLUTION DE L'EQUATION INTRINSEQUE
        AOLD = 0.
        do while(abs(ASOL-AOLD).GE.1.E-6)
            AOLD = ASOL
            ASOL = a0/(1+exp(-DIST/AOLD))
        enddo
        write(101,'(2E13.5)') DIST, ASOL

!! INITIALISATION DU SYSTEME H2
        rL(3) = -DIST/2
        R1    = 0.
        R2    = 0.
        PHI0  = 4/exp(DIST)

!! BOUCLE SUR LE NOMBRE DE CONFIGURATION POUR CHAQUE DISTANCE
        do INDX = 1,NMBR
!! BOUCLE SUR LE NOMBRE DE CORRELATION
            do NCRL = 1,CORL
                do


!! CALCUL D'UNE NOUVELLE CONFIGURATION
!! GENERATION DES NOUVELLES POSITIONS
                    call random_number(VRND)
                    THT1 = VRND(1)*2._r*pi
                    THT2 = VRND(2)*2._r*pi
                    PHI1 = VRND(3)*pi
                    PHI2 = VRND(4)*pi

                    R1NEW = R1 + PASM*[cos(THT1)*sin(PHI1),sin(THT1)*sin(PHI1),cos(PHI1)]
                    R2NEW = R2 + PASM*[cos(THT2)*sin(PHI2),sin(THT2)*sin(PHI2),cos(PHI2)]

!! CALCUL DES NOUVELLES DISTANCES
                    R1R = norm2(R1NEW+rL)
                    R1L = norm2(R1NEW-rL)
                    R2R = norm2(R2NEW+rL)
                    R2L = norm2(R2NEW-rL)
                    R12 = norm2(R1NEW-R2NEW)

!! CALCUL DE LA NOUVELLE FONCTION D'ONDE
                    FR12 = exp(R12/(ALPH*(1+BETA*R12)))
                    PHI1 = exp(-R1L/ASOL) + exp(-R1R/ASOL)
                    PHI2 = exp(-R2L/ASOL) + exp(-R2R/ASOL)
                    PHI = PHI1*PHI2*FR12

!! ACCEPTATION OU REJET DE LA NOUVELLE FONCTION D'ONDE
                    if (VRND(5).LT.(PHI/PHI0)**2) exit
                enddo
                R1   = R1NEW
                R2   = R2NEW
                PHI0 = PHI
            enddo

!! ECRITURE DANS LE FICHIER ssi S = 3 A
            if(abs(DIST-3).LT.STEP/2) write(102,'(6E13.5)') R1, R2


!! CALCUL DE L'ENERGIE DE LA NOUVELLE CONFIGURATION NON IGNOREE
!! CALCUL DE L'ENERGIE POTENTIELLE
            V    = -(e**2)/(4*pi*EPS0*1.E-10)*(1/R1L + 1/R2L + 1/R1R + 1/R2R - 1/R12)

!! CALCUL DES LAPLACIENS
            LPCf = exp(R12/(ALPH*(1+BETA*R12)))*(R12+2*ALPH*(1+BETA*R12))/(ALPH**2*R12*(1+BETA*R12)**4)
            LPH1 = -(exp(-R1L/ASOL)*(2/R1L-1/ASOL) + exp(-R1R/ASOL)*(2/R1R-1/ASOL))/ASOL
            LPH2 = -(exp(-R2L/ASOL)*(2/R2L-1/ASOL) + exp(-R2R/ASOL)*(2/R2R-1/ASOL))/ASOL

!! CALCUL DES GRADIENTS
            GRD1 = -sum((exp(-R1L/ASOL)*(R1-rL)/R1L+exp(-R1R/ASOL)*(R1+rL)/R1R)*(R1-R2)) &
                 *  exp(R12/(ALPH*(1.+BETA*R12)))/(R12*ASOL*ALPH*(1.+BETA*R12)**2)
            GRD2 = -sum((exp(-R2L/ASOL)*(R2-rL)/R2L+exp(-R2R/ASOL)*(R2+rL)/R2R)*(R1-R2)) &
                 *  exp(R12/(ALPH*(1.+BETA*R12)))/(R12*ASOL*ALPH*(1.+BETA*R12)**2)

!! CALCUL DU LAPLACIEN TOTAL
            LPHT = PHI2*FR12*LPH1 + 2*PHI1*PHI2*LPCf + 2*PHI2*GRD1 + PHI1*FR12*LPH2 + 2*PHI1*GRD2

!! CALCUL DE L'ENERGIE CINETIQUE ET DE L'ENERGIE TOTALE
            T    = -(HBAR*1.E10)**2*LPHT/(2*m)
            EPSN(INDX) = T/PHI + V
        enddo

!! CALCUL FINAL DE L'ENERGIE POUR UNE DISTANCE INTER-ATOMIQUE ET ECRITURE DANS LE FIHCIER
        U = (e**2/(4*pi*EPS0*1.E-10*DIST) + sum(EPSN)/NMBR)/e
        write(103,'(2E13.5)') DIST, U
    enddo

    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SOUS-ROUTINE INIT_RANDOM_SEED
    subroutine INIT_RANDOM_SEED()
      integer              :: INDX, N, CLCK
      integer, allocatable :: SEED(:)

        call random_seed(size = N)
        call system_clock(COUNT = CLCK)

        SEED = CLCK + 37*(/(INDX, INDX=0, N-1)/)
        call random_seed(PUT = SEED)

        deallocate(SEED)
    end subroutine

end program HYDROGENE