module turb_enhance
    use variable_precision, only: wp

    implicit none

    character(len=*), parameter, private :: ModuleName='TURB_ENHANCE'

    private

    public autoconversion_enhance, accretion_enhance, selfcollection_enhance
contains

    subroutine autoconversion_enhance(k_turb, turb_edr, mu, lam)

        implicit none

        !*******************************************************************************
        ! Autoconversion enhancement from Seifert and Beheng (2001, Atmos. Res.)       *
        !*******************************************************************************

        ! Onishi kernel (of 29 July 2015)
        real(wp), parameter :: kc1_a1 = 3.985e-03_wp
        real(wp), parameter :: kc1_a2 = 6.210e-03_wp
        real(wp), parameter :: kc1_a3 = 1.331e+00_wp
        real(wp), parameter :: kc1_b1 = 1.381e+01_wp
        real(wp), parameter :: kc1_b2 = 9.980e+00_wp
        real(wp), parameter :: kc1_b3 = 5.018e-01_wp
        real(wp), parameter :: kc1_c1 = 6.325e+00_wp
        real(wp), parameter :: kc1_c2 = -9.238e-01_wp
        real(wp), parameter :: kc1_c3 = -1.528e-01_wp
        real(wp), parameter :: kc1_bet = 2.026e-03_wp

        real(wp), intent(in) :: turb_edr, mu, lam
        real(wp), intent(inout) :: k_turb
        real(wp) :: kc_alf,kc_rad,kc_sig,kc_bet,prey,Re,diss,nu_c,D_mean

        character(len=*), parameter :: RoutineName='AUTOCONVERSION_ENHANCE'

        nu_c = 0.0_wp

        kc_alf = ( kc1_a1 + kc1_a2 * nu_c )/ ( 1.0_wp + kc1_a3 * nu_c )
        kc_rad = ( kc1_b1 + kc1_b2 * nu_c )/ ( 1.0_wp + kc1_b3 * nu_c )
        kc_sig = ( kc1_c1 + kc1_c2 * nu_c )/ ( 1.0_wp + kc1_c3 * nu_c )
        kc_bet = kc1_bet
        prey   = -0.125_wp

        D_mean = (1+mu)/lam
        diss = turb_edr * 1.0E4_wp
        Re = MAX(10000.0_wp * EXP ( 1.0_wp/6.0_wp*LOG(diss/100.0_wp) ), 2000.0_wp)
        k_turb = 1.0_wp + diss* EXP(prey*LOG(Re)) * (kc_bet + kc_alf * &
                EXP( -1.0_wp* ((((D_mean/2.0_wp)*1.e6_wp-kc_rad)/kc_sig)**2) ))

    end subroutine autoconversion_enhance

    subroutine accretion_enhance(k_turb, turb_edr, mu, lam)

        implicit none

        !*******************************************************************************
        ! Autoconversion enhancement from Seifert and Beheng (2001, Atmos. Res.)       *
        !*******************************************************************************

        ! Onishi kernel (of 29 July 2015)
        real(wp), parameter :: D_min = 8.0e-05_wp

        real(wp), intent(in) :: turb_edr, mu, lam
        real(wp), intent(inout) :: k_turb
        real(wp) :: diss,D_mean

        D_mean = (1+mu)/lam
        diss = turb_edr * 1.0E4_wp
        k_turb = 1.0_wp + 0.8e-3_wp * diss * (D_min/D_mean)**(2.0_wp)

    end subroutine accretion_enhance

    subroutine selfcollection_enhance(k_turb, turb_edr, mu, lam)

        implicit none

        !*******************************************************************************
        ! Autoconversion enhancement from Seifert and Beheng (2001, Atmos. Res.)       *
        !*******************************************************************************

        ! Onishi kernel (of 29 July 2015)
        real(wp), parameter :: D_min = 8.0e-05_wp

        real(wp), intent(in) :: turb_edr, mu, lam
        real(wp), intent(inout) :: k_turb
        real(wp) :: diss,D_mean

        D_mean = (1+mu)/lam
        diss = turb_edr * 1.0E4_wp
        k_turb = 1.0_wp + 0.8e-3_wp * diss * (D_min/D_mean)**(2.0_wp)

    end subroutine selfcollection_enhance

end module turb_enhance