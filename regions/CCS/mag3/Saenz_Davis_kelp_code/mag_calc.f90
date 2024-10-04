
!Notes on fortran building (gfortran):
!f2py -c mag_kinds_mod.f90 mag_parameters_mod.f90 mag_calc.f90 -m mag_calc_fortran --debug --f90flags="-m64 -g -c -O3 -fopenmp -fdec-math -fopenmp -march=native -funroll-loops" -lgomp
!
! MODULE mag_calc_mod
!
!   use iso_c_binding, only: c_int, c_float, c_double
!   use mag_kinds_mod
!   use mag_parameters_mod
!   implicit none
!
!   PUBLIC
!
! CONTAINS
!
!
!   SUBROUTINE c_mag_calc(lat,lon,sst,par,swh,cmag,no3,nflux,mask, &
!                       params, &
!                       Q,B,d_B,d_Q,Growth2,d_Be,d_Bm,d_Ns,harv,GRate,B_N)
!     real(c_float),   dimension(nx1,ny1), intent(in) :: lat,lon,sst,par,swh,cmag,no3,nflux
!     integer(c_int),  dimension(nx1,ny1), intent(in) :: mask
!     real(c_float),   dimension(npar), intent(in) :: params
!     real(c_float),   dimension(nx1,ny1), intent(inout) :: Q,B,d_B,d_Q,Growth2,d_Be,d_Bm,d_Ns,harv,GRate,B_N
!     CALL mag_calc(lat,lon,sst,par,swh,cmag,no3,nflux,mask, &
!                   params, &
!                   Q,B,d_B,d_Q,Growth2,d_Be,d_Bm,d_Ns,harv,GRate,B_N)
!   END SUBROUTINE c_mag_calc

  SUBROUTINE mag_calc(lat,lon,sst,par,chl,swh,mwp,cmag,no3,nflux, &
                      do_harvest,seed_now, mask, GD_count, B_0000, n_harv, t_harv, &
                      params, &
                      Q,B,d_B,d_Q,Growth2,d_Be,d_Bm,d_Ns,harv,GRate,B_N,Gave,Dave, &
                      min_lim,gQout,gTout,gEout,gHout, &
                      turn_of_day)

   use mag_kinds_mod
   use mag_parameters_mod
   use OMP_LIB, only : omp_set_num_threads

    implicit none

    ! input/output variables
    !----------------------------------------------------
    integer, parameter :: nx1 = 2160  ! apparently can't pull this from module :\
    integer, parameter :: ny1 = 4320
    integer, parameter :: npar1 = 37

    real(kind=r8), INTENT(IN) :: &
    lat(nx1,ny1), &
    lon(nx1,ny1), &
    sst(nx1,ny1), &
    par(nx1,ny1), &
    chl(nx1,ny1), &
    swh(nx1,ny1), &
    mwp(nx1,ny1), &
    cmag(nx1,ny1), &
    no3(nx1,ny1), &
    nflux(nx1,ny1)

    integer(kind=i4), INTENT(IN) :: &
    do_harvest(nx1,ny1)

    integer(kind=i4), INTENT(INOUT) :: &
    seed_now(nx1,ny1), &
    mask(nx1,ny1), &
    GD_count(nx1,ny1), &
    n_harv(nx1,ny1), &
    t_harv(nx1,ny1), &
    min_lim(nx1,ny1)

    !integer(kind=i4), INTENT(IN) :: &
    !ny,ny,npar

    real(kind=r8), INTENT(IN) :: &
      params(npar1)

    real(kind=r8), INTENT(INOUT) :: &
    Q(nx1,ny1), &
    B(nx1,ny1), &
    B_0000(nx1,ny1), &
    d_B(nx1,ny1), &
    d_Q(nx1,ny1), &
    Growth2(nx1,ny1), &
    d_Be(nx1,ny1), &
    d_Bm(nx1,ny1), &
    d_Ns(nx1,ny1), &
    harv(nx1,ny1), &
    GRate(nx1,ny1), &
    B_N(nx1,ny1), &
    Gave(nx1,ny1), &
    Dave(nx1,ny1), &
    gQout(nx1,ny1), &
    gTout(nx1,ny1), &
    gEout(nx1,ny1), &
    gHout(nx1,ny1)

    logical, INTENT(IN) :: &
    turn_of_day

    ! shared internal variables
    !----------------------------------------------------
    real(kind=r8)    ::   KsNO3,PARs,PARc

    ! private internal variables
    !----------------------------------------------------
    logical :: bad_forcing
    integer(kind=i4) :: i,j,thread
    real(kind=r8) :: &
    day_h,         &
    lambda,        &
    vQ,            &
    NO3_u,         &
    vNuTw_NO3,     &
    Uptake_NO3,     &
    UptakeN,     &
    gQ,     &
    gT,     &
    gE,     &
    gH,     &
    Growth,     &
    WP,     &
    M_Wave,     &
    M,     &
    f_harvest,    &
    par_watts,    &
    atten,    &
    chlmin,    &
    A,    &
    Bnew,    &
    N_new,    &
    dNs,    &
    harv1,  &
    Qmin,   &
    Qmax,   &
    mlim,   &
    dE,     &
    d_Q_growth, &
    N, B_calc, N_calc, G_calc, dBdt, dNdt, Upc1, Upc2, Q_new

    ! whole model domain calculations for current timestep
    ! ------------------------------------------------------------
    KsNO3 = params(mp_spp_Ks_NO3)
    PARs = params(mp_spp_PARs)
    PARc = params(mp_spp_PARc)
    Qmin = params(mp_spp_Qmin)
    Qmax = params(mp_spp_Qmax)

  CALL omp_set_num_threads(10)

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED) &
  !$OMP PRIVATE(i,j,thread,day_h,lambda,vQ,NO3_u,vNuTw_NO3,Uptake_NO3,UptakeN,gQ,gT,gE,gH,Growth,WP,M_Wave, &
  !$OMP M,f_harvest,par_watts,A,Bnew,N_new,dNs,harv1)
  !$OMP DO SCHEDULE(DYNAMIC,chunk_size)

    ! ----------------------------------------------------------------
   !   do i=1708,1710
    !do j=1001,1001
    !do j=425,425
    do j=1,ny1
      !do i=1001,1001
      !do i=375,375
      do i=1,nx1
        !print *,"Dude",i,j
      if (mask(i,j) > 0) then

        if (seed_now(i,j) > 0) then
            B(i,j) = params(mp_spp_seed)
            Q(i,j) = Qmin + no3(i,j)*(Qmax-Qmin)/35.0_r8
            t_harv(i,j) = 0
        endif

        ! Nutrient Uptake
        ! ----------------------------------------------------------
        lambda = lambda_NO3(cmag(i,j),mwp(i,j),params(mp_spp_CD),params(mp_spp_Vmax), &
                            params(mp_spp_Ks_NO3),no3(i,j))
        !print *,j,lambda
        ! Quota-limited uptake: maximum uptake when Q is minimum and
        ! approaches zero as Q increases towards maximum; Possible that Q
        ! is greater than Qmax. Set any negative values to zero.
        if (params(mp_N_uptake_type) == 0) then
			vQ = (Qmax-Q(i,j))/(Qmax-Qmin)
			vQ = max(c0,vQ)
			vQ = min(c1,vQ)

		! new vQ formulation that permits high uptake even with high Q, so there is not negative feedback
		! between high growth rate/high Q usage and reduced growth rate
		! The reasoning behind this uptake curve is that seaweed would not reduce uptake under high-growth,
		! high-nutrient conditions, just because stores are full - it's really a time stepping issue. High Q->
		! leads to nutrient limitation in one time-step.  If at moderate Q, seaweeds are also prevented from
		! uptaking short-pulsed nutrients, like closer to the scale of the timestep.  Basically the linear
		! function above does not allow response to changing nutrients at rates that seaweeds are
		! capable of.
		elseif (params(mp_N_uptake_type) == 1) then
			vQ = 1.0 - 1.0/(1.0+(max(0.0,Qmax-Q(i,j)))*50.0/(Qmax-Qmin))
			vQ = max(c0,vQ)
			vQ = min(c1,vQ)
		endif

        ! Below is what we call "Uptake Factor." It varies betwen 0
        ! and 1 and includes kinetically limited uptake and
        ! mass-transfer-limited uptake (oscillatory + uni-directional flow)
        NO3_u = no3(i,j)*1000.0_r8
        vNuTw_NO3 = NO3_u / (KsNO3 * ((NO3_u/KsNO3)  + 0.5_r8 * (lambda+sqrt(lambda**2 + c4 * (NO3_u/KsNO3)))))
        vNuTw_NO3 = max(c0,vNuTw_NO3)
        vNuTw_NO3 = min(c1,vNuTw_NO3)

        ! Uptake Rate [mg N/g(dry)/d]
        ! Nutrient Uptake Rate = Max Uptake * v[Ci,u,Tw] * vQ
        ! converted from umol N/m2/d -> mg N/g(dry)/d by 14.0067 / 1e3
        Uptake_NO3 = params(mp_spp_Vmax) * vNuTw_NO3 * vQ ! [umol/m2/d]
        Uptake_NO3 = Uptake_NO3 * mw_n / 1.e3_r8
        UptakeN = Uptake_NO3
        UptakeN = UptakeN / params(mp_spp_dry_sa)
        !print *,lambda,vQ,vNuTw_NO3,UptakeN,Q(i,j)

        ! Growth
        ! ----------------------------------------------------------
        ! Growth, nitrogen movement from Ns to Nf = umax*gQ*gT*gE*gH; [per day]
        ! Output:
        !   Growth, [h-1]
        !   gQ, quota-limited growth
        !       from Wheeler and North 1980 Fig. 2
        !   gT, temperature-limited growth
        !       piecewise approach taken from Broch and Slagstad 2012 (for sugar
        !       kelp) and optimized for Macrocystis pyrifera
        !   gE, light-limited growth
        !       from Dean and Jacobsen 1984
        !   gH, carrying capacity-limited growth

        ! nutrient (quota-based) limitation
        !gQ = (Q(i,j) - params(mp_spp_Qmin)) / Q(i,j) ! Droop equation
        gQ = (Q(i,j) - Qmin) / Q(i,j) * Qmax/(Qmax-Qmin) ! Droop scaled from 0-1
        !gQ = (Q(i,j) - Qmin) / (Qmax-Qmin) ! Freider et al.

		if (params(mp_Q_lim_type) == 0) then
			gQ = (Q(i,j) - Qmin) / Q(i,j) * Qmax/(Qmax-Qmin) ! Droop scaled from 0-1
		elseif (params(mp_Q_lim_type) == 1) then
			gQ = (Q(i,j) - Qmin) / (Qmax-Qmin) ! Freider et al.
		elseif (params(mp_Q_lim_type) == 2) then
			gQ = 1.0 - exp(-11.0 * (Q(i,j) - Qmin) / (Qmax-Qmin)) ! more permissive than Droop scaled
		endif

        gQ = max(c0,gQ)
        gQ = min(c1,gQ)

        ! temperature limitation
        gT = temp_lim(sst(i,j),params(mp_spp_Topt1),params(mp_spp_K1), &
                      params(mp_spp_Topt2),params(mp_spp_K2))

        ! light limitation
        par_watts = par(i,j) ! *2.515376387217542_r8 <-- no longer internally converting par - should be in W/m2
        ! attentuation according to MARBL
        chlmin = max(0.02_r8,chl(i,j))
        chlmin = min(30.0_r8,chlmin)
        if (chlmin < 0.13224_r8) then
            atten = -0.000919_r8*(chlmin**0.3536_r8) ! 1/cm
        else
            atten = -0.001131_r8*(chlmin**0.4562_r8) ! 1/cm
        endif
        par_watts = par_watts * exp(atten*params(mp_farm_depth)*100.0_r8)
        if (par_watts < PARc) then
          gE = c0
        elseif (par_watts > PARs) then
          gE = c1
        else
          gE = (par_watts-PARc)/(PARs-PARc)*exp(-(par_watts-PARc)/(PARs-PARc)+c1)
        endif

        ! consider daylength if timestep is > 1/2 a day
        if (params(mp_dt_mag) > 0.51) then
            day_h = daylength(lat(i,j),lon(i,j),c0,params(mp_dte)) ! daylength in h
            gE = min(day_h,gE*day_h)
        endif

        ! Carrying capacity
        ! ----------------------------------------------------------
          ! gH -> density-limited growth (ranges from the max growth rate to the death rate, excluding wave mortality)
          ! This expression follows Xiao et al (2019 and ignores wave mortality when
        ! thinking about the death rate

        Bnew = B(i,j) !*params(mp_spp_line_sep) ! converting from g/m2 to g/m

        !A = params(mp_spp_kcap_rate)/(params(mp_spp_kcap)**(-1.44_r8))
        !gH = A*Bnew**(-1.44_r8)
        A = params(mp_spp_kcap_rate)/(params(mp_spp_kcap)**(-0.75_r8))
        gH = A*Bnew**(-0.75_r8)

		if (params(mp_growth_lim_type) == 0) then
			! interacting limitation terms growth model
			Growth =  min(gH,params(mp_spp_Gmax_cap)) * gT * gE * gQ
		else
			! independent light and nutrient limitation growth model
			Growth =  min(gH,params(mp_spp_Gmax_cap)) * gT * min(gE,gQ)
		endif

        !print *,day_h,Growth,gQ,gT ,gE , gH

        gQout(i,j) = gQ
        gTout(i,j) = gT
        gEout(i,j) = gE
        gHout(i,j) = gH

        !           1, 2, 3, 4
        mlim = min(gQ,gT,gE,gH)
        if (gT == mlim) then  ! depending on floating point type, this may not work in fortran
            min_lim(i,j) = 2
        elseif (gE == mlim) then
            min_lim(i,j) = 3
        elseif (gQ == mlim) then
            min_lim(i,j) = 1
        elseif (gH == mlim) then
            min_lim(i,j) = 4
        endif


        ! Mortality
        ! ----------------------------------------------------------
        ! d_wave = frond loss due to waves; dependent on Hs, significant
        ! wave height [m]; Rodrigues et al. 2013 demonstrates linear relationship
        ! between Hs and frond loss rate in Macrocystis [d-1] (continuous)
        ! Duarte and Ferreira (1993) find a linear relationship between wave power and mortality in red seaweed.
        if (params(mp_breakage_type) == breakage_Duarte_Ferreira) then
          !WP = rho.*g.^2/(64*pi)*swh.^2.*Tw
          WP = 1025._r8*9.8_r8**2 / (64._r8*pi)*swh(i,j)**2 * mwp(i,j) /1.e3_r8 ! [kW]
          ! [Duarte and Ferreira (1993), in daily percentage]
          M_wave = 2.3_r8*1e-4_r8*WP + 2.2_r8*1e-3_r8 * params(mp_wave_mort_factor)
        else !if (params(mp_breakage_type) == breakage_Rodrigues) then
          M_wave  = params(mp_spp_death) * swh(i,j) * params(mp_wave_mort_factor)
        endif

        ! M = M_wave + general Mortality rate; [d-1]
        M = params(mp_spp_death) + M_wave

        ! limit growth to nflux
        ! ----------------------------------------------------------
        ! Comparing the nitrogen taken up to the amount of nitrogen fluxed in per m2
        dNs = UptakeN * params(mp_dt_mag) * B(i,j)   ! [mg-N/ m2]
        if (params(mp_N_flux_limit) > c0) then
            N_new = nflux(i,j) * mw_n * 864.0  ! /100 * 86400 * 14.006 [microM/cm2/s]->[mmol-N/m2/day] -> [mg-N/m2/day]

        if ((dNs > N_new) .and. (N_new > c0)) then
              dNs = N_new
            elseif (N_new < c0) then
              dNs = c0
            endif
        endif
        UptakeN = dNs/params(mp_dt_mag)/B(i,j) ! [mg-N/g-dry/day] - Redefining Uptake -- is this needed?
        ! Uptake2 = dNs # [mg-N] --> saving for post-processing

        ! find exudation
        dE = params(mp_spp_E) * (Q(i,j)-Qmin) * params(mp_dt_mag)

        ! growth via solver
        ! ----------------------------------------------------------

        if (vQ <= c0) then
            Upc1 = c0
        else
            Upc1 = UptakeN/vQ
        endif
        Upc2 = 50.0_r8/(Qmax-Qmin)
        N = Q(i,j)*B(i,j)

        call growth_forward_euler(N, B(i,j), params(mp_spp_Gmax_cap), A,  &
                    params(mp_spp_kcap), Qmin, Qmax, params(mp_spp_E), M, &
                    gE, gT, params(mp_growth_lim_type), Upc1, Upc2, params(mp_dt_mag), &
                    B_calc, N_calc, G_calc, dBdt, dNdt)


        ! Output terms
        ! ----------------------------------------------------------
        d_B(i,j) = dBdt
        Q_new = (N+dNdt) / (B(i,j)+dBdt)
        d_Q(i,j) = Q_new - Q(i,j)
        Q(i,j) = Q_new

        GRate(i,j) = G_calc
        Growth2(i,j) = Growth2(i,j) + G_calc * B(i,j) * params(mp_dt_mag) !--- this will be inaccurate for non-forward euler growth solvers?

        d_Be(i,j) = d_Be(i,j) + (N_calc - Qmin*B_calc)*params(mp_spp_E)*params(mp_dt_mag) !--- this will be inaccurate for non-forward euler growth solvers?
        d_Bm(i,j) = d_Bm(i,j) + B_calc * M * params(mp_dt_mag) !--- this will be inaccurate for non-forward euler growth solvers?
        d_Ns(i,j) = d_Ns(i,j) + UptakeN * params(mp_dt_mag) * B(i,j) ! mg N per m2 --- this will be inaccurate for non-forward euler growth solvers?
        B_N(i,j) = 1.0e3_r8/Q(i,j)


        ! legacy code ------------------------------------------------------------------
!         d_B(i,j) = Growth * B(i,j) * params(mp_dt_mag) - M * B(i,j) * params(mp_dt_mag)
!         d_Q_growth = Q(i,j) * (c1/(c1 + Growth * params(mp_dt_mag)) - c1)
!         d_Q(i,j) = UptakeN * params(mp_dt_mag) + d_Q_growth - dE
!
!         !print *,j,Q(i,j),d_Q(i,j),B(i,j),d_B(i,j),Growth,gQ,gT,gE,gH,UptakeN
!         !print *,j,par(i,j),day_h,gE,d_Q(i,j),UptakeN,lambda,NO3_u,vNuTw_NO3
!
!         GRate(i,j) = Growth
!         Growth2(i,j) = Growth2(i,j) + Growth * B(i,j) * params(mp_dt_mag)
!         d_Be(i,j) = d_Be(i,j) + dE * B(i,j)
!         d_Bm(i,j) = d_Bm(i,j) + B(i,j) * M * params(mp_dt_mag)
!         d_Ns(i,j) = d_Ns(i,j) + UptakeN * params(mp_dt_mag) * B(i,j) ! mg N per m2
!         B_N(i,j) = 1.0e3_r8/Q(i,j)

        ! end legacy code --------------------------------------------------------------

        ! Update State Variables
        ! ----------------------------------------------------------
        B(i,j) = B(i,j) + d_B(i,j)
        Q(i,j) = Q(i,j) + d_Q(i,j)


        ! Harvest
        ! ----------------------------------------------------------

        ! Older fixed harvest code
        !if (do_harvest(i,j) == 1) then
        !    if (params(mp_harvest_type) == 0) then
        !      ! harvest to seed weight
        !      harv1 = max(c0,(B(i,j)-params(mp_spp_seed)/params(mp_spp_line_sep))/B(i,j))
        !    else ! params[mp_harvest_type] == 1:
        !      ! fractionally harvest, but not if below seed weight
        !      if (B(i,j) < params(mp_spp_seed)) then
        !          harv1 = c0
        !      else
        !          harv1 = params(mp_harvest_f)
        !      endif
        !    endif
        !
        !    harv1 = harv1 * B(i,j)  ! biomass harvested
        !    B(i,j) = B(i,j) - harv1
        !    harv(i,j) = harv(i,j) + harv1
        !endif

        ! increment growth/death running averages
        if (seed_now(i,j) > 0) then
           Gave(i,j) = Growth ! growth [1/timestep]
           Dave(i,j) = M    ! death [1/timestep]
        endif
        Gave(i,j) = Gave(i,j) + (Growth-Gave(i,j)) / (params(mp_harvest_avg_period) / params(mp_dt_mag) )
        Dave(i,j) = Dave(i,j) + (M-Dave(i,j)) / (params(mp_harvest_avg_period) / params(mp_dt_mag) )

        seed_now(i,j) = 0

        if (turn_of_day) then

            ! counter for if death rate exceeds growth rate
            if (M > Growth) then
                GD_count(i,j) = GD_count(i,j) + 1
            else
                GD_count(i,j) = 0
            endif

            !print *, do_harvest(i,j)
            !print *, n_harv(i,j),t_harv(i,j),params(mp_harvest_nmax)
            !print *, B(i,j),params(mp_harvest_kg)*params(mp_spp_line_sep)*1.0e3_r8
            !print *, GD_count(i,j),params(mp_harvest_avg_period)

            ! check for, and perform harvest


            harv1 = c0
            if (params(mp_harvest_schedule) == 0) then
                ! fixed harvest
                if (do_harvest(i,j) == 1) then
                    if (params(mp_harvest_type) == 0) then
                        ! harvest to seed weight
                        if (B(i,j) > c0) then
                            harv1 = max(c0,(B(i,j)-params(mp_spp_seed)/params(mp_spp_line_sep))/B(i,j))
                        endif
                        n_harv(i,j) = n_harv(i,j) + 1
                        t_harv(i,j) = t_harv(i,j) + 1

                    else ! params[mp_harvest_type] == 1:
                        ! fractionally harvest, but not if below mp_harvest_kg
                        !if (B(i,j) >= params(mp_harvest_kg)*params(mp_spp_line_sep)*1.0e3_r8) then
                        !    harv1 = c0
                        !else
                            harv1 = params(mp_harvest_f)
                            n_harv(i,j) = n_harv(i,j) + 1
                            t_harv(i,j) = t_harv(i,j) + 1
                        !endif

                    endif
                endif
            else
                !  conditional harvest
                if (do_harvest(i,j) == 2) then
                    if (params(mp_harvest_type) == 1) then
                        ! we are done for the season, if not already harvested
                        harv1 = params(mp_harvest_f)
                    else
                        harv1 = 0.99
                    endif
                    mask(i,j) = 0
                    n_harv(i,j) = n_harv(i,j) + 1
                    t_harv(i,j) = t_harv(i,j) + 1
                elseif (do_harvest(i,j) == 1) then

                    if (params(mp_harvest_schedule) == 1) then

                        ! within harvest span - check for declining growth
                        !if Gave[i,j]/Dave[i,j] < 1.0:  ! test of declining biomass over time
                        if (GD_count(i,j) > params(mp_harvest_avg_period)) then ! test of negative growth over time
                            if (params(mp_harvest_type) == 1) then
                                ! we are done for the season, if not already harvested
                                harv1 = params(mp_harvest_f)
                            else
                                harv1 = 0.99
                                !harv1 = max(c0,(B(i,j)-params(mp_spp_seed)/params(mp_spp_line_sep))/B(i,j))
                            endif
                            mask(i,j) = 0
                            n_harv(i,j) = n_harv(i,j) + 1
                            t_harv(i,j) = t_harv(i,j) + 1
                        endif

                        if (harv1 < 1.0e-4_r8) then
                            if (t_harv(i,j) < (params(mp_harvest_nmax)-1)) then   !!!!!!!!!!!!!!!!! add -1
                                ! if not already harvested, check to see if conditions are OK for incremental harvest
                                if (B(i,j) >= params(mp_harvest_kg)*params(mp_spp_line_sep)*1.0e3_r8) then
                                    harv1 = params(mp_harvest_f)
                                    n_harv(i,j) = n_harv(i,j) + 1
                                    t_harv(i,j) = t_harv(i,j) + 1
                                endif
                            endif
                        endif

                    else ! period harvest -- params[mp_harvest_schedule] == 2:

                        ! do harvest according to "type"
                        if (params(mp_harvest_type) == 0) then
                            ! harvest to seed weight
                            harv1 = max(c0,(B(i,j)-params(mp_spp_seed)/params(mp_spp_line_sep))/B(i,j))
                            n_harv(i,j) = n_harv(i,j) + 1
                            t_harv(i,j) = t_harv(i,j) + 1
                        else ! params[mp_harvest_type] == 1:
                            ! fractionally harvest, but not if below mp_harvest_kg
                            if (B(i,j) > params(mp_harvest_kg)*params(mp_spp_line_sep)*1.0e3_r8) then
                                harv1 = params(mp_harvest_f)
                                n_harv(i,j) = n_harv(i,j) + 1
                                t_harv(i,j) = t_harv(i,j) + 1
                            endif
                        endif

                    endif

                endif
            endif
        endif
        !if (harv1 > c0) then  ! reset mean biomass growth/death records so they don't swing around?
        !   Gave(i,j) = Growth ! growth [1/timestep]
        !   Dave(i,j) = M    ! death [1/timestep]
        !endif
        !print *, harv1, n_harv(i,j), t_harv(i,j)

        harv1 = harv1 * B(i,j)  ! biomass harvested
        B(i,j) = B(i,j) - harv1
        harv(i,j) = harv(i,j) + harv1

        ! ## sanity checks, at end of calc to prevent output of bad data
        if (mask(i,j) > 0) then
            bad_forcing = .false.
            ! FORTRAN test for NaN - not equal to itself!
            if (sst(i,j) /= sst(i,j)) then
              bad_forcing = .true.
            endif
            if (par(i,j) /= par(i,j)) then
              bad_forcing = .true.
            endif
            if (no3(i,j) /= no3(i,j)) then
              bad_forcing = .true.
            endif
            if (swh(i,j) /= swh(i,j)) then
              bad_forcing = .true.
            endif
            if (cmag(i,j) /= cmag(i,j)) then
              bad_forcing = .true.
            endif
            if (nflux(i,j) /= nflux(i,j)) then
              bad_forcing = .true.
            endif
            if (B(i,j) /= B(i,j)) then
              bad_forcing = .true.
            endif
            if (Q(i,j) /= Q(i,j)) then
              bad_forcing = .true.
            endif
            ! not sure what to set things too - can't use nan in fortran?  maybe using ieee stuff
            if (bad_forcing) then
              mask(i,j) = 0
              B(i,j) = -999.9
              Q(i,j) = -999.9
              harv(i,j) = -999.9
              d_B(i,j) = -999.9
              d_Q(i,j) = -999.9
              GRate(i,j) = -999.9
              Growth2(i,j) = -999.9
              d_Be(i,j) = -999.9
              d_Bm(i,j) = -999.9
              d_Ns(i,j) = -999.9
              B_N(i,j) = -999.9
              n_harv(i,j) = 0
              Gave(i,j) = -999.9
              Dave(i,j) = -999.9
              min_lim(i,j) = 0
              gQout(i,j) = -999.9
              gTout(i,j) = -999.9
              gEout(i,j) = -999.9
              gHout(i,j) = -999.9
            endif

        endif

      endif
      enddo
    enddo

  !$OMP END DO
  !$OMP END PARALLEL


  end SUBROUTINE mag_calc


  PURE SUBROUTINE growth_forward_euler(N, B, Gmax, A, Kslope, Qmin, Qmax, E, D, &
                    Llim, Tlim, growth_lim_type, Upc1, Upc2, dt, &
                    B_calc, N_calc, G_calc, dBdt, dNdt)
!
!     Generates the forward euler method (delta) biomass at time t+1, incorporating
!     growth rate dependencies on Q (modified Droop version) and crowding, of the
!     form growth_rate_max = A*Biomass**(Kslope).
!
!     Parameters
!     ----------
!     N : step initial Nutrient mass (mg / m-2; tracer, what we are solving for)
!     B : step initial biomass (g / m-2; tracer, what we are solving for)
!     Gmax : maximum growth rate (1/day)
!     A : crowding constant, dependent on species parameters
!     Kslope : crowding constant, dependent on species parameters
!     Qmin : minimum species Nutrient storage (mg Nut/g biomass)
!     Qmax : maximum Nutrient storage (mg Nut/g biomass)
!     E : exudation rate (1/day)
!     D : death rate (1/day)
!     Llim : light limitation, considered constant [0-1]
!     Tlim : temperature limitation, considered constant  [0-1]
!     growth_lim_type : growth rate calc switch
!     Upc1: Uptake constant 1 (environment+mechalis menton limited)
!     Upc2: Uptake contant 2 (50./(Qmax-Qmin), part of the modified Droop Q-based uptake limtation formulation
!     dt : time step (days)
!
!     Returns
!     -------
!     B_calc : best estimate of biomass used for step calc
!     N_calc : best estimate of nutrient mass used for step calc
!     G_mid = best estimate of biological growth rate used for step calc
!     dBdt = biomass change at t+dt
!     dNdt = Q change at t+dt
!
    use mag_kinds_mod, only : r8

    real(kind=r8), INTENT(IN) :: &
    N, B, Gmax, A, Kslope, Qmin, Qmax, E, D, &
    Llim, Tlim, growth_lim_type, Upc1, Upc2, dt

    real(kind=r8), INTENT(OUT) :: &
    B_calc, N_calc, G_calc, dBdt, dNdt

    ! internal variables
    !----------------------------------------------------
    real(kind=r8) :: &
    Gmax_crowding, Q, Qlim

    Gmax_crowding = A * B**Kslope
    Gmax_crowding = max(0.0_r8,Gmax_crowding)
    Gmax_crowding = min(1.0_r8,Gmax_crowding)

    ! Calculate G and Qlim based on current values
    Q = N/B
    Qlim = (Q - Qmin) / Q * Qmax/(Qmax-Qmin)
    if (growth_lim_type == 0) then
        G_calc = min(Gmax,Gmax_crowding) * Tlim * Llim * Qlim
    else
        G_calc = min(Gmax,Gmax_crowding) * Tlim * min(Llim,Qlim)
    endif

    ! debug terms
    !UptakeN = Upc1*(1.0 - 1.0/(1.0 + (Qmax - N/B) * Upc2))
    !dQdt = Q * (1.0 / (1.0 + G * dt) - 1.0) + (UptakeN - (Q - Qmin) * E) * dt
    !alt_dQdt = (N+dNdt)/(B+dBdt)-Q #dQdt

    dBdt = B * (G_calc - D) * dt
    dNdt = (B*Upc1*(1.0_r8 - 1.0_r8/(1.0_r8 + (Qmax - N/B) * Upc2)) - N*D - (N - Qmin*B)*E)*dt

    ! return B,Q,G in addition to dBdt and dQdt for calculation of some output terms
    B_calc = B
    N_calc = N

end SUBROUTINE growth_forward_euler

!end MODULE mag_calc_mod




