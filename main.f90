program gwgen

    use csv_file,      only : csv_write
    use parametersmod, only : sp, tfreeze
    use weathergenmod, only : metvars_in, metvars_out, weathergen, init_weathergen, &
                              rmsmooth
    use randomdistmod, only : ran_seed
    use geohashmod,    only: geohash

    implicit none

    character(200) :: infile
    character(200) :: output
    character(200) :: dailyfile

    integer, parameter :: n = 1
    integer, parameter :: n_curr = n + 1
    integer, parameter :: n_tot = n * 2 + 1
    character(12), parameter :: NOSTATION = '__NOSTATION'

    integer :: i, d, rmin, rd, lmin, ld, curr_start, curr_end, end_counter = n
    integer :: year(n_tot) = 0     ! month of previous, current, next and the month after the next month
    integer :: month(n_tot) = 0    ! year of previous, current, next and the month after the next month

    character(12) :: stationid(n_tot) = NOSTATION  ! previous, current and next station

    real(sp) :: lon(n_tot) = 0  ! longitude of the station of previous, current, next and the month after the next month
    real(sp) :: lat(n_tot) = 0  ! latitude of the station of previous, current, next and the month after the next month

    real(sp) :: mprec(n_tot) = 0    ! precipitation amount of current month
    integer  :: mwet(n_tot) = 0      ! number of wet days of current month
    real(sp) :: mtmin(n_tot) = 0  ! min temperture of previous, current, next and the month after the next month
    real(sp) :: mtmax(n_tot) = 0  ! min temperture of previous, current, next and the month after the next month
    real(sp) :: mcloud(n_tot) = 0 ! cloudiness of previous, current, next and the month after the next month
    real(sp) :: mwind(n_tot) = 0 ! wind speed of previous, current, next and the month after the next month

    real(sp), target :: tmin_sm(n_tot * 31) = 0   ! smoothed daily values of min temperature
    real(sp), target :: tmax_sm(n_tot * 31) = 0   ! smoothed daily values of max temperature
    real(sp), target :: cloud_sm(n_tot * 31) = 0  ! smoothed daily values of cloudiness
    real(sp), target :: wind_sm(n_tot * 31) = 0   ! smoothed daily values of wind speed

    real(sp) :: bcond_tmin(2) = 0      ! boundary conditions of min temp for smoothing
    real(sp) :: bcond_tmax(2) = 0      ! boundary conditions of max temp for smoothing
    real(sp) :: bcond_cloud(2) = 0     ! boundary conditions of cloud for smoothing
    real(sp) :: bcond_wind(2) = 0      ! boundary conditions of wind speed for smoothing

    integer :: i_consecutives(n_tot - 1)

    ! pointers to the tmin, tmax, cloud and wind speed values of the current month
    real(sp), pointer :: mtmin_curr(:), mtmax_curr(:), mcloud_curr(:), mwind_curr(:)

    integer  :: mwetd_sim, wet_day
    real(sp) :: mprec_sim

    integer :: pdaydiff

    real(sp) :: precdiff

    real(sp) :: tmindiff, tmin_acc

    type(metvars_in)  :: met_in
    type(metvars_out) :: met_out
    type(metvars_out) :: met_out_save

    integer :: ndm(n_tot) = 0

    real(sp) :: prec_t

    integer :: i_count, i_linecount = n_curr + 1

    type(metvars_out), dimension(31) :: month_met

    ! -----------------------------------------------------------------------------
    ! ------------------- Defaults for the namelist parameters --------------------
    ! -----------------------------------------------------------------------------
    integer :: seed = -30000 ! seed for the random number generator. Has no effect if use_geohash is .true.
    ! logical to determines whether the geohashmod module shall be used for
    ! defining the random seed. If .true., the first three columns must be
    ! stationid, longitude, latitude. Otherwise longitude and latitude should be
    ! skipped
    logical :: use_geohash = .true.
    logical :: ldebug = .false.  ! switch to output additional informations
    logical :: lreset = .true. ! switch to reset the residuals every month

    ! -----------------------------------------------------------------------------
    ! ------------------- END. Defaults for the namelist parameters ---------------
    ! -----------------------------------------------------------------------------

    namelist / main_ctl / seed, use_geohash, ldebug, lreset

    open(101,file='weathergen.nml',status='old')
    read(101, main_ctl)
    call init_weathergen(101)
    close(101)

    !initialize random state


    call getarg(1,infile)
    call getarg(2,output)

    open(10,file=infile,status='old')
    ! consume header
    read(10,*)

    ! read in the first n months and calculate the wet days


    if (use_geohash) then
        do i=n_curr,n_tot-1

            read(10,*)stationid(i),lon(i),lat(i),year(i), month(i),mtmin(i),mtmax(i), &
                   mcloud(i),mwind(i),mprec(i),mwet(i)
            ndm(i) = ndaymonth(year(i), month(i))

        end do
    else
        do i=n_curr,n_tot-1

            read(10,*)stationid(i),year(i), month(i),mtmin(i),mtmax(i), &
                       mcloud(i),mwind(i),mprec(i),mwet(i)
            ndm(i) = ndaymonth(year(i), month(i))

        end do
        call ran_seed(seed,met_in%rndst)
    endif

    ! open the output file
    dailyfile = trim(output)

    open(30,file=dailyfile,status='unknown')
    write(30,'(a)', advance='no') "id,year,month,day,tmin,tmax,mean_cloud,wind,prcp,wet_day"
    if (ldebug) then
        write(30,'(a)', advance='no') ",tmin_bias,tmin_mn,tmin_sd,wind_bias,wind_intercept_bias,"
        write(30,'(a)', advance='no') "wind_mn,wind_sd,resid_tmin,resid_tmax,resid_cloud,resid_wind,"
        write(30, '(a)') "unorm_tmin,unorm_tmax,unorm_cloud,unorm_wind"
    else
        write(30, '(a)') ''
    endif

    do  !read the input file until the end

        ! nullify pointers
        nullify(mtmin_curr)
        nullify(mtmax_curr)
        nullify(mcloud_curr)
        nullify(mwind_curr)

        i_linecount = i_linecount + 1

        !read in one month of summary weather station data from a text file
        if (use_geohash) then
            read(10,*,end=99)stationid(n_tot),lon(n_tot),lat(n_tot),year(n_tot),month(n_tot), &
                           mtmin(n_tot),mtmax(n_tot),mcloud(n_tot),mwind(n_tot),mprec(n_tot), &
                           mwet(n_tot)
            ndm(n_tot) = ndaymonth(year(n_tot),month(n_tot))
            if (stationid(n_curr) /= stationid(n)) then
                call ran_seed(geohash(lon(n_curr), lat(n)),met_in%rndst)
            end if
        else
            read(10,*,end=99)stationid(n_tot),lon(n_tot),lat(n_tot),year(n_tot),month(n_tot), &
                           mtmin(n_tot),mtmax(n_tot),mcloud(n_tot),mwind(n_tot),mprec(n_tot), &
                           mwet(n_tot)
            ndm(n_tot) = ndaymonth(year(n_tot),month(n_tot))
        endif

        goto 110

        ! ------ this part is skipped if we are not already at the end of the file -------------
99 continue
            end_counter = end_counter - 1
            stationid(n_tot) = NOSTATION
        ! ---------------------------------------------------------------------------------------
110 continue

        ! ------ check for invalid values
        if (mprec(n_tot) < 0.0 .and. mprec(n_tot) > -0.1) then
            mprec(n_tot) = 0.0
        elseif (mprec(n_tot) < -0.1) then
            write(0,*) "Invalid precipitation value", mprec(n_tot), "mm/d at line ", i_linecount
            stop 1
        endif
        if (mtmin(n_tot) + tfreeze < 0.0) then
            write(0,*) "Invalid minimum temperature value", mtmin(n_tot), "degC at line ", i_linecount
            stop 1
        elseif (mtmax(n_tot) + tfreeze < 0.0) then
            write(0,*) "Invalid maximum temperature value", mtmax(n_tot), "degC at line ", i_linecount
            stop 1
        elseif (mcloud(n_tot) < 0.0 .or. mcloud(n_tot) > 1.0) then
            write(0,*) "Invalid cloud fraction ", mcloud(n_tot), "at line ", i_linecount
            stop 1
        elseif (mwind(n_tot) < 0.0) then
            write(0,*) "Invalid wind speed ", mwind(n_tot), "m/s at line ", i_linecount
            stop 1
        endif

        i_consecutives(:) = are_consecutive_months(stationid,year,month)

        if (any(i_consecutives(n:1:-1) == 0)) then
            lmin = n_curr - minloc(i_consecutives(n:1:-1), 1) + 1
            ld = max(1, SUM(ndm(:lmin - 1)))
            bcond_tmin(1) = mtmin(lmin - 1)
            bcond_tmax(1) = mtmax(lmin - 1)
            bcond_cloud(1) = mcloud(lmin - 1)
            bcond_wind(1) = mwind(lmin - 1)
        else
            lmin = 1
            ld = 1
        end if

        if (any(i_consecutives(n+1:) == 0)) then
            rmin = n + minloc(i_consecutives(n+1:), 1)
            bcond_tmin(2) = mtmin(rmin)
            bcond_tmax(2) = mtmax(rmin)
            bcond_cloud(2) = mcloud(rmin)
            bcond_wind(2) = mwind(rmin)
        else
            rmin = n_tot
            bcond_tmin(2) = mtmin(n_tot)
            bcond_tmax(2) = mtmax(n_tot)
            bcond_cloud(2) = mcloud(n_tot)
            bcond_wind(2) = mwind(n_tot)
        end if
        rd = SUM(ndm(:rmin))

        call rmsmooth(mtmin(lmin:rmin), ndm(lmin:rmin), bcond_tmin, tmin_sm(ld:rd))
        call rmsmooth(mtmax(lmin:rmin), ndm(lmin:rmin), bcond_tmax, tmax_sm(ld:rd))
        call rmsmooth(mcloud(lmin:rmin), ndm(lmin:rmin), bcond_cloud, cloud_sm(ld:rd))
        call rmsmooth(mwind(lmin:rmin), ndm(lmin:rmin), bcond_wind, wind_sm(ld:rd))

        curr_start = sum(ndm(:n))
        curr_end = curr_start + ndm(n+1) - 1

        mtmin_curr => tmin_sm(curr_start:curr_end)
        mtmax_curr => tmax_sm(curr_start:curr_end)
        mcloud_curr => cloud_sm(curr_start:curr_end)
        mwind_curr => wind_sm(curr_start:curr_end)

        met_in%prec = mprec(n_curr)
        met_in%wetd = real(mwet(n_curr))
        met_in%wetf = real(mwet(n_curr)) / real(ndm(n_curr))

        !initialize weather residuals and other variables that carry over from one day to the next
        !these and the random state below should be reset once per station

        if (i_consecutives(n_curr - 1) == 0) then
            met_out%pday(1) = .false.
            met_out%pday(2) = .false.
            met_out%resid = 0.
            i_count = 1
        else if (lreset) then
            met_out%resid = 0.
            i_count = 1
        else
            i_count = 2
        end if
        prec_t = max(0.5,0.05 * mprec(n_curr))  !set quality threshold for preciptation amount

        ! save the current state of the residuals and pday
        met_out_save = met_out

        do

            mwetd_sim = 0
            mprec_sim = 0.
            tmin_acc = 0.0

            !start day loop

            do d = 1,ndm(n_curr)

                met_in%tmin  = mtmin_curr(d)
                met_in%tmax  = mtmax_curr(d)
                met_in%cldf  = real(mcloud_curr(d))
                met_in%wind  = real(mwind_curr(d))
                met_in%pday  = met_out_save%pday
                met_in%resid = met_out_save%resid

                call weathergen(met_in,met_out)

                met_in%rndst = met_out%rndst

                month_met(d) = met_out

                if (met_out%prec > 0.) then
                    mwetd_sim = mwetd_sim + 1
                    mprec_sim = mprec_sim + met_out%prec
                end if
                tmin_acc = tmin_acc + met_out%tmin

            end do

            !end of month
            ! NOTE: THIS HAS TO BE ENABLED IF NECESSARY
            tmindiff = 1.0  !abs(mtmin(n_curr) - tmin_acc / ndm(n_curr))

            ! Reset met_out_save after initialization
            if (i_count == 1) then
                if (i_consecutives(n_curr - 1) == 0) then
                    met_out_save = met_out
                else
                    met_out_save%resid = met_out%resid
                end if
            endif

            if (mprec(n_curr) <= 0.1 .and. tmindiff < 2.5) then

                pdaydiff = 0
                precdiff = 0.

                exit

            else if (i_count >= 2) then  !enforce at least two times over the month to get initial values ok

                pdaydiff = abs(mwet(n_curr) - mwetd_sim)
                precdiff = abs(mprec(n_curr) - mprec_sim)

                ! restrict simulated total monthly precip to +/-5% or 0.5 mm of observed value
                if (pdaydiff <= 1 .and. precdiff <= prec_t .and. tmindiff < 2.5) then
                    exit

                else if (i_count > 10000000) then
                    write (*,*) "No good solution found after 10000000 iterations at line", i_linecount
                    stop 1

                end if

            end if

            i_count = i_count + 1

        end do

        !write out final results for this station/year/month combo and

        do d = 1,ndm(n_curr)
            if (month_met(d)%prec > 0) then
                wet_day = 1
            else
                wet_day = 0
            end if
            call csv_write(30,stationid(n_curr),advance=.false.)
            call csv_write(30,year(n_curr),advance=.false.)
            call csv_write(30,month(n_curr),advance=.false.)
            call csv_write(30,d,advance=.false.)
            call csv_write(30,month_met(d)%tmin,advance=.false.)
            call csv_write(30,month_met(d)%tmax,advance=.false.)
            call csv_write(30,month_met(d)%cldf,advance=.false.)
            call csv_write(30,month_met(d)%wind,advance=.false.)
            call csv_write(30,month_met(d)%prec,advance=.false.)
            call csv_write(30,wet_day,advance=.not. ldebug)
            if (ldebug) then
                call csv_write(30,month_met(d)%tmin_bias,advance=.false.)
                call csv_write(30,month_met(d)%tmin_mn,advance=.false.)
                call csv_write(30,month_met(d)%tmin_sd,advance=.false.)
                call csv_write(30,month_met(d)%wind_bias,advance=.false.)
                call csv_write(30,month_met(d)%wind_intercept_bias,advance=.false.)
                call csv_write(30,month_met(d)%wind_mn,advance=.false.)
                call csv_write(30,month_met(d)%wind_sd,advance=.false.)
                call csv_write(30,month_met(d)%resid(1),advance=.false.)
                call csv_write(30,month_met(d)%resid(2),advance=.false.)
                call csv_write(30,month_met(d)%resid(3),advance=.false.)
                call csv_write(30,month_met(d)%resid(4),advance=.false.)
                call csv_write(30,month_met(d)%unorm(1),advance=.false.)
                call csv_write(30,month_met(d)%unorm(2),advance=.false.)
                call csv_write(30,month_met(d)%unorm(3),advance=.false.)
                call csv_write(30,month_met(d)%unorm(4),advance=.true.)
            end if
        end do

        ! set boundary conditions for next timestep
        ! left boundary condition: the end of the last month
        ! right boundary condition: we try to use the month after the next two months
        ! (see code above) but if that does not work, we use here the end of the
        ! month after the next month
        bcond_tmin(:) = (/ tmin_sm(ld), tmin_sm(rd) /)
        bcond_tmax(:) = (/ tmax_sm(ld), tmax_sm(rd) /)
        bcond_cloud(:) = (/ cloud_sm(ld), cloud_sm(rd) /)
        bcond_wind(:) = (/ wind_sm(ld), wind_sm(rd) /)

        stationid(1:n_tot-1) = stationid(2:n_tot)
        year(1:n_tot-1) = year(2:n_tot)
        month(1:n_tot-1) = month(2:n_tot)
        mprec(1:n_tot-1) = mprec(2:n_tot)
        mwet(1:n_tot-1) = mwet(2:n_tot)
        mtmin(1:n_tot-1) = mtmin(2:n_tot)
        mtmax(1:n_tot-1) = mtmax(2:n_tot)
        mcloud(1:n_tot-1) = mcloud(2:n_tot)
        mwind(1:n_tot-1) = mwind(2:n_tot)
        ndm(1:n_tot-1) = ndm(2:n_tot)
        lat(1:n_tot-1) = lat(2:n_tot)
        lon(1:n_tot-1) = lat(2:n_tot)

        !loop to next station/year/month
        if (end_counter == 0) then
            exit
        end if

    end do

    close(10)
    close(30)

    contains

    !-------------------------------------

    function are_consecutive_months(stationid, year, month) result(r)
        ! Check whether the given entries belong to consecutive months
        !
        ! This function checks whether the given `months` are consecutive. The
        ! return is a logical array with length ``size(stationid) - 1``.

        ! The first value is True if the second month is a consecutive month
        ! of the first month, the second value is True if the third month is
        ! a consecutive month of the second one, and so on. Months are only
        ! seen as consecutive if the `stationid` of the corresponding months
        ! are the same

        character(12), intent(in) :: stationid(:) ! the station ids of the entries
        integer, intent(in)       :: year(:) ! the years of the entries
        integer, intent(in)       :: month(:) ! the months of the entries
        integer                   :: i, r(size(stationid) - 1)

        do i=1,size(r)
            if (stationid(i) /= stationid(i + 1)) then
                r(i) = 0
            else
                if (month(i) == 12) then
                    r(i) = l2i(((year(i + 1) == year(i) + 1) .and. (month(i+1) == 1)))
                else
                    r(i) = l2i(((year(i + 1) == year(i)) .and. (month(i+1) == month(i) + 1)))
                end if
            end if
        end do

    end function are_consecutive_months

    !-------------------------------------

    integer function l2i(l)
        ! convenience function to convert a logical to an integer in [0, 1]

        logical :: l ! the logical to convert

        if (l) then
            l2i = 1
        else
            l2i = 0
        end if

    end function l2i

    !-------------------------------------

    integer function ndaymonth(year,month)
        ! get the number of days in a month given the year in calendar years AD and the month number

        integer, intent(in) :: year ! the year of the month
        integer, intent(in) :: month ! the month

        integer, parameter, dimension(12) :: std_year = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
        integer, parameter, dimension(12) :: leapyear = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

        if (mod(year,400) == 0) then

            ndaymonth = leapyear(month)

        else if (mod(year,100) == 0) then

            ndaymonth = std_year(month)

        else if (mod(year,4) == 0) then

            ndaymonth = leapyear(month)

        else

            ndaymonth = std_year(month)

        end if

    end function ndaymonth

!-------------------------------------

end program gwgen
