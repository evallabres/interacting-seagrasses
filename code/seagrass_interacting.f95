
 program model

!*************************
!     N species 
!*************************

    INTEGER :: N
    INTEGER :: tf, ti, iseed, e1, e2, ncells, nsx, nsy, naloop, nak, nrunt, ncx, ncy
    REAL :: Lx, Ly, rs, dtotal, timestep, nw, nrunrand, p
    REAL, PARAMETER :: pi = 3.1415927
    INTEGER :: nmorthi, nadeath1, nadeath2, nmarka, nmarkb
    REAL :: dn, nubdens, timep, mur0int, mur0f, mur0i, nmorth, nmorthrand, secit
    INTEGER :: na, itc, s,sp, partant, dpart, tm, picint, sm, q, qf, m, seedf

    REAL, ALLOCATABLE :: u(:,:,:), r(:,:,:), phi(:,:), npartcel(:,:,:), rpart(:,:,:), rprop(:)
    REAL, ALLOCATABLE :: ui(:,:,:), ri(:,:,:), phii(:,:), npartceli(:,:,:), rparti(:,:,:)
    REAL, ALLOCATABLE :: dtotalt(:,:,:,:), davg(:, :), phibgauss(:), rhogauss(:), dsat(:,:,:), dsatseed(:,:)
    INTEGER, ALLOCATABLE :: partinf(:,:,:), apexinf(:,:,:), partinfi(:,:,:), apexinfi(:,:,:)

    
    INTEGER, ALLOCATABLE ::  inipart(:), inin(:), inirangeix(:), inirangefx(:), inirangeiy(:), inirangefy(:), nnorm(:)
    INTEGER, ALLOCATABLE :: npartv(:,:,:,:), navg(:,:)
    REAL, ALLOCATABLE :: rho0(:), rhoe(:), rhomin(:), rhomax(:), nu0(:), nub0(:)
    REAL, ALLOCATABLE :: mur0(:), dt(:), phib0(:), phibe(:), phibmax(:), phibmin(:), phiw(:), rhomaxf(:), rhoint(:)
    REAL, ALLOCATABLE :: a(:,:), morty(:),  dmax(:),  dnorm(:), dini(:), alpha(:), nrunR(:)

    INTEGER, ALLOCATABLE :: ndeath(:), nram(:), nexcl(:), ndeathpart(:)
    INTEGER, ALLOCATABLE ::  nf(:), part(:), napex(:), npart(:), nmax(:), tci(:), nfant(:), npartant(:), nrun(:)
   

!*********
!   READ PARAMETERS
!*********

    open(unit=1, file='param.txt', status='old')
    read(1,*) N
    close(1)

    ALLOCATE(inin(N), inirangeix(N), inirangefx(N), inirangeiy(N), inirangefy(N), nnorm(N),nmax(N))
    ALLOCATE(rho0(N), rhoe(N), rhomin(N), rhomax(N), nu0(N), nub0(N))
    ALLOCATE(mur0(N), dt(N), phib0(N), phibe(N), phibmax(N), phibmin(N), phiw(N), rhomaxf(N), rhoint(N))
    ALLOCATE(a(N,N), morty(N), dmax(N),  dnorm(N), dini(N), alpha(N))

    open(unit=1, file='param.txt', status='old')

    read(1,*) N

    read(1,*) ti
    read(1,*) tf
    read(1,*) iseed
    read(1,*) seedf
    read(1,*) e1
    read(1,*) e2

    read(1,*) ncx
    read(1,*) ncy
    read(1,*) rs
    read(1,*) tm

    read(1,*) sm
    read(1,*) mur0i
    read(1,*) mur0f
    read(1,*) qf

    do s = 1, N
        read(1,*) rho0(s), rhoe(s), rhomin(s), rhomax(s)
        read(1,*) phib0(s), phibe(s), phibmin(s), phibmax(s)
        read(1,*) nu0(s)
        read(1,*) nub0(s)
        read(1,*) mur0(s)
        read(1,*) phiw(s)
        read(1,*) dmax(s)
        read(1,*) dnorm(s)
        read(1,*) dini(s)
        read(1,*) alpha(s)
        do sp = 1, N
            read(1,*) a(s,sp)
        enddo 
        read(1,*) morty(s)

        read(1,*) inirangeix(s)
        read(1,*) inirangefx(s)
        read(1,*) inirangeiy(s)
        read(1,*) inirangefy(s)
    enddo

    close(1)


!********** Calculate maximum number of particles and cells, and time

    nsx = int(ncx/2)
    nsy = int(ncy/2)

    Lx = real(ncx)*rs
    Ly = real(ncy)*rs

    naf = 10**e1-1
    np = 10**e2-1

    nmax = int(dmax*(rs/100)**2)

    nnorm = int(dnorm*(rs/100)**2)

    inin = int(int(dini)*((rs)/100)**2)

    do s = 1, N
        if (inin(s) == 0) then
            inin(s) = 1
        endif
    enddo
    
!******** time
!This is the time interval that I need for each apex, to develop a new shoot
    do s= 1,N
        dt(s)=(rho0(s)/nu0(s))
    enddo
!I choose the time step of the code, as the characteristic time of one of the species. So all apexs of one species, and the other, just a percentatge. 
    timestep = dt(1)

!*********
!   ALLOCATE
!*********

    
    ALLOCATE(r(N, naf, 2), u(N, naf, 2), rpart(N, np, 2), npartcel(N,-nsx:nsx,-nsy:nsy))
    ALLOCATE(ri(N, naf, 2), ui(N, naf, 2), rparti(N, np, 2), npartceli(N,-nsx:nsx,-nsy:nsy))
    ALLOCATE(phi(N,naf), phii(N,naf), dsat(N,qf,seedf), dsatseed(N,qf))
    ALLOCATE(npartv(N, tf, qf, seedf), navg(N, tf), rprop(2))
    ALLOCATE(partinf(N, np, 5), apexinf(N, 0:naf, 1), partinfi(N, np, 5), apexinfi(N, 0:naf, 1))
    ALLOCATE(phibgauss(1), rhogauss(naf), dtotalt(N, tf, qf, seedf), davg(N, tf))
    ALLOCATE(ndeath(N), nram(N), nexcl(N), ndeathpart(N),inipart(N), )
    ALLOCATE(nf(N), part(N), napex(N), npart(N), tci(N), nfant(N), npartant(N), nrun(N), nrunR(N))


!partinf has the information for every particle of partinf(s, na, it, apex1, apex2, celx, cely, celdx, celdy). It goes from 0:np particles because the zero-th particle is the one the corresponds to the apexes that have not ramified yet. Its a comodin for all of them

!*********
! INITIALIZE
!*********

    open(3, file = 'npart.txt')
    open(35, file = 'npart.avg.txt')
    open(4, file = 'density.txt')
    open(45, file = 'density.avg.txt')

    open(10, file = 'hyst.txt')
    open(20, file = 'pics.txt')
    open(30, file = 'time.txt')

    open(100, file = 'colordensity.1.txt')
    open(200, file = 'colordensity.2.txt')
    open(300, file = 'colordensity.3.txt')
    open(400, file = 'colordensity.4.txt')
    open(500, file = 'colordensity.5.txt')
    open(600, file = 'colordensity.6.txt')


!********************
!REALISATION LOOP
!********************
    dtotalt = 0
    npartv = 0

 loopseed: do m = 1, seedf
                iseed = iseed + 1
                call ran_ini(iseed)

    inipart = 0
    ri = 0
    phii = 0
    ui = 0
    npartceli = 0
    partinfi = 0
    apexinfi = 0

    mur0int = (mur0f - mur0i)/qf

 loopi:  do s = 1, N
!***************Negative x-cells
    do iceldx = inirangeix(s), -1
!************** Negative y-cells
                do iceldy = inirangeiy(s), -1
                    do k = 1, inin(s)
                        inipart(s) = inipart(s) + 1
                        rand = ran_u()
                        rand2 = ran_u()
                        rand5 = ran_u()
!**********************
                        ri(s, inipart(s), 1) = (iceldx+1)*rs - rand*rs
                        ri(s, inipart(s), 2) =(iceldy+1)*rs - rand2*rs
!*********** Random ramification unit vector
                        phii(s, inipart(s)) = 360 * rand5
                        ui(s, inipart(s), 1) = sin(phii(s, inipart(s))*(pi/180))
                        ui(s, inipart(s), 2) = cos(phii(s, inipart(s))*(pi/180))
!***********************
                        rparti(s, inipart(s), :) = ri(s, inipart(s), :)
                        partinfi(s, inipart(s),:) = (/  tci(s), inipart(s), 0, iceldx, iceldy/)
                        npartceli(s, iceldx, iceldy) = npartceli(s, iceldx, iceldy) + 1
                        apexinfi(s, inipart(s), 1) = inipart(s)
                    enddo
                enddo
!************** Positive y-cells
                do iceldy = 1, inirangefy(s)
                    do k = 1, inin(s)
                        inipart(s) = inipart(s) + 1
                        rand = ran_u()
                        rand2 = ran_u()
                        rand5 = ran_u()
!**********************
                        ri(s, inipart(s), 1) = (iceldx+1)*rs - rand*rs
                        ri(s, inipart(s), 2) = (iceldy-1)*rs + rand2*rs
!*********** Random ramification unit vector
                        phii(s, inipart(s)) = 360 * rand5
                        ui(s, inipart(s), 1) = sin(phii(s, inipart(s))*(pi/180))
                        ui(s, inipart(s), 2) = cos(phii(s, inipart(s))*(pi/180))
!***********************
                        rparti(s, inipart(s), :) = ri(s, inipart(s), :)
                        partinfi(s, inipart(s),:) = (/  tci(s), inipart(s), 0, iceldx, iceldy/)
                        npartceli(s, iceldx, iceldy) = npartceli(s, iceldx, iceldy) + 1
                        apexinfi(s, inipart(s), 1) = inipart(s)
                    enddo
                enddo
            enddo
!***************Positive x-cells
    do iceldx = 1, inirangefx(s)
!************** Negative y-cells
                do iceldy = inirangeiy(s), -1
                    do k = 1, inin(s)
                        inipart(s) = inipart(s) + 1
                        rand = ran_u()
                        rand2 = ran_u()
                        rand5 = ran_u()
!**********************
                        ri(s, inipart(s), 1) = (iceldx-1)*rs + rand*rs
                        ri(s, inipart(s), 2) =(iceldy+1)*rs - rand2*rs
!*********** Random ramification unit vector
                        phii(s, inipart(s)) = 360 * rand5
                        ui(s, inipart(s), 1) = sin(phii(s, inipart(s))*(pi/180))
                        ui(s, inipart(s), 2) = cos(phii(s, inipart(s))*(pi/180))
!***********************
                        rparti(s, inipart(s), :) = ri(s, inipart(s), :)
                        partinfi(s, inipart(s),:) = (/  tci(s), inipart(s), 0, iceldx, iceldy/)
                        npartceli(s, iceldx, iceldy) = npartceli(s, iceldx, iceldy) + 1
                        apexinfi(s, inipart(s), 1) = inipart(s)
                    enddo
                enddo
!************** Positive y-cells
                do iceldy = 1, inirangefy(s)
                    do k = 1, inin(s)
                        inipart(s) = inipart(s) + 1
                        rand = ran_u()
                        rand2 = ran_u()
                        rand5 = ran_u()
!**********************
                        ri(s, inipart(s), 1) = (iceldx-1)*rs + rand*rs
                        ri(s, inipart(s), 2) = (iceldy-1)*rs + rand2*rs
!*********** Random ramification unit vector
                        phii(s, inipart(s)) = 360 * rand5
                        ui(s, inipart(s), 1) = sin(phii(s, inipart(s))*(pi/180))
                        ui(s, inipart(s), 2) = cos(phii(s, inipart(s))*(pi/180))
!***********************
                        rparti(s, inipart(s), :) = ri(s, inipart(s), :)
                        partinfi(s, inipart(s),:) = (/  tci(s), inipart(s), 0, iceldx, iceldy/)
                        npartceli(s, iceldx, iceldy) = npartceli(s, iceldx, iceldy) + 1
                        apexinfi(s, inipart(s), 1) = inipart(s)
                    enddo
                enddo
            enddo
        enddo loopi

!************************
!INTRINSIC BRANCHING LOOP
!************************
       
    loopQ: do q = 1, qf

        mur0(sm) = mur0i + q*mur0int

        print *, 'iseed', iseed, mur0(sm)
        
!**************
!   INITIALIZE PARTICLES
!**************
    
        part = inipart
        npart = inipart
        nf = inipart
        npartant = inipart
        nfant = inipart
        napex = inipart
        nrun = inipart
        nrunt = 0
        do s=1,N
            nrunR(s) = timestep/dt(s)*real(nfant(s))
            nrun(s) = int(nrunR(s))
            nrunrand = nrunR(s) - real(nrun(s))
            rand = ran_u()
            if (nrunrand > rand) then
               nrun(s) = nrun(s) + 1
            endif
            nrunt = nrunt + nrun(s)
        enddo
        p = nrun(1)/real(nrunt)
        nram = 0
        ndeathpart = 0
        ndeath = 0
        nexcl = 0

        secit = 0.
        r = 0.
        u = 0.
        npartcel = 0
        rpart = 0
        partinf = 0
        apexinf = 0

        r = ri
        u = ui
        npartcel = npartceli
        rpart = rparti
        partinf = partinfi
        apexinf = apexinfi

 
!************************
!MAIN TIME LOOP
!************************
 loopt: do itc = 2, tf

!***********************
!PREPARING THE NEXT ITERATION
!***********************

            nfant =  napex
            npartant = npart
            nrunt = 0
            do s=1,N
                nrunR(s) = timestep/dt(s)*real(nfant(s))
                nrun(s) = int(nrunR(s))
                nrunrand = nrunR(s) - real(nrun(s))
                rand = ran_u()
                if (nrunrand > rand) then
                    nrun(s) = nrun(s) + 1
                endif
                nrunt = nrunt + nrun(s)
            enddo
            p = nrun(1)/real(nrunt)

            !write(30, *) nrunR(1), nrun(1), nrunR(2), nrun(2), nrunt

            !ndeath = 0
            !ndeathpart = 0
            !nram = 0
            !nexcl = 0

!************************
!APEX LOOP
!************************
            

 loopap:   do naloop = 1, nrunt

                ndeath = 0
                ndeathpart = 0
                nram = 0
                nexcl = 0

                rand = ran_u()
                if (rand <= p) then
                    s=1
                else
                    s=2
                endif

                        na = i_ran(nf(s))
                        !if (s==2) then
                            !write(30, *) na
                        !endif
                        rprop(:) = r(s, na, :) + rho0(s) * u(s, na, :)
              
!************************
!PERIODIC BOUNDARY CONDITIONS
!************************
                        if (rprop(1) >= Lx/2) then
                            rprop(1) = rprop(1) - Lx
                        endif
                        if (rprop(1) < -Lx/2) then
                            rprop(1) = rprop(1) + Lx
                        endif
                        if (rprop(2) >= Ly/2) then
                            rprop(2) = rprop(2) - Ly
                        endif
                        if (rprop(2) < -Ly/2) then
                            rprop(2) = rprop(2) + Ly
                        endif

!************************
!EXCLUSION PRINCIPLE in density form
!************************
                        if  (rprop(1)>=0) then
                            icelx = int((rprop(1)+rs)/rs)
                        elseif (rprop(1) > -Lx/2) then
                            icelx = int((rprop(1)-rs)/rs)
                        else
                            icelx = -nsx
                        endif
                        if  (rprop(2)>=0) then
                            icely = int((rprop(2)+rs)/rs)
                        elseif (rprop(2) > -Ly/2) then
                            icely = int((rprop(2)-rs)/rs)
                        else
                            icely = -nsy
                        endif
                     
!**************
                        nw = 0 
                        do sp = 1, N
                           nw = nw + a(s, sp)*npartcel(sp, icelx, icely)
                        enddo 
                        if (nw >= nmax(s)) then
!************* Excluded
                            nexcl(s) = nexcl(s) + 1
                        !if (s==2) then
                            !write(30, *) 'nexcl', nexcl(s)
                        !endif
                            goto 40
                        else
!************** The same particle must have associated 2 apexes, because it can have just ramified. I have to update the list properly, and erase this apex from the particle
                            
                            partant = apexinf(s, na, 1)
                            !if (s==2) then
                                !write(30, *) 'partant', partant
                            !endif
                            if (partinf(s,partant,2) == na) then
                                partinf(s, partant, 2) = 0
                            elseif (partinf(s,partant,3) == na) then
                                partinf(s, partant, 3) = 0
                            endif
!**************
                            r(s, na, :) = rprop(:)
                            npartcel(s, icelx, icely) = npartcel(s, icelx, icely) + 1
                            part(s) = part(s) + 1
                            npart(s) = npart(s) + 1
                            rpart(s, part(s), :) = r(s, na, :)
                            partinf(s, part(s), :) = (/itc, na, 0, icelx, icely/)
                            apexinf(s,na,1) = part(s)
                        endif


!************************
!RHIZOME BRANCHING
!************************
                            dn = 0
                            do sp = 1, N
                                dn = dn + a(s, sp)*npartcel(sp, icelx, icely)/(nnorm(s))
                            enddo 
                            nubdens = nub0(s)-alpha(s)/8 + alpha(s)*dn*(1 - dn)
                            ! Multiply for dt, the time that the apex does not elongate.  
                            ! This is the time that an apex at the same position has to ramify. 
                            prob = nubdens*dt(s)
                            rand = ran_u()
                            if (prob > rand .AND. r(s,na,1) /= 0 .AND. r(s,na,2) /= 0) then
                                call ran_gv(phibgauss, 1)
                                phibgauss = phib0(s) + phibe(s) * phibgauss
                                    if (phibgauss(1)>phibmax(s)) then
                                        phibgauss(1)=phibmax(s)
                                    elseif (phibgauss(1)<phibmin(s)) then
                                    phibgauss(1)=phibmin(s)
                                    endif
                                !Apex counting
                                nf(s) = nf(s) + 1
                                napex(s) = napex(s) + 1
                                nram(s) = nram(s) + 1
                                !Label of the ramification in the shoot
                                partinf(s, part(s), 3) = nf(s)
                                apexinf(s, nf(s), 1) = part(s)
                                !if (s==2) then
                                    !write(30, *) 'naram', nf(s)
                                !endif
                                rand2 = ran_u()
                                if (rand2 >= 0.5) then
                                    phi(s, nf(s)) = phi(s, na) + phibgauss(1)
                                else
                                    phi(s, nf(s)) = phi(s, na) - phibgauss(1)
                                endif
                                u(s, nf(s), 1) = sin(phi(s, nf(s))*(pi/180))
                                u(s, nf(s), 2) = cos(phi(s, nf(s))*(pi/180))
                                r(s, nf(s), :) = r(s, na, :)
                                endif

                                !do i= 1, part(1)
                                    !print *, i, partinf(1,i,2), partinf(1,i,3)
                                !enddo
                                !print *, '-----------'
                                !do i= 1, nf(1)
                                    !print *, i, apexinf(1,i,1)
                                !enddo
                                !print *, '-----------'
!***************
!RHIZOME WIGGLE
!***************
!(I give half of the probability to non-wiggles, and the other two quarts to wiggle to the left and to the right. Is it correct?)
                                rand2 = ran_u()
                                if (0.85 > rand2 .and. rand2 >= 0.5) then
                                    rand = ran_u()
                                    phi(s, na) = phi(s, na) + rand*phiw(s)
                                    u(s, na, 1) = sin(phi(s, na)*(pi/180))
                                    u(s, na, 2) = cos(phi(s, na)*(pi/180))
                                elseif (0.5 > rand2 .and. rand2 >= 0.15) then
                                    rand = ran_u()
                                    phi(s, na) = phi(s, na) - rand*phiw(s)
                                    u(s, na, 1) = sin(phi(s, na)*(pi/180))
                                    u(s, na, 2) = cos(phi(s, na)*(pi/180))
                                endif
!************************
!SHOOT MORTALITY
!************************
                        !if (itc > tm .and. s == 1) then
                        !    mur0(s) = mur0(s) + timestep*morty(s)
                        !endif
                        ! The number of shoot that will be dead at each time step
    40                  probmort =  npartant(s)*mur0(s)*timestep
                        ! I divide them by the number of apex runs I do in the time step. 
                        !nt = 0
                        !do sp = 1,N
                        !    nt = nt + nrun(sp)
                        !enddo
                        nmorth = probmort/real(nrun(s))
                        nmorthi = int(nmorth)
                        nmorthrand = nmorth - nmorthi
                        rand = ran_u()
                        if (nmorthrand > rand) then
                            nmorthi = nmorthi + 1
                        endif
                        !print *, mur0(s), mur0(s)*timestep, npart(s), napex(:), nt, nmorth, nmorthi
                        
           loopmort:    do k = 1, nmorthi
                             ! We chose particle we kill at random
                            dpart = i_ran(part(s))
                            !Find cell of the particle, and other info
                            icelx = partinf(s,dpart,4)
                            icely = partinf(s,dpart,5)
                            nadeath1 = partinf(s,dpart, 2)
                            nadeath2 = partinf(s,dpart, 3)
                            nmarka = 0
                            nmarkb = 0
        !********************** Update the list
                                rpart(s, dpart,:) = rpart(s, part(s),:)
                                partinf(s, dpart,:) = partinf(s, part(s),:)
                                rpart(s, part(s),:) = 0
                                partinf(s, part(s),:) = 0
                                npartcel(s, icelx, icely) = npartcel(s, icelx, icely) - 1
                                npart(s) = npart(s) - 1
                                part(s) = part(s)-1
                                ndeathpart(s) = ndeathpart(s) + 1
!********************** !Update the apex label of the moved particle. If it does not carries apexes, nak=0, and I will be updating apexinf(0).
                                !if (s==2) then
                                    !write(30, *) 'dpart', dpart
                                !endif
                                nak = partinf(s, dpart, 2)
                                apexinf(s, nak, 1) = dpart
                                nak2 = partinf(s, dpart, 3)
                                apexinf(s, nak2, 1) = dpart
!******************  The killed apexes might be the last in the list.
                            if (nadeath1 /= 0 .or. nadeath2 /= 0) then
                                if (nadeath1 /= 0 .and. nadeath1 == nf(s)) then
                                    u(s,nadeath1, :) = 0
                                    r(s,nadeath1, :) = 0
                                    phi(s,nadeath1) = 0
                                    apexinf(s,nadeath1,:) = 0
                                    ndeath(s) = ndeath(s) + 1
                                    napex(s) = napex(s) - 1
                                    nf(s) = nf(s) - 1
                                    nmarka = 1
                                endif
                                if (nadeath2 /= 0 .and. nadeath2 == nf(s)) then
                                    u(s,nadeath2, :) = 0
                                    r(s,nadeath2, :) = 0
                                    phi(s,nadeath2) = 0
                                    apexinf(s,nadeath2,:) = 0
                                    ndeath(s) = ndeath(s) + 1
                                    napex(s) = napex(s) - 1
                                    nf(s) = nf(s) - 1
                                    nmarkb = 1
                                endif
                                if (nadeath1 /= 0 .and. nadeath1 == nf(s)) then
                                    u(s,nadeath1, :) = 0
                                    r(s,nadeath1, :) = 0
                                    phi(s,nadeath1) = 0
                                    apexinf(s,nadeath1,:) = 0
                                    ndeath(s) = ndeath(s) + 1
                                    napex(s) = napex(s) - 1
                                    nf(s) = nf(s) - 1
                                    nmarka = 1
                                    !goto 8
                                endif
                                if (nadeath1 /= 0 .and. nmarka == 0) then
                                    npk = apexinf(s, nf(s), 1)
                                    if (partinf(s,npk,2) == nf(s)) then
                                        partinf(s, npk, 2) = nadeath1
                                    elseif (partinf(s,npk,3) == nf(s)) then
                                        partinf(s, npk, 3) = nadeath1
                                    endif
                                    u(s,nadeath1, :) = u(s,nf(s), :)
                                    r(s,nadeath1, :) = r(s,nf(s), :)
                                    phi(s,nadeath1) = phi(s,nf(s))
                                    apexinf(s,nadeath1,:) = apexinf(s,nf(s),:)
                                    u(s,nf(s), :) = 0
                                    r(s,nf(s), :) = 0
                                    phi(s,nf(s)) = 0
                                    apexinf(s,nf(s),:) = 0
                                    ndeath(s) = ndeath(s) + 1
                                    napex(s) = napex(s) - 1
                                    nf(s) = nf(s) - 1
                                endif
                                if (nadeath2 /= 0 .and. nmarkb == 0) then
                                    npk = apexinf(s, nf(s), 1)
                                    if (partinf(s,npk,2) == nf(s)) then
                                        partinf(s, npk, 2) = nadeath2
                                    elseif (partinf(s,npk,3) == nf(s)) then
                                        partinf(s, npk, 3) = nadeath2
                                    endif
                                    u(s,nadeath2, :) = u(s,nf(s), :)
                                    r(s,nadeath2, :) = r(s,nf(s), :)
                                    phi(s,nadeath2) = phi(s,nf(s))
                                    apexinf(s,nadeath2,:) = apexinf(s,nf(s),:)
                                    u(s,nf(s), :) = 0
                                    r(s,nf(s), :) = 0
                                    phi(s,nf(s)) = 0
                                    apexinf(s,nf(s),:) = 0
                                    ndeath(s) = ndeath(s) + 1
                                    napex(s) = napex(s) - 1
                                    nf(s) = nf(s) - 1
                                endif
                            endif
     !8
                        enddo loopmort

                        if (napex(s)<=0) then
                            print *, 'the bush', s , 'has no more apexes at time', itc*dt(s), iseed
                            if (s==1) then
                                p = 0.
                            else
                                p = 1
                            endif
                            !The leftover particles are not killed by the simulation,
                            ! But with this we make them invisible for the other plant,
                            ! and their density is zero
                            npartcel(s,:,:) = 0
                            partinf(s,:,:) = 0
                            npart(s) = 0
                        endif
        

!************************
!END APEX LOOP
!************************

    !if (s==2) then
        !do i= 1, part(2)
            !write(30, *) i, partinf(2,i,2), partinf(2,i,3)
        !enddo
        !write(30, *) '-----------'
        !do i= 1, nf(2)
            !write(30, *) i, apexinf(2,i,1)
        !enddo
        !write(30, *) '-----------'
        !do s = 1, N
            !write(30, *) 'itc=', itc
            !write(30, *) s, npart(s), ndeathpart(s)
           !write(30, *) s , napex(s), ndeath(s), nram(s), nexcl(s), nrun(s), p
        !enddo
        
        !write(30, *) '----------------'
    !endif


        enddo loopap


!***********************
!        AVERAGED DENSITY AND PARTICLE NUMBER
!***********************
                

        secit = secit + timestep


!***************** Density

        npartv(:, itc, q, m) = npart(:)
        write(3,*) secit, npartv(:, itc, q, m)

!***************** Number of particles
        do s = 1, N
                ncells = 0
                dtotal = 0
!***************** Find the furthest cells
                nsmaxx = maxval(partinf(s,:,4))
                nsminx = minval(partinf(s,:,4))
                nsmaxy = maxval(partinf(s,:,5))
                nsminy = minval(partinf(s,:,5))
                do iceldx = nsminx, nsmaxx
                    do iceldy = nsminy, nsmaxy
                        if (npartcel(s, iceldx, iceldy)/= 0) then
                            dtotal = dtotal + npartcel(s, iceldx, iceldy)/(rs)**2
                            ncells = ncells + 1
                        endif
                    enddo
                enddo
                if (ncells /= 0) then
                    dtotalt(s, itc, q, m) = dtotal/ncells
                    write(s, *) secit, dtotalt(s, itc, q, m)*10000
                else
                    dtotalt(s, itc, q, m) = 0
                    write(s, *) secit, dtotalt(s, itc, q, m)*10000
                endif     

                dsat(s,q,m) = dtotalt(s, itc, q, m)*10000
         enddo


         write(4, *) secit, dtotalt(:, itc, q, m)*10000




!************************
!PICTURES 
!************************

    picint = int(tf/6)

 ld: if ( itc == picint .or. itc == 2*picint .or. itc == 3*picint .or. &
                itc ==  4*picint .or. itc == 5*picint .or. itc == 6*picint) then

 la: if (N==2) then

    do iceldx = -nsx, -1
        do iceldy = -nsy, -1
            if (npartcel(1, iceldx, iceldy)/= 0 .or. npartcel(2, iceldx, iceldy)/= 0) then
                write(itc/picint*100,*)   (-rs/2 + rs*(iceldx + 1))/100, (-rs/2 + rs*(iceldy + 1))/100, &
                            (npartcel(2,iceldx, iceldy)- npartcel(1, iceldx, iceldy))/rs**2*10000, &
                            npartcel(2,iceldx, iceldy), npartcel(1, iceldx, iceldy)
            endif
        enddo
        do iceldy = 1,nsy
            if (npartcel(1, iceldx, iceldy)/= 0 .or. npartcel(2, iceldx, iceldy)/= 0) then
                write(itc/picint*100,*)  ( -rs/2 + rs*(iceldx + 1))/100, (rs/2 + rs*(iceldy - 1))/100, &
                            (npartcel(2,iceldx, iceldy)- npartcel(1, iceldx, iceldy))/rs**2*10000, &
                            npartcel(2,iceldx, iceldy), npartcel(1, iceldx, iceldy)
            endif
        enddo
        write(itc/picint*100,*) 
    enddo
    do iceldx = 1, nsx
        do iceldy = -nsy, -1
            if (npartcel(1, iceldx, iceldy)/= 0 .or. npartcel(2, iceldx, iceldy)/= 0) then
                write(itc/picint*100,*)   (rs/2 + rs*(iceldx - 1))/100, (-rs/2 + rs*(iceldy + 1))/100, &
                            (npartcel(2,iceldx, iceldy)- npartcel(1, iceldx, iceldy))/rs**2*10000, &
                            npartcel(2,iceldx, iceldy), npartcel(1, iceldx, iceldy)
            endif
        enddo
        do iceldy = 1,nsy
            if (npartcel(1, iceldx, iceldy)/= 0 .or. npartcel(2, iceldx, iceldy)/= 0) then
                write(itc/picint*100,*)   (rs/2 + rs*(iceldx - 1))/100, (rs/2 + rs*(iceldy - 1))/100, &
                            (npartcel(2,iceldx, iceldy)- npartcel(1, iceldx, iceldy))/rs**2*10000, &
                            npartcel(2,iceldx, iceldy), npartcel(1, iceldx, iceldy)
                         
            endif
        enddo
        write(itc/picint*100,*) 
    enddo
    
        if (itc == picint) then
            call system('gnuplot -p colordensity.1.gnu')
        elseif (itc == 2*picint) then
            call system('gnuplot -p colordensity.2.gnu')
        elseif (itc == 3*picint) then
            call system('gnuplot -p colordensity.3.gnu')
        elseif (itc == 4*picint) then
            call system('gnuplot -p colordensity.4.gnu')
        elseif (itc == 5*picint) then
            call system('gnuplot -p colordensity.5.gnu')
        elseif (itc == 6*picint) then
            call system('gnuplot -p colordensity.6.gnu')
        endif
        
        write(20,*) timestep*itc, 0

    endif la

    endif ld

!************************
!END MAIN TIME LOOP
!************************

     !do i= 1, part(2)
     !    print *, i, partinf(2,i,2), partinf(2,i,3)
     !enddo
     !print *, '-----------'
     !do i= 1, nf(2)
     !    print *, i, apexinf(2,i,1)
     !enddo
     !print *, '----------------'
     !print *, '----------------'
     !do s = 1, N
     !       print *, itc
     !       print *, s, npart(s), ndeathpart(s)
     !       print *, s , napex(s), ndeath(s), nram(s), nexcl(s), nrun(s), p
     !enddo
     !print *, '---------------- ' , iseed
    
     if (npart(2) < napex(2)) then
        print *, iseed, npart(2),  napex(2), 'problem'
     endif

    enddo loopt

    write(3, *)
    write(4, *)

!************************
!END MORTALITY LOOP
!************************

 enddo loopQ


!************************
!END SEED LOOP
!************************

 enddo loopseed

    dsatseed = 0
    do q = 1, qf
        do m = 1, seedf
            do s = 1, N
                dsatseed(s, q)  = dsatseed(s, q) + dsat(s,q,m)
            enddo
        enddo
        write(10,*) mur0i + mur0int*q, dsatseed(:, q)/seedf !, dsatseed(2, q)/seedf
    enddo


    do q = 1, qf
        davg = 0
        navg = 0
        secit = 0
        do it = 2, tf
            secit = secit + dt(1)
            do m = 1, seedf
                do s = 1, N
                    navg(s, it)  = navg(s, it) + npartv(s, it, q, m)
                    davg(s, it)  = davg(s, it) + dtotalt(s, it, q, m)
                enddo
            enddo
            write(35, *) mur0i + mur0int*q, secit, navg(:, it)/seedf
            write(45, *) mur0i + mur0int*q, secit, davg(:, it)*10000/seedf
        enddo
            write(35, *)
            write(45, *)
    enddo

    call system('gnuplot -p density.gnu')
    call system('gnuplot -p density.avg.gnu')
    call system('gnuplot -p npart.gnu')
    call system('gnuplot -p npart.avg.gnu')
    call system('gnuplot -p hyst.gnu')
    call system('gnuplot -p densrel.gnu')
    call system('gnuplot -p densrel.2.gnu')
    


    CALL CPU_TIME(timep)
    write(30, *) 'time model', timep/60, 'minutes.'

    end program model
    


!!**********************

        subroutine ran_ini(iseed0)
            double precision dseed
                parameter(ip=1279)
            parameter(np=14)
            parameter(nbit=31)
            parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
            integer ix(ip)
            dimension g(0:m)

            data c0,c1,c2/2.515517,0.802853,0.010328/
            data d1,d2,d3/1.432788,0.189269,0.001308/
        
            common /ixx/ ix
            common /icc/ ic
            common /gg/ g
        
            dseed=iseed0
                do i=1,ip
                ix(i)=0
                do j=0,nbit-1
                if(rand_xx(dseed).lt.0.5) ix(i)=ibset(ix(i),j)
                enddo
                enddo
            ic=0
        
            pi=4.0d0*datan(1.0d0)
            do i=m/2,m
            p=1.0-real(i+1)/(m+2)
            t=sqrt(-2.0*log(p))
            x=t-(c0+t*(c1+c2*t))/(1.0+t*(d1+t*(d2+t*d3)))
            g(i)=x
            g(m-i)=-x
            enddo

            u2th=1.0-real(m+2)/m*sqrt(2.0/pi)*g(m)*exp(-g(m)*g(m)/2)
            u2th=nn1*sqrt(u2th)
            do i=0,m
                g(i)=g(i)/u2th
            enddo

            return
            end

                subroutine ran_read(iunit)
                parameter(ip=1279)
            parameter(np=14)
            parameter(m=2**np)
                integer ix(ip)
            dimension g(0:m)
                common /ixx/ ix
                common /icc/ ic
            common /gg/ g
                read (iunit,*) ic
                read (iunit,*) (ix(i),i=1,ip)
            read (iunit,*) (g(i),i=0,m)
                return
                end

                subroutine ran_write(iunit)
                parameter(ip=1279)
            parameter(np=14)
            parameter(m=2**np)
                integer ix(ip)
            dimension g(0:m)
            common /ixx/ ix
            common /icc/ ic
            common /gg/ g
                write (iunit,*) ic
                write (iunit,*) (ix(i),i=1,ip)
            write (iunit,*) (g(i),i=0,m)
                return
                end

                function i_ran(n)

                parameter(ip=1279)
            parameter(iq=418)
            parameter(is=ip-iq)
            integer ix(ip)
            common /ixx/ ix
            common /icc/ ic
            ic=ic+1
                if(ic.gt.ip) ic=1
            if(ic.gt.iq) then
                ix(ic)=ieor(ix(ic),ix(ic-iq))
            else
                    ix(ic)=ieor(ix(ic),ix(ic+is))
                endif
                i_ran=ix(ic)
                if (n.gt.0) i_ran=mod(i_ran,n)+1
                return
            end

            function ran_u()
                parameter(ip=1279)
            parameter(iq=418)
            parameter(is=ip-iq)
                parameter (rmax=2147483647.0)
            integer ix(ip)
            common /ixx/ ix
            common /icc/ ic
            ic=ic+1
                if(ic.gt.ip) ic=1
            if(ic.gt.iq) then
                ix(ic)=ieor(ix(ic),ix(ic-iq))
            else
                    ix(ic)=ieor(ix(ic),ix(ic+is))
                endif
                ran_u=real(ix(ic))/rmax
            return
            end


            function ran_g()
                parameter(ip=1279)
            parameter(iq=418)
            parameter(np=14)
            parameter(nbit=31)
            parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
            parameter(is=ip-iq)

            integer ix(ip)
            dimension g(0:m)

            common /ixx/ ix
            common /icc/ ic
            common /gg/ g

            ic=ic+1
                if(ic.gt.ip) ic=1
            if(ic.gt.iq) then
                ix(ic)=ieor(ix(ic),ix(ic-iq))
            else
                    ix(ic)=ieor(ix(ic),ix(ic+is))
                endif
            i=ishft(ix(ic),-np1)
            i2=iand(ix(ic),nn)
            ran_g=i2*g(i+1)+(nn1-i2)*g(i)
            return
            end


            subroutine ran_bm(x1,x2)
                parameter(ip=1279)
            parameter(iq=418)
            parameter(is=ip-iq)
                parameter (rmax=2147483647.0)
            integer ix(ip)
            common /ixx/ ix
            common /icc/ ic
            data pi2 /6.2831853072/
            ic=ic+1
                if (ic.gt.ip) ic=1
            if(ic.gt.iq) then
                ix(ic)=ieor(ix(ic),ix(ic-iq))
            else
                    ix(ic)=ieor(ix(ic),ix(ic+is))
                endif
            u=pi2*real(ix(ic))/rmax
            ic=ic+1
                if(ic.gt.ip) ic=1
            if(ic.gt.iq) then
                ix(ic)=ieor(ix(ic),ix(ic-iq))
            else
                    ix(ic)=ieor(ix(ic),ix(ic+is))
                endif
            v=(real(ix(ic))+0.5)/rmax
            v=sqrt(-2.0*log(v))
            x1=cos(u)*v
            x2=sin(u)*v
            return
            end


            subroutine ran_gv(u,n)
            parameter(ip=1279)
            parameter(iq=418)
            parameter(np=14)
            parameter(nbit=31)
            parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
            parameter(is=ip-iq)
            dimension g(0:m)
            dimension u(n)
            dimension ix(ip)
            common /gg/ g
            common /ixx/ ix
            common /icc/ic

                n1=0
        10      if (ic.lt.iq) then
                  kmax=min(n-n1,iq-ic)
              do 99 k=1,kmax
              ic=ic+1
              ix(ic)=ieor(ix(ic),ix(ic+is))
              i=ishft(ix(ic),-np1)
              i2=iand(ix(ic),nn)
              u(n1+k)=i2*g(i+1)+(nn1-i2)*g(i)
        99        continue
                else
                  kmax=min(n-n1,ip-ic)
              do 98 k=1,kmax
              ic=ic+1
              ix(ic)=ieor(ix(ic),ix(ic-iq))
              i=ishft(ix(ic),-np1)
              i2=iand(ix(ic),nn)
              u(n1+k)=i2*g(i+1)+(nn1-i2)*g(i)
        98        enddo
                endif
                if(ic.ge.ip) ic=0
                n1=n1+kmax
                if (n1.lt.n) goto 10
                
            return
            end


              function rand_xx(dseed)
              double precision a,c,xm,rm,dseed
              parameter (xm=2.d0**32,rm=1.d0/xm,a=69069.d0,c=1.d0)
              dseed=mod(dseed*a+c,xm)
              rand_xx=dseed*rm
              return
              end

!********************************
