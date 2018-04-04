program rod 

  use mcnp_random 
  implicit none

  integer, parameter :: I8 = selected_int_kind(18)
  integer, parameter :: R8 = selected_real_kind(15,307)
  integer, parameter :: n0 = 1000
  integer, parameter :: source_n=0
  integer, parameter :: n_max = (n0*1.4+source_n)
  integer, parameter :: nrun = 100
  integer, parameter :: nrun_throw = 50
  integer, parameter :: run_act=nrun-nrun_throw
  integer, parameter :: sites_bin=25
  real(R8), parameter :: length =4.0d+0
  integer, parameter :: bins= 101
  real(R8), parameter :: binwidth=length/(bins) 
  real(R8), parameter :: loc0 = length/2.0d+0
  real(R8), parameter :: nubar = 2.414d+0
  real(R8), parameter :: error = 1.0/n0*1.0
  real(R8), parameter :: source_pos=0d+0
  real(R8), parameter :: bias=0.0
  integer  :: i, absorb, leak, j, knt
  real(R8) :: locn, rn, rx, den, Navg, wght, kinf, L2,D, kpath,kpathold, xpos_ratio, kact
  integer :: n, run, coll, mult, k, source
  real(R8) :: sigc, sigs,sigf, siga, sigt, omega, PI, Clength, keff, kpathrat, nsig, navgn, pj, Hs, Hsold
  real(R8), dimension (1:bins) ::xpos, xposnew
  real(R8), dimension (1:bins) :: flux
  real(R8), dimension (1:8) :: fis_dist
  real(R8), dimension (1:((bins)/sites_bin)) :: shannon
  real(R8)  ::   wtot, wcum, prob, posf, poss, start, finish, distance
  real(R8), dimension (1:nrun-nrun_throw) :: nmax

!////////////////////////////////////////////////////////////////////////////////////   
  
   type neutron
       real(R8):: number_n
       real(R8):: nposition
       real(R8):: angle_n
       real(R8):: move_n
   end type neutron

type(neutron) fbank(1:int(n_max)), fbanknew(1:int(n_max)),keep(1:int(n_max))

!//////////////////////////////////////////////////////////////////////////////////////
 
    den=19.1d+0           !g/cc
    Navg=6.022d+23       !avagodros number
    wght=238.0289d+0       !molecular weight
    sigc= 3.4d-24       !cm^2
    sigf= 4.19d-24
    sigs= 10.0d-24
    sigc=den*Navg*sigc/(wght)
    sigf=den*Navg*sigf/(wght)
    sigs=den*Navg*sigs/(wght)
    siga= sigc+sigf
    sigt= siga+sigs
    fis_dist=(/ 0.0317223, 0.2034294, 0.5396285, 0.843598, 0.9705439, 0.9972232, 0.9998554, 1.00/)

!-------------------------------Actual keff---------------------------------------

    !Assumes a bare slab
    kinf=nubar*sigf/siga
    D=1.0d+0/(3.00*sigt)!*(1.0d+0-2.00/(3.00*int(wght))))
    L2=D/siga
    PI=4.D0*DATAN(1.D0)
    Clength=PI*(L2/(kinf-1d+0))**0.5d+0-2.0d+0*2.1312d+0*D
    kact=(nubar*sigf)/(siga+D*(PI/(length+2.1312*2.0*D))**2.0d+0)
    
!////////////////////   Initializes components   ////////////////////////////////////  

    n=n0
    do i=1,int(n_max)
        fbank(i)%number_n= 0.0
        fbanknew(i)%number_n=0.0
        fbank(i)%nposition=loc0
        fbanknew(i)%nposition=fbank(i)%nposition
        fbank(i)%angle_n=2.0
        fbanknew(i)%angle_n=fbank(i)%angle_n

    end do
    do i=1,int(n0+source_n)

        fbank(i)%number_n= 1.0
        fbanknew(i)%number_n=1.0
    end do
    
    run=0
    write(*,'(/)')
    write(*,'(2x,A,F4.1,A,F4.1,A)'), "The Bias is ",(1+bias)/0.02,"% Forward and  ",dabs(bias-1)/0.02,"% Backward"
!////////////////  The beginning of the successive runs   ///////////////////////////////

    call RN_init_problem( 1234567_I8, 1 )

do while (run<=nrun)                                 !ensures the number of iterations occur
    xpos_ratio=20.0d+0
    kpathrat=10.5d+0
    kpathold = 200d+0
    Hsold=10.0
    call cpu_time(start)
    
    do while (dabs(kpathrat)>error .or. xpos_ratio>error)  !ensures keff and flux is converged before the next run
      leak= 0; absorb=0;
      mult=0; coll=0;distance=0.0;
      i=1
      xpos(1:bins)=0.0;
      source=source_n
      call RN_init_particle( int(i,I8) )

      do while(sum(fbank%number_n)>=1.0)
        if (i>n_max) i=1;                               !index remornalization
                    
        if (fbank(i)%number_n<1.00) then                !No particle present
            i=i+1
        else                                            !particle exist and moves
            fbank(i)%move_n=-log(rang())/sigt;          !how far the particle moves
            if ((fbank(i)%angle_n)==2.0) then
                omega=(2.00*rang()-1.00)
                fbank(i)%angle_n=omega                  !direction of movement
!                
            end if
            distance=distance+(fbank(i)%move_n)
!=========================   Movement   ===========================================

            locn=(fbank(i)%nposition+fbank(i)%move_n*fbank(i)%angle_n)/binwidth+1.0
            poss=(fbank(i)%nposition/binwidth)+1.0
                        
!==================== Leakage out left or right ======================================

            if (length<(locn-1.0)*binwidth) then              ! right
                leak=leak+1
                fbank(i)%number_n=fbank(i)%number_n-1.0
                fbanknew(i)%number_n=fbanknew(i)%number_n-1.0
                fbanknew(i)%angle_n=2.0
                fbank(i)%angle_n=2.0
                xpos(ceiling(poss):bins)=xpos(ceiling(poss):bins)+binwidth
                xpos(floor(poss))=xpos(floor(poss))+binwidth*(ceiling(poss)-poss)
                distance=distance+length-(locn-1.0)*binwidth
                
            else if((locn-1.0)*binwidth<0.0) then               ! left
                leak=leak+1
                fbank(i)%number_n=fbank(i)%number_n-1.0
                fbanknew(i)%number_n=fbanknew(i)%number_n-1.0
                fbanknew(i)%angle_n=2.0
                fbank(i)%angle_n=2.0
                xpos(1:floor(poss))=xpos(1:floor(poss))+binwidth
                xpos(ceiling(poss))=xpos(ceiling(poss))+binwidth*(poss-floor(poss))
                distance=distance-(locn-1.0)*binwidth+(poss-1.0)*binwidth 
            else

!-------------------------- Flux math -----------------------------

                if (fbank(i)%angle_n>0.0) then
                    xpos(ceiling(poss):floor(locn))=xpos(ceiling(poss):floor(locn))+binwidth
                    xpos(ceiling(locn))=xpos(ceiling(locn))+binwidth*(locn-floor(locn))
                    xpos(floor(poss))=xpos(floor(poss))+binwidth*(ceiling(poss)-poss)
                else
                    xpos(ceiling(locn):floor(poss))=xpos(ceiling(locn):floor(poss))+binwidth
                    xpos(floor(locn))=xpos(floor(locn))+binwidth*(ceiling(locn)-locn)
                    xpos(ceiling(poss))=xpos(ceiling(poss))+binwidth*(poss-floor(poss)) 
                end if
            
!----------------------Determine the type of reaction-------------------------------
            
                rn=rang();
                if (rn<sigc/sigt) then               !capture
                    absorb =absorb+1
                    fbank(i)%number_n=fbank(i)%number_n-1.0
                    fbanknew(i)%number_n=fbanknew(i)%number_n-1.0
                    fbanknew(i)%angle_n=2.0
                    fbank(i)%angle_n=2.0

                else if (rn<(sigc+sigf)/sigt) then   !fission
                    rx=rang();
                    absorb=absorb+1;
                    k=1
                    do while (rx>fis_dist(k))
                        k=k+1
                    end do
                    mult=mult+k  
                    do while(k>1.0)
                        j=1
                        do while (fbanknew(j)%number_n>0.0)
                            j=j+1
                        end do   
                        fbanknew(j)%number_n=1.0
                        fbanknew(j)%nposition=fbank(i)%nposition+fbank(i)%move_n*fbank(i)%angle_n
                        omega=(2.00*rang()-1.00+bias)
                        if(dabs(omega)>1.0) omega=sign(1.0d+0,omega)
                        fbanknew(j)%angle_n=omega  
                        k=k-1
                    end do 
                    fbank(i)%number_n=fbank(i)%number_n-1.0
                    fbanknew(i)%number_n=fbanknew(i)%number_n-1.0
                    fbanknew(i)%angle_n=2.0
                    fbank(i)%angle_n=2.0                               
                        
                else                                  !scatter
                        
                    fbanknew(i)%nposition=fbank(i)%nposition+fbank(i)%move_n*fbank(i)%angle_n 
                    fbank(i)%nposition=fbanknew(i)%nposition 
                    fbanknew(i)%angle_n=(2.0d+0*rang()-1.0d+0)
                    fbank(i)%angle_n=fbanknew(i)%angle_n
                    coll=coll+1
                    
                end if

!---------------------------------------------------------------------------------------

            end if             !movement
        end if                 !particle present
      end do                   !sum(fbank)=0

!-------------------------------Calculated keff----------------------------------------

       kpath=(nubar*sigf/sigt*(coll+absorb)/n+kinf*absorb/n)/2.0d+0!nubar*sigf*distance/n
       kpathrat=(kpath/kpathold-1.0d+0)!*(kpath/kpathold)**run
       kpathold=kpath 
!print *, kpathrat 
      
!--------------------------Fission Bank Renormalization----------------------------------

         wtot=sum(fbanknew%number_n)
         wcum=0.0
         k=0
!print *, wtot
         do j=1,int(n_max)
             prob=fbanknew(j)%number_n*real(n-k)/(wtot-wcum)
             wcum=wcum + fbanknew(j)%number_n
             knt =prob + rang()
             do i=1,int(knt)
                 k=k+1;
                 keep(k)=fbanknew(j)
                 
             end do
         end do

         n=sum(keep%number_n)


!------------------------ Renormalizing the fission bank ------------------------------

         k=1
         do while (k<int(n_max))
             fbank(k)%number_n=0.0;fbank(k)%nposition=0.0;fbank(k)%angle_n=0.0;fbank(k)%move_n=0.0;
             fbanknew(k)%number_n=0.0;fbanknew(k)%nposition=0.0;fbanknew(k)%angle_n=0.0;fbanknew(k)%move_n=0.0;
             k=k+1
         end do

         k=1
         do while(sum(keep%number_n)>0) 
             j=1
             do while (keep(j)%number_n<1)
                 j=j+1
             end do
             fbank(k)=keep(j)
         fbanknew(k)=fbank(k)
             k=k+1
             keep(j)%number_n=0.0
         end do

         k=1 
!         do while(k<n)
!             if (fbank(k)%nposition==fbank(k+1)%nposition .and. fbank(k)%angle_n==fbank(k+1)%angle_n) then
!                fbank(k)%angle_n=2.0
!                fbanknew(k)%angle_n=2.0
!             end if
!             k=k+1
!         end do 
         
!-----------------------------Shannon Entropy------------------------------------------

    do i=1,int(bins/sites_bin)
      pj=0
      do j=int(1+sites_bin*(i-1)),int(1+sites_bin*i)
        pj=pj+ xpos(j)
      end do
      shannon(i)=pj/sum(xpos)*log(pj/sum(xpos))/log(2.0)
    end do
    Hs=-sum(shannon)
    xpos_ratio=dabs(Hs/Hsold-1.0)
    Hsold=Hs
!print *, dabs(kpathrat), dabs(kpath), run!1.00+real(mult-absorb-leak)/n, n

!----------------------------- Add source----------------------------------------------
        k=1
        do while (source>0.0)
           fbank(n+k)%number_n=1.0
           fbanknew(n+k)%number_n=1.0
           fbank(n+k)%angle_n=2.0
           fbanknew(n+k)%angle_n=2.0
           fbanknew(n+k)%nposition=source_pos
           fbanknew(n+k)%nposition=source_pos
           k=k+1
           source=source-1
        end do

!-------------------------------------------------------------------------------------
    end do          !Keff


    if (run>nrun_throw) then
      xposnew(1:bins)=(xpos(1:bins)/maxval(xpos)+xposnew(1:bins)*(run-nrun_throw))/(run+1.00-nrun_throw)
      keff=(keff*(run-nrun_throw)+kpath)/(run+1.00-nrun_throw)
      nmax(run-nrun_throw)= keff
      
    else
      xposnew(1:bins)=xpos(1:bins)/maxval(xpos)
      keff=kpath
    end if

if (run-nrun_throw==0) print *,  "------------------Start Count----------------------------"
if (mod(run,int(nrun/5))==0) write(*,'(10x,A,I3,5x,A,F6.4)') 'run ',run,'keff=', keff

    run= run+1
    call cpu_time(finish)
!print *, finish-start
!print *,  'Time =',finish-start

    end do    !runs

!---------------------------------------------  d error  ----------------------------------------------------------
   navgn=sum(nmax)/(run-nrun_throw)
   nsig=(1.00/(run-nrun_throw-1.00)*1/(run-nrun_throw)*sum((nmax(1:nrun-nrun_throw)**2.00-navgn**2.00)))
   nsig=sqrt(abs(nsig))


!-----------------------------------flux calc---------------------------------------

i=1
flux(1:bins)=0
posf=0
do while(i<=bins)
  flux(i)=(dsin(PI*(posf+2.1312*D)/(length+2.0*2.1312*D))+dsin(PI*(binwidth+posf+2.1312*D)/(length+2.0*2.1312*D)))/2.0
  posf=posf+binwidth
  i=i+1
end do


!---------------------------------Print stuff--------------------------------------

print *, "kact= ", kact
print *, "keff= ", keff, "sig",nsig
print *, "dkeff[$]=", abs(keff-kact)/.0065
 print *, maxloc(xposnew), maxloc(flux) 
write(*,'(/)')


 OPEN(UNIT=12, FILE="rodn.txt", ACTION="write", STATUS="replace")
 write(12,*), xposnew/maxval(xposnew)
 write(12,*), flux 
 close(12)

end program rod
