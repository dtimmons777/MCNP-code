program rod 

  use mcnp_random 
  implicit none

  integer, parameter :: I8 = selected_int_kind(18)
  integer, parameter :: R8 = selected_real_kind(15,307)
  integer, parameter :: n0 = 1000
  integer, parameter :: nrun = 100
  integer, parameter :: nrun_throw = 20
  integer, parameter :: run_act=nrun-nrun_throw
  integer, parameter :: sites_bin=10
  real(R8), parameter :: length =5.0d+0
  integer, parameter :: bins= 101
  real(R8), parameter :: binwidth=length/(bins-1) 
  real(R8), parameter :: loc0 = 2d+0
  real(R8), parameter :: nubar = 2.414d+0
  real(R8), parameter :: error = .001
  integer, parameter :: source=100
  real(R8), parameter :: source_pos=2d+0
  integer  :: i, absorb, leak,j,ii, jj,knt
  real(R8) :: total, locn, rn, direction, rx, den, Navg, wght, kinf, L2,D, kpath,kpathold, pos, xpos_ratio, kact
  integer :: bin0, n,run, particle, nleft, coll, inscat, mult,k
  real(R8) :: sigc, sigs,sigf, siga, sigt, move, omega, PI, Clength, keff, kpathrat, errork,nmax,nmin,nsig,navgn, pj, Hs, Hsold
  real(R8), dimension (1:bins) ::xpos, xposnew, flux, xposmid 
  real(R8), dimension (1:8) :: fis_dist
  real(R8), dimension (1:((bins)/sites_bin)) :: shannon
  real(R8)  ::   wtot, wcum, prob
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
   type neutron
       real(R8):: number_n
       real(R8):: nposition
       real(R8):: angle_n
       real(R8):: move_n
   end type neutron

type(neutron) fbank(1:int(n0*1.5)), fbanknew(1:int(n0*1.5)),keep(1:int(n0*1.5))

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
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
!-----------------------------------------------------Actual keff-------------------------------------------------------
!Assumes a bare slab
    kinf=nubar*sigf/siga
    D=1.0d+0/(sigt*3.00)!*(1.0d+0-2.00/(3.00*nint(wght))))
    L2=D/siga
    PI=4.D0*DATAN(1.D0)
    Clength=PI*(L2/(kinf-1d+0))**0.5d+0-2.0d+0*2.1312d+0*D
    kact=(nubar*sigf)/(siga+D*(PI/(length+2.1312*2.0*D))**2.0d+0)
    
!////////////////////////////////////////////   Initializes components   /////////////////////////////////////////////////////  

    bin0= loc0/binwidth
    n=n0
    locn= loc0
    do i=1,int(n0*1.5)
        fbank(i)%number_n= 0.0
        fbanknew(i)%number_n=0.0
        fbank(i)%nposition= 0.0
        fbanknew(i)%nposition=fbank(i)%nposition
        fbank(i)%angle_n=0
        fbanknew(i)%angle_n=fbank(i)%angle_n


    end do
    do i=1,int(n0)
        fbank(i)%number_n= 1.0
        fbanknew(i)%number_n=1.0
    end do
    
    run=0
    write(*,'(/)')


!////////////////////////////////////////  The beginning of the successive runs   ///////////////////////////////////////////

    call RN_init_problem( 1234567_I8, 1 )

do while (run<nrun)                                 !ensures the number of iterations occur
    xpos(1:bins)=0.0
    xpos_ratio=20.0d+0
    kpathrat=10.5d+0
    kpathold = 200d+0
    Hsold=10.0
	
    do while (abs(kpathrat)>error )!.or. xpos_ratio>error)                !ensures keff and flux is converged before the next run
      leak= 0; absorb=0
      mult=0; coll=0
      i=1
      call RN_init_particle( int(i,I8) )


      do while(sum(fbank%number_n)>=1)
        if (i>n) i=1;                   !index remornalization
            
        
        if (fbank(i)%number_n<1.00) then                        !No particle present
            i=i+1
        else                                        !particle exist and moves
            fbank(i)%move_n=-log(rang())/sigt;                 !how far the particle moves
            if (abs(fbank(i)%angle_n)==0.0) then
                omega=sign(1.0d+0,2.00*rang()-1.00) 
                fbank(i)%angle_n=(2.0d+0*rang()-1.0d+0)                        !direction of movement
            end if


!==========================================   Movement   ===================================================================






!========================================= Leakage out left or right ========================================================

            if (length<fbank(i)%nposition+fbank(i)%move_n*fbank(i)%angle_n) then              ! right
                leak=leak+1
                fbank(i)%number_n=fbank(i)%number_n-1.0
                fbanknew(i)%number_n=fbanknew(i)%number_n-1.0
                fbanknew(i)%angle_n=0.0
                
            else if(fbank(i)%nposition+fbank(i)%move_n*fbank(i)%angle_n<0) then               ! left
                leak=leak+1
                fbank(i)%number_n=fbank(i)%number_n-1.0
                fbanknew(i)%number_n=fbanknew(i)%number_n-1.0
                fbanknew(i)%angle_n=0.0
                 
            else

!---------------------------------------Determine the type of reaction-----------------------------------------------

                rn=rang();
                if (rn<sigc/sigt) then               !capture
                    absorb =absorb+1
                    fbank(i)%number_n=fbank(i)%number_n-1.0
                    fbanknew(i)%number_n=fbanknew(i)%number_n-1.0
                     fbanknew(i)%angle_n=0.0

                else if (rn<(sigc+sigf)/sigt) then   !fission
                    rx=rang();
                    absorb=absorb+1;
                    k=0
                    do while (rx>fis_dist(k))
                        k=k+1
                    end do
                    mult=mult+k  
                    do while(k>1)
                        j=1
                        do while (fbanknew(j)%number_n>0.0)
                            j=j+1
                        end do   
                        fbanknew(j)%number_n=1.0
                        fbanknew(j)%nposition=fbank(i)%nposition+fbank(i)%move_n*fbank(i)%angle_n
                        fbanknew(j)%angle_n=(2.0d+0*rang()-1.0d+0)
                        k=k-1
                    end do 
                    fbank(i)%number_n=fbank(i)%number_n-1.0
                    fbanknew(i)%number_n=fbanknew(i)%number_n-1.0
                    fbanknew(i)%angle_n=0.0
                                                   
                        
                else                                  !scatter
                        
                    fbanknew(i)%nposition=fbank(i)%nposition+fbank(i)%move_n*fbank(i)%angle_n 
                    fbank(i)%nposition=fbank(i)%nposition+fbank(i)%move_n*fbank(i)%angle_n 
                    fbanknew(i)%angle_n=(2.0d+0*rang()-1.0d+0)
                    fbank(i)%angle_n=(2.0d+0*rang()-1.0d+0)
                    coll=coll+1
                    
                end if

!--------------------------------------------------------------------------------------------------------------------------

            end if             !movement
        end if                 !particle present
      end do                   !sum(fbank)=0

!-------------------------------------------------Calculated keff---------------------------------------------------------      

       kpath=sum(fbanknew%number_n)/n
       kpathrat=(kpath/kpathold-1.0)*(kpath/kpathold)**run
       kpathold=kpath 

!-------------------------------------------------Fission Bank Renormalization---------------------------------------------

         wtot=sum(fbanknew%number_n)
         wcum=0.0
         k=0
!print *, wtot
         do j=1,int(n0*1.5)
             prob=fbanknew(j)%number_n*real(n-k)/(wtot-wcum)
             wcum=wcum + fbanknew(j)%number_n
             knt =prob + rang()
             do i=1,int(knt)
                 k=k+1;
                 keep(k)=fbanknew(j)
                 
             end do
         end do

         n=sum(keep%number_n)



!------------------------------------------------------------------------------------------------------------------------
         k=1
         do while (k<int(n0*1.5))
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

!print *, 1.00+real(mult-absorb-leak)/n, n

    end do          !Keff


    if (run>=nrun_throw) then
     ! xposnew(0:bins-1)=(xpos(0:bins-1)/maxval(xpos)+xposnew(0:bins-1)*(run+1.00-nrun_throw))/(run+2.00-nrun_throw)
      keff=(keff*(run+1.00-nrun_throw)+kpath)/(run+2.00-nrun_throw)
      !nmax(run-nrun_throw)= kpath
      
    else
     ! xposnew(0:bins-1)=xpos(0:bins-1)/maxval(xpos)
      keff=kpath
    end if

if (run-nrun_throw==0) print *,  "------------------Start Count----------------------------"
if (mod(run,int(nrun/5))==0) write(*,'(10x,A,I3,5x,A,F6.4)') 'run ',run,'keff=', keff

    run= run+1


    end do    !runs






print *, "kact= ", kact, Clength
print *, "keff= ", keff, "dkeff[$]=", abs(keff-kact)/.0065
write(*,'(/)')





end program rod
