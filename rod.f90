program rod 

  use mcnp_random 
  implicit none

  integer, parameter :: I8 = selected_int_kind(18)
  integer, parameter :: R8 = selected_real_kind(15,307)
  real(R8), parameter :: n0 = 100000d+0
  integer, parameter :: nrun = 200
  integer, parameter :: nrun_throw = 50
  integer, parameter :: run_act=nrun-nrun_throw
  integer, parameter :: sites_bin=10
  real(R8), parameter :: length =4.0d+0
  integer, parameter :: bins= 101
  real(R8), parameter :: binwidth=length/(bins-1) 
  real(R8), parameter :: loc0 = 2d+0
  real(R8), parameter :: nubar = 2.414d+0
  real(R8), parameter :: error = .001
  integer, parameter :: source=100
  real(R8), parameter :: source_pos=2d+0
  integer  :: i, absorb, leak,j,ii
  real(R8) :: total, locn, rn, direction, rx, den, Navg, wght, kinf, L2,D, kpath,kpathold, pos, xpos_ratio, kact, lmove
  integer :: bin0, n,run, particle, nleft, coll, inscat, mult,k
  real(R8) :: sigc, sigs,sigf, siga, sigt, move, omega, PI, Clength, keff, kpathrat, errork,nsig,navgn, pj, Hs, Hsold
  real(R8), dimension (0:bins-1) ::xpos, xposnew, flux, xposmid,fbank,fbanknew, wgh
  real(R8), dimension (0:7) :: fis_dist
  real(R8), dimension (0:((bins-1)/sites_bin)) :: shannon
  real(R8) :: start, finish, wtot, wcum, prob, knt
  real(R8), dimension (0:nrun-nrun_throw) :: nmax

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
 
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
    D=1.0d+0/(3.0*sigt)!(sigs*3.00*(1.0d+0-2.00/(3.00*int(wght))))
    L2=D/siga
    PI=4.D0*DATAN(1.D0)
    Clength=PI*(L2/(kinf-1d+0))**0.5d+0-2.0d+0*2.1312d+0*D
    kact=(nubar*sigf)/(siga+D*(PI/(length+2.1312*2.0*D))**2.0d+0)
    
!////////////////////////////////////////////   Initializes components   /////////////////////////////////////////////////////  

  bin0= loc0/binwidth
  n=n0
  locn= loc0
  fbank(0:bins-1)= floor(real(n0/bins))
  fbanknew(0:bins-1)=fbank(0:bins-1)
  xposmid(0:bins-1)=0
  run=0
  nmax(0:nrun-nrun_throw)=0

print *, "kact=", kact, "Crit Length=", Clength 

!////////////////////////////////////////  The beginning of the successive runs   ////////////////////////////////////////////

call RN_init_problem( 1234567_I8, 1 )

do while (run<nrun)                                 !ensures the number of iterations occur
    xpos(0:bins-1)=0.0
    xpos_ratio=20.0d+0
    kpathrat=100.5d+0
    kpathold = 200d+0
    Hsold=100.0
    call cpu_time(start)
	
    do while (abs(kpathrat)>error .or. xpos_ratio>error)                !ensures keff and flux is converged before the next run
      leak= 0; absorb=0
      mult=0; coll=0
      xpos(0:bins-1)=0.0
      i=maxval(maxloc(fbank))
      call RN_init_particle( int(i,I8) )
      do while(sum(fbank)>0)
        if (i>bins-1) then                          !index remornalization
            i=0
        end if

        if (fbank(i)<1) then                        !No particle present
            i=i+1
        else                                        !particle exist and moves
            locn=i*binwidth
            move=-log(rang())/sigt;                 !how far the particle moves
            omega=sign(1.0d+0,2.00*rang()-1.00) 
            move=move*(2.0d+0*rang()-1d+0)                        !direction of movement

!========================================= Leakage out left or right ========================================================

            if (length<locn+move) then              ! right
                leak=leak+1
                fbank(i)=fbank(i)-1
                fbanknew(i)=fbanknew(i)-1
                xpos(i:bins-1)=xpos(i:bins-1)+binwidth
            else if(locn+move<0) then               ! left
                leak=leak+1
                fbank(i)=fbank(i)-1
                fbanknew(i)=fbanknew(i)-1
                xpos(0:i)=xpos(0:i)+binwidth 
            else

! ================================  track the particles movement for flux calculation  ========================================
 
                lmove=(locn+move)/binwidth
                if (move>0) then
                    xpos(i:floor(lmove))=xpos(i:floor(lmove))+binwidth
                    ii=-1
                    if (lmove-floor(lmove)>rang()) then
                        xpos(ceiling(lmove))= xpos(ceiling(lmove))+binwidth
                        ii=1
                    end if
                else if (move==0) then
                    xpos(i)=xpos(i)
                    ii=0
                else
                    xpos(ceiling(lmove):i)=xpos(ceiling(lmove):i)+binwidth
                    ii=1
                    if (lmove-floor(lmove)<rang()) then
                        xpos(floor(lmove))= xpos(floor(lmove))+binwidth
                        ii=-1
                    end if
                end if

!---------------------------------------Determine the type of reaction-----------------------------------------------

                rn=rang();
                if (rn<sigc/sigt) then               !capture
                    absorb =absorb+1
                    fbank(i)=fbank(i)-1
                    fbanknew(i)=fbanknew(i)-1

                else if (rn<(sigc+sigf)/sigt) then   !fission
                    rx=rang();
                    absorb=absorb+1;
                    k=0
                    do while (rx>fis_dist(k))
                        k=k+1
                    end do
                        fbank(i)=fbank(i)-1
                        fbanknew(i)=fbanknew(i)-1
                        if (ii>0) then   
                            fbanknew(ceiling(lmove))=fbanknew(ceiling(lmove))+k  
                        else if (ii<0) then   
                            fbanknew(floor(lmove))=fbanknew(floor(lmove))+k  
                        else
                            fbanknew(i)=fbanknew(i)+k 
                        end if 
                        mult=mult+k             
                        
                else                                  !scatter
                        if (ii>0) then   
                            fbanknew(ceiling(lmove))=fbanknew(ceiling(lmove))+1
                            fbank(ceiling(lmove))=fbank(ceiling(lmove))+1  
                        else if (ii<0) then   
                            fbanknew(floor(lmove))=fbanknew(floor(lmove))+1
                            fbank(floor(lmove))=fbank(floor(lmove))+1  
                        else
                            fbanknew(i)=fbanknew(i)+1
                            fbank(i)=fbank(i)+1 
                        end if 
                    fbank(i)=fbank(i)-1
                    fbanknew(i)=fbanknew(i)-1
                    coll=coll+1 
                end if
!--------------------------------------------------------------------------------------------------------------------------
            end if             !movement
        end if                 !particle present
      end do                   !sum(fbank)=0

!-------------------------------------------------Calculated keff---------------------------------------------------------      

       kpath=sum(fbanknew)/n
       kpathrat=(kpath/kpathold-1.0)*(kpath/kpathold)**run
       kpathold=kpath 
       

!----------------------------------------This Renormalizes the fission bank-----------------------------------------------

  fbanknew=(fbanknew/kpath)  
  do i=0,bins-1
    if (fbanknew(i)-floor(fbanknew(i))>rang()) then
      fbanknew(i)=ceiling(fbanknew(i))
    else
      fbanknew(i)=floor(fbanknew(i))
    end if
  end do
  n=sum(fbanknew)
    fbank(0:bins-1)=fbanknew(0:bins-1)

!print *, n
	
!-------------------------------------------Shannon Entropy--------------------------------------------------------------
	do i=0,int(bins-1)/sites_bin,1
	  pj=0
	  do j=int(sites_bin*i),int(sites_bin-1+sites_bin*i)
	    pj=pj+xpos(j)/sum(xpos)*log(xpos(j)/sum(xpos))/log(2.0)
	  end do
	  shannon(i)=pj
	end do
	Hs=-sum(shannon)
	xpos_ratio=abs(Hs/Hsold-1.00)
	Hsold=Hs


end do                         !Keff converged
!print *, xpos_ratio
!---------------------------------------xpos/keff averaging --------------------------------------------- 
    
    if (run>=nrun_throw) then
      xposnew(0:bins-1)=(xpos(0:bins-1)/maxval(xpos)+xposnew(0:bins-1)*(run+1.00-nrun_throw))/(run+2.00-nrun_throw)
      keff=(keff*(run+1.00-nrun_throw)+kpath)/(run+2.00-nrun_throw)
      nmax(run-nrun_throw)= keff
      
    else
      xposnew(0:bins-1)=xpos(0:bins-1)/maxval(xpos)
      keff=kpath
    end if

!print *, "keff=", keff, " for run", run+1
  run=run+1
  call cpu_time(finish)
!  print *, finish-start  
end do                                  !ends the number of runs

!---------------------------------------------  d error  ----------------------------------------------------------
   navgn=sum(nmax)/(run-nrun_throw-1)
   nsig=(1.00/(run-nrun_throw-2.00)*1/(run-nrun_throw-1.00)*sum((nmax(0:nrun-nrun_throw-1)**2.00-navgn**2.00)))
   nsig=sqrt(abs(nsig))

!---------------------------------------------------flux calc--------------------------------------------------------
i=1
flux(0:bins-1)=0.0
pos=0.0
flux(0)=dsin(PI*(pos+2.1312*D)/(length+2.00*2.1312*D))
do while(i<bins)
    pos=pos+binwidth
    flux(i)=dsin(PI*(pos+2.1312*D)/(length+2.00*2.1312*D))
    i=i+1
end do

!----------------------------------------------------output-------------------------------------------------------------
 print *, "keff=", keff, "dkeff[$]=",abs(kact-keff)/.0065
 print *, maxloc(xposnew), maxloc(flux), "   Calc_dkeff[$]=", nsig*keff/.0065
 OPEN(UNIT=12, FILE="rod.txt", ACTION="write", STATUS="replace")
 write(12,*), xposnew/maxval(xposnew)
 write(12,*), flux 
 close(12)

end program rod
