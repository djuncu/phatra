module subs

  contains

!----------------------------------------------------------------  
!laplacian------------------
!!!!!!!!!CHANGE IN v0.9.6: BOUNDARY CONDITIONS TOP AND BOTTOM
!!!!! AND ADDED BCTYPE
!-----------------
	function laplace(a,d,lptype,topvalue,bottomvalue)
	implicit none
	real,intent(in) 		    :: d
	real,optional,intent(in)	    :: topvalue,bottomvalue
	real,dimension(:,:),intent(in) 	    :: a
	real,dimension(size(a,1),size(a,2)) :: laplace
	integer 			    :: i,j
	character(len=*)		    :: lptype
	
    
	! boundary conditions top and bottom 0:
	! a(top+0.5*d,:)=-a(top-0.5*d) 
	! a(-0.5*d,:) = -a(0.5*d)


	select case (trim(lptype))

				
	case('lp_vy') ! Boundary conditions for vy: top/bottom vy = 0	
	

		do j=1,size(a,2)
		do i=1,size(a,1)
						
			if(i==1.and.j==1) then 
			laplace(i,j) = (a(i+1,j) + a(i,size(a,2)) + a(i,j+1) - 4.*a(i,j)) / d**2.
			elseif(i==1.and.j==size(a,2)) then
			laplace(i,j) = (a(i+1,j) + a(i,j-1) + a(i,1) - 4.*a(i,j)) / d**2.
			elseif(i==1) then
			laplace (i,j) = (a(i+1,j) + a(i,j-1) + a(i,j+1) - 4.*a(i,j)) / d**2.  
			elseif(i==size(a,1).and.j==1) then
			laplace(i,j) = (a(i-1,j) + a(i,size(a,2)) + a(i,j+1) - 4.*a(i,j)) / d**2.   
			elseif(i==size(a,1).and.j==size(a,2)) then
			laplace(i,j) = (a(i-1,j) + a(i,j-1) + a(i,1) - 4.*a(i,j)) / d**2.	
			elseif(i==size(a,1)) then
			laplace(i,j) = (a(i-1,j) + a(i,j-1) + a(i,j+1) - 4.*a(i,j)) / d**2.
			elseif(j==1) then
			laplace(i,j) = (a(i-1,j) + a(i+1,j) + a(i,size(a,2)) + a(i,j+1) - 4.*a(i,j)) / d**2.
			elseif(j==size(a,2)) then
			laplace(i,j) = (a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,1) - 4.*a(i,j)) / d**2.
			else
			laplace(i,j) = (a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1) - 4.*a(i,j)) / d**2.   
			endif	
				
		enddo
		enddo

	case('lp_T') ! temperature
	
				
		do j=1,size(a,2)
		do i=1,size(a,1)

			if(present(topvalue)) then
			if(i==1.and.j==1) then 
			laplace(i,j) = (a(i+1,j) + a(i,size(a,2)) + a(i,j+1) + 2.*bottomvalue - 5.*a(i,j)) / d**2.
			elseif(i==1.and.j==size(a,2)) then
			laplace(i,j) = (a(i+1,j) + a(i,j-1) + a(i,1) + 2.*bottomvalue - 5.*a(i,j)) / d**2.
			elseif(i==1) then
			laplace(i,j) = (a(i+1,j) + a(i,j-1) + a(i,j+1) + 2.*bottomvalue - 5.*a(i,j)) / d**2.  
			elseif(i==size(a,1).and.j==1) then
			laplace(i,j) = (a(i-1,j) + a(i,size(a,2)) + a(i,j+1) + 2.*topvalue - 5.*a(i,j)) / d**2.   
			elseif(i==size(a,1).and.j==size(a,2)) then
			laplace(i,j) = (a(i-1,j) + a(i,j-1) + a(i,1) + 2.*topvalue - 5.*a(i,j)) / d**2.	
			elseif(i==size(a,1)) then
			laplace(i,j) = (a(i-1,j) + a(i,j-1) + a(i,j+1) + 2.*topvalue - 5.*a(i,j)) / d**2.
			elseif(j==1) then
			laplace(i,j)=(a(i-1,j)+a(i+1,j)+a(i,size(a,2))+a(i,j+1)-4.*a(i,j)) / d**2.
			elseif(j==size(a,2)) then
			laplace(i,j)=(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,1)-4.*a(i,j)) / d**2.
			else
			laplace(i,j)=(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4.*a(i,j)) / d**2.   
			endif
				
			else
			print*, 'specify top and bottom values in laplace'
			endif 
		
		enddo
		enddo

	case('lp_vx')	!boundary conditions for vx			
		
		
		do j=1,size(a,2)
		do i=1,size(a,1)

			if(present(topvalue)) then
			if(i==1.and.j==1) then 
			laplace(i,j) = (a(i+1,j) + a(i,size(a,2)-1) + a(i,j+1) + 2.*bottomvalue - 5.*a(i,j)) / d**2.
			elseif(i==1.and.j==size(a,2)) then
			laplace(i,j) = (a(i+1,j) + a(i,j-1) + a(i,2) + 2.*bottomvalue - 5.*a(i,j)) / d**2.
			elseif(i==1) then
			laplace(i,j) = (a(i+1,j) + a(i,j-1) + a(i,j+1) + 2.*bottomvalue - 5.*a(i,j)) / d**2.  
			elseif(i==size(a,1).and.j==1) then
			laplace(i,j)=(a(i-1,j) + a(i,size(a,2)-1) + a(i,j+1) + 2.*topvalue - 5.*a(i,j)) / d**2.   
			elseif(i==size(a,1).and.j==size(a,2)) then
			laplace(i,j) = (a(i-1,j) + a(i,j-1) + a(i,2) + 2.*topvalue - 5.*a(i,j)) / d**2.	
			elseif(i==size(a,1)) then
			laplace(i,j) = (a(i-1,j) + a(i,j-1) + a(i,j+1) + 2.*topvalue - 5.*a(i,j)) / d**2.
			elseif(j==1) then
			laplace(i,j)=(a(i-1,j)+a(i+1,j)+a(i,size(a,2)-1)+a(i,j+1)-4.*a(i,j)) / d**2.
			elseif(j==size(a,2)) then
			laplace(i,j)=(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,2)-4.*a(i,j)) / d**2.
			else
			laplace(i,j)=(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4.*a(i,j)) / d**2.        
			endif
				
			else
			print*, 'specify top and bottom values in laplace'
			endif 
		
		enddo
		enddo

	case DEFAULT
	print*,'wrong input in laplace'
	stop
	end select

		
 
	end function laplace

!----------------------------------------------------------------	
!----------------------------------------------------------------
!first derivative in x-direction  (1st order)
!---------------------------------------
	function ddx(a,d,dtype)
	implicit none
	real,intent(in)		        :: d
	real,dimension(:,:),intent(in)  :: a
	real,dimension(:,:),allocatable :: ddx
	integer			        :: i,j
	character(len=*), intent(in)   :: dtype	
	
	select case (trim(dtype))

	case ('center_to_face') ! differences between cell centers, at cell faces
	
		 !  (P)---|dP/dx|---(P)
		 ! ||-face; ()-center	

		allocate(ddx(size(a,1),size(a,2)+1))

 		do j=1,size(a,2)+1
		do i=1,size(a,1)
			
			if(j==1.or.j==size(a,2)+1) then
			ddx(i,j) = (a(i,1) - a(i,size(a,2))) / d   
			else
			ddx(i,j) = (a(i,j) - a(i,j-1)) / d    
			endif
		
		enddo
		enddo

	case ('face_to_center') ! differences between cell faces, at cell centers
	
		 ! |vx|---(dvx/dx)---|vx|
		 ! ||-face; ()-center

		allocate(ddx(size(a,1),size(a,2)-1))

		do j=1,size(a,2)-1
		do i=1,size(a,1)
			
			if(j==size(a,2)-1) then
			ddx(i,j) = (a(i,1) - a(i,j)) / d   
			else
			ddx(i,j) = (a(i,j+1) - a(i,j)) / d    
			endif
		
		enddo
		enddo

	case ('fc_to_fc') ! differences between locations of same type: face/face, center/center
	

		allocate(ddx(size(a,1),size(a,2)))

		do j=1,size(a,2)
		do i=1,size(a,1)

			if(j.eq.1.or.j.eq.size(a,2)) then
			ddx(i,j) = (a(i,2) - a(i,size(a,2)-1)) / (2. * d)
			else
			ddx(i,j) = (a(i,j+1) - a(i,j-1)) / (2. * d)
			endif
		
		enddo
		enddo	

	end select	

	end function ddx


!----------------------------------------------------------------	
!----------------------------------------------------------------
!first derivative in y-direction, 1st order
!---------------------------------------
	function ddy(a,d,dtype,top,bottom)
	implicit none
	real,intent(in)			:: d
	real,optional,intent(in)	:: top,bottom
	real,dimension(:,:),intent(in)  :: a
	real,dimension(:,:),allocatable :: ddy
	integer				:: i,j
	character(len=*), intent(in)   :: dtype
  
	select case (trim(dtype))

	
	case ('center_to_face') ! differences between cell centers, at cell faces
	
		! -> y-axis (P)---|dP/dy|---(P)
		! ||-face; ()-center

		allocate(ddy(size(a,1)+1,size(a,2)))
		
				
 		do j=1,size(a,2)
		do i=1,size(a,1)+1
			
			if(i==1.or.i==size(a,1)+1) then
			ddy(i,j) = 0.
			else
			ddy(i,j) = (a(i,j) - a(i-1,j)) / d     
			endif
			
		enddo
		enddo

				

		

	case ('face_to_center') ! differences between cell faces, at cell centers  
	
		! -> y-axis |vy|---(dvy/dy)---|vy|
		! ||-face; ()-center

		allocate(ddy(size(a,1)-1,size(a,2)))

		do j=1,size(a,2)
		do i=1,size(a,1)-1
			ddy(i,j) = (a(i+1,j) - a(i,j)) / d     
		enddo
		enddo

	case ('fc_to_fc') ! differences between locations of same type: face/face, center/center
	
		 !!! CHECK BOUNDARY CONDITIONS

		allocate(ddy(size(a,1),size(a,2)))

		do j=1,size(a,2)
		do i=1,size(a,1)

			if(i.eq.1.or.i.eq.2) then
			ddy(i,j) = (a(i+1,j) - bottom) / (2. * d)
			elseif(i.eq.size(a,1).or.i.eq.size(a,1)-1) then
			ddy(i,j) = (top - a(i-1,j)) / (2. * d)
			else
			ddy(i,j) = (a(i+1,j) - a(i-1,j)) / (2. * d)
			endif
		
		enddo
		enddo	

	end select

	end function ddy



!-------------------------------------------------
! interpolation between the different locations on the staggered grid
! bilinear from face to face,  linear from center to face
!-------------------------------------------
subroutine interpolate(Fin,Fout,direction,top,bottom)
	implicit none
	real,intent(in),dimension(:,:) :: Fin
	real,intent(out),dimension(:,:),allocatable :: Fout
	real,optional,intent(in) :: top,bottom
	integer::i,j,nxa,nya
	character(len=*),intent(in) :: direction
	
	nya=size(Fout,1)
	nxa=size(Fout,2) 	
	
	select case (trim(direction))

	case ('vy_to_xpos')  !interpolate from vy to vx positions
	

		allocate(Fout(size(Fin,1)-1,size(Fin,2)+1))

		do j=1,nxa
		do i=1,nya
	
			if(j.eq.1.or.j.eq.nxa) then
				Fout(i,j) = ( Fin(i,1) + Fin(i+1,1) + Fin(i,nxa-1) + Fin(i+1,nxa-1) ) / 4.
			else
				Fout(i,j) = ( Fin(i,j-1) + Fin(i+1,j-1) + Fin(i,j) + Fin(i+1,j) ) / 4.	
			endif
	
		enddo
		enddo	


	case ('vx_to_ypos') !interpolate from vx to vy positions
	

		allocate(Fout(size(Fin,1)+1,size(Fin,2)-1))

		do j=1,nxa
		do i=1,nya

			if(i.eq.1) then
				Fout(i,j) = bottom
			elseif(i.eq.nya) then	
				Fout(i,j) = top	
			else
				Fout(i,j) = ( Fin(i,j) + Fin(i,j+1) + Fin(i-1,j+1) + Fin(i-1,j) ) / 4.	
			endif

		enddo
		enddo	

	case ('center_to_ypos') ! interpolate from cell-centers to vy positions
	
		allocate(Fout(size(Fin,1)+1,size(Fin,2)))		

		do j=1,nxa
		do i=1,nya
			
			if(i.eq.1) then
				Fout(i,j) = bottom
			elseif(i.eq.nya) then
				Fout(i,j) = top
			else
				Fout(i,j) = ( Fin(i,j) + Fin(i-1,j) ) / 2.
			endif 
		enddo
		enddo	

	case ('center_to_corner') !interpolate from cell-centers to cell-corners
	

		allocate(Fout(size(Fin,1)+1,size(Fin,2)+1))

		nya=size(Fout,1); nxa=size(Fout,2) 

		
		do j=1,nxa
		do i=1,nya

		if(i.eq.1.and.j.eq.1) then
			Fout(i,j) =  Fin(i,j) + Fin(i,nxa-1)
		elseif(i.eq.1.and.j.eq.nxa) then
			Fout(i,j) =  Fin(i,1) + Fin(i,j-1)
		elseif(i.eq.nya.and.j.eq.nxa) then
			Fout(i,j) = Fin(i-1,j-1) + Fin(i-1,1)
		elseif(i.eq.nya.and.j.eq.1) then
			Fout(i,j) = Fin(i-1,nxa-1) + Fin(i-1,j)
		elseif(i.eq.nya) then
			Fout(i,j) = Fin(i-1,j-1) + Fin(i-1,j)
		elseif(i.eq.1) then
			Fout(i,j) = Fin(i,j) + Fin(i,j-1)
		elseif(j.eq.1.or.j.eq.nxa) then
			Fout(i,j) = (Fin(i-1,nxa-1) + Fin(i-1,1) + Fin(i,1) + Fin(i,nxa-1)) / 4.
		else
			Fout(i,j) = (Fin(i-1,j-1) + Fin(i-1,j) + Fin(i,j) + Fin(i,j-1)) / 4.
		endif

		enddo
		enddo

	
	end select			
		
	
    	end subroutine interpolate

!---------------------------------------------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------
!advection-of vx--v*grad(F)
!3rd order accurate upwind scheme----------------------
!!!!!CHANGE IN v0.9.6 TOP AND BOTTOM BOUNDARY CONDITIONS
!---------------------------------
	subroutine advection(d,vx,vy,nxa,nya,F,vgrad,atype,topvalue,bottomvalue)
	implicit none
	real,intent(in)		      	    :: d
	real,intent(in),optional	    :: topvalue,bottomvalue	
	integer,intent(in)	     	    :: nxa,nya
	real,intent(in),dimension(:,:)      :: vx,vy,F
	real,dimension(nya,nxa) 	    :: xcomp,ycomp 
	real,dimension(:,:),allocatable,intent(out) :: vgrad
	integer				    :: i,j
	character(len=*),intent(in)	    :: atype
 	
	! boundary conditions top and bottom are no-slip: 
	! vx(top+0.5*d,:)=-vx(top-0.5*d) 
	! vx(-0.5*d,:) = -vx(0.5*d)

	select case (trim(atype))

	
	case('advect_vx') ! advect vx
	

	allocate(vgrad(nya,nxa))

	do j=1,nxa
	do i=1,nya

		
		if(vx(i,j).ge.0) then
		
		if(j.eq.1) then
		xcomp(i,j)=vx(i,j)*(2.*F(i,j+1)+3.*F(i,j)-6.*F(i,nxa-1)+F(i,nxa-2)) 
		elseif(j.eq.2) then
		xcomp(i,j)=vx(i,j)*(2.*F(i,j+1)+3.*F(i,j)-6.*F(i,j-1)+F(i,nxa-1))  
		elseif(j.eq.nxa) then
		xcomp(i,j)=vx(i,j)*(2.*F(i,2)+3.*F(i,j)-6.*F(i,j-1)+F(i,j-2))
		else
		xcomp(i,j)=vx(i,j)*(2.*F(i,j+1)+3.*F(i,j)-6.*F(i,j-1)+F(i,j-2))
		endif
			
		else
		
		if(j.eq.nxa-1) then
		xcomp(i,j)=vx(i,j)*(-1.*F(i,2)+6.*F(i,j+1)-3.*F(i,j)-2.*F(i,j-1))
		elseif(j.eq.nxa) then
		xcomp(i,j)=vx(i,j)*(-1.*F(i,3)+6.*F(i,2)-3.*F(i,j)-2.*F(i,j-1))
		elseif(j.eq.1) then
		xcomp(i,j)=vx(i,j)*(-1.*F(i,j+2)+6.*F(i,j+1)-3.*F(i,j)-2.*F(i,nxa-1))
		else
		xcomp(i,j)=vx(i,j)*(-1.*F(i,j+2)+6.*F(i,j+1)-3.*F(i,j)-2.*F(i,j-1))
		endif
		
		endif

		if(vy(i,j).ge.0) then

		if(i.eq.2) then
		ycomp(i,j)=vy(i,j)*(2.*F(i+1,j)+3.*F(i,j)-6.*F(i-1,j) - F(i-1,j))
		elseif(i.eq.1) then
		ycomp(i,j)=vy(i,j)*(2.*F(i+1,j) + 3.*F(i,j) + 6.*F(i,j) - 3.*F(i,j))
		elseif(i.eq.nya) then
		ycomp(i,j)=vy(i,j)*(-2.*F(i,j) + 3.*F(i,j) - 6.*F(i-1,j) + F(i-2,j))
		else
		ycomp(i,j)=vy(i,j)*(2.*F(i+1,j) + 3.*F(i,j) - 6.*F(i-1,j) + F(i-2,j))
		endif

		else
		
		if(i.eq.nya-1) then
		ycomp(i,j) = vy(i,j)*(1.*F(i+1,j) + 6.*F(i+1,j) - 3.*F(i,j) - 2.*F(i-1,j))
		elseif(i.eq.nya) then
		ycomp(i,j) = vy(i,j)*(3.*F(i,j) - 6.*F(i,j) - 3.*F(i,j) - 2.*F(i-1,j))
		elseif(i.eq.1) then
		ycomp(i,j)=vy(i,j)*(-1.*F(i+2,j) + 6.*F(i+1,j) - 3.*F(i,j) + 2.*F(i,j))
		else
		ycomp(i,j)=vy(i,j)*(-1.*F(i+2,j) + 6.*F(i+1,j) - 3.*F(i,j) - 2.*F(i-1,j))
		endif

		endif

	enddo
	enddo
	
	
	case('advect_vy') ! advect vy
	

	allocate(vgrad(nya,nxa))

	do j=1,nxa
	do i=1,nya

		if(vx(i,j).ge.0) then
		
		if(j.eq.1) then
		xcomp(i,j)=vx(i,j)*(2.*F(i,j+1)+3.*F(i,j)-6.*F(i,nxa)+F(i,nxa-1)) 
		elseif(j.eq.2) then
		xcomp(i,j)=vx(i,j)*(2.*F(i,j+1)+3.*F(i,j)-6.*F(i,j-1)+F(i,nxa))  
		elseif(j.eq.nxa) then
		xcomp(i,j)=vx(i,j)*(2.*F(i,1)+3.*F(i,j)-6.*F(i,j-1)+F(i,j-2))
		else
		xcomp(i,j)=vx(i,j)*(2.*F(i,j+1)+3.*F(i,j)-6.*F(i,j-1)+F(i,j-2))
		endif
			
		else
		
		if(j.eq.nxa-1) then
		xcomp(i,j)=vx(i,j)*(-1.*F(i,1)+6.*F(i,j+1)-3.*F(i,j)-2.*F(i,j-1))
		elseif(j.eq.nxa) then
		xcomp(i,j)=vx(i,j)*(-1.*F(i,2)+6.*F(i,1)-3.*F(i,j)-2.*F(i,j-1))
		elseif(j.eq.1) then
		xcomp(i,j)=vx(i,j)*(-1.*F(i,j+2)+6.*F(i,j+1)-3.*F(i,j)-2.*F(i,nxa))
		else
		xcomp(i,j)=vx(i,j)*(-1.*F(i,j+2)+6.*F(i,j+1)-3.*F(i,j)-2.*F(i,j-1))
		endif
		
		endif

		if(vy(i,j).ge.0) then

		if(i.eq.2) then
		ycomp(i,j) = vy(i,j)*(2.*F(i+1,j) + 3.*F(i,j) - 5.*F(i-1,j))
		elseif(i.eq.1) then
		ycomp(i,j) = vy(i,j)*(2.*F(i+1,j) - 2.*F(i,j))
		elseif(i.eq.nya) then
		ycomp(i,j) = vy(i,j)*(5.*F(i,j) - 6.*F(i-1,j) + F(i-2,j))
		else
		ycomp(i,j) = vy(i,j)*(2.*F(i+1,j) + 3.*F(i,j) - 6.*F(i-1,j) + F(i-2,j))
		endif

		else
		
		if(i.eq.nya-1) then
		ycomp(i,j)=vy(i,j)*(5.*F(i+1,j)-3.*F(i,j)-2.*F(i-1,j))
		elseif(i.eq.nya) then
		ycomp(i,j)=vy(i,j)*(2.*F(i,j)-2.*F(i-1,j))
		elseif(i.eq.1) then
		ycomp(i,j)=vy(i,j)*(-1.*F(i+2,j)+6.*F(i+1,j)-5.*F(i,j))
		else
		ycomp(i,j)=vy(i,j)*(-1.*F(i+2,j)+6.*F(i+1,j)-3.*F(i,j)-2.*F(i-1,j))
		endif

		endif

	enddo
	enddo

	end select

	vgrad = (xcomp+ycomp) / (6.*d)
    	end subroutine advection  
!----------------------------------------------------------------
!----------------------------------------------------------------
!donor cell upwind scheme------
!---------------------------------
	subroutine donorcell(d,dt,vx,vy,nx,ny,F,donor)
	implicit none
	real,intent(in)::d,dt
	integer,intent(in)::nx,ny
	real,intent(in),dimension(:,:)::F,vx,vy
	real,dimension(ny,nx),intent(out)::donor
	real,dimension(ny,nx+1)::fluxh
	real,dimension(ny+1,nx)::fluxv
	integer::i,j
 
		
	!cell flux computation	

	!horizontal flux
	forall(j=2:nx+1,i=1:ny,vx(i,j).ge.0) fluxh(i,j)=F(i,j-1)*vx(i,j)
	forall(j=1:nx,i=1:ny,vx(i,j).lt.0) fluxh(i,j)=F(i,j)*vx(i,j)
	
	
	forall(i=1:ny,vx(i,1).ge.0) fluxh(i,1)=F(i,nx)*vx(i,1) !left boundary
	forall(i=1:ny,vx(i,nx+1).lt.0) fluxh(i,nx+1)=F(i,1)*vx(i,nx+1) !right boundary
	
	!vertical flux
	forall(i=2:ny,j=1:nx,vy(i,j).ge.0) fluxv(i,j)=F(i-1,j)*vy(i,j)
	forall(i=2:ny,j=1:nx,vy(i,j).lt.0) fluxv(i,j)=F(i,j)*vy(i,j)

	fluxv(1,:)=0.	!bottom boundary
	fluxv(ny+1,:)=0. !top boundary
    	
	
	forall(j=1:nx,i=1:ny) donor(i,j)=(dt/d)*(fluxh(i,j)+fluxv(i,j)-fluxh(i,j+1)-fluxv(i+1,j))
	

	end subroutine donorcell  
!--------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
! diffusive term with variable material properties
!----------------------------------------------------

	subroutine diffusion(a,b,p,q,d,dfield,dtype,vxtop,vxbot)
	implicit none
	real,intent(in) 	      	 :: d
	real,optional,intent(in) 	 :: vxtop,vxbot ! boundary values for vx 
	real,dimension(:,:),intent(in) 	 :: a,b,p,q ! a,b: primary/secondary velocity, p,q: viscosity at cell centres/corners
	real,dimension(:,:),allocatable,intent(out) :: dfield
	real 				 :: visc_top,visc_bot ! boundary values for viscosity
	integer 			 :: i,j,nxa,nya,nyb,nxb
	character(len=*),intent(in)    	 :: dtype

	nya = size(a,1); nxa = size(a,2)
	nyb = size(b,1); nxb = size(b,2)

	allocate(dfield(nya,nxa))

	visc_top = q(size(q,1),1)
	visc_bot = q(1,1)

	select case (trim(dtype)) ! choose if diffusion for vx or vy equation

	
	! boundary conditions top and bottom are no-slip: 
	! vx(top+0.5*d,:)=-vx(top-0.5*d) 
	! vx(-0.5*d,:) = -vx(0.5*d)
	

	case ('diffuse_vx') ! diffusion for vx equation
	

	do j=1,nxa
	do i=1,nya

	if(i.eq.1.and.j.eq.1) then
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i,j+1) - a(i,j)) - p(i,nxb) * (a(i,j) - a(i,nxa-1))) & 
			     + (1./d**2.) * (q(i+1,j) * (a(i+1,j) - a(i,j) + b(i+1,j) - b(i+1,nxb)) &
					  - q(i,j) * (a(i,j) + a(i,j) + b(i,j) - b(i,nxb)))
	elseif(i.eq.1.and.j.eq.nxa) then
		dfield(i,j) = (2./d**2.) * (p(i,1) * (a(i,2) - a(i,j)) - p(i,j-1) * (a(i,j) - a(i,j-1))) & 
			     + (1./d**2.) * (q(i+1,j) * (a(i+1,j) - a(i,j) + b(i+1,1) - b(i+1,j-1)) &
					  - q(i,j) * (a(i,j) - a(i,j) + b(i,1) - b(i,j-1)))
	elseif(i.eq.nya.and.j.eq.nxa) then
		dfield(i,j) = (2./d**2.) * (p(i,1) * (a(i,2) - a(i,j)) - p(i,j-1) * (a(i,j) - a(i,j-1))) & 
			     + (1./d**2.) * (q(i+1,j) * (-1.*a(i,j) - a(i,j) + b(i+1,1) - b(i+1,j-1)) &
					  - q(i,j) * (a(i,j) - a(i-1,j) + b(i,1) - b(i,j-1)))
	elseif(i.eq.nya.and.j.eq.1) then
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i,j+1) - a(i,j)) - p(i,nxb) * (a(i,j) - a(i,nxa-1))) & 
			     + (1./d**2.) * (q(i+1,j) * (-1.*a(i,j) - a(i,j) + b(i+1,j) - b(i+1,nxb)) &
					  - q(i,j) * (a(i,j) - a(i-1,j) + b(i,j) - b(i,nxb)))
	elseif(i.eq.1) then
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i,j+1) - a(i,j)) - p(i,j-1) * (a(i,j) - a(i,j-1))) &
			     + (1./d**2.) * (q(i+1,j) * (a(i+1,j) - a(i,j) + b(i+1,j) - b(i+1,j-1)) &
					  - q(i,j) * (a(i,j) + a(i,j) + b(i,j) - b(i,j-1)))
	elseif(i.eq.nya) then		
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i,j+1) - a(i,j)) - p(i,j-1) * (a(i,j) - a(i,j-1))) &
			     + (1./d**2.) * (q(i+1,j) * (a(i,j) - a(i,j) + b(i+1,j) - b(i+1,j-1)) &
					  - q(i,j) * (a(i,j) - a(i-1,j) + b(i,j) - b(i,j-1)))
	elseif(j.eq.1) then
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i,j+1) - a(i,j)) - p(i,nxb) * (a(i,j) - a(i,nxa-1))) &
			     + (1./d**2.) * (q(i+1,j) * (a(i+1,j) - a(i,j) + b(i+1,j) - b(i+1,nxb)) &
					  - q(i,j) * (a(i,j) - a(i-1,j) + b(i,j) - b(i,nxb)))
	elseif(j.eq.nxa) then
		dfield(i,j) = (2./d**2.) * (p(i,1) * (a(i,2) - a(i,j)) - p(i,j-1) * (a(i,j) - a(i,j-1))) &
			     +(1./d**2.) * (q(i+1,j) * (a(i+1,j) - a(i,j) + b(i+1,1) - b(i+1,j-1)) &
					  - q(i,j) * (a(i,j) - a(i-1,j) + b(i,1) - b(i,j-1)))
	else
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i,j+1) - a(i,j)) - p(i,j-1) * (a(i,j) - a(i,j-1))) &
			     + (1./d**2.) * (q(i+1,j) * (a(i+1,j) - a(i,j) + b(i+1,j) - b(i+1,j-1)) &
					  - q(i,j) * (a(i,j) - a(i-1,j) + b(i,j) - b(i,j-1)))
	endif

	enddo
	enddo

	case ('diffuse_vy') ! diffusion for vy equation

	
	do j=1,nxa
	do i=1,nya

	if(i.eq.1.and.j.eq.1) then
		dfield(i,j) = (2./d**2.) * (p(i,j)*(a(i+1,j) - a(i,j)) - visc_bot*(a(i,j) - 0.)) &
			     + (1./d**2.) * (q(i,j+1)*(b(i,j+1) + b(i,j+1) + a(i,j+1) - a(i,j)) &
					  - q(i,j)*(b(i,j) + b(i,j) + a(i,j) - a(i,nxa)))
	elseif(i.eq.1.and.j.eq.nxa) then
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i+1,j) - a(i,j)) - visc_bot * (a(i,j) - 0.)) &
			     + (1./d**2.) * (q(i,j+1) * (b(i,j+1) + b(i,j+1) + a(i,1) - a(i,j)) &
					  - q(i,j) * (b(i,j) + b(i,j) + a(i,j) - a(i,j-1)))
	elseif(i.eq.nya.and.j.eq.nxa) then
		dfield(i,j) = (2./d**2.) * (visc_top * (0. - a(i,j)) - p(i-1,j) * (a(i,j) - a(i-1,j))) &
			     + (1./d**2.) * (q(i,j+1) * (-1.* b(i-1,j+1) - b(i-1,j+1) + a(i,1) - a(i,j)) &
					  - q(i,j) * (-1.* b(i-1,j) - b(i-1,j) + a(i,j) - a(i,j-1)))
	elseif(i.eq.nya.and.j.eq.1) then
		dfield(i,j) = (2./d**2.) * (visc_top * (0. - a(i,j)) - p(i-1,j) * (a(i,j) - a(i-1,j))) &
			     + (1./d**2.) * (q(i,j+1) * (-1.* b(i-1,j+1) - b(i-1,j+1) + a(i,j+1) - a(i,j)) &
					  - q(i,j) * (-1.* b(i-1,j) - b(i-1,j) + a(i,j) - a(i,nxa)))
	elseif(i.eq.1) then
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i+1,j) - a(i,j)) - visc_bot * (a(i,j) - 0.)) &
			     + (1./d**2.) * (q(i,j+1) * (b(i,j+1) + b(i,j+1) + a(i,j+1) - a(i,j)) &
					  - q(i,j) * (b(i,j) + b(i,j) + a(i,j) - a(i,j-1)))
	elseif(i.eq.nya) then
		dfield(i,j) = (2./d**2.) * (visc_top * (0. - a(i,j)) - p(i-1,j) * (a(i,j) - a(i-1,j))) &
			     + (1./d**2.) * (q(i,j+1) * (-1.* b(i-1,j+1) - b(i-1,j+1) + a(i,j+1) - a(i,j)) &
					  - q(i,j) * (-1.* b(i-1,j) - b(i-1,j) + a(i,j) - a(i,j-1)))
	elseif(j.eq.1) then
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i+1,j) - a(i,j)) - p(i-1,j) * (a(i,j) - a(i-1,j))) &
			     + (1./d**2.) * (q(i,j+1) * (b(i,j+1) - b(i-1,j+1) + a(i,j+1) - a(i,j)) &
					  - q(i,j) * (b(i,j) - b(i-1,j) + a(i,j) - a(i,nxa)))
	elseif(j.eq.nxa) then
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i+1,j) - a(i,j)) - p(i-1,j) * (a(i,j) - a(i-1,j))) &
			     + (1./d**2.) * (q(i,j+1) * (b(i,j+1) - b(i-1,j+1) + a(i,1) - a(i,j)) &
					  - q(i,j) * (b(i,j) - b(i-1,j) + a(i,j) - a(i,j-1)))
	else
		dfield(i,j) = (2./d**2.) * (p(i,j) * (a(i+1,j) - a(i,j)) - p(i-1,j) * (a(i,j) - a(i-1,j))) &
			     + (1./d**2.) * (q(i,j+1) * (b(i,j+1) - b(i-1,j+1) + a(i,j+1) - a(i,j)) &
					  - q(i,j) * (b(i,j) - b(i-1,j) + a(i,j) - a(i,j-1)))
	endif
	
	enddo
	enddo
			
	end select

	end subroutine diffusion
!---------------------------------------------------
!freeze the flow
!------------------
	

	subroutine freeze(fl,vx,vy)
	implicit none
	real,dimension(:,:),intent(in)	  :: fl
	real,dimension(:,:),intent(inout) :: vx,vy 
	integer 			  :: i,j,nx,ny

	nx = size(fl,2)
	ny = size(fl,1)

	
	!set velocity to 0 in cells where liquid fraction < 0.5
	! if this holds as well for one of their neighbouring cells (or if it is adjacent to a wall)
	do j=1,nx
	do i=1,ny
		if(fl(i,j).lt.0.5) then
			if(i.eq.1.or.i.eq.ny) then
				vx(i,j)   = 0.
				vx(i,j+1) = 0.	
				vy(i,j)   = 0.
				vy(i+1,j) = 0.

				if(j.eq.1) then
					vx(i,nx+1) = 0.	
				elseif(j.eq.nx) then	
					vx(i,1)    = 0.					
				endif

			elseif(j.eq.1) then

				if(fl(i+1,j).lt.0.5.or.fl(i+1,j+1).lt.0.5.or.fl(i,j+1).lt.0.5.or.fl(i-1,j+1).lt.0.5.or. &
				fl(i-1,j).lt.0.5.or.fl(i-1,nx).lt.0.5.or.fl(i,nx).lt.0.5.or.fl(i+1,nx).lt.0.5) then
		
				vx(i,1)    = 0.
				vx(i,2)    = 0.	
				vy(i,1)    = 0.
				vy(i+1,1)  = 0.
				vx(i,nx+1) = 0.
				
				endif
			
			elseif(j.eq.nx) then
				
				if(fl(i+1,j).lt.0.5.or.fl(i+1,1).lt.0.5.or.fl(i,1).lt.0.5.or.fl(i-1,1).lt.0.5.or. &
				fl(i-1,j).lt.0.5.or.fl(i-1,j-1).lt.0.5.or.fl(i,j-1).lt.0.5.or.fl(i+1,j-1).lt.0.5) then

				vx(i,nx)   = 0.
				vx(i,nx+1) = 0.	
				vy(i,nx)   = 0.
				vy(i+1,nx) = 0.
				vx(i,1)    = 0.				
	
				endif

			else
				if(fl(i+1,j).lt.0.5.or.fl(i+1,j+1).lt.0.5.or.fl(i,j+1).lt.0.5.or.fl(i-1,j+1).lt.0.5.or. &
				fl(i-1,j).lt.0.5.or.fl(i-1,j-1).lt.0.5.or.fl(i,j-1).lt.0.5.or.fl(i+1,j-1).lt.0.5) then

				vx(i,j)   = 0.
				vx(i,j+1) = 0.	
				vy(i,j)   = 0.
				vy(i+1,j) = 0.
				
				endif
			endif
		endif
	enddo
	enddo

	end subroutine freeze
!----------------------------------------------------------------
!----------------------------------------------------------------
!poisson iteration------------
!-laplacian(a)=b--------------------------

	function iteration_2DPoisson(a,b,d,alpha)
	implicit none
	real, intent(in)::d,alpha
	real,dimension(:,:),intent(inout)::a
	real,dimension(:,:),intent(in)::b
	real,dimension(:,:),allocatable::residue
	real::iteration_2DPoisson
	integer::i,j,nxa,nya
  
    	nya=size(a,1); nxa=size(a,2)  
	allocate(residue(nya,nxa))
	 
   	
	call residue_2DPoisson(a,b,d,residue)

   	!compute rms of the Residues
   	iteration_2DPoisson = sqrt(sum(residue*residue)/real(nxa*nya))   
  	 
  	!update u-field
  	 a = a + 0.25 * alpha * d**2. * residue
     
 	 end function iteration_2DPoisson

!residues------------
!---------------------------
	
	  subroutine residue_2DPoisson(a,b,d,res_f)
	  implicit none
	  real, intent(in)::d
	  real,dimension(:,:),intent(in)::a,b
	  real,dimension(:,:),intent(out)::res_f
	  integer::i,j,nxa,nya
  
	nya=size(res_f,1); nxa=size(res_f,2) 	

		!compute residues
		do j=1,nxa
		do i=1,nya
			if(i==1.and.j==1) then
			res_f(i,j)=((a(i+1,j)+a(i,nxa)+a(i,j+1)-3.*a(i,j))/d**2.)-b(i,j)
			elseif(i==1.and.j==nxa) then
			res_f(i,j)=((a(i+1,j)+a(i,j-1)+a(i,1)-3.*a(i,j))/d**2.)-b(i,j)
			elseif(i==1) then
			res_f(i,j)=((a(i+1,j)+a(i,j-1)+a(i,j+1)-3.*a(i,j))/d**2.)-b(i,j)
			elseif(i==nya.and.j==1) then
			res_f(i,j)=((a(i-1,j)+a(i,nxa)+a(i,j+1)-3.*a(i,j))/d**2.)-b(i,j)
			elseif(i==nya.and.j==nxa) then
			res_f(i,j)=((a(i-1,j)+a(i,j-1)+a(i,1)-3.*a(i,j))/d**2.)-b(i,j)
			elseif(i==nya) then
			res_f(i,j)=((a(i-1,j)+a(i,j-1)+a(i,j+1)-3.*a(i,j))/d**2.)-b(i,j)
			elseif(j==1) then
			res_f(i,j)=((a(i-1,j)+a(i+1,j)+a(i,nxa)+a(i,j+1)-4.*a(i,j))/d**2.)-b(i,j)
			elseif(j==nxa) then
			res_f(i,j)=((a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,1)-4.*a(i,j))/d**2.)-b(i,j)
			else
			res_f(i,j)=((a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4.*a(i,j))/d**2.)-b(i,j)
			endif
		enddo
		enddo
 	 end subroutine residue_2DPoisson
!-----------------------
!-----------------------
!GRID-RESTRICTION-------	
	
 	subroutine restrict(a,b) !a:fine grid, b:coarse grid
	   implicit none
	   real,dimension(:,:) :: a,b
	   integer :: i,j
   
   
	   do j = 1,size(b,2)
	   do i = 1,size(b,1)
   
		   b(i,j) = ( a(2*i-1,2*j-1) + a(2*i-1,2*j) + a(2*i,2*j-1) + a(2*i,2*j) ) / 4.
   
	   enddo
	   enddo
   
	 end subroutine restrict
!-----------------------
!-----------------------
!GRID-PROLONGATION-------   

	   subroutine prolongate(b,a) !a:fine grid, b:coarse grid
	    implicit none
	     real,dimension(:,:)::a,b
	     integer::i,j
   
   
	   !prolongate
	   do j=1,size(b,2)
	   do i=1,size(b,1)
   	
		    a(2*i-1,2*j-1) = b(i,j) 
   		    a(2*i-1,2*j)   = b(i,j)
		    a(2*i,2*j-1)   = b(i,j)
  		    a(2*i,2*j)     = b(i,j)
	   enddo
	   enddo
    
       
	  !forall(i=1:size(a,1)) a(i,1)=a(i,size(a,2))

   	end subroutine prolongate
!-----------------------
!-----------------------
!Vcycle-(written by Paul Tackley)------   

  recursive function Vcycle_2DPoisson(u_f,rhs,h) result (resV)
    implicit none
    real resV
    real,intent(inout):: u_f(:,:)  ! arguments
    real,intent(in)   :: rhs(:,:),h
    integer         :: nx,ny,nxc,nyc, i,j  ! local variables
    real,allocatable:: res_c(:,:),corr_c(:,:),res_f(:,:),corr_f(:,:)
    real            :: alpha=0.7, res_rms

    nx=size(u_f,2); ny=size(u_f,1)  ! must be power of 2 plus 1
    nxc=nx/2.; nyc=ny/2.  ! coarse grid size

    if (min(nx,ny)>4) then  ! not the coarsest level

       allocate(res_f(ny,nx),corr_f(ny,nx), &
            corr_c(nyc,nxc),res_c(nyc,nxc))

       !---------- take 2 iterations on the fine grid--------------
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha) 
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)

       !---------- restrict the residue to the coarse grid --------
       call residue_2DPoisson(u_f,rhs,h,res_f) 
       call restrict(res_f,res_c)

       !---------- solve for the coarse grid correction -----------
       corr_c = 0.  
       res_rms = Vcycle_2DPoisson(corr_c,res_c,2*h) ! *RECURSIVE CALL*

       !---- prolongate (interpolate) the correction to the fine grid 
       call prolongate(corr_c,corr_f)

       !---------- correct the fine-grid solution -----------------
       u_f = u_f - corr_f  

       !---------- two more smoothing iterations on the fine grid---
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)

       deallocate(res_f,corr_f,res_c,corr_c)

    else  

       !----- coarsest level (ny=4): iterate to get 'exact' solution
		
       do i = 1,100
          res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       end do

    end if

    resV = res_rms   ! returns the rms. residue

  end function Vcycle_2DPoisson
!--------------------------------------
!--------------------------------------

	SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
          END SUBROUTINE


!--------------------
!rounding!!
!!-------------------

	subroutine customround(a)
        
        implicit none
        real, dimension(:,:),intent(inout) :: a
        integer,dimension(size(a,1),size(a,2))  :: i
	integer :: k

	do k=4,2,-1
	
        i(:,:) = ANINT(a(:,:)*(10.**k)) !Multiply by 10^k and round to int
        a(:,:)  = i(:,:) / dble(10.**k)     !Floating point divide by 10^k

  	enddo

        end subroutine customround

!----------------------------------------------------------------
!----------------------------------------------------------------
 end module subs
!----------------------------------------------------------------
!----------------------------------------------------------------


!----------------------------------------------------------------
!----------------------------------------------------------------
!------------------------------------------------------------------------
!-Advection & phase transition------using pressure-velocity formulation-- staggered grid
!---------------------------------------------------------------------

program phatra
use subs
	
	implicit none
	integer :: i,j,nx,ny,fieldsaver,iteration,timestep,totalsaver,datasaver
	real,dimension(:,:),allocatable :: T,vgradvx,vgradvy,dTdx,k,phi,p,rhs,vx,vy,donor,pcorr
	real,dimension(:,:),allocatable :: vxguess,vyguess,vy_xpos,vx_ypos,T_ypos,force,buoyancy,diffusevx,diffusevy
	real,dimension(:,:),allocatable :: nu_corner,nu
	real :: d,rhs_rms,alpha,res_rms,dt,dt_dif,dt_adv,time,Re,Pr,runtime,eps,Tm,St,start,finish,lambda,Ra
	real :: Ttop,Tbot, vxbot, vxtop,phi_int,rn,starttime,forcefactor,beta,nondim_factor,T_init
	real :: vy_int,kinetic,savefields,savedata,a_dif_vel,a_adv_vel,a_dif_visc,a_adv_visc,v_ini
	character(len=9)  ::  freezetype,flowtype,output_type
	real*8,parameter ::pi = 3.14159265358979323846
	character(len=40) :: filename, fd2 = '(2X,F20.7))' 
	character(len=10) :: fd1
	character(len=40) :: fd	
	character(len=4)  :: savenumber
	logical :: loaddata,first_step

	namelist /inputs/ nx,ny,alpha,runtime,a_dif_vel,a_adv_vel,a_dif_visc,a_adv_visc,beta, & 	
				  Re,Pr,eps,Tm,St,Ra,loaddata,forcefactor, &
			          freezetype,flowtype,T_init,savefields,savedata,v_ini,output_type

	open(1,file='parameters.txt',status='old')
	read(1,inputs)
	close(1)
	
	call cpu_time(start)

	!set format descriptor for saving files
	write(fd1,'(I4)') nx
	fd = '('//trim(adjustl(fd1))//fd2
		
	
	allocate(T(ny,nx))
	allocate(P(ny,nx))
	allocate(Pcorr(ny,nx))
	allocate(phi(ny,nx))
    	allocate(k(ny,nx))
	allocate(rhs(ny,nx))
	allocate(donor(ny,nx))

	allocate(T_ypos(ny+1,nx))
	allocate(buoyancy(ny+1,nx))
	allocate(diffusevy(ny+1,nx))
	allocate(vx_ypos(ny+1,nx))
	allocate(vyguess(ny+1,nx))	
	allocate(vy(ny+1,nx))
	allocate(vgradvy(ny+1,nx))

	allocate(vx(ny,nx+1))
	allocate(diffusevx(ny,nx+1))
	allocate(force(ny,nx+1))
	allocate(vy_xpos(ny,nx+1))
	allocate(vxguess(ny,nx+1))
	allocate(vgradvx(ny,nx+1))
	
	select case(trim(freezetype))
	case('VISCOSITY')
		allocate(nu(ny,nx))
		allocate(nu_corner(ny+1,nx+1))
	end select



	!----grid----
	!------------
	!--------vy----------
	! -      -	-
	! -	 -	-	
	!-vx-----P,T*------vx---------------*also phi and k
	! -      -	-
	! -	 -	-	
	!--------vy---------
	
	d = 1. / real(ny)
	
	!Top and bottom boundary values
	Ttop = 0.
	Tbot = 1.
	vxbot = 0.
	vxtop = 0.
		

	!do i=1,ny
	!force(i,:) = (1.-4.*((i/real(ny))-0.5)**2.)*forcefactor
	!enddo
	
	force=forcefactor

	select case(trim(flowtype))
	case('RAYLEIGH')
		nondim_factor = Pr
	case('CHANNEL')
		nondim_factor = 1./Re
	end select

	!-------------------------------------------
	! Initializations
	!-------------------------------------

	if(loaddata) then

	!load input data from file
	!--------

		!load temperature data
		filename = 'dump/final/T_final.dat'
		open(11,file=filename)
		do i=1,ny
		read(11,*) (T(i,j), j=1,nx)
		enddo
		close(11)

		!load velocity data
		filename = 'dump/final/vx_final.dat'
		open(12,file=filename)
		do i=1,ny
		read(12,*) (vx(i,j), j=1,nx)
		enddo
		close(12)

		filename = 'dump/final/vy_final.dat'
		open(13,file=filename)
		do i=1,ny
		read(13,*) (vy(i,j), j=1,nx)
		enddo
		close(13)

		!load pressure
		filename = 'dump/final/P_final.dat'
		open(14,file=filename)
		do i=1,ny
		read(14,*) (P(i,j), j=1,nx)
		enddo
		close(14)
		
		!load time
		filename = 'dump/final/information.dat'	
		open(unit=15,form='formatted',status='old',action='read',position='rewind',file=filename)
		read(15,*) time,totalsaver,timestep
		close(15)	

	
	else

		time       	= 0. 	
		totalsaver      = 1
		timestep   	= 0
		
	
		! Pressure
		!---------------------------------
		P(:,:)     	= 1.
			
		! Velocity
		!--------------------------	
		!call init_random_seed() !resets the random fields for each compilation
	
		!random velocity field
		call random_number(vx)	
		call random_number(vy)
	
		vx = (vx - 0.5)*v_ini
		vy = (vy - 0.5)*v_ini

		! Temperature 
		!----------------------
		
		T(:,:) = T_init

		
	endif	

	


	!----------------------------------
	!velocity BC
	vx(1,:)    = (2.*vxbot + vx(2,:))  / 3.
	vx(ny,:)   = (2.*vxtop + vx(ny-1,:)) / 3.
	vx(:,1)    = vx(:,nx+1) 
	vy(ny+1,:) = 0.
	vy(1,:)    = 0.
	
  	!boundary conditions for temperature
	T(1,:)  = (2.*Tbot + T(2,:))  / 3.
	T(ny,:) = (2.*Ttop + T(ny-1,:))  / 3.
	

	phi = 0.5 * ( 1. - tanh( ( Tm-T ) / eps ) ) ! liquid fraction
	phi_int = sum(phi)/real(nx*ny)	
	call interpolate(vx,vx_ypos,'vx_to_ypos',vxtop,vxbot)


!----------------------	
!time loop------------
!----------------------

fieldsaver = 1
datasaver  = 1
starttime = time
first_step = .TRUE.

do while (time < (starttime + runtime))


	phi = 0.5 * ( 1. - tanh( ( Tm-T ) / eps ) ) ! liquid fraction
	k = 1. / (1. + ( 1. / St ) / ( 2. * ( cosh( (Tm-T) / eps ) )**2. )) !apparent conductivity
	
	select case(trim(freezetype))
	case('VISCOSITY')
		nu = 1. * exp(beta*(1. - phi))  ! viscosity
	end select
	

	select case(trim(freezetype))
	case ('VELOCITY')
		call freeze(phi,vx,vy)

		! compute timestep
		dt_adv = a_adv_vel * min( d / maxval(abs(vx(:,:))) , d / maxval(abs(vy(:,:))) )
		dt_dif = a_dif_vel * d**2. / max(Pr,maxval(k))
		dt = min(dt_dif,dt_adv)
	
		! compute the terms for momentum diffusion 
		! -----------------------------------------------------------
		diffusevx = laplace(vx,d,'lp_vx',vxtop,vxbot)
		diffusevy = laplace(vy,d,'lp_vy')		

		

	case('VISCOSITY')

		! compute timstep
		dt_adv = a_adv_visc * min( d / maxval(abs(vx(:,:))) , d / maxval(abs(vy(:,:))) )
		dt_dif = a_dif_visc * d**2. / max(Pr*maxval(nu),maxval(k))
		dt = min(dt_dif,dt_adv)


		! compute the terms for momentum diffusion 
		! -----------------------------------------------------------
		call interpolate(nu,nu_corner,'center_to_corner')

		call diffusion(vx,vy,nu,nu_corner,d,diffusevx,'diffuse_vx',vxtop,vxbot)
		call diffusion(vy,vx,nu,nu_corner,d,diffusevy,'diffuse_vy',vxtop,vxbot)

	
	end select


	! compute the other terms needed to solve momentum and energy equation
	!---------------------------------------------------------------
	
	call interpolate(vy,vy_xpos,'vy_to_xpos')
	call interpolate(vx,vx_ypos,'vx_to_ypos',vxtop,vxbot)
	call interpolate(T,T_ypos,'center_to_ypos',Ttop,Tbot)
	
	call advection(d,vx,vy_xpos,nx+1,ny,vx,vgradvx,'advect_vx',vxtop,vxbot) 
	call advection(d,vx_ypos,vy,nx,ny+1,vy,vgradvy,'advect_vy') 

	select case(trim(flowtype))
	case('RAYLEIGH')
		buoyancy = Pr * Ra * T_ypos
	case('CHANNEL')
		buoyancy = 0.
	end select

	call donorcell(d,dt,vx,vy,nx,ny,T,donor) ! call donorcell before updating velocities (is this correct?)


	!save integrated phi,vy, kinetic energy and effective rayleigh number at every timestep
	!------------------------------------------------------------
	phi_int = sum(phi)/real(nx*ny)
	kinetic = sum(0.5*(vx_ypos+vy)*(vx_ypos+vy))/real(nx*(ny+1))

        if(savedata.eq.0.or.runtime/(time-starttime+1e-15).lt.savedata/real(datasaver)) then
	if(first_step) then
		open(16,file='data.dat',status='replace',position='append')
	else
		open(16,file='data.dat',status='old',position='append')
        endif

	write(16,'(F14.7,F14.7,F14.7)') time,phi_int,kinetic	
	close(16)
	
	datasaver = datasaver + 1
        first_step = .FALSE.
	endif
	!------------------------------------
	!compute the pressure
	!by guessing velocities using P from the previous timestep
	!-----------------------------------
	
	!1. guess velocities ----------------------
	!------------------------------------------

	vxguess = vx + dt * ( nondim_factor * diffusevx - vgradvx -  ddx(P,d,'center_to_face') + force )
	vyguess = vy + dt * ( nondim_factor * diffusevy - vgradvy -  ddy(P,d,'center_to_face') + buoyancy )

		
	!------------------------------------
	!boundary conditions 
	vxguess(:,1)    = vxguess(:,nx+1)
	vxguess(1,:)    = (2.*vxbot + vx(2,:))  / 3.
	vxguess(ny,:)   = (2.*vxtop + vx(ny-1,:)) / 3.	
	vyguess(ny+1,:) = 0.
	vyguess(1,:)    = 0.

	select case (trim(freezetype))
	case('VELOCITY')
		call freeze(phi,vxguess,vyguess)
	end select

	
	!2. laplacian(P_corr)= div(vguess)/dt
	!(P_new = P_old + P_corr)
	!-------------------------------------

	!compute div(vguess)/dt as the rhs of pressure poisson equation
	rhs = ( ddx(vxguess,d,'face_to_center') + ddy(vyguess,d,'face_to_center') ) / (dt)
	
	!compute residues of the rhs field
  	rhs_rms = sqrt( sum(rhs*rhs) / real(nx*ny) )
	res_rms = rhs_rms	
	
	Pcorr(:,:) = 0.
	iteration = 0

	!---------------------------
	!-solve poisson's equation for P_corr (iteratively, using the multigrid solver)
	!--------------------------------------------------

 
 	do while(res_rms/rhs_rms.ge.1.e-3)
		res_rms = Vcycle_2DPoisson(Pcorr,rhs,d)
		!print*, res_rms/rhs_rms
		
		if(res_rms/rhs_rms.eq.res_rms/rhs_rms+1.0) then
		print*, 'scheiss in vcycle'
		stop
		endif

		if(iteration.eq.100000) then		
		print*, 'no convergence'
		stop
		endif		

		iteration=iteration+1
	enddo
     	
	
	!update the pressure field
	P = P + Pcorr
	
	
	!------------------------------------
	!velocity
	!-----------------------------------
	
	!compute new velocity field
	vx = vx + dt * ( nondim_factor * diffusevx - vgradvx -  ddx(P,d,'center_to_face') + force )
	vy = vy + dt * ( nondim_factor * diffusevy - vgradvy -  ddy(P,d,'center_to_face') + buoyancy )

		

	!------------------------------------
	!boundary conditions 
	vx(:,1)    = vx(:,nx+1) 
	vx(1,:)    = (2.*vxbot + vx(2,:))  / 3.
	vx(ny,:)   = (2.*vxtop + vx(ny-1,:)) / 3.
	vy(ny+1,:) = 0.
	vy(1,:)    = 0.

	
	!------------------------------------
	!Temperature
	!-----------------------------------
	
	!compute new temperature field
	T = T + dt * nondim_factor * k *( laplace(T,d,'lp_T',Ttop,Tbot) / Pr ) + donor 

                 
      	!------------------------------------
	!boundary conditions  (have to be changed at some point)
	T(1,:)  = (2.*Tbot + T(2,:))  / 3.
	T(ny,:) = (2.*Ttop + T(ny-1,:))  / 3.
	
   	!save data
	!----------
	if(runtime/(time-starttime+1e-15).lt.savefields/real(fieldsaver)) then
	
		!naming of files
		if (totalsaver.le.9) then
		write(savenumber,'(A3,I1)') '000', totalsaver
	
		elseif(totalsaver.le.99) then
		write(savenumber,'(A2,I2)') '00', totalsaver

		elseif(totalsaver.le.999) then
		write(savenumber,'(A1,I3)') '0', totalsaver

		elseif(totalsaver.le.99) then
		write(savenumber,'(I4)')  totalsaver

		endif
	
		!save temperature data
		filename = 'dump/temperature'//savenumber//'.dat'
		open(3,file=filename)
		do i=1,ny
		write(3,fd) (T(i,j), j=1,nx)
		enddo
		close(3)

		!save velocity data
		filename = 'dump/vx'//savenumber//'.dat'
		open(4,file=filename)
		do i=1,ny
		write(4,fd) (vx(i,j), j=1,nx)
		enddo
		close(4)

		filename = 'dump/vy'//savenumber//'.dat'
		open(5,file=filename)
		do i=1,ny
		write(5,fd) (vy(i,j), j=1,nx)
		enddo
		close(5)

		!save pressure
		filename = 'dump/pressure'//savenumber//'.dat'
		open(7,file=filename)
		do i=1,ny
		write(7,fd) (P(i,j), j=1,nx)
		enddo
		close(7)
		
		!save liquid fraction
		filename = 'dump/phi'//savenumber//'.dat'
		open(9,file=filename)
		do i=1,ny
		write(9,fd) (phi(i,j), j=1,nx)
		enddo
		close(9)
		
		!save k
		!filename = 'dump/ca'//savenumber//'.dat'
		!open(10,file=filename)
		!do i=1,ny
		!write(10,fd) (k(i,j), j=1,nx)
		!enddo
		!close(10)

		
		!at every save-step: output to screen is file number and the accordant simulation time
		write(*,*)  'File',fieldsaver,'saved at time:', time
	
		!save counter
		fieldsaver      = fieldsaver+1
		totalsaver = totalsaver+1
		
	endif


		
	
	!timestep
	time=time+dt
	timestep=timestep+1
		
	select case(trim(output_type))
	case('LONG')
		print*, 'time:',time
		print*, 'dt_dif:', dt_dif
		print*, 'dt_adv:', dt_adv
		print*, '-----------'
	end select
	
enddo

	
!----------------------
!--save final data----
!---------------------

		call interpolate(vx,vx_ypos,'vx_to_ypos',vxtop,vxbot)

		!save integrated phi,vy, kinetic energy and effective rayleigh number at every timestep
		phi_int = sum(phi)/real(nx*ny)
		kinetic = sum(0.5*(vx_ypos+vy)*(vx_ypos+vy))/real(nx*(ny+1))

		open(16,file='data.dat',status='old',position='append')
		write(16,'(F14.7,F14.7,F14.7)') time,phi_int,kinetic	
		close(16)
		

		select case (trim(freezetype))
		case('VELOCITY')
			call freeze(phi,vxguess,vyguess)
		end select		

				!save temperature data
		filename = 'dump/final/T_final.dat'
		open(11,file=filename,status='unknown')
		do i=1,ny
		write(11,fd) (T(i,j), j=1,nx)
		enddo
		close(11)

		!save velocity data
		filename = 'dump/final/vx_final.dat'
		open(12,file=filename)
		do i=1,ny
		write(12,fd) (vx(i,j), j=1,nx)
		enddo
		close(12)

		filename = 'dump/final/vy_final.dat'
		open(13,file=filename)
		do i=1,ny
		write(13,fd) (vy(i,j), j=1,nx)
		enddo
		close(13)

		!save pressure
		filename = 'dump/final/P_final.dat'
		open(14,file=filename)
		do i=1,ny
		write(14,fd) (P(i,j), j=1,nx)
		enddo
		close(14)
		
		!save time
		filename = 'dump/final/information.dat'	
		open(15,file=filename)
		write(15,'(F14.7,I4,I9)') time,totalsaver,timestep
		close(15)

      call cpu_time(finish)
              print '("Time = ",f12.3," seconds.")',finish-start

end program phatra
