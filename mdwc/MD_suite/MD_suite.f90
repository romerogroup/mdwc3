module npt_MD_suite

implicit none

contains

	subroutine cross_product(a,b,crossp)
	 !a very simple implementation of the cross product
		  implicit none
		  real(8),intent(in) ::a(3),b(3)
		  real(8),intent(out)::crossp(3)

		  crossp(1)=a(2)*b(3)-a(3)*b(2)
		  crossp(2)=a(3)*b(1)-a(1)*b(3)
		  crossp(3)=a(1)*b(2)-a(2)*b(1)
		  return
	end subroutine cross_product

	subroutine invertmat(mat,matinv,n)
		implicit none
	    integer,intent(in) :: n
		real(8),intent(in) :: mat(n,n)
		real(8),intent(out)   :: matinv(n,n)
		real(8)               :: a(n,n),div
		integer               :: IPIV(n), INFO
		!Here only for a 3*3 matrix
		a=mat
		div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+&
		&a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
		div=1.d0/div
		matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
		matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
		matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
		matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
		matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
		matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
		matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
		matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
		matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
		return
	end subroutine invertmat

	subroutine rotation(rotmat,angle,axe)
		!This subroutine will calculate the rotational matrix rotmat for a
		!3-dim vector around an axis 'axe' by the angle 'angle'.
		implicit none
		real(8),INTENT(IN) :: angle
		real(8),INTENT(IN) :: axe(3)
		real(8):: rotator(3,3)
		real(8):: rotmat(3,3)

		!Define Rotation Matrix
		rotator(1,1)=dcos(angle)+(axe(1)**2)*(1.d0-dcos(angle))
		rotator(1,2)=axe(1)*axe(2)*(1.d0-dcos(angle))-axe(3)*dsin(angle)
		rotator(1,3)=axe(1)*axe(3)*(1.d0-dcos(angle))+axe(2)*dsin(angle)

		rotator(2,1)=axe(2)*axe(1)*(1.d0-dcos(angle))+axe(3)*dsin(angle)
		rotator(2,2)=dcos(angle)+(axe(2)**2)*(1.d0-dcos(angle))
		rotator(2,3)=axe(2)*axe(3)*(1.d0-dcos(angle))-axe(1)*dsin(angle)

		rotator(3,1)=axe(3)*axe(1)*(1.d0-dcos(angle))-axe(2)*dsin(angle)
		rotator(3,2)=axe(3)*axe(2)*(1.d0-dcos(angle))+axe(1)*dsin(angle)
		rotator(3,3)=dcos(angle)+(axe(3)**2)*(1.d0-dcos(angle))
		rotmat(:,:)=rotator(:,:)

		!do i=1,3
		!   vector2(i)=rotator(i,1)*vector(1)+rotator(i,2)*vector(2)+rotator(i,3)*vector(3)
		!enddo
		!vector(:)=vector2(:)
	end subroutine rotation

	subroutine gausdist(nat,vel_in)
		!generates 3*nat random numbers distributed according to  exp(-.5*vel_in**2)
		implicit none
		integer, intent(in) :: nat
		real(8), dimension(3,nat), intent(out) :: vel_in
		!temporal storage
		real:: s1,s2
		real(8):: t1,t2,tt
		! On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
		real(8),parameter:: eps=1.d-8
		real(8),dimension(3*nat)::  vxyz
		integer:: i,j
		do i=1,3*nat-1,2
			call random_number(s1)
			t1=eps+(1.d0-2.d0*eps)*dble(s1)
			call random_number(s2)
			t2=dble(s2)
			tt=sqrt(-2.d0*log(t1))
			vxyz(i)=tt*cos(6.28318530717958648d0*t2)
			vxyz(i+1)=tt*sin(6.28318530717958648d0*t2)
		enddo
		call random_number(s1)
		t1=eps+(1.d0-2.d0*eps)*dble(s1)
		call random_number(s2)
		t2=dble(s2)
		tt=sqrt(-2.d0*log(t1))
		vxyz(3*nat)=tt*cos(6.28318530717958648d0*t2)
		!call elim_moment(nat,vel_in,amass)
		do j=1, nat
			do i=1, 3
				vel_in(i,j)= vxyz(3*(j-1) + i)
			end do
		end do
		return
	end subroutine gausdist

    subroutine gausdist_cell(vel_lat_in)
		! generates 3*3 random numbers distributed according to  exp(-.5*vxyz**2) for the cell vectors
		implicit none
		!real(8), dimension(3,3), intent(in) :: latvec_in
		real(8), dimension(3,3), intent(out) :: vel_lat_in
		!real(8), dimension(3,3), intent(out) :: vel_lat_in
		integer:: i,j
		real:: s1,s2
		real(8) :: t1,t2,tt
		! On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
		real(8),parameter:: eps=1.d-8
		real(8)::  vlat(9)

		do i=1,3*3-1,2
			call random_number(s1)
			t1=eps+(1.d0-2.d0*eps)*dble(s1)
			call random_number(s2)
			t2=dble(s2)
			tt=sqrt(-2.d0*log(t1))
			vlat(i)=tt*cos(6.28318530717958648d0*t2)
			vlat(i+1)=tt*sin(6.28318530717958648d0*t2)
		enddo
		call random_number(s1)
		t1=eps+(1.d0-2.d0*eps)*dble(s1)
		call random_number(s2)
		t2=dble(s2)
		tt=sqrt(-2.d0*log(t1))
		vlat(3*3)=tt*cos(6.28318530717958648d0*t2)
		do j=1, 3
			do i=1, 3
				vel_lat_in(i,j)= vlat(3*(j-1) + i)
			end do
		end do
		!call elim_torque_cell(latvec_in,vel_lat_in)
		return
    end subroutine gausdist_cell
    
    subroutine init_vel_atoms(amu,temp, nat, vel_out)!, vel_lat_in_bug)
		implicit none
		!*********************************************************************************************
		integer, intent(in) :: nat
		!integer, intent(in) :: nsoften
		real(8), intent(in) :: temp
		!real(8), intent(in) :: bmass
		real(8), intent(in) :: amu(nat)
		!real(8):: amass(nat)
		real(8),dimension(3,nat), intent(out) :: vel_out
		real(8),dimension(3,nat) :: vel_to_scale
		real(8),dimension(3,nat) :: vel_scaled
		real(8),dimension(3) :: p_total
		integer:: iat,i_dim
		real(8), parameter :: Ha_eV=27.21138386d0 ! 1 Hartree, in eV
		real(8), parameter :: kb_HaK=8.617343d-5/27.21138386d0 ! Boltzmann constant in Ha/K
		real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
		real(8), parameter :: temp_fac_lat=1.d-1 !This percentage of the temperature that should be given to the lattice 
		real(8):: rescale_vel, E_kin

		!Get random Gaussian distributed atomic velocities
		call gausdist(nat,vel_to_scale)
		!scale volocities to the right temperature
		E_kin= 0.d0
		do iat=1,nat
			do i_dim=1,3
				E_kin=E_kin+(vel_to_scale(i_dim,iat)*vel_to_scale(i_dim,iat)*amu(iat))/2.d0
			end do
		end do
		!Now rescale the velocities to give the exact temperature
		rescale_vel=sqrt(3.d0*(nat-1.d0)*kb_HaK*temp/(2.d0*E_kin))
		vel_scaled(:,:)=vel_to_scale(:,:)*rescale_vel
		!get rid of center of mass momentum
		p_total= 0.d0
		do iat=1, nat
			p_total= p_total + amu(iat)*vel_scaled(:,iat)
		end do
		p_total= p_total/nat
		do iat=1, nat
			vel_scaled(:,iat) = vel_scaled(:,iat) - p_total(:)/amu(iat)
		end do
		vel_out= vel_scaled
	end subroutine init_vel_atoms
	
   !subroutine init_vel_lattice(bmass,temp, h_in, h_norm, h_dot_norm, h_dot_out)!, vel_lat_in_bug)
   subroutine init_vel_lattice(bmass,temp, h_in, h_dot_out)!
		implicit none
		!*********************************************************************************************

		real(8), intent(in) :: temp
		real(8), intent(in) :: bmass
		real(8),dimension(3,3), intent(in) :: h_in
		!real(8),dimension(3), intent(out) :: h_norm
		!real(8),dimension(3), intent(out) :: h_dot_norm
		real(8),dimension(3) :: h_norm
		real(8),dimension(3) :: h_dot_norm
		real(8),dimension(3,3), intent(out) :: h_dot_out
		real(8),dimension(3,3) :: h_dot
		real(8),dimension(3,3) :: vel_scaled
		integer:: ih,i_dim
		real(8), parameter :: Ha_eV=27.21138386d0 ! 1 Hartree, in eV
		real(8), parameter :: kb_HaK=8.617343d-5/27.21138386d0 ! Boltzmann constant in Ha/K
		real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
		real(8), parameter :: temp_fac_lat=1.d-1 !This percentage of the temperature that should be given to the lattice 
		real(8):: rescale_vel, cell_kin

		!Get random Gaussian distributed cell velocities
		call gausdist_cell(h_dot)
		cell_kin= sum(bmass*(h_dot*h_dot))/2.d0

		!Now rescale the velocities to give the exact temperature
		rescale_vel=sqrt(3.d0*9.d0*kb_HaK*temp/(2.d0*cell_kin))
		h_dot(:,:)=h_dot(:,:)*rescale_vel
		!Now get some cell velocity
		!call gausdist_cell(latvec_in,vel_lat_in)
		!get rid of center of mass momentum
		!h_dot_out= h_dot
		do ih=1, 3
			h_dot_norm(ih)= sqrt(sum(h_dot(:,ih)*h_dot(:,ih)))
			h_norm(ih)= sqrt(sum(h_in(:,ih)*h_in(:,ih)))
			h_dot_out(:,ih) = (h_dot_norm(ih)/h_norm(ih))*h_in(:,ih)
			!h_dot(:,ih) = (h_dot_norm/1.d0)*h_in(:,ih)
		end do
		h_dot_out= h_dot
	end subroutine init_vel_lattice	

	subroutine get_afx_t(fcart_in, h_inv, amu, nat, afx_t)
		implicit none
		integer, intent(in) :: nat
		real(8), intent(in) :: fcart_in(3,nat)
		real(8),dimension(3,3):: h_inv !inverse matrix of h
		real(8), intent(in) :: amu(nat)
		real(8), intent(out) :: afx_t(3,nat)
		integer :: iat
		do iat=1, nat
			afx_t(:, iat)= matmul(h_inv,fcart_in(:,iat))/amu(iat)
		enddo
		return
	end subroutine get_afx_t

	subroutine get_G_tot(h, h_dot, g_tot)
		implicit none
		real(8), intent(in):: h(3,3)
		real(8), intent(in):: h_dot(3,3)
		real(8),dimension(3,3):: h_trans
		real(8),dimension(3,3):: h_dot_trans
		real(8),dimension(3,3):: g
		real(8),dimension(3,3):: gdot
		real(8),dimension(3,3):: ginv
		real(8), intent(out):: g_tot(3,3)
		integer :: i
		h_trans= 0.d0
		h_dot_trans= 0.d0
		do i=1,3
			h_trans(:,i)=h(i,:)
			h_dot_trans(:,i)=h_dot(i,:)
		enddo
		gdot=matmul(h_dot_trans,h)+matmul(h_trans,h_dot)
		g=matmul(h_trans,h)
		call invertmat(g,ginv,3)
		g_tot=matmul(ginv,gdot)
		return	
	end subroutine get_G_tot
	
	subroutine get_sigma(strten_in, sigma)
		implicit none
		real(8), intent(in) :: strten_in(6)
		real(8), intent(out):: sigma(3,3)
		!ATTENTION: THIS IS WITH TAKING INTO ACCOUNT THE STRESS TENSOR GIVEN IN INPUT
		!strten_in in hartree/borh^3
		sigma(1,1)=strten_in(1)
		sigma(2,2)=strten_in(2)
		sigma(3,3)=strten_in(3)
		sigma(1,2)=strten_in(6)
		sigma(2,1)=strten_in(6)
		sigma(1,3)=strten_in(5)
		sigma(3,1)=strten_in(5)
		sigma(2,3)=strten_in(4)
		sigma(3,2)=strten_in(4)
		return
	end subroutine get_sigma
	
	subroutine get_volume(h, vol)
		implicit none
		real(8), intent(in) :: h(3,3)
		real(8), intent(out) :: vol
		vol= h(1,1)*h(2,2)*h(3,3)-h(1,1)*h(2,3)*h(3,2)-h(1,2)*h(2,1)*h(3,3)+&
		  h(1,2)*h(2,3)*h(3,1)+h(1,3)*h(2,1)*h(3,2)-h(1,3)*h(2,2)*h(3,1)
		return
	end subroutine get_volume
	
	subroutine get_v_mat(amu, v_in, vol, nat, v_mat)
		implicit none
		integer, intent(in) :: nat
		real(8), intent(in) :: v_in(3,nat)
		real(8), intent(in) :: amu(nat)
		real(8), intent(in) :: vol
		real(8), intent(out) :: v_mat(3,3)
		integer :: iat, i, j
		v_mat=0.d0
		do iat=1,nat
			do i=1,3
	   			do j=1,3
			  		v_mat(i,j)= v_mat(i,j)+amu(iat)*v_in(i,iat)*v_in(j,iat)/vol
		   		enddo
			enddo
		enddo
		return
	end subroutine get_v_mat
	
	subroutine get_eta(h, eta)
		implicit none
		real(8), intent(in) :: h(3,3)
		real(8), intent(out):: eta(3,3)
		real(8), dimension(3) :: crossp
		eta= 0.d0
		crossp= 0.d0
		!calculation of eta, eta= [b_coss_c, c_cross_a, a_cross_b]
		call cross_product(h(:,2),h(:,3),crossp); eta(:,1)=crossp
		call cross_product(h(:,3),h(:,1),crossp); eta(:,2)=crossp
		call cross_product(h(:,1),h(:,2),crossp); eta(:,3)=crossp
		return
	end subroutine get_eta

	subroutine get_asv(v_in, amu, nat, asv) 
	!velocity dependent part of the aceleration over the 
	!thermostat degree of freedom
		implicit none
		integer, intent(in) :: nat
		real(8), intent(in) :: v_in(3,nat)
		real(8) :: temp(nat)
		real(8), intent(in) :: amu(nat)
		real(8), intent(out) :: asv
        temp= sum(v_in*v_in, dim=1)
        asv= sum(amu*temp)
		return
	end subroutine get_asv	
	
	
!***********************************************************************
!BEGIN OF BOND CONSTRAIN
!***********************************************************************	
	subroutine sigma_bond(r1, r2, d, sb)
		implicit none
		real(8), intent(in) :: r1(3)
		real(8), intent(in) :: r2(3)
		real(8), intent(in) :: d
		real(8), intent(out) :: sb !value of sigma bond
		real(8) :: r12(3)
		r12= r2 - r1
		sb= dot_product(r12,r12)
		sb= d**2.d0 - sb
		return
	end subroutine sigma_bond
	
	subroutine constrain_bond(r1, r2, sb)
		implicit none
		real(8), intent(in) :: r1(3)
		real(8), intent(in) :: r2(3)
		real(8), intent(out) :: sb !value of sigma bond
		real(8) :: r12(3)
		r12= r2 - r1
		sb= dot_product(r12,r12)
		return
	end subroutine constrain_bond
	
	subroutine sigma_atom_fix(r1, r2, d, sb)
		implicit none
		real(8), intent(in) :: r1(3)
		real(8), intent(in) :: r2(3)
		real(8), intent(in) :: d
		real(8), intent(out) :: sb !value of sigma bond
		real(8) :: r12(3)
		!write(*,*) 'sigma_atom_fix_start'
		!write(*,*) 'r1 ', r1
		!write(*,*) 'r2 ', r2
		r12= r2 - r1
		sb= dot_product(r12,r12)
		!write(*,*) 'sb ', sb
		sb= d**2.d0 - sb
		!write(*,*) 'sigma_atom_fix_end'
		return
	end subroutine sigma_atom_fix
	subroutine grad_sigma_bond(r1,r2, gsb)
		implicit none
		real(8), intent(in) :: r1(3)
		real(8), intent(in) :: r2(3)
		real(8), intent(out) :: gsb(2,3) !value of gradients of sigma bond
		gsb(1,:) = 2.d0*(r2-r1) !gradient respect r1
		gsb(2,:) = 2.d0*(r1-r2) !gradient respect r2
	end subroutine grad_sigma_bond
	subroutine grad_atom_fix(r1,r2, gsb)
		implicit none
		real(8), intent(in) :: r1(3)
		real(8), intent(in) :: r2(3)
		real(8), intent(out) :: gsb(3) !value of gradients of sigma bond
		gsb(:) = 2.d0*(r2-r1) !gradient respect r1
		!gsb(2,:) = 2.d0*(r1-r2) !gradient respect r2
	end subroutine grad_atom_fix
	subroutine lambda_bond(sb, gsb, gsb_o, dt, lb)
		implicit none
		real(8), intent(in) :: sb
		real(8), intent(in) :: gsb(2,3)!value of gradients of sigma bond
		real(8), intent(in) :: gsb_o(2,3)!value of gradients of sigma bond old
		real(8), intent(in) :: dt
		real(8), intent(out) :: lb !lambda value
		real(8) :: temp
		temp= (dt**2.d0)*sum(gsb*gsb_o)
		lb= sb/temp
	end subroutine lambda_bond
	subroutine lambda_atom_fix(sb, gsb, gsb_o, dt, lb)
		implicit none
		real(8), intent(in) :: sb
		real(8), intent(in) :: gsb(3)!value of gradients of sigma bond
		real(8), intent(in) :: gsb_o(3)!value of gradients of sigma bond old
		real(8), intent(in) :: dt
		real(8), intent(out) :: lb !lambda value
		real(8) :: temp
		temp= (dt**2.d0)*dot_product(gsb,gsb_o)!sum(gsb*gsb_o)
		lb= sb/temp
	end subroutine lambda_atom_fix
	subroutine lambda_bond_test(sb, gsb, gsb_o, dt, lb)
		implicit none
		real(8), intent(in) :: sb
		real(8), intent(in) :: gsb(2,3)!value of gradients of sigma bond
		real(8), intent(in) :: gsb_o(2,3)!value of gradients of sigma bond old
		real(8), intent(in) :: dt
		real(8), intent(out) :: lb !lambda value
		real(8) :: temp
		temp= (dt**2.d0)*sum(gsb*gsb_o)
		lb= sb/temp
	end subroutine lambda_bond_test

	subroutine apply_bond_constrain(r1,r2,d,dt, r1_out, r2_out)
		implicit none
		real(8), intent(in) :: r1(3)
		real(8), intent(in) :: r2(3)
		real(8), intent(in) :: d
		real(8), intent(in) :: dt
		real(8), intent(out) :: r1_out(3)
		real(8), intent(out) :: r2_out(3)
		real(8) :: sb !value of sigma bond
		real(8) :: lb !lambda value
		real(8) :: gsb(2,3)!value of gradients of sigma bond
		real(8) :: gsb_o(2,3)!value of gradients of sigma bond old
		real(8) :: r1_l(3)
		real(8) :: r2_l(3)
		r1_l= r1
		r2_l= r2
		call grad_sigma_bond(r1_l, r2_l, gsb_o)
		call sigma_bond(r1_l, r2_l, d, sb)
		!write(*,*) 'sb_0', sb
		do while (abs(sb) > 1.0d-8)
			call grad_sigma_bond(r1_l, r2_l, gsb)
			call lambda_bond(sb, gsb, gsb_o, dt, lb)
			r1_l= r1_l - (dt**2.d0)*lb*gsb(1,:)
			r2_l= r2_l - (dt**2.d0)*lb*gsb(2,:)
			call sigma_bond(r1_l, r2_l, d, sb)
		end do
		r1_out= r1_l
		r2_out= r2_l
		return
	end subroutine apply_bond_constrain
	
	subroutine get_force_bond_constrain(r1,r2,gsb_o,d,dt, F1_out, F2_out)
		implicit none
		real(8), intent(in) :: r1(3)
		real(8), intent(in) :: r2(3)
		real(8), intent(in) :: gsb_o(2,3)!value of gradients of sigma bond old
		real(8), intent(in) :: d
		real(8), intent(in) :: dt
		real(8), intent(out) :: F1_out(3)
		real(8), intent(out) :: F2_out(3)
		real(8) :: sb !value of sigma bond
		real(8) :: lb !lambda value
		real(8) :: gsb(2,3)!value of gradients of sigma bond
		call grad_sigma_bond(r1, r2, gsb)
		call sigma_bond(r1, r2, d, sb)
		call lambda_bond(sb, gsb, gsb_o, dt, lb)
		F1_out= -1.d0*(dt*dt)*lb*gsb(1,:)
		F2_out= -1.d0*(dt*dt)*lb*gsb(2,:)
		return
	end subroutine get_force_bond_constrain

	subroutine get_force_atom_fix_constrain(r1,r2,gsaf_o,dt, F1_out)
		implicit none
		real(8), intent(in) :: r1(3)
		real(8), intent(in) :: r2(3)
		real(8), intent(in) :: gsaf_o(3)!value of gradients of sigma atom fix
		real(8), intent(in) :: dt
		real(8), intent(out) :: F1_out(3)
		real(8) :: saf !value of sigma atom fix
		real(8) :: laf !lambda atom fix
		real(8) :: gsaf(3)!value of gradients of sigma atom fix
		real(8) :: d
		d= 0.d0
		call grad_atom_fix(r1, r2, gsaf)
		call sigma_atom_fix(r1, r2, d, saf)
		call lambda_atom_fix(saf, gsaf, gsaf_o, dt, laf)
		F1_out= -1.d0*(dt*dt)*laf*gsaf(:)
		return
	end subroutine get_force_atom_fix_constrain

	subroutine get_force_bond_constrain_test(r1,r2,gsb_o,d,dt, lb_out, F1_out, F2_out, gsb_out)
		implicit none
		real(8), intent(in) :: r1(3)
		real(8), intent(in) :: r2(3)
		real(8), intent(in) :: gsb_o(2,3)!value of gradients of sigma bond old
		real(8), intent(in) :: d
		real(8), intent(in) :: dt
		real(8), intent(out) :: lb_out
		real(8), intent(out) :: F1_out(3)
		real(8), intent(out) :: F2_out(3)
		real(8), intent(out) :: gsb_out(2,3)
		real(8) :: sb !value of sigma bond
		real(8) :: lb !lambda value
		real(8) :: gsb(2,3)!value of gradients of sigma bond
		call grad_sigma_bond(r1, r2, gsb)
		call sigma_bond(r1, r2, d, sb)
		call lambda_bond(sb, gsb, gsb_o, dt, lb)
		lb_out= lb
		F1_out= -1.d0*(dt*dt)*lb*gsb(1,:)
		F2_out= -1.d0*(dt*dt)*lb*gsb(2,:)
		gsb_out= gsb 
		return
	end subroutine get_force_bond_constrain_test
!***********************************************************************
!END OF BOND CONSTRAIN
!***********************************************************************
!***********************************************************************
!BEGIN OF CELL PARAMETER CONSTRAIN
!***********************************************************************
	subroutine sigma_cell_para(h, a, scp)
		implicit none
		real(8), intent(in) :: h(3)!cell vector to apply constrain
		real(8), intent(in) :: a !value of cell parameter
		real(8), intent(out) :: scp !value of sigma cell parameter
		scp= dot_product(h,h)
		scp= a**2.d0 - scp
		return
	end subroutine sigma_cell_para

	subroutine constrain_cell_para(h, scp)
		implicit none
		real(8), intent(in) :: h(3)!cell vector to apply constrain
		real(8), intent(out) :: scp !value of sigma cell parameter
		scp= dot_product(h,h)
		return
	end subroutine constrain_cell_para
	
	subroutine grad_sigma_cell_para(h, gscp)
		implicit none
		real(8), intent(in) :: h(3)
		real(8), intent(out) :: gscp(3) !value of gradients of sigma bond
		gscp(:) = -2.d0*h(:) !gradient respect r1
	end subroutine grad_sigma_cell_para
	subroutine lambda_cell_para(scp, gscp, gscp_o, dt, lcp)
		implicit none
		real(8), intent(in) :: scp
		real(8), intent(in) :: gscp(3)!value of gradients of sigma bond
		real(8), intent(in) :: gscp_o(3)!value of gradients of sigma bond old
		real(8), intent(in) :: dt
		real(8), intent(out) :: lcp !lambda value
		real(8) :: temp
		temp= (dt**2.d0)*sum(gscp*gscp_o)
		lcp= scp/temp
	end subroutine lambda_cell_para

	subroutine apply_cell_para_constrain(h,a,dt, h_out)
		implicit none
		real(8), intent(in) :: h(3)
		real(8), intent(in) :: a
		real(8), intent(in) :: dt
		real(8), intent(out) :: h_out(3)
		real(8) :: scp !value of sigma cell parameter
		real(8) :: lcp !lambda value cell parameter
		real(8) :: gscp(3)!value of gradients of sigma bond
		real(8) :: gscp_o(3)!value of gradients of sigma bond old
		real(8) :: h_l(3)
		h_l= h
		call grad_sigma_cell_para(h_l, gscp_o)
		call sigma_cell_para(h_l, a, scp)
		!write(*,*) 'sb_0', sb
		do while (abs(scp) > 1.0d-8)
			call grad_sigma_cell_para(h_l, gscp)
			call lambda_cell_para(scp, gscp, gscp_o, dt, lcp)
			h_l= h_l - (dt**2.d0)*lcp*gscp(:)
			call sigma_cell_para(h_l, a, scp)
		end do
		h_out= h_l
		return
	end subroutine apply_cell_para_constrain
	
	subroutine get_force_cell_parameter_constrain(h,gscp_o,a,dt, F_h_out)
		implicit none
		real(8), intent(in) :: h(3)
		real(8), intent(in) :: gscp_o(3)
		real(8), intent(in) :: a
		real(8), intent(in) :: dt
		real(8), intent(out) :: F_h_out(3)
		real(8) :: scp !value of sigma cell parameter
		real(8) :: lcp !lambda value cell parameter
		real(8) :: gscp(3)!value of gradients of sigma bond
		!real(8) :: gscp_o(3)!value of gradients of sigma bond old
		real(8) :: h_l(3)
		h_l= h
		call grad_sigma_cell_para(h_l, gscp)
		call sigma_cell_para(h_l, a, scp)
		call lambda_cell_para(scp, gscp, gscp_o, dt, lcp)
		F_h_out= -1.d0*(dt**2.d0)*lcp*gscp(:)
		return
	end subroutine get_force_cell_parameter_constrain	
!********************************************	
!***********************************************************************
!END OF CELL PARAMETER CONSTRAIN
!***********************************************************************
!***********************************************************************
!BEGIN OF volume CONSTRAIN
!***********************************************************************
	subroutine sigma_vol(h, v0, sv)
		implicit none
		real(8), intent(in) :: h(3,3)!cell vector to apply constrain
		real(8), intent(in) :: v0 !value of cell parameter
		real(8), intent(out) :: sv !value of sigma cell parameter
		real(8) :: vol
		call get_volume(h, vol)
		sv= vol - v0
		return
	end subroutine sigma_vol
	
	subroutine grad_sigma_vol(h, gsv)
		implicit none
		real(8), intent(in) :: h(3,3)
		real(8), intent(out) :: gsv(3,3) !value of gradients of sigma bond
		real(8) :: crossp(3)
		call cross_product(h(:,2),h(:,3),crossp); gsv(:,1)=crossp
		call cross_product(h(:,3),h(:,1),crossp); gsv(:,2)=crossp
		call cross_product(h(:,1),h(:,2),crossp); gsv(:,3)=crossp
	end subroutine grad_sigma_vol
	
	subroutine lambda_vol(sv, gsv, gsv_o, dt, lv)
		implicit none
		real(8), intent(in) :: sv!value of sigma cell parameter
		real(8), intent(in) :: gsv(3,3)!value of gradients of sigma volume
		real(8), intent(in) :: gsv_o(3,3)!value of gradients of sigma volume old
		real(8), intent(in) :: dt
		real(8), intent(out) :: lv !lambda value
		real(8) :: temp
		temp= (dt**2.d0)*sum(gsv*gsv_o)
		lv= sv/temp
	end subroutine lambda_vol

	subroutine apply_vol_constrain(h,v0,dt, h_out)
		implicit none
		real(8), intent(in) :: h(3,3)
		real(8), intent(in) :: v0
		real(8), intent(in) :: dt
		real(8), intent(out) :: h_out(3,3)
		real(8) :: sv !value of sigma volume
		real(8) :: lv !lambda volume
		real(8) :: gsv(3,3)!value of gradients of sigma volume
		real(8) :: gsv_o(3,3)!value of gradients of sigma volume old
		real(8) :: h_l(3,3)
		h_l= h
		call grad_sigma_vol(h_l, gsv_o)
		call sigma_vol(h_l, v0, sv)
		!write(*,*) 'sb_0', sb
		do while (abs(sv) > 1.0d-8)
			call grad_sigma_vol(h_l, gsv)
			call lambda_vol(sv, gsv, gsv_o, dt, lv)
			h_l= h_l - (dt**2.d0)*lv*gsv
			call sigma_vol(h_l, v0, sv)
		end do
		h_out= h_l
		return
	end subroutine apply_vol_constrain
	
	subroutine get_force_vol_constrain(h,gsv_o,v0,dt, F_h_out)
		implicit none
		real(8), intent(in) :: h(3,3)
		real(8), intent(in) :: gsv_o(3,3)
		real(8), intent(in) :: v0
		real(8), intent(in) :: dt
		real(8), intent(out) :: F_h_out(3,3)
		real(8) :: sv !value of sigma volume
		real(8) :: lv !lambda volume
		real(8) :: gsv(3,3)!value of gradients of sigma volume
		!real(8) :: gsv_o(3,3)!value of gradients of sigma volume old
		real(8) :: h_l(3,3)
		h_l= h
		call grad_sigma_vol(h_l, gsv)
		call sigma_vol(h_l, v0, sv)
		call lambda_vol(sv, gsv, gsv_o, dt, lv)
		F_h_out= -1.d0*(dt*dt)*lv*gsv
		return
	end subroutine get_force_vol_constrain	
!***********************************************************************
!END OF volume CONSTRAIN
!***********************************************************************	
!***********************************************************************
!BGIN OF ANGLE COS CONSTRAIN BETWEEN PARTICLES
!***********************************************************************
	subroutine cos_ijk(rij, rik, cos_out)
		implicit none
		real(8), intent(in) :: rij(3)
		real(8), intent(in) :: rik(3)
		real(8), intent(out) :: cos_out
		real(8) :: rij_norm !temporal variables
		real(8) :: rik_norm
		rij_norm= dot_product(rij,rij)**0.5d0
		rik_norm= dot_product(rik,rik)**0.5d0
		cos_out=(dot_product(rij,rik))/(rij_norm*rik_norm)
		return
	end subroutine cos_ijk
	
	subroutine sigma_cos(rij, rik, cos_n, sc)
		implicit none
		real(8), intent(in) :: rij(3)
		real(8), intent(in) :: rik(3)
		real(8), intent(in) :: cos_n
		real(8), intent(out) :: sc !value of sigma cos
		real(8) :: cosijk
		
		call cos_ijk(rij, rik, cosijk)
		sc= cos_n - cosijk
		return
	end subroutine sigma_cos

	subroutine constrain_cos(rij, rik, sc)
		implicit none
		real(8), intent(in) :: rij(3)
		real(8), intent(in) :: rik(3)
		real(8), intent(out) :: sc !value of sigma cos
		real(8) :: cosijk		
		call cos_ijk(rij, rik, cosijk)
		sc= cosijk
		return
	end subroutine constrain_cos
	
	subroutine Dcos_1(rij, rik, Dcos1)
		implicit none
		real(8), intent(in) :: rij(3)
		real(8), intent(in) :: rik(3)
		real(8), intent(out) :: Dcos1(3)
		real(8) :: rij_norm !temporal variables
		real(8) :: rik_norm
		real(8) :: cosijk
		rij_norm= dot_product(rij,rij)**0.5d0
		rik_norm= dot_product(rik,rik)**0.5d0
		call cos_ijk(rij, rik, cosijk)
		Dcos1= rij/(rij_norm*rik_norm) - rik*cosijk/(rik_norm**2.d0)
		!Dcos1= rij/(rij_norm*rik_norm) 
		return
	end subroutine Dcos_1	

	subroutine Dcos_2(rij, rik, Dcos2)
		implicit none
		real(8), intent(in) :: rij(3)
		real(8), intent(in) :: rik(3)
		real(8), intent(out) :: Dcos2(3)
		real(8) :: rij_norm !temporal variables
		real(8) :: rik_norm
		real(8) :: cosijk
		rij_norm= dot_product(rij,rij)**0.5d0
		rik_norm= dot_product(rik,rik)**0.5d0
		call cos_ijk(rij, rik, cosijk)
		Dcos2= rik/(rij_norm*rik_norm) - rij*cosijk/(rij_norm**2.d0)
		!Dcos2= rik/(rij_norm*rik_norm)
		return
	end subroutine Dcos_2

	subroutine grad_sigma_cos(rij, rik, gsc)
		implicit none
		real(8), intent(in) :: rij(3)
		real(8), intent(in) :: rik(3)
		real(8), intent(out) :: gsc(3,3)
		real(8) :: Dcos1(3) !temporal variables
		real(8) :: Dcos2(3)
		call Dcos_1(rij, rik, Dcos1)
		call Dcos_2(rij, rik, Dcos2)
		gsc(1,:) = -1.d0*(Dcos1 + Dcos2)
		gsc(2,:) = Dcos2
		gsc(3,:) = Dcos1
		return
	end subroutine grad_sigma_cos

	subroutine lambda_cos(sc, gsc, gsc_o, dt, lc)
		implicit none
		real(8), intent(in) :: sc !value of sigma cos
		real(8), intent(in) :: gsc(3,3)!value of gradients of sigma cos
		real(8), intent(in) :: gsc_o(3,3)!value of gradients of sigma cos old
		real(8), intent(in) :: dt
		real(8), intent(out) :: lc !lambda value
		real(8) :: temp
		temp= (dt**2.d0)*sum(gsc*gsc_o)
		lc= sc/temp
	end subroutine lambda_cos

	subroutine apply_angle_constrain(ri,rj,rk,cos_n,dt, ri_out, rj_out, rk_out)
		implicit none
		real(8), intent(in) :: ri(3)
		real(8), intent(in) :: rj(3)
		real(8), intent(in) :: rk(3)
		real(8), intent(in) :: cos_n
		real(8), intent(in) :: dt
		real(8), intent(out) :: ri_out(3)
		real(8), intent(out) :: rj_out(3)
		real(8), intent(out) :: rk_out(3)
		real(8) :: sc !value of sigma cos
		real(8) :: lc !lambda value
		real(8) :: gsc(3,3)!value of gradients of sigma cos
		real(8) :: gsc_o(3,3)!value of gradients of sigma cos old
		real(8) :: ri_l(3)
		real(8) :: rj_l(3)
		real(8) :: rk_l(3)
		real(8) :: rij(3)
		real(8) :: rik(3)
		real(8) :: rij_norm0 !initial length
		real(8) :: rik_norm0 !initial length
		real(8) :: rij_norm
		real(8) :: rik_norm
		ri_l= ri
		rj_l= rj
		rk_l= rk
		rij= rj_l - ri_l
		rik= rk_l - ri_l
		rij_norm0= dot_product(rij,rij)**0.5d0
		rik_norm0= dot_product(rik,rik)**0.5d0
		rij= rij/rij_norm0
		rik= rik/rik_norm0
		call grad_sigma_cos(rij, rik, gsc_o)
		call sigma_cos(rij, rik, cos_n, sc)
		!write(*,*) 'sb_0', sb
		do while (abs(sc) > 1.0d-8)
			call grad_sigma_cos(rij, rik, gsc)
			call lambda_cos(sc, gsc, gsc_o, dt, lc)
			!ri_l= ri_l - 0.d0*(dt**2.d0)*lc*gsc(1,:)
			rj_l= rj_l + (dt**2.d0)*lc*gsc(2,:)
			rk_l= rk_l + (dt**2.d0)*lc*gsc(3,:)
			rij= rj_l - ri_l
			rik= rk_l - ri_l
			rij_norm= dot_product(rij,rij)**0.5d0
			rik_norm= dot_product(rik,rik)**0.5d0
			rij= rij/rij_norm
			rik= rik/rik_norm
			call sigma_cos(rij, rik, cos_n, sc)
		end do
		ri_out= ri_l
		rj_out= rj_l
		rk_out= rk_l
		return
	end subroutine apply_angle_constrain
	
	subroutine get_force_angle_constrain(ri,rj,rk,gsc_o,cos_n,dt, Fi_out, Fj_out, Fk_out)
		implicit none
		real(8), intent(in) :: ri(3)
		real(8), intent(in) :: rj(3)
		real(8), intent(in) :: rk(3)
		real(8), intent(in) :: gsc_o(3,3)!value of gradients of sigma cos old
		real(8), intent(in) :: cos_n
		real(8), intent(in) :: dt
		real(8), intent(out) :: Fi_out(3)
		real(8), intent(out) :: Fj_out(3)
		real(8), intent(out) :: Fk_out(3)
		real(8) :: sc !value of sigma cos
		real(8) :: lc !lambda value
		real(8) :: gsc(3,3)!value of gradients of sigma cos
		!real(8) :: gsc_o(3,3)!value of gradients of sigma cos old
		real(8) :: rij(3)
		real(8) :: rik(3)
		real(8) :: rij_norm0 !initial length
		real(8) :: rik_norm0 !initial length
		real(8) :: rij_norm
		real(8) :: rik_norm
		rij= rj - ri
		rik= rk - ri
		rij_norm0= dot_product(rij,rij)**0.5d0
		rik_norm0= dot_product(rik,rik)**0.5d0
		rij= rij/rij_norm0
		rik= rik/rik_norm0
		call grad_sigma_cos(rij, rik, gsc)
		call sigma_cos(rij, rik, cos_n, sc)
		call lambda_cos(sc, gsc, gsc_o, dt, lc)
		Fi_out= -1.d0*0.d0*(dt**2.d0)*lc*gsc(1,:)
		Fj_out= (dt**2.d0)*lc*gsc(2,:)
		Fk_out= (dt**2.d0)*lc*gsc(3,:)
		return
	end subroutine get_force_angle_constrain

	subroutine grad_sigma_cos_cell(rij, rik, gsc)
		implicit none
		real(8), intent(in) :: rij(3)
		real(8), intent(in) :: rik(3)
		real(8), intent(out) :: gsc(3,3)
		real(8) :: Dcos1(3) !temporal variables
		real(8) :: Dcos2(3)
		call Dcos_1(rij, rik, Dcos1)
		call Dcos_2(rij, rik, Dcos2)
		gsc(1,:) = 0.d0
		gsc(2,:) = Dcos2
		gsc(3,:) = Dcos1
		return
	end subroutine grad_sigma_cos_cell

	subroutine apply_cell_angle_constrain(rj,rk,cos_n,dt, rj_out, rk_out)
		implicit none
		!real(8), intent(in) :: ri(3)
		real(8), intent(in) :: rj(3)
		real(8), intent(in) :: rk(3)
		real(8), intent(in) :: cos_n
		real(8), intent(in) :: dt
		!real(8), intent(out) :: ri_out(3)
		real(8), intent(out) :: rj_out(3)
		real(8), intent(out) :: rk_out(3)
		real(8) :: ri(3)
		real(8) :: sc !value of sigma cos
		real(8) :: lc !lambda value
		real(8) :: gsc(3,3)!value of gradients of sigma cos
		real(8) :: gsc_o(3,3)!value of gradients of sigma cos old
		real(8) :: ri_l(3)
		real(8) :: rj_l(3)
		real(8) :: rk_l(3)
		real(8) :: rij(3)
		real(8) :: rik(3)
		real(8) :: rij_norm0 !initial length
		real(8) :: rik_norm0 !initial length
		real(8) :: rij_norm
		real(8) :: rik_norm
		ri= 0.d0
		ri_l= ri
		rj_l= rj
		rk_l= rk
		rij= rj_l - ri_l
		rik= rk_l - ri_l
		rij_norm0= dot_product(rij,rij)**0.5d0
		rik_norm0= dot_product(rik,rik)**0.5d0
		rij= rij/rij_norm0
		rik= rik/rik_norm0
		call grad_sigma_cos_cell(rij, rik, gsc_o)
		call sigma_cos(rij, rik, cos_n, sc)
		!write(*,*) 'sb_0', sb
		do while (abs(sc) > 1.0d-8)
			call grad_sigma_cos_cell(rij, rik, gsc)
			call lambda_cos(sc, gsc, gsc_o, dt, lc)
			!ri_l= ri_l - 0.d0*(dt**2.d0)*lc*gsc(1,:)
			rj_l= rj_l + (dt**2.d0)*lc*gsc(2,:)
			rk_l= rk_l + (dt**2.d0)*lc*gsc(3,:)
			rij= rj_l - ri_l
			rik= rk_l - ri_l
			rij_norm= dot_product(rij,rij)**0.5d0
			rik_norm= dot_product(rik,rik)**0.5d0
			rij= rij/rij_norm
			rik= rik/rik_norm
			call sigma_cos(rij, rik, cos_n, sc)
		end do
		!ri_out= ri_l
		rj_out= rj_l
		rk_out= rk_l
		return
	end subroutine apply_cell_angle_constrain

	subroutine get_force_cell_angle_constrain(rj,rk,gsc_o,cos_n,dt, F_j_out, F_k_out)
		implicit none
		!real(8), intent(in) :: ri(3)
		real(8), intent(in) :: rj(3)
		real(8), intent(in) :: rk(3)
		real(8), intent(in) :: gsc_o(3,3)
		real(8), intent(in) :: cos_n
		real(8), intent(in) :: dt
		!real(8), intent(out) :: ri_out(3)
		real(8), intent(out) :: F_j_out(3)
		real(8), intent(out) :: F_k_out(3)
		real(8) :: ri(3)
		real(8) :: sc !value of sigma cos
		real(8) :: lc !lambda value
		real(8) :: gsc(3,3)!value of gradients of sigma cos
		!real(8) :: gsc_o(3,3)!value of gradients of sigma cos old
		real(8) :: ri_l(3)
		real(8) :: rj_l(3)
		real(8) :: rk_l(3)
		real(8) :: rij(3)
		real(8) :: rik(3)
		real(8) :: rij_norm0 !initial length
		real(8) :: rik_norm0 !initial length
		real(8) :: rij_norm
		real(8) :: rik_norm
		ri= 0.d0
		ri_l= ri
		rj_l= rj
		rk_l= rk
		rij= rj_l - ri_l
		rik= rk_l - ri_l
		rij_norm0= dot_product(rij,rij)**0.5d0
		rik_norm0= dot_product(rik,rik)**0.5d0
		rij= rij/rij_norm0
		rik= rik/rik_norm0
		call grad_sigma_cos_cell(rij, rik, gsc)
		call sigma_cos(rij, rik, cos_n, sc)
		call lambda_cos(sc, gsc, gsc_o, dt, lc)
		!F_i_l= ri_l - 0.d0*(dt**2.d0)*lc*gsc(1,:)
		F_j_out= (dt*dt)*lc*gsc(2,:)
		F_k_out= (dt*dt)*lc*gsc(3,:)
		return
	end subroutine get_force_cell_angle_constrain

	subroutine get_r_t(h, x, nat, r_out)
	implicit none
	integer, intent(in) :: nat
	real(8), intent(in) :: h(3,3)
	real(8), intent(in) :: x(3,nat) !reduced coordinates
	real(8), intent(out) :: r_out(3,nat)
	r_out= matmul(h,x)
	return
	end subroutine get_r_t
	
	subroutine get_x_t(h, r, nat, x_out)
	implicit none
	integer, intent(in) :: nat
	real(8), intent(in) :: h(3,3)
	real(8), intent(in) :: r(3,nat) !reduced coordinates
	real(8), intent(out) :: x_out(3,nat)
	real(8) :: h_inv(3,3)
	call invertmat(h,h_inv,3)
	x_out= matmul(h_inv,r)
	return
	end subroutine get_x_t
	
	subroutine get_x_dot(h, v_in, nat, x_dot_out)
	implicit none
	integer, intent(in) :: nat
	real(8), intent(in) :: h(3,3)
	real(8), intent(in) :: v_in(3,nat) !reduced coordinates
	real(8), intent(out) :: x_dot_out(3,nat)
	real(8) :: h_inv(3,3)
	call invertmat(h,h_inv,3)
	x_dot_out= matmul(h_inv,v_in)
	return
	end subroutine get_x_dot

	subroutine md_npt_no_constrains(latvec_in,xred_in, x_t_dot_in,fcart_in,strten_in,&
						&vel_in,vel_lat_in, target_pressure_habohr, &
						&amu, Qmass, bmass, dtion_md, temp, s_in, s_in_dot, &
						&correc_steps,&
						&md_steps,&
						&nat,&
						&s_out, s_out_dot, h_out, h_dot_out, x_out, x_dot_out, v_out)
	!**********************************************************************************************
	! MOLECULAR DYNAMICS, USING THE LAGRANGIAN METHOD DEVELOPED BY PARRINELLO AND RAHMAN IN:
	![1]M. Parrinello and A. Rahman, J. Appl. Phys. J. Chem. Phys. 521, 2384 (1981).
	!THE CODE HERE, FOLLOWS THE DERIVATIONS AS PRESENTED BY AMSLER AND GOEDECKER IN:
	![2]M. Amsler and S. Goedecker, J. Chem. Phys. 133, (2010).
	!DESCRIPTION OF THE INPUT:
	!latvec_in: 3x3 matrix with the lattice vectors(a,b,c) in bohr, in the code latvec_in -> h, h(:,1)= a, h(:,2)= b
	!xred_in: matrix with the reduced positions of the atoms in the primitive lattice defined 
	!by h, in the code xred_in -> s with dimension (3, number of atoms)
	!fcart_in: matrix with the forces over every atom in hartree/bohr with dimension (3, number of atoms)
	!fcart_in is transformed into fs which is the force in reduced coordinates in hartree with dimension (3, number of atoms)
	!strten_in: the strain tensor in hartree/bohr^3, in the code strten_in-> sigma
	!vel_in: initial velocity of the atoms in bohr/second, in the code vel_in-> v
	!vel_lat_in: initial velocity of the lattice in the code vel_lat_in-> h_dot
	!g= (h^T)*h

		!use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat,ntime_md,&
		!		       &char_type,bmass,mdmin,dtion_md,strfact,units,usewf_md
		!use defs_basis
		!use utilities
		implicit none
		integer, intent(in) :: nat
		integer, intent(in) :: correc_steps
		integer, intent(in) :: md_steps
		real(8), intent(in) :: latvec_in(3,3)
		real(8), intent(in) :: xred_in(3,nat)
		real(8), intent(in) :: x_t_dot_in(3,nat)
		real(8), intent(in) :: fcart_in(3,nat)
		real(8), intent(in) :: strten_in(6)
		real(8), intent(in) :: vel_in(3,nat)
		real(8), intent(in) :: vel_lat_in(3,3)		
		real(8), intent(in) :: target_pressure_habohr
		real(8), intent(in) :: amu(nat)
		real(8), intent(in) :: Qmass !mass of the thermostat
		real(8), intent(in) :: bmass
		real(8), intent(in) :: dtion_md
		real(8), intent(in) :: temp
		real(8), intent(in) :: s_in !thermostat variable at time t
		real(8), intent(in) :: s_in_dot !thermostat variable at time t+ dt
		real(8), intent(out) :: s_out !thermostat variable at time t
		real(8), intent(out) :: s_out_dot !thermostat variable at time t+ dt
		!real(8), intent(out) :: as_out
		!real(8), intent(out) :: asvt_out
		real(8),dimension(3,3), intent(out):: h_out
		real(8),dimension(3,3), intent(out):: h_dot_out
		!real(8),dimension(3,3), intent(out) :: ah_out
		real(8),dimension(3,nat), intent(out):: x_out!reduced coordinates of atoms
		real(8),dimension(3,nat), intent(out):: x_dot_out !reduced velocity of atoms
		real(8),dimension(3,nat), intent(out):: v_out
		!real(8),dimension(3,nat), intent(out):: ax_out
		!*********************************************************************
		!Variables for my MD part
		real(8), parameter :: Ha_eV=27.21138386d0 ! 1 Hartree, in eV
		real(8), parameter :: kb_HaK=8.617343d-5/27.21138386d0 ! Boltzmann constant in Ha/K
		real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
		real(8):: amass(nat)
		real(8):: pressure_md(3,3)
		real(8):: vol_t
		real(8):: vol_t_dt
		!integer :: nat= 2
		real(8),dimension(3,nat):: afx_t !reduced aceleration over atom
		!real(8),dimension(3,nat):: vpospred
		real(8),dimension(3,nat):: x_t !xred_in position of the atoms in reduced coordinates
		real(8),dimension(3,nat):: x_t_dt !xred_in position of the atoms in reduced coordinates
		real(8),dimension(3,nat):: r_t !xred_in position of the atoms in reduced coordinates
		real(8) :: r1_out(3)
		real(8) :: r2_out(3)
		real(8) :: r3_out(3)
		real(8),dimension(3,nat):: x_t_dot !v in reduced coordinates
		real(8),dimension(3,nat):: x_c_dot !v in reduced coordinates
		real(8),dimension(3,nat):: v_t !vel_in
		real(8),dimension(3,nat):: v_c !vel_in
		real(8),dimension(3,nat):: ax_t !vel_in
		real(8),dimension(3,nat):: ax_c !vel_in
		real(8),dimension(3,3):: h_t !matix with the lattice vectors
		real(8),dimension(3,3):: h_t_dt !matix with the lattice vectors
		real(8),dimension(3,3):: h_t_inv! trasnpose matrix of h
		real(8),dimension(3,3):: h_c_dot !derivative of h, matrix with the lattice velocity vel_lat_in
		real(8),dimension(3,3):: h_t_dot !derivative of h, matrix with the lattice velocity vel_lat_in
		real(8),dimension(3,3):: v_t_mat
		real(8),dimension(3,3):: v_c_mat
		real(8),dimension(3,3):: g_t_tot
		real(8),dimension(3,3):: g_c_tot
		real(8),dimension(3,3):: ah_t
		real(8),dimension(3,3):: ah_c
		real(8),dimension(3,3):: unitmat
		real(8),dimension(3,3):: f_h !force for the lattice equation of motion
		real(8),dimension(3,3):: eta_t
		real(8),dimension(3,3):: eta_t_dt
		real(8),dimension(3,3):: sigma
		real(8):: crossp(3)
		real(8) :: s_t !thermostat variable at time t
		real(8) :: s_t_dt !thermostat variable at time t+ dt
		real(8) :: s_t_dot !time derivative of thermostat variable at time t
		real(8) :: s_c_dot !time derivative of thermostat variable at correction steps
		real(8) :: asvt ! aceleration dependent on the velocity over the thermostat at time t
		real(8) :: asvt1 ! aceleration dependent on the velocity over the thermostat at time t
		real(8) :: asvc ! aceleration dependent on the velocity over the thermostat at correction step
		real(8) :: as_t ! aceleration over the thermostat at time t
		real(8) :: as_c ! aceleration over the thermostat at correction steps
		real(8):: dt
		real(8):: ndf
		!real(8):: temp1
		!real(8):: temp2
		integer:: md_i,correc_step_i,  i, j, iat, a, b, c
		!Assign masses to each atom (for MD)
		x_t=xred_in !matrix with reduce coordinates
		dt=dtion_md !define time step
		h_t=latvec_in !latice matrix h=[a, b, c] in bohr
		h_t_dot= vel_lat_in
		call invertmat(h_t,h_t_inv,3)
		v_t=vel_in !copy velocity matrix [3, number_of_atoms]
		x_t_dot= x_t_dot_in!matmul(h_t_inv,v_t)
		s_t= s_in
		s_t_dot= s_in_dot
		ndf= 3.d0*nat
		unitmat=0.d0
		do i=1,3
			unitmat(i,i)=1.d0
		enddo
		pressure_md(:,:)= target_pressure_habohr*unitmat(:,:) !pressure matrix hartree/bohr^3
		call get_afx_t(fcart_in, h_t_inv, amu, nat, afx_t)
        call get_sigma(strten_in, sigma)
        call get_volume(h_t, vol_t)
        call get_eta(h_t, eta_t)
		call get_G_tot(h_t, h_t_dot, g_t_tot)
		call get_v_mat(amu, v_t, vol_t, nat, v_t_mat)  
		call get_asv(v_t, amu, nat, asvt)
        do md_i= 1, md_steps		
			ax_t= afx_t/(s_t**2.d0) - 2.d0*(s_t_dot/s_t)*x_t_dot - matmul(g_t_tot,x_t_dot)
			as_t= (asvt/s_t - (ndf+1.d0)*kb_HaK*temp/s_t)/Qmass
			ah_t= v_t_mat - sigma - pressure_md
			ah_t= matmul(ah_t,eta_t)/bmass
			
			x_t_dt= x_t + dt*x_t_dot + 0.5d0*dt*dt*ax_t
			s_t_dt= s_t + dt*s_t_dot + 0.5d0*dt*dt*as_t
			h_t_dt= h_t + dt*h_t_dot + 0.5d0*dt*dt*ah_t

			call get_volume(h_t_dt, vol_t_dt)
			call get_eta(h_t_dt, eta_t_dt)

			!prediction correction
			!prediction, calculate input for correction data
			!temp1= sum(x_t_dot*x_t_dot)
			!write(*,*) 'x_t_dot 1', temp1
			x_c_dot= x_t_dot + dt*ax_t
			s_c_dot= s_t_dot + dt*as_t
			h_c_dot= h_t_dot + dt*ah_t
			call get_G_tot(h_t_dt, h_c_dot, g_c_tot)
			v_c= s_t_dt*(matmul(h_t_dt, x_c_dot))
			call get_v_mat(amu, v_c, vol_t_dt, nat, v_c_mat) 			 
			call get_asv(v_c, amu, nat, asvc)
			
			!correction steps
			do correc_step_i=1, correc_steps
				ax_c= afx_t/(s_t_dt**2.d0) - 2.d0*(s_c_dot/s_t_dt)*x_c_dot &
					& - matmul(g_c_tot,x_c_dot)
					
				as_c= (asvc/s_t_dt - (ndf+1.d0)*kb_HaK*temp/s_t_dt)/Qmass
					
				ah_c = v_c_mat - sigma - pressure_md
				ah_c= matmul(ah_c,eta_t_dt)/bmass
				x_c_dot= x_t_dot + 0.5d0*dt*(ax_c + ax_t)
				s_c_dot= s_t_dot + 0.5d0*dt*(as_c + as_t)
				h_c_dot= h_t_dot + 0.5d0*dt*(ah_c + ah_t)
				call get_G_tot(h_t_dt, h_c_dot, g_c_tot)
				v_c= s_t_dt*(matmul(h_t_dt, x_c_dot))
				call get_v_mat(amu, v_c, vol_t_dt, nat, v_c_mat)
				call get_asv(v_c, amu, nat, asvc) 
			end do !correction
			x_t= x_t_dt
			x_t_dot= x_c_dot
			h_t= h_t_dt
			h_t_dot= h_c_dot
			call invertmat(h_t,h_t_inv,3)
			s_t= s_t_dt
			s_t_dot= s_c_dot
			v_t= s_t*(matmul(h_t, x_t_dot))
			call get_afx_t(fcart_in, h_t_inv, amu, nat, afx_t)
			!call get_sigma(strten_in, sigma)
			call get_volume(h_t, vol_t)
			call get_eta(h_t, eta_t)
			call get_G_tot(h_t, h_t_dot, g_t_tot)
			call get_v_mat(amu, v_t, vol_t, nat, v_t_mat)  
			call get_asv(v_t, amu, nat, asvt)
		end do !md
		s_out= s_t
		s_out_dot= s_t_dot
		!as_out= as_t
		!asvt_out= asvt
		h_out= h_t
		h_dot_out= h_t_dot
		!ah_out= ah_t
		x_out= x_t
		x_dot_out= x_t_dot
		!ax_out= ax_t
		v_out= v_t
		return
	end subroutine md_npt_no_constrains

	subroutine cell_constrains_volume(h_in,dt, cell_para_valu, cell_angl_valu,&
						&volu_valu,cell_para_const,cell_angl_const,&
						&bool_cell_para_cons,bool_cell_angl_cons,bool_volu,&
						&nat,numb_cell_angl_cons, numb_cell_para_cons,volu_cons,&
						&h_out)
		implicit none
		integer, intent(in) :: nat, numb_cell_angl_cons, numb_cell_para_cons,volu_cons
		integer, intent(in) :: bool_cell_para_cons
		integer, intent(in) :: bool_cell_angl_cons
		integer, intent(in) :: bool_volu
		integer, intent(in) :: cell_para_const(numb_cell_para_cons,1)
		integer, intent(in) :: cell_angl_const(numb_cell_angl_cons,2)
		real(8), intent(in) :: cell_para_valu(numb_cell_para_cons)
		real(8), intent(in) :: cell_angl_valu(numb_cell_angl_cons)
		real(8), intent(in) :: volu_valu(volu_cons)
		real(8), intent(in) :: h_in(3,3)
		real(8), intent(in) :: dt
		real(8), intent(out) :: h_out(3,3)
		real(8) :: h_l(3,3)
		real(8) :: Force(3,3)
		real(8) :: F_a(3)
		real(8) :: F_b(3)
		real(8) :: F_c(3)
		real(8) :: F_h(3,3)
		real(8) :: agsc_o(numb_cell_angl_cons,3,3)!arra gradient of sigma cos
		real(8) :: agscp_o(numb_cell_para_cons,3)!arra gradient of sigma cell para
		real(8) :: gsv_o(3,3) !gradient sigma volume
		real(8) :: gsc(3,3)!arra gradient of sigma cos
		real(8) :: gscp(3)!arra gradient of sigma bond
		real(8) :: gsv(3,3)!value of gradients of sigma volume
		integer :: sigma_cell_para_state(numb_cell_para_cons)
		integer :: sigma_angle_cell_state(numb_cell_angl_cons)
		integer :: sigma_volu_state
		!real(8) :: rij(3)
		!real(8) :: rik(3)
		real(8) :: sc !value of sigma cos
		real(8) :: scp !value of sigma cell parameter
		real(8) :: sv !value of sigma volume
		real(8) :: lv !lambda volum
		real(8) :: volu
		integer :: i, j, a, b, c
		h_l= h_in
		if (((bool_cell_angl_cons == 1) .AND. (bool_cell_para_cons == 1))&
				& .AND. (bool_volu == 1)) then
			Force= 0.d0
			sigma_cell_para_state= 0
			sigma_angle_cell_state = 0
			sigma_volu_state = 0
			agsc_o = 0.d0
			agscp_o = 0.d0
			call grad_sigma_vol(h_l, gsv_o)
			call sigma_vol(h_l, volu_valu(1), sv)
			!write(*,*) 'sv' ,sv
			if (abs(sv) < 1.0d-5) then
				sigma_volu_state = 1
			end if
			do i=1, numb_cell_para_cons
				a= cell_para_const(i,1)
				call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
				call grad_sigma_cell_para(h_l(:,a), gscp)
				agscp_o(i,:)= gscp(:)
				if (abs(scp) < 1.0d-5) then
					sigma_cell_para_state(i) = 1
				end if
			end do!numb_cell_para_cons
			
			do i=1, numb_cell_angl_cons
				a= cell_angl_const(i,1)
				b= cell_angl_const(i,2)
				call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
				call grad_sigma_cos(h_l(:,a), h_l(:,b), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_angle_cell_state(i) = 1
				end if
			end do!numb_cell_angl_cons
			j=0
			do while (((product(sigma_cell_para_state) == 0) .OR. &
						&(product(sigma_angle_cell_state) == 0)) .OR. &
						& (sigma_volu_state== 0))
!********************************************
!Force calculation
				Force= 0.d0
				F_h= 0.d0
				j= j +1
				!write(*,*)'j ', j
				!write(*,*)'product(sigma_cell_para_state) ', product(sigma_cell_para_state)
				!write(*,*)'product(sigma_angle_cell_state) ', product(sigma_angle_cell_state)
				!write(*,*)'sigma_volu_state ', sigma_volu_state
				call get_force_vol_constrain(h_l,gsv_o,volu_valu(1),dt, F_h)
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_cell_parameter_constrain(h_l(:,a),&
						&agscp_o(i,:),cell_para_valu(i),dt, F_a)
					Force(:,a)= Force(:,a) + F_a
				end do!numb_cell_para_cons				
				do i=1,  numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call  get_force_cell_angle_constrain(h_l(:,a), h_l(:,b),&
						&agsc_o(i,:,:),cell_angl_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_cell_para_cons
!End Force calculation
!********************************************
				h_l= h_l + Force + F_h
!********************************************
!re calculate sigma states
				call sigma_vol(h_l, volu_valu(1), sv)
				!write(*,*) 'sv' ,sv
				if (abs(sv) < 1.0d-5) then
					sigma_volu_state = 1
				else
					sigma_volu_state = 0
				end if
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
					if (abs(scp) < 1.0d-5) then
						sigma_cell_para_state(i) = 1
					else
						sigma_cell_para_state(i) = 0
					end if
				end do!numb_cell_para_cons
				
				do i=1, numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_angle_cell_state(i) = 1
					else
						sigma_angle_cell_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! if bond angle constrain
		
		if (((bool_cell_angl_cons == 0) .AND. (bool_cell_para_cons == 1)) &
			& .AND. (bool_volu == 1)) then
			Force= 0.d0
			sigma_cell_para_state= 0
			sigma_volu_state = 0
			agscp_o = 0.d0
			call grad_sigma_vol(h_l, gsv_o)
			call sigma_vol(h_l, volu_valu(1), sv)
			if (abs(sv) < 1.0d-5) then
				sigma_volu_state = 1
			end if
			do i=1, numb_cell_para_cons
				a= cell_para_const(i,1)
				call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
				call grad_sigma_cell_para(h_l(:,a), gscp)
				agscp_o(i,:)= gscp(:)
				if (abs(scp) < 1.0d-5) then
					sigma_cell_para_state(i) = 1
				end if
			end do!numb_cell_para_cons
			
			do while ((product(sigma_cell_para_state) == 0) .OR. &
						& (sigma_volu_state== 0))
!********************************************
!Force calculation
				Force= 0.d0
				F_h= 0.d0
				call get_force_vol_constrain(h_l,gsv_o,volu_valu(1),dt, F_h)
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_cell_parameter_constrain(h_l(:,a),&
						&agscp_o(i,:),cell_para_valu(i),dt, F_a)
					Force(:,a)= Force(:,a) + F_a
				end do!numb_cell_para_cons				
!End Force calculation
!********************************************
				h_l= h_l + Force + F_h
!********************************************
!re calculate sigma states
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
					if (abs(scp) < 1.0d-5) then
						sigma_cell_para_state(i) = 1
					else
						sigma_cell_para_state(i) = 0
					end if
				end do!numb_cell_para_cons
				call sigma_vol(h_l, volu_valu(1), sv)
				if (abs(sv) < 1.0d-5) then
					sigma_volu_state = 1
				else
					sigma_volu_state = 0
				end if
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! cell para constrain
		
		if (((bool_cell_angl_cons == 1) .AND. (bool_cell_para_cons == 0)) &
			& .AND. (bool_volu == 1)) then
			Force= 0.d0
			sigma_angle_cell_state = 0
			sigma_volu_state = 0
			agsc_o = 0.d0
			call grad_sigma_vol(h_l, gsv_o)
			call sigma_vol(h_l, volu_valu(1), sv)
			if (abs(sv) < 1.0d-5) then
				sigma_volu_state = 1
			end if
			do i=1, numb_cell_angl_cons
				a= cell_angl_const(i,1)
				b= cell_angl_const(i,2)
				call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
				call grad_sigma_cos(h_l(:,a), h_l(:,b), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_angle_cell_state(i) = 1
				end if
			end do!numb_cell_angl_cons
			
			do while (product(sigma_angle_cell_state) == 0 .OR. &
						& (sigma_volu_state== 0))
!********************************************
!Force calculation
				Force= 0.d0	
				F_h= 0.d0
				call get_force_vol_constrain(h_l,gsv_o,volu_valu(1),dt, F_h)		
				do i=1,  numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call  get_force_cell_angle_constrain(h_l(:,a), h_l(:,b),&
						&agsc_o(i,:,:),cell_angl_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_cell_para_cons
!End Force calculation
!********************************************
				h_l= h_l + Force + F_h
!********************************************
!re calculate sigma states
				do i=1, numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_angle_cell_state(i) = 1
					else
						sigma_angle_cell_state(i) = 0
					end if
				end do!numb_angl_cons
				call sigma_vol(h_l, volu_valu(1), sv)
				if (abs(sv) < 1.0d-5) then
					sigma_volu_state = 1
				else
					sigma_volu_state = 0
				end if
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! if cell_angl constraine

		if ((bool_cell_angl_cons == 1) .AND. (bool_cell_para_cons == 1)&
			& .AND. (bool_volu == 0)) then
			Force= 0.d0
			sigma_cell_para_state= 0
			sigma_angle_cell_state = 0
			agsc_o = 0.d0
			agscp_o = 0.d0
			do i=1, numb_cell_para_cons
				a= cell_para_const(i,1)
				call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
				call grad_sigma_cell_para(h_l(:,a), gscp)
				agscp_o(i,:)= gscp(:)
				if (abs(scp) < 1.0d-5) then
					sigma_cell_para_state(i) = 1
				end if
			end do!numb_cell_para_cons
			
			do i=1, numb_cell_angl_cons
				a= cell_angl_const(i,1)
				b= cell_angl_const(i,2)
				call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
				call grad_sigma_cos(h_l(:,a), h_l(:,b), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_angle_cell_state(i) = 1
				end if
			end do!numb_cell_angl_cons
			
			do while ((product(sigma_cell_para_state) == 0) .OR. (product(sigma_angle_cell_state) == 0))
!********************************************
!Force calculation
				Force= 0.d0
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_cell_parameter_constrain(h_l(:,a),&
						&agscp_o(i,:),cell_para_valu(i),dt, F_a)
					Force(:,a)= Force(:,a) + F_a
				end do!numb_cell_para_cons				
				do i=1,  numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call  get_force_cell_angle_constrain(h_l(:,a), h_l(:,b),&
						&agsc_o(i,:,:),cell_angl_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_cell_para_cons
!End Force calculation
!********************************************
				h_l= h_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
					if (abs(scp) < 1.0d-5) then
						sigma_cell_para_state(i) = 1
					else
						sigma_cell_para_state(i) = 0
					end if
				end do!numb_cell_para_cons
				
				do i=1, numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_angle_cell_state(i) = 1
					else
						sigma_angle_cell_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! if bond angle constrain
		
		if ((bool_cell_angl_cons == 0) .AND. (bool_cell_para_cons == 1) &
			& .AND. (bool_volu == 0)) then
			Force= 0.d0
			sigma_cell_para_state= 0
			agscp_o = 0.d0
			do i=1, numb_cell_para_cons
				a= cell_para_const(i,1)
				call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
				call grad_sigma_cell_para(h_l(:,a), gscp)
				agscp_o(i,:)= gscp(:)
				if (abs(scp) < 1.0d-5) then
					sigma_cell_para_state(i) = 1
				end if
			end do!numb_cell_para_cons
			
			do while ((product(sigma_cell_para_state) == 0))
!********************************************
!Force calculation
				Force= 0.d0
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_cell_parameter_constrain(h_l(:,a),&
						&agscp_o(i,:),cell_para_valu(i),dt, F_a)
					Force(:,a)= Force(:,a) + F_a
				end do!numb_cell_para_cons				
!End Force calculation
!********************************************
				h_l= h_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
					if (abs(scp) < 1.0d-5) then
						sigma_cell_para_state(i) = 1
					else
						sigma_cell_para_state(i) = 0
					end if
				end do!numb_cell_para_cons

!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! cell para constrain
		
		if ((bool_cell_angl_cons == 1) .AND. (bool_cell_para_cons == 0)&
			& .AND. (bool_volu == 0)) then
			Force= 0.d0
			sigma_angle_cell_state = 0
			agsc_o = 0.d0
			do i=1, numb_cell_angl_cons
				a= cell_angl_const(i,1)
				b= cell_angl_const(i,2)
				call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
				call grad_sigma_cos(h_l(:,a), h_l(:,b), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_angle_cell_state(i) = 1
				end if
			end do!numb_cell_angl_cons
			
			do while (product(sigma_angle_cell_state) == 0)
!********************************************
!Force calculation
				Force= 0.d0			
				do i=1,  numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call  get_force_cell_angle_constrain(h_l(:,a), h_l(:,b),&
						&agsc_o(i,:,:),cell_angl_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_cell_para_cons
!End Force calculation
!********************************************
				h_l= h_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_angle_cell_state(i) = 1
					else
						sigma_angle_cell_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! if cell_angl constraine

		if ((bool_cell_angl_cons == 0) .AND. (bool_cell_para_cons == 0) &
			& .AND. (bool_volu == 1)) then
			Force= 0.d0
			sigma_volu_state = 0
			call grad_sigma_vol(h_l, gsv_o)
			call sigma_vol(h_l, volu_valu(1), sv)
			!write(*,*) 'sv' ,sv
			if (abs(sv) < 1.0d-5) then
				sigma_volu_state = 1
			end if
			
			do while (sigma_volu_state== 0)
!********************************************
!Force calculation
				Force= 0.d0
				F_h= 0.d0
				call get_force_vol_constrain(h_l,gsv_o,volu_valu(1),dt, F_h)
				
!End Force calculation
!********************************************
				h_l= h_l + Force + F_h
!********************************************
!re calculate sigma states
				call sigma_vol(h_l, volu_valu(1), sv)
				!!write(*,*) 'sv' ,sv
				if (abs(sv) < 1.0d-5) then
					sigma_volu_state = 1
				else
					sigma_volu_state = 0
				end if
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! cell para constrain


		if ((bool_cell_angl_cons == 0) .AND. (bool_cell_para_cons == 0) &
			& .AND. (bool_volu == 0)) then
			Force= 0.d0				
!End Force calculation
!********************************************
				h_l= h_l + Force 
!********************************************
		end if! cell para constrain

		h_out= h_l
		return
	end subroutine cell_constrains_volume


	subroutine cell_constrains(h_in,dt, cell_para_valu, cell_angl_valu,&
						&cell_para_const,cell_angl_const,&
						&bool_cell_para_cons,bool_cell_angl_cons,&
						&nat,numb_cell_angl_cons, numb_cell_para_cons,&
						&h_out)
		implicit none
		integer, intent(in) :: nat, numb_cell_angl_cons, numb_cell_para_cons
		integer, intent(in) :: bool_cell_para_cons
		integer, intent(in) :: bool_cell_angl_cons
		integer, intent(in) :: cell_para_const(numb_cell_para_cons,1)
		integer, intent(in) :: cell_angl_const(numb_cell_angl_cons,2)
		real(8), intent(in) :: cell_para_valu(numb_cell_para_cons)
		real(8), intent(in) :: cell_angl_valu(numb_cell_angl_cons)
		real(8), intent(in) :: h_in(3,3)
		real(8), intent(in) :: dt
		real(8), intent(out) :: h_out(3,3)
		real(8) :: h_l(3,3)
		real(8) :: Force(3,3)
		real(8) :: F_a(3)
		real(8) :: F_b(3)
		real(8) :: F_c(3)
		real(8) :: agsc_o(numb_cell_angl_cons,3,3)!arra gradient of sigma cos
		real(8) :: agscp_o(numb_cell_para_cons,3)!arra gradient of sigma cell para
		real(8) :: gsc(3,3)!arra gradient of sigma cos
		real(8) :: gscp(3)!arra gradient of sigma bond
		integer :: sigma_cell_para_state(numb_cell_para_cons)
		integer :: sigma_angle_cell_state(numb_cell_angl_cons)
		!real(8) :: rij(3)
		!real(8) :: rik(3)
		real(8) :: sc !value of sigma cos
		real(8) :: scp !value of sigma cell parameter
		integer :: i, j, a, b, c
		h_l= h_in
		if ((bool_cell_angl_cons == 1) .AND. (bool_cell_para_cons == 1)) then
			Force= 0.d0
			sigma_cell_para_state= 0
			sigma_angle_cell_state = 0
			agsc_o = 0.d0
			agscp_o = 0.d0
			do i=1, numb_cell_para_cons
				a= cell_para_const(i,1)
				call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
				call grad_sigma_cell_para(h_l(:,a), gscp)
				agscp_o(i,:)= gscp(:)
				if (abs(scp) < 1.0d-5) then
					sigma_cell_para_state(i) = 1
				end if
			end do!numb_cell_para_cons
			
			do i=1, numb_cell_angl_cons
				a= cell_angl_const(i,1)
				b= cell_angl_const(i,2)
				call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
				call grad_sigma_cos(h_l(:,a), h_l(:,b), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_angle_cell_state(i) = 1
				end if
			end do!numb_cell_angl_cons
			j=0
			do while ((product(sigma_cell_para_state) == 0) .OR. (product(sigma_angle_cell_state) == 0))
!********************************************
!Force calculation
				j= j +1
				!!write(*,*)'j ', j
				!!write(*,*)'product(sigma_cell_para_state) ', product(sigma_cell_para_state)
				!!write(*,*)'product(sigma_angle_cell_state) ', product(sigma_angle_cell_state)
				Force= 0.d0
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_cell_parameter_constrain(h_l(:,a),&
						&agscp_o(i,:),cell_para_valu(i),dt, F_a)
					Force(:,a)= Force(:,a) + F_a
				end do!numb_cell_para_cons				
				do i=1,  numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call  get_force_cell_angle_constrain(h_l(:,a), h_l(:,b),&
						&agsc_o(i,:,:),cell_angl_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_cell_para_cons
!End Force calculation
!********************************************
				h_l= h_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
					!write(*,*)'scp ', scp, ' i ', i
					if (abs(scp) < 1.0d-5) then
						sigma_cell_para_state(i) = 1
					else
						sigma_cell_para_state(i) = 0
					end if
				end do!numb_cell_para_cons
				
				do i=1, numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
					!write(*,*)'sc ', sc, ' i ', i
					if (abs(sc) < 1.0d-5) then
						sigma_angle_cell_state(i) = 1
					else
						sigma_angle_cell_state(i) = 0
					end if
				end do!numb_angl_cons
				!write(*,*)'j ', j
				!write(*,*)'product(sigma_cell_para_state) ', product(sigma_cell_para_state)
				!write(*,*) (product(sigma_cell_para_state) .EQ. 0)
				!write(*,*)'product(sigma_angle_cell_state) ', product(sigma_angle_cell_state)
				!write(*,*) (product(sigma_angle_cell_state) .EQ. 0)
				!write(*,*) ((product(sigma_cell_para_state) .EQ. 0) .AND. (product(sigma_angle_cell_state) .EQ. 0))
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! if bond angle constrain
		
		if ((bool_cell_angl_cons == 0) .AND. (bool_cell_para_cons == 1)) then
			Force= 0.d0
			sigma_cell_para_state= 0
			agscp_o = 0.d0
			do i=1, numb_cell_para_cons
				a= cell_para_const(i,1)
				call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
				call grad_sigma_cell_para(h_l(:,a), gscp)
				agscp_o(i,:)= gscp(:)
				if (abs(scp) < 1.0d-5) then
					sigma_cell_para_state(i) = 1
				end if
			end do!numb_cell_para_cons
			
			do while ((product(sigma_cell_para_state) == 0))
!********************************************
!Force calculation
				Force= 0.d0
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_cell_parameter_constrain(h_l(:,a),&
						&agscp_o(i,:),cell_para_valu(i),dt, F_a)
					Force(:,a)= Force(:,a) + F_a
				end do!numb_cell_para_cons				
!End Force calculation
!********************************************
				h_l= h_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_cell_para_cons
					a= cell_para_const(i,1)
					call sigma_cell_para(h_l(:,a), cell_para_valu(i), scp)
					if (abs(scp) < 1.0d-5) then
						sigma_cell_para_state(i) = 1
					else
						sigma_cell_para_state(i) = 0
					end if
				end do!numb_cell_para_cons

!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! cell para constrain
		
		if ((bool_cell_angl_cons == 1) .AND. (bool_cell_para_cons == 0)) then
			Force= 0.d0
			sigma_angle_cell_state = 0
			agsc_o = 0.d0
			do i=1, numb_cell_angl_cons
				a= cell_angl_const(i,1)
				b= cell_angl_const(i,2)
				call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
				call grad_sigma_cos(h_l(:,a), h_l(:,b), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_angle_cell_state(i) = 1
				end if
			end do!numb_cell_angl_cons
			
			do while (product(sigma_angle_cell_state) == 0)
!********************************************
!Force calculation
				Force= 0.d0			
				do i=1,  numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call  get_force_cell_angle_constrain(h_l(:,a), h_l(:,b),&
						&agsc_o(i,:,:),cell_angl_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_cell_para_cons
!End Force calculation
!********************************************
				h_l= h_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_cell_angl_cons
					a= cell_angl_const(i,1)
					b= cell_angl_const(i,2)
					call sigma_cos(h_l(:,a), h_l(:,b), cell_angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_angle_cell_state(i) = 1
					else
						sigma_angle_cell_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! if cell_angl constraine
		h_out= h_l
		return
	end subroutine cell_constrains

	subroutine atomic_constrains(r_in,dt, bond_valu, angl_valu,&
						&bond_const,angl_const,&
						&bool_bond_cons,bool_angl_cons,&
						&nat,numb_angl_cons, numb_bond_cons,&
						&r_out)
		implicit none
		integer, intent(in) :: nat, numb_angl_cons, numb_bond_cons
		integer, intent(in) :: bool_bond_cons
		integer, intent(in) :: bool_angl_cons
		integer, intent(in) :: bond_const(numb_bond_cons,2)
		integer, intent(in) :: angl_const(numb_angl_cons,3)
		real(8), intent(in) :: bond_valu(numb_bond_cons)
		real(8), intent(in) :: angl_valu(numb_angl_cons)
		real(8), intent(in) :: r_in(3,nat)
		real(8), intent(in) :: dt
		real(8), intent(out) :: r_out(3,nat)
		real(8) :: r_l(3,nat)
		real(8) :: Force(3,nat)
		real(8) :: F_a(3)
		real(8) :: F_b(3)
		real(8) :: F_c(3)
		real(8) :: agsc_o(numb_angl_cons,3,3)!arra gradient of sigma cos
		real(8) :: agsb_o(numb_bond_cons,2,3)!arra gradient of sigma bond
		real(8) :: gsc(3,3)!arra gradient of sigma cos
		real(8) :: gsb(2,3)!arra gradient of sigma bond
		integer :: sigma_bond_state(numb_bond_cons)
		integer :: sigma_cos_state(numb_angl_cons)
		!real(8) :: rij(3)
		!real(8) :: rik(3)
		real(8) :: sc !value of sigma cos
		real(8) :: sb !value of sigma bond
		integer :: i, j, a, b, c
		r_l= r_in
		if ((bool_angl_cons == 1) .AND. (bool_bond_cons == 1)) then
			Force= 0.d0
			sigma_bond_state= 0
			sigma_cos_state = 0
			agsc_o = 0.d0
			agsb_o = 0.d0
			do i=1, numb_bond_cons
				a= bond_const(i,1)
				b= bond_const(i,2)
				call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
				call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
				agsb_o(i,:,:)= gsb(:,:)
				if (abs(sb) < 1.0d-5) then
					sigma_bond_state(i) = 1
				end if
			end do!numb_bond_cons
			
			do i=1, numb_angl_cons
				a= angl_const(i,1)
				b= angl_const(i,2)
				c= angl_const(i,3)
				call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
				call grad_sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_cos_state(i) = 1
				end if
			end do!numb_angl_cons
			
			do while ((product(sigma_bond_state) == 0) .OR. (product(sigma_cos_state) == 0))
!********************************************
!Force calculation
				Force= 0.d0
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_bond_constrain(r_l(:,a), r_l(:,b),&
						&agsb_o(i,:,:),bond_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_bond_cons				
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call get_force_angle_constrain(r_l(:,a),r_l(:,b),&
						&r_l(:,c),agsc_o(i,:,:),angl_valu(i),dt,&
						& F_a, F_b, F_c)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
					Force(:,c)= Force(:,c) + F_c
				end do!numb_angl_cons
!End Force calculation
!********************************************
				r_l= r_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
					if (abs(sb) < 1.0d-5) then
						sigma_bond_state(i) = 1
					else
						sigma_bond_state(i) = 0
					end if
				end do!numb_bond_cons
				
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_cos_state(i) = 1
					else
						sigma_cos_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! if bond angle constrain
		
		if ((bool_angl_cons == 0) .AND. (bool_bond_cons == 1)) then
			Force= 0.d0
			sigma_bond_state= 0
			agsb_o = 0.d0
			do i=1, numb_bond_cons
				a= bond_const(i,1)
				b= bond_const(i,2)
				call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
				call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
				agsb_o(i,:,:)= gsb(:,:)
				if (abs(sb) < 1.0d-5) then
					sigma_bond_state(i) = 1
				end if
			end do
!*********************************************
			do while (product(sigma_bond_state) == 0)
!********************************************
!Force calculation
				Force= 0.d0
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_bond_constrain(r_l(:,a), r_l(:,b),&
						&agsb_o(i,:,:),bond_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_bond_cons				
!End Force calculation
!********************************************
				r_l= r_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
					if (abs(sb) < 1.0d-5) then
						sigma_bond_state(i) = 1
					else
						sigma_bond_state(i) = 0
					end if
				end do!numb_bond_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if !only bond
		
		if ((bool_angl_cons == 1) .AND. (bool_bond_cons == 0)) then
			Force= 0.d0
			sigma_cos_state = 0
			agsc_o = 0.d0
			do i=1, numb_angl_cons
				a= angl_const(i,1)
				b= angl_const(i,2)
				c= angl_const(i,3)
				call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
				call grad_sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_cos_state(i) = 1
				end if
			end do
!********************************************
			do while (product(sigma_cos_state) == 0)
!********************************************
!Force calculation
				Force= 0.d0
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call get_force_angle_constrain(r_l(:,a),r_l(:,b),&
						&r_l(:,c),agsc_o(i,:,:),angl_valu(i),dt,&
						& F_a, F_b, F_c)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
					Force(:,c)= Force(:,c) + F_c
				end do!numb_angl_cons
!End Force calculation
!********************************************
				r_l= r_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_cos_state(i) = 1
					else
						sigma_cos_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if !only angle
		if ((bool_angl_cons == 0) .AND. (bool_bond_cons == 0)) then
			Force= 0.d0
!********************************************
				r_l= r_l + Force
!********************************************
!re calculate sigma states
		end if !no costrain		
		
		r_out= r_l
		return
	end subroutine atomic_constrains
	
	subroutine atomic_constrains_atom_fix(r_in,dt, bond_valu, angl_valu,&
						&atom_fix_valu,atom_fix_cord,&
						&bond_const,angl_const,atom_fix_const,&
						&bool_bond_cons,bool_angl_cons,bool_atom_fix_cons,&
						&nat,numb_angl_cons, numb_bond_cons,numb_atom_fix_cons,&
						&r_out)
		implicit none
		integer, intent(in) :: nat, numb_angl_cons, numb_bond_cons,numb_atom_fix_cons
		integer, intent(in) :: bool_bond_cons
		integer, intent(in) :: bool_angl_cons
		integer, intent(in) :: bool_atom_fix_cons
		integer, intent(in) :: bond_const(numb_bond_cons,2)
		integer, intent(in) :: angl_const(numb_angl_cons,3)
		integer, intent(in) :: atom_fix_const(numb_atom_fix_cons,1)
		real(8), intent(in) :: bond_valu(numb_bond_cons)
		real(8), intent(in) :: angl_valu(numb_angl_cons)
		real(8), intent(in) :: atom_fix_valu(numb_atom_fix_cons,3)
		real(8), intent(in) :: atom_fix_cord(numb_atom_fix_cons,3)
		real(8), intent(in) :: r_in(3,nat)
		real(8), intent(in) :: dt
		real(8), intent(out) :: r_out(3,nat)
		real(8) :: r_l(3,nat)
		real(8) :: Force(3,nat)
		real(8) :: Force_atom(3,nat)
		real(8) :: Force_bond(3,nat)
		real(8) :: F_a(3)
		real(8) :: F_b(3)
		real(8) :: F_c(3)
		real(8) :: agsc_o(numb_angl_cons,3,3)!arra gradient of sigma cos
		real(8) :: agsb_o(numb_bond_cons,2,3)!arra gradient of sigma bond
		real(8) :: agaf_o(numb_atom_fix_cons,3)
		real(8) :: gsc(3,3)!arra gradient of sigma cos
		real(8) :: gsb(2,3)!arra gradient of sigma bond
		real(8) :: gaf(3)
		real(8) :: x(3)! mask from atom_fix_cord
		integer :: sigma_bond_state(numb_bond_cons)
		integer :: sigma_cos_state(numb_angl_cons)
		integer :: sigma_atom_fix_state(numb_atom_fix_cons)
		!real(8) :: rij(3)
		!real(8) :: rik(3)
		real(8) :: sc !value of sigma cos
		real(8) :: sb !value of sigma bond
		real(8) :: saf !value of sigma atom fix
		integer :: i, j, a, b, c
		r_l= r_in
!beging of bool_atom_fix_cons == 1
		if ((bool_angl_cons == 1) .AND. (bool_bond_cons == 1) &
			&.AND. (bool_atom_fix_cons == 1)) then
			Force= 0.d0
			sigma_bond_state= 0
			sigma_cos_state = 0
			agsc_o = 0.d0
			agsb_o = 0.d0
			sigma_atom_fix_state= 0
			agaf_o= 0.d0
			do i=1, numb_atom_fix_cons
				a= atom_fix_const(i,1)
				x= atom_fix_cord(i,:)
				call sigma_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
									&0.d0, saf)
				call grad_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
									&gaf)
				agaf_o(i,:)= gaf
				if (abs(saf) < 1.0d-5) then
					sigma_atom_fix_state(i) = 1
				end if
			end do!numb_atom_fix_cons
			do i=1, numb_bond_cons
				a= bond_const(i,1)
				b= bond_const(i,2)
				call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
				call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
				agsb_o(i,:,:)= gsb(:,:)
				if (abs(sb) < 1.0d-5) then
					sigma_bond_state(i) = 1
				end if
			end do!numb_bond_cons
			
			do i=1, numb_angl_cons
				a= angl_const(i,1)
				b= angl_const(i,2)
				c= angl_const(i,3)
				call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
				call grad_sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_cos_state(i) = 1
				end if
			end do!numb_angl_cons
			
			do while ((product(sigma_bond_state) == 0) .OR. &
						&(product(sigma_cos_state) == 0) .OR. &
						& (product(sigma_atom_fix_state) == 0))
!********************************************
!Force calculation
				Force= 0.d0
				Force_atom= 0.d0
				!Force_bond= 0.d0
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					call get_force_atom_fix_constrain(x*r_l(:,a),&
							&x*atom_fix_valu(i,:),agaf_o(i,:),dt, F_a)
					Force_atom(:,a)= Force_atom(:,a) + F_a
				end do!numb_atom_fix_cons				
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_bond_constrain(r_l(:,a), r_l(:,b),&
						&agsb_o(i,:,:),bond_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_bond_cons				
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call get_force_angle_constrain(r_l(:,a),r_l(:,b),&
						&r_l(:,c),agsc_o(i,:,:),angl_valu(i),dt,&
						& F_a, F_b, F_c)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
					Force(:,c)= Force(:,c) + F_c
				end do!numb_angl_cons
				!fixing the particles that must be fixed
				!it means that the bond force don't move partiles that
				!should be fixed
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					Force(:,a)= (1.d0 - x)*Force(:,a)
				end do!numb_atom_fix_cons
!End Force calculation
!********************************************
				r_l= r_l + Force + Force_atom
!********************************************
!re calculate sigma states
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					call sigma_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
										&0.d0, saf)
					!write(*,*) 'saf  ', saf
					if (abs(saf) < 1.0d-5) then
						sigma_atom_fix_state(i) = 1
					else
						sigma_atom_fix_state(i) = 0
					end if
				end do!numb_atom_fix_cons
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
					!write(*,*) 'i  ',i, ' sb ',sb
					if (abs(sb) < 1.0d-5) then
						sigma_bond_state(i) = 1
					else
						sigma_bond_state(i) = 0
					end if
				end do!numb_bond_cons
				
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
					!write(*,*) 'sc  ', sc
					if (abs(sc) < 1.0d-5) then
						sigma_cos_state(i) = 1
					else
						sigma_cos_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! if bond angle  atom constrain
		
		if ((bool_angl_cons == 0) .AND. (bool_bond_cons == 1) &
			&.AND. (bool_atom_fix_cons == 1)) then
			Force= 0.d0
			sigma_bond_state= 0
			agsb_o = 0.d0
			sigma_atom_fix_state= 0
			agaf_o= 0.d0
			do i=1, numb_atom_fix_cons
				a= atom_fix_const(i,1)
				x= atom_fix_cord(i,:)
				call sigma_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
									&0.d0, saf)
				call grad_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
									&gaf)
				agaf_o(i,:)= gaf
				if (abs(saf) < 1.0d-5) then
					sigma_atom_fix_state(i) = 1
				end if
			end do!numb_atom_fix_cons
			do i=1, numb_bond_cons
				a= bond_const(i,1)
				b= bond_const(i,2)
				call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
				call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
				agsb_o(i,:,:)= gsb(:,:)
				if (abs(sb) < 1.0d-5) then
					sigma_bond_state(i) = 1
				end if
			end do
!*********************************************
			do while ((product(sigma_bond_state) == 0) .OR. &
						&(product(sigma_atom_fix_state) == 0))
!********************************************
!Force calculation
				Force= 0.d0
				Force_atom= 0.d0
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					call get_force_atom_fix_constrain(x*r_l(:,a),&
							&x*atom_fix_valu(i,:),agaf_o(i,:),dt, F_a)
					Force_atom(:,a)= Force_atom(:,a) + F_a
				end do!numb_atom_fix_cons	
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_bond_constrain(r_l(:,a), r_l(:,b),&
						&agsb_o(i,:,:),bond_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_bond_cons
				!fixing the particles that must be fixed
				!it means that the bond force don't move partiles that
				!should be fixed
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					Force(:,a)= (1.d0 - x)*Force(:,a)
				end do!numb_atom_fix_cons				
!End Force calculation
!********************************************
				r_l= r_l + Force + Force_atom
!********************************************
!re calculate sigma states
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					call sigma_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
										&0.d0, saf)
					if (abs(saf) < 1.0d-5) then
						sigma_atom_fix_state(i) = 1
					else
						sigma_atom_fix_state(i) = 0
					end if
				end do!numb_atom_fix_cons
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
					if (abs(sb) < 1.0d-5) then
						sigma_bond_state(i) = 1
					else
						sigma_bond_state(i) = 0
					end if
				end do!numb_bond_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if !only bond
		
		if ((bool_angl_cons == 1) .AND. (bool_bond_cons == 0) &
			&.AND. (bool_atom_fix_cons == 1)) then
			Force= 0.d0
			sigma_cos_state = 0
			agsc_o = 0.d0
			sigma_atom_fix_state= 0
			agaf_o= 0.d0
			do i=1, numb_atom_fix_cons
				a= atom_fix_const(i,1)
				x= atom_fix_cord(i,:)
				call sigma_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
									&0.d0, saf)
				call grad_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
									&gaf)
				agaf_o(i,:)= gaf
				if (abs(saf) < 1.0d-5) then
					sigma_atom_fix_state(i) = 1
				end if
			end do!numb_atom_fix_cons
			do i=1, numb_angl_cons
				a= angl_const(i,1)
				b= angl_const(i,2)
				c= angl_const(i,3)
				call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
				call grad_sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_cos_state(i) = 1
				end if
			end do
!********************************************
			do while ((product(sigma_cos_state) == 0) .OR. &
						& (product(sigma_atom_fix_state) == 0) )
!********************************************
!Force calculation
				Force= 0.d0
				Force_atom= 0.d0
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					call get_force_atom_fix_constrain(x*r_l(:,a),&
							&x*atom_fix_valu(i,:),agaf_o(i,:),dt, F_a)
					Force_atom(:,a)= Force_atom(:,a) + F_a
				end do!numb_atom_fix_cons
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call get_force_angle_constrain(r_l(:,a),r_l(:,b),&
						&r_l(:,c),agsc_o(i,:,:),angl_valu(i),dt,&
						& F_a, F_b, F_c)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
					Force(:,c)= Force(:,c) + F_c
				end do!numb_angl_cons
				!fixing the particles that must be fixed
				!it means that the bond force don't move partiles that
				!should be fixed
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					Force(:,a)= (1.d0 - x)*Force(:,a)
				end do!numb_atom_fix_cons
!End Force calculation
!********************************************
				r_l= r_l + Force + Force_atom
!********************************************
!re calculate sigma states
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					call sigma_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
										&0.d0, saf)
					if (abs(saf) < 1.0d-5) then
						sigma_atom_fix_state(i) = 1
					else
						sigma_atom_fix_state(i) = 0
					end if
				end do!numb_atom_fix_cons
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_cos_state(i) = 1
					else
						sigma_cos_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if !only angle
		if ((bool_angl_cons == 0) .AND. (bool_bond_cons == 0) &
			&.AND. (bool_atom_fix_cons == 1)) then
			!write(*,*) 'here'
			Force= 0.d0
			sigma_atom_fix_state= 0
			agaf_o= 0.d0
			do i=1, numb_atom_fix_cons
				a= atom_fix_const(i,1)
				x= atom_fix_cord(i,:)
				!write(*,*) ' a ',a, ' x ', x
				call sigma_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
									&0.d0, saf)
				!write(*,*)' a ',a, ' saf ', saf
				call grad_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
									&gaf)
				!write(*,*)' a ',a, ' gaf ', gaf
				agaf_o(i,:)= gaf
				if (abs(saf) < 1.0d-6) then
					sigma_atom_fix_state(i) = 1
				end if
			end do!numb_atom_fix_cons
			do while ((product(sigma_atom_fix_state) == 0) )
!Force calculation
				Force= 0.d0
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					call get_force_atom_fix_constrain(x*r_l(:,a),&
							&x*atom_fix_valu(i,:),agaf_o(i,:),dt, F_a)
					Force(:,a)= Force(:,a) + F_a
					!write(*,*) ' a ',a, ' F_a ', F_a
					!write(*,*) ' a ',a, ' x ', x
				end do!numb_atom_fix_cons				
!********************************************
				r_l= r_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_atom_fix_cons
					a= atom_fix_const(i,1)
					x= atom_fix_cord(i,:)
					call sigma_atom_fix(x*r_l(:,a), x*atom_fix_valu(i,:),&
										&0.d0, saf)
					if (abs(saf) < 1.0d-6) then
						sigma_atom_fix_state(i) = 1
					else
						sigma_atom_fix_state(i) = 0
					end if
				end do!numb_atom_fix_cons
			end do!while over sigmas
		end if !only atom_fix_const		
!end of bool_atom_fix_cons = 1

		if ((bool_angl_cons == 1) .AND. (bool_bond_cons == 1)&
			&.AND. (bool_atom_fix_cons == 0)) then
			Force= 0.d0
			sigma_bond_state= 0
			sigma_cos_state = 0
			agsc_o = 0.d0
			agsb_o = 0.d0
			do i=1, numb_bond_cons
				a= bond_const(i,1)
				b= bond_const(i,2)
				call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
				call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
				agsb_o(i,:,:)= gsb(:,:)
				if (abs(sb) < 1.0d-5) then
					sigma_bond_state(i) = 1
				end if
			end do!numb_bond_cons
			
			do i=1, numb_angl_cons
				a= angl_const(i,1)
				b= angl_const(i,2)
				c= angl_const(i,3)
				call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
				call grad_sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_cos_state(i) = 1
				end if
			end do!numb_angl_cons
			
			do while ((product(sigma_bond_state) == 0) .OR. (product(sigma_cos_state) == 0))
!********************************************
!Force calculation
				Force= 0.d0
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_bond_constrain(r_l(:,a), r_l(:,b),&
						&agsb_o(i,:,:),bond_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_bond_cons				
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call get_force_angle_constrain(r_l(:,a),r_l(:,b),&
						&r_l(:,c),agsc_o(i,:,:),angl_valu(i),dt,&
						& F_a, F_b, F_c)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
					Force(:,c)= Force(:,c) + F_c
				end do!numb_angl_cons
!End Force calculation
!********************************************
				r_l= r_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
					if (abs(sb) < 1.0d-5) then
						sigma_bond_state(i) = 1
					else
						sigma_bond_state(i) = 0
					end if
				end do!numb_bond_cons
				
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_cos_state(i) = 1
					else
						sigma_cos_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if! if bond angle constrain
		
		if ((bool_angl_cons == 0) .AND. (bool_bond_cons == 1)&
			&.AND. (bool_atom_fix_cons == 0)) then
			Force= 0.d0
			sigma_bond_state= 0
			agsb_o = 0.d0
			do i=1, numb_bond_cons
				a= bond_const(i,1)
				b= bond_const(i,2)
				call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
				call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
				agsb_o(i,:,:)= gsb(:,:)
				if (abs(sb) < 1.0d-5) then
					sigma_bond_state(i) = 1
				end if
			end do
!*********************************************
			do while (product(sigma_bond_state) == 0)
!********************************************
!Force calculation
				Force= 0.d0
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					!call grad_sigma_bond(r_l(:,a), r_l(:,b), gsb)
					call get_force_bond_constrain(r_l(:,a), r_l(:,b),&
						&agsb_o(i,:,:),bond_valu(i),dt, F_a, F_b)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
				end do!numb_bond_cons				
!End Force calculation
!********************************************
				r_l= r_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_bond_cons
					a= bond_const(i,1)
					b= bond_const(i,2)
					call sigma_bond(r_l(:,a), r_l(:,b), bond_valu(i), sb)
					if (abs(sb) < 1.0d-5) then
						sigma_bond_state(i) = 1
					else
						sigma_bond_state(i) = 0
					end if
				end do!numb_bond_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if !only bond
		
		if ((bool_angl_cons == 1) .AND. (bool_bond_cons == 0)&
			&.AND. (bool_atom_fix_cons == 0)) then
			Force= 0.d0
			sigma_cos_state = 0
			agsc_o = 0.d0
			do i=1, numb_angl_cons
				a= angl_const(i,1)
				b= angl_const(i,2)
				c= angl_const(i,3)
				call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
				call grad_sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), gsc)
				agsc_o(i,:,:)= gsc(:,:)
				if (abs(sc) < 1.0d-5) then
					sigma_cos_state(i) = 1
				end if
			end do
!********************************************
			do while (product(sigma_cos_state) == 0)
!********************************************
!Force calculation
				Force= 0.d0
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call get_force_angle_constrain(r_l(:,a),r_l(:,b),&
						&r_l(:,c),agsc_o(i,:,:),angl_valu(i),dt,&
						& F_a, F_b, F_c)
					Force(:,a)= Force(:,a) + F_a
					Force(:,b)= Force(:,b) + F_b
					Force(:,c)= Force(:,c) + F_c
				end do!numb_angl_cons
!End Force calculation
!********************************************
				r_l= r_l + Force
!********************************************
!re calculate sigma states
				do i=1, numb_angl_cons
					a= angl_const(i,1)
					b= angl_const(i,2)
					c= angl_const(i,3)
					call sigma_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), angl_valu(i), sc)
					if (abs(sc) < 1.0d-5) then
						sigma_cos_state(i) = 1
					else
						sigma_cos_state(i) = 0
					end if
				end do!numb_angl_cons
!End re calculate sigma states
!********************************************
			end do !while over sigmas
		end if !only angle
		if ((bool_angl_cons == 0) .AND. (bool_bond_cons == 0)&
			&.AND. (bool_atom_fix_cons == 0)) then
			Force= 0.d0
!********************************************
				r_l= r_l + Force
!********************************************
!re calculate sigma states
		end if !no costrain	
		
		r_out= r_l
		return
	end subroutine atomic_constrains_atom_fix	



	subroutine atomic_constrains_evaluation(r_in,dt, bond_valu, angl_valu,&
						&atom_fix_valu,atom_fix_cord,&
						&bond_const,angl_const,atom_fix_const,&
						&bool_bond_cons,bool_angl_cons,bool_atom_fix_cons,&
						&nat,numb_angl_cons, numb_bond_cons,numb_atom_fix_cons,&
						&value_bond_constrain_out, value_cos_constrain_out)
		implicit none
		integer, intent(in) :: nat, numb_angl_cons, numb_bond_cons,numb_atom_fix_cons
		integer, intent(in) :: bool_bond_cons
		integer, intent(in) :: bool_angl_cons
		integer, intent(in) :: bool_atom_fix_cons
		integer, intent(in) :: bond_const(numb_bond_cons,2)
		integer, intent(in) :: angl_const(numb_angl_cons,3)
		integer, intent(in) :: atom_fix_const(numb_atom_fix_cons,1)
		real(8), intent(in) :: bond_valu(numb_bond_cons)
		real(8), intent(in) :: angl_valu(numb_angl_cons)
		real(8), intent(in) :: atom_fix_valu(numb_atom_fix_cons,3)
		real(8), intent(in) :: atom_fix_cord(numb_atom_fix_cons,3)
		real(8), intent(in) :: r_in(3,nat)
		real(8), intent(in) :: dt
		real(8), intent(out) :: value_bond_constrain_out(numb_bond_cons)
		real(8), intent(out) :: value_cos_constrain_out(numb_angl_cons)
		real(8) :: r_l(3,nat)
		real(8) :: Force(3,nat)
		real(8) :: Force_atom(3,nat)
		real(8) :: Force_bond(3,nat)
		real(8) :: F_a(3)
		real(8) :: F_b(3)
		real(8) :: F_c(3)
		real(8) :: agsc_o(numb_angl_cons,3,3)!arra gradient of sigma cos
		real(8) :: agsb_o(numb_bond_cons,2,3)!arra gradient of sigma bond
		real(8) :: agaf_o(numb_atom_fix_cons,3)
		real(8) :: gsc(3,3)!arra gradient of sigma cos
		real(8) :: gsb(2,3)!arra gradient of sigma bond
		real(8) :: gaf(3)
		real(8) :: x(3)! mask from atom_fix_cord
		integer :: sigma_bond_state(numb_bond_cons)
		integer :: sigma_cos_state(numb_angl_cons)
		integer :: sigma_atom_fix_state(numb_atom_fix_cons)
		!real(8) :: rij(3)
		!real(8) :: rik(3)
		real(8) :: sc !value of sigma cos
		real(8) :: sb !value of sigma bond
		real(8) :: saf !value of sigma atom fix
		integer :: i, j, a, b, c
		r_l= r_in
		value_bond_constrain_out= 0.d0
		value_cos_constrain_out= 0.d0
		if ((bool_angl_cons == 1) .AND. (bool_bond_cons == 1)) then
			do i=1, numb_bond_cons
				a= bond_const(i,1)
				b= bond_const(i,2)
				call constrain_bond(r_l(:,a), r_l(:,b), sb)
				value_bond_constrain_out(i)= sb
			end do!numb_bond_cons
			
			do i=1, numb_angl_cons
				a= angl_const(i,1)
				b= angl_const(i,2)
				c= angl_const(i,3)
				call constrain_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), sc)
				value_cos_constrain_out(i)= sc
			end do!numb_angl_cons
		end if! if bond angle constrain
		
		if ((bool_angl_cons == 0) .AND. (bool_bond_cons == 1)) then
			do i=1, numb_bond_cons
				a= bond_const(i,1)
				b= bond_const(i,2)
				call constrain_bond(r_l(:,a), r_l(:,b), sb)
				value_bond_constrain_out(i)= sb
			end do!numb_bond_consgmas
		end if !only bond
		
		if ((bool_angl_cons == 1) .AND. (bool_bond_cons == 0)) then
			do i=1, numb_angl_cons
				a= angl_const(i,1)
				b= angl_const(i,2)
				c= angl_const(i,3)
				call constrain_cos(r_l(:,b) - r_l(:,a),r_l(:,c) - r_l(:,a), sc)
				value_cos_constrain_out(i)= sc
			end do!numb_angl_cons
		end if !only angle
		return
	end subroutine atomic_constrains_evaluation
	
	subroutine get_af_t(fcart_in, amu, nat, af_t)
		implicit none
		integer, intent(in) :: nat
		real(8), intent(in) :: fcart_in(3,nat)
		real(8), intent(in) :: amu(nat)
		real(8), intent(out) :: af_t(3,nat)
		integer :: iat
		do iat=1, nat
			af_t(:, iat)= fcart_in(:,iat)/amu(iat)
		enddo
		return
	end subroutine get_af_t
	
	subroutine md_npt_constrains(latvec_in,xred_in, x_t_dot_in,fcart_in,strten_in,&
						&vel_in,vel_lat_in, bond_valu, angl_valu, cell_para_valu,&
						&cell_angl_valu, volu_valu, atom_fix_valu, atom_fix_cord,&
						&target_pressure_habohr, &
						&amu, Qmass, bmass, dtion_md, temp, s_in, s_in_dot, &
						&bond_const,angl_const, cell_para_const,&
						&cell_angl_const,atom_fix_const, correc_steps,&
						&md_steps,bool_bond_cons,bool_angl_cons,bool_cell_para_cons,&
						&bool_cell_angl_cons,bool_volu,bool_atom_fix_cons,&
						&volu_cons, nat, numb_cell_angl_cons,numb_cell_para_cons,&
						&numb_angl_cons, numb_bond_cons,numb_atom_fix_cons,&
						&s_out, s_out_dot, pressure_out, volu_out,&
						&bond_constrain_out, cos_constrain_out, h_out,&
						& h_dot_out, x_out, x_dot_out, v_out)
	!**********************************************************************************************
	! MOLECULAR DYNAMICS, USING THE LAGRANGIAN METHOD DEVELOPED BY PARRINELLO AND RAHMAN IN:
	![1]M. Parrinello and A. Rahman, J. Appl. Phys. J. Chem. Phys. 521, 2384 (1981).
	!THE CODE HERE, FOLLOWS THE DERIVATIONS AS PRESENTED BY AMSLER AND GOEDECKER IN:
	![2]M. Amsler and S. Goedecker, J. Chem. Phys. 133, (2010).
	!DESCRIPTION OF THE INPUT:
	!latvec_in: 3x3 matrix with the lattice vectors(a,b,c) in bohr, in the code latvec_in -> h, h(:,1)= a, h(:,2)= b
	!xred_in: matrix with the reduced positions of the atoms in the primitive lattice defined 
	!by h, in the code xred_in -> s with dimension (3, number of atoms)
	!fcart_in: matrix with the forces over every atom in hartree/bohr with dimension (3, number of atoms)
	!fcart_in is transformed into fs which is the force in reduced coordinates in hartree with dimension (3, number of atoms)
	!strten_in: the strain tensor in hartree/bohr^3, in the code strten_in-> sigma
	!vel_in: initial velocity of the atoms in bohr/second, in the code vel_in-> v
	!vel_lat_in: initial velocity of the lattice in the code vel_lat_in-> h_dot
	!g= (h^T)*h

		!use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat,ntime_md,&
		!		       &char_type,bmass,mdmin,dtion_md,strfact,units,usewf_md
		!use defs_basis
		!use utilities
		implicit none
		integer, intent(in) :: volu_cons, nat, numb_cell_angl_cons,numb_cell_para_cons,&
								&numb_angl_cons, numb_bond_cons,numb_atom_fix_cons
		integer, intent(in) :: correc_steps
		integer, intent(in) :: md_steps
		integer, intent(in) :: bool_bond_cons
		integer, intent(in) :: bool_angl_cons
		integer, intent(in) :: bool_cell_para_cons
		integer, intent(in) :: bool_cell_angl_cons
		integer, intent(in) :: bool_volu
		integer, intent(in) :: bool_atom_fix_cons
		integer, intent(in) :: bond_const(numb_bond_cons,2)
		integer, intent(in) :: angl_const(numb_angl_cons,3)
		integer, intent(in) :: cell_para_const(numb_cell_para_cons,1)
		integer, intent(in) :: cell_angl_const(numb_cell_angl_cons,2)
		integer, intent(in) :: atom_fix_const(numb_atom_fix_cons,1)
		real(8), intent(in) :: latvec_in(3,3)
		real(8), intent(in) :: xred_in(3,nat)
		real(8), intent(in) :: x_t_dot_in(3,nat)
		real(8), intent(in) :: fcart_in(3,nat)
		real(8), intent(in) :: strten_in(6)
		real(8), intent(in) :: vel_in(3,nat)
		real(8), intent(in) :: vel_lat_in(3,3)
		real(8), intent(in) :: bond_valu(numb_bond_cons)
		real(8), intent(in) :: angl_valu(numb_angl_cons)
		real(8), intent(in) :: cell_para_valu(numb_cell_para_cons)
		real(8), intent(in) :: cell_angl_valu(numb_cell_angl_cons)
		real(8), intent(in) :: volu_valu(volu_cons)
		real(8), intent(in) :: atom_fix_valu(numb_atom_fix_cons,3)
		real(8), intent(in) :: atom_fix_cord(numb_atom_fix_cons,3)
		real(8), intent(in) :: target_pressure_habohr
		real(8), intent(in) :: amu(nat)
		real(8), intent(in) :: Qmass !mass of the thermostat
		real(8), intent(in) :: bmass
		real(8), intent(in) :: dtion_md
		real(8), intent(in) :: temp
		real(8), intent(in) :: s_in !thermostat variable at time t
		real(8), intent(in) :: s_in_dot !thermostat variable at time t+ dt
		real(8), intent(out) :: s_out !thermostat variable at time t
		real(8), intent(out) :: s_out_dot !thermostat variable at time t+ dt
		!real(8), intent(out) :: as_out
		!real(8), intent(out) :: asvt_out
		real(8),dimension(3,3), intent(out):: h_out
		real(8),dimension(3,3), intent(out):: h_dot_out
		!real(8),dimension(3,3), intent(out) :: ah_out
		real(8),dimension(3,nat), intent(out):: x_out!reduced coordinates of atoms
		real(8),dimension(3,nat), intent(out):: x_dot_out !reduced velocity of atoms
		real(8),dimension(3,nat), intent(out):: v_out
		real(8),dimension(md_steps), intent(out) :: pressure_out
		real(8),dimension(md_steps), intent(out) :: volu_out
		real(8),dimension(md_steps, numb_bond_cons), intent(out) :: bond_constrain_out
		real(8),dimension(md_steps, numb_angl_cons), intent(out) :: cos_constrain_out
		!real(8),dimension(3,nat), intent(out):: ax_out
		!*********************************************************************
		!Variables for my MD part
		real(8), parameter :: Ha_eV=27.21138386d0 ! 1 Hartree, in eV
		real(8), parameter :: kb_HaK=8.617343d-5/27.21138386d0 ! Boltzmann constant in Ha/K
		real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
		real(8):: amass(nat)
		real(8):: pressure_md(3,3)
		real(8):: vol_t
		real(8):: vol_t_dt
		!integer :: nat= 2
		real(8),dimension(3,nat):: afx_t !reduced aceleration over atom
		!real(8),dimension(3,nat):: vpospred
		real(8),dimension(3,nat):: x_t !xred_in position of the atoms in reduced coordinates
		real(8),dimension(3,nat):: x_t_dt !xred_in position of the atoms in reduced coordinates
		real(8),dimension(3,nat):: r_t !xred_in position of the atoms in reduced coordinates
		real(8),dimension(3,nat):: r_t_dt
		real(8) :: r1_out(3)
		real(8) :: r2_out(3)
		real(8) :: r3_out(3)
		real(8),dimension(3,nat):: x_t_dot !v in reduced coordinates
		real(8),dimension(3,nat):: x_c_dot !v in reduced coordinates
		real(8),dimension(3,nat):: v_t !vel_in
		real(8),dimension(3,nat):: v_c !vel_in
		real(8),dimension(3,nat):: ax_t !vel_in
		real(8),dimension(3,nat):: ax_c !vel_in
		real(8),dimension(3,3):: h_t !matix with the lattice vectors
		real(8),dimension(3,3):: h_t_dt !matix with the lattice vectors
		real(8),dimension(3,3):: h_t_inv! trasnpose matrix of h
		real(8),dimension(3,3):: h_c_dot !derivative of h, matrix with the lattice velocity vel_lat_in
		real(8),dimension(3,3):: h_t_dot !derivative of h, matrix with the lattice velocity vel_lat_in
		real(8),dimension(3,3):: v_t_mat
		real(8),dimension(3,3):: v_c_mat
		real(8),dimension(3,3):: g_t_tot
		real(8),dimension(3,3):: g_c_tot
		real(8),dimension(3,3):: ah_t
		real(8),dimension(3,3):: ah_c
		real(8),dimension(3,3):: unitmat
		real(8),dimension(3,3):: f_h !force for the lattice equation of motion
		real(8),dimension(3,3):: eta_t
		real(8),dimension(3,3):: eta_t_dt
		real(8),dimension(3,3):: sigma
		real(8) :: value_bond_constrain(numb_bond_cons)
		real(8) :: value_cos_constrain(numb_angl_cons)
		real(8):: crossp(3)
		real(8) :: s_t !thermostat variable at time t
		real(8) :: s_t_dt !thermostat variable at time t+ dt
		real(8) :: s_t_dot !time derivative of thermostat variable at time t
		real(8) :: s_c_dot !time derivative of thermostat variable at correction steps
		real(8) :: asvt ! aceleration dependent on the velocity over the thermostat at time t
		real(8) :: asvt1 ! aceleration dependent on the velocity over the thermostat at time t
		real(8) :: asvc ! aceleration dependent on the velocity over the thermostat at correction step
		real(8) :: as_t ! aceleration over the thermostat at time t
		real(8) :: as_c ! aceleration over the thermostat at correction steps
		real(8):: dt
		real(8):: ndf
		!real(8):: temp1
		!real(8):: temp2
		integer:: md_i,correc_step_i,  i, j, iat, a, b, c
		!Assign masses to each atom (for MD)
		x_t=xred_in !matrix with reduce coordinates
		dt=dtion_md !define time step
		h_t=latvec_in !latice matrix h=[a, b, c] in bohr
		h_t_dot= vel_lat_in
		call invertmat(h_t,h_t_inv,3)
		v_t=vel_in !copy velocity matrix [3, number_of_atoms]
		x_t_dot= x_t_dot_in!matmul(h_t_inv,v_t)
		s_t= s_in
		s_t_dot= s_in_dot
		ndf= 3.d0*nat
		
		unitmat=0.d0
		do i=1,3
			unitmat(i,i)=1.d0
		enddo
		pressure_md(:,:)= target_pressure_habohr*unitmat(:,:) !pressure matrix hartree/bohr^3
		call get_afx_t(fcart_in, h_t_inv, amu, nat, afx_t)
        call get_sigma(strten_in, sigma)
        call get_volume(h_t, vol_t)
        call get_eta(h_t, eta_t)
		call get_G_tot(h_t, h_t_dot, g_t_tot)
		call get_v_mat(amu, v_t, vol_t, nat, v_t_mat)  
		call get_asv(v_t, amu, nat, asvt)
        do md_i= 1, md_steps
			!write(*,*) 'md_steps  ',md_i
			ax_t= afx_t/(s_t**2.d0) - 2.d0*(s_t_dot/s_t)*x_t_dot - matmul(g_t_tot,x_t_dot)
			as_t= (asvt/s_t - (ndf+1.d0)*kb_HaK*temp/s_t)/Qmass
			ah_t= v_t_mat - sigma - pressure_md
			ah_t= matmul(ah_t,eta_t)/bmass
			
			x_t_dt= x_t + dt*x_t_dot + 0.5d0*dt*dt*ax_t
			s_t_dt= s_t + dt*s_t_dot + 0.5d0*dt*dt*as_t
			h_t_dt= h_t + dt*h_t_dot + 0.5d0*dt*dt*ah_t

			call get_volume(h_t_dt, vol_t_dt)
			call get_eta(h_t_dt, eta_t_dt)

			!prediction correction
			!prediction, calculate input for correction data
			!temp1= sum(x_t_dot*x_t_dot)
			!!write(*,*) 'x_t_dot 1', temp1
			x_c_dot= x_t_dot + dt*ax_t
			s_c_dot= s_t_dot + dt*as_t
			h_c_dot= h_t_dot + dt*ah_t
			call get_G_tot(h_t_dt, h_c_dot, g_c_tot)
			v_c= s_t_dt*(matmul(h_t_dt, x_c_dot))
			call get_v_mat(amu, v_c, vol_t_dt, nat, v_c_mat) 			 
			call get_asv(v_c, amu, nat, asvc)
			
			!correction steps
			do correc_step_i=1, correc_steps
				ax_c= afx_t/(s_t_dt**2.d0) - 2.d0*(s_c_dot/s_t_dt)*x_c_dot &
					& - matmul(g_c_tot,x_c_dot)
					
				as_c= (asvc/s_t_dt - (ndf+1.d0)*kb_HaK*temp/s_t_dt)/Qmass
					
				ah_c = v_c_mat - sigma - pressure_md
				ah_c= matmul(ah_c,eta_t_dt)/bmass
				x_c_dot= x_t_dot + 0.5d0*dt*(ax_c + ax_t)
				s_c_dot= s_t_dot + 0.5d0*dt*(as_c + as_t)
				h_c_dot= h_t_dot + 0.5d0*dt*(ah_c + ah_t)
				call get_G_tot(h_t_dt, h_c_dot, g_c_tot)
				v_c= s_t_dt*(matmul(h_t_dt, x_c_dot))
				call get_v_mat(amu, v_c, vol_t_dt, nat, v_c_mat)
				call get_asv(v_c, amu, nat, asvc) 
			end do !correction
			!begin constrain section
			!apply all the constrains over the cell
			call cell_constrains_volume(h_t_dt,dt, cell_para_valu, cell_angl_valu,&
						&volu_valu,cell_para_const,cell_angl_const,&
						&bool_cell_para_cons,bool_cell_angl_cons,bool_volu,&
						&nat,numb_cell_angl_cons, numb_cell_para_cons,volu_cons,&
						&h_t)
			!get atomic positions in catesian coordinates
			call get_r_t(h_t, x_t_dt, nat, r_t_dt)
			!apply all the atomic constrains
			!call atomic_constrains(r_t_dt,dt, bond_valu, angl_valu,&
			!			&bond_const,angl_const,&
			!			&bool_bond_cons,bool_angl_cons,&
			!			&nat,numb_angl_cons, numb_bond_cons,&
			!			&r_t)
			call atomic_constrains_atom_fix(r_t_dt,dt, bond_valu, angl_valu,atom_fix_valu,atom_fix_cord,&
						&bond_const,angl_const,atom_fix_const,&
						&bool_bond_cons,bool_angl_cons,bool_atom_fix_cons,&
						&nat,numb_angl_cons, numb_bond_cons,numb_atom_fix_cons,&
						&r_t)
			call atomic_constrains_evaluation(r_t,dt, bond_valu, angl_valu,&
						&atom_fix_valu,atom_fix_cord,&
						&bond_const,angl_const,atom_fix_const,&
						&bool_bond_cons,bool_angl_cons,bool_atom_fix_cons,&
						&nat,numb_angl_cons, numb_bond_cons,numb_atom_fix_cons,&
						&value_bond_constrain, value_cos_constrain)
			!get new reduced coordinates			
			call get_x_t(h_t, r_t, nat, x_t)
			!end constrain section
			!x_t= x_t_dt
			x_t_dot= x_c_dot
			h_t_dot= h_c_dot
			call invertmat(h_t,h_t_inv,3)
			s_t= s_t_dt
			s_t_dot= s_c_dot
			v_t= s_t*(matmul(h_t, x_t_dot))
			call get_afx_t(fcart_in, h_t_inv, amu, nat, afx_t)
			!call get_sigma(strten_in, sigma)
			call get_volume(h_t, vol_t)
			call get_eta(h_t, eta_t)
			call get_G_tot(h_t, h_t_dot, g_t_tot)
			call get_v_mat(amu, v_t, vol_t, nat, v_t_mat) !test get_v_mat(amu, v_t, vol_t, nat, v_t_mat)  
			pressure_out(md_i)= v_t_mat(1,1) + v_t_mat(2,2) + v_t_mat(3,3)
			call get_asv(v_t, amu, nat, asvt)
			volu_out(md_i)= vol_t
			bond_constrain_out(md_i,:)= value_bond_constrain
			cos_constrain_out(md_i,:)= value_cos_constrain
		end do !md
		s_out= s_t
		s_out_dot= s_t_dot
		!as_out= as_t
		!asvt_out= asvt
		h_out= h_t
		h_dot_out= h_t_dot
		!ah_out= ah_t
		x_out= x_t
		x_dot_out= x_t_dot
		!ax_out= ax_t
		v_out= v_t
		return
	end subroutine md_npt_constrains
	
	subroutine md_nvt_constrains(latvec_in,xred_in,fcart_in,&
						&vel_in, bond_valu, angl_valu,&
						& atom_fix_valu, atom_fix_cord,&
						&amu, Qmass, dtion_md, temp, s_in, s_in_dot, &
						&bond_const,angl_const,&
						&atom_fix_const, correc_steps,&
						&md_steps,bool_bond_cons,bool_angl_cons,&
						&bool_atom_fix_cons,&
						&nat,&
						&numb_angl_cons, numb_bond_cons,numb_atom_fix_cons,&
						&s_out, s_out_dot, x_out, v_out)
	!**********************************************************************************************
	! MOLECULAR DYNAMICS, USING THE LAGRANGIAN METHOD DEVELOPED BY PARRINELLO AND RAHMAN IN:
	![1]M. Parrinello and A. Rahman, J. Appl. Phys. J. Chem. Phys. 521, 2384 (1981).
	!THE CODE HERE, FOLLOWS THE DERIVATIONS AS PRESENTED BY AMSLER AND GOEDECKER IN:
	![2]M. Amsler and S. Goedecker, J. Chem. Phys. 133, (2010).
	!DESCRIPTION OF THE INPUT:
	!latvec_in: 3x3 matrix with the lattice vectors(a,b,c) in bohr, in the code latvec_in -> h, h(:,1)= a, h(:,2)= b
	!xred_in: matrix with the reduced positions of the atoms in the primitive lattice defined 
	!by h, in the code xred_in -> s with dimension (3, number of atoms)
	!fcart_in: matrix with the forces over every atom in hartree/bohr with dimension (3, number of atoms)
	!fcart_in is transformed into fs which is the force in reduced coordinates in hartree with dimension (3, number of atoms)
	!strten_in: the strain tensor in hartree/bohr^3, in the code strten_in-> sigma
	!vel_in: initial velocity of the atoms in bohr/second, in the code vel_in-> v
	!vel_lat_in: initial velocity of the lattice in the code vel_lat_in-> h_dot
	!g= (h^T)*h

		!use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat,ntime_md,&
		!		       &char_type,bmass,mdmin,dtion_md,strfact,units,usewf_md
		!use defs_basis
		!use utilities
		implicit none
		integer, intent(in) :: nat,&
								&numb_angl_cons, numb_bond_cons,numb_atom_fix_cons
		integer, intent(in) :: correc_steps !steps for correction alght
		integer, intent(in) :: md_steps !steps for md integration
		integer, intent(in) :: bool_bond_cons
		integer, intent(in) :: bool_angl_cons
		integer, intent(in) :: bool_atom_fix_cons
		integer, intent(in) :: bond_const(numb_bond_cons,2)
		integer, intent(in) :: angl_const(numb_angl_cons,3)
		integer, intent(in) :: atom_fix_const(numb_atom_fix_cons,1)
		real(8), intent(in) :: latvec_in(3,3)
		real(8), intent(in) :: xred_in(3,nat)
		real(8), intent(in) :: fcart_in(3,nat)
		real(8), intent(in) :: vel_in(3,nat)
		real(8), intent(in) :: bond_valu(numb_bond_cons)
		real(8), intent(in) :: angl_valu(numb_angl_cons)
		real(8), intent(in) :: atom_fix_valu(numb_atom_fix_cons,3)
		real(8), intent(in) :: atom_fix_cord(numb_atom_fix_cons,3)
		real(8), intent(in) :: amu(nat)
		real(8), intent(in) :: Qmass !mass of the thermostat
		real(8), intent(in) :: dtion_md
		real(8), intent(in) :: temp
		real(8), intent(in) :: s_in !thermostat variable at time t
		real(8), intent(in) :: s_in_dot !thermostat variable at time t+ dt
		real(8),dimension(3,nat), intent(out):: x_out !position at time t
		real(8),dimension(3,nat), intent(out):: v_out !position at time t
		real(8), intent(out) :: s_out !thermostat variable at time t
		real(8), intent(out) :: s_out_dot !thermostat variable at time t+ dt
		!*********************************************************************
		!Variables for my MD part
		!real(8), parameter :: Ha_eV=27.21138386 ! 1 Hartree, in eV
		real(8) :: kb_HaK! Boltzmann constant in Ha/K
		!real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
		real(8):: m(nat) !mass of atoms
		real(8):: vol
		!integer :: nat= 2
		real(8),dimension(3,nat):: r_t !position at time t
		real(8),dimension(3,nat):: r_t_dt !position at time t+ dt
		real(8) :: s_t !thermostat variable at time t
		real(8) :: s_t_dt !thermostat variable at time t+ dt
		real(8) :: s_t_dot !time derivative of thermostat variable at time t
		real(8),dimension(3,nat):: af_t !aceleration force at time t F/m
		real(8),dimension(3,nat):: af_t_dt !aceleration force at time t + dt
		real(8),dimension(3,nat):: ar_t !atomic aceleration at time t
		real(8),dimension(3,nat):: ar_c !atomic aceleration at corrections steps
		real(8) :: asvt ! aceleration dependent on the velocity over the thermostat at time t
		real(8) :: asvc ! aceleration dependent on the velocity over the thermostat at correction step
		real(8) :: as_t ! aceleration over the thermostat at time t
		real(8) :: as_c ! aceleration over the thermostat at correction steps
		real(8),dimension(3,nat):: v_t !velocity at time t
		real(8),dimension(3,nat):: v_c !velocity at correction steps
		real(8) :: s_c_dot !time derivative of thermostat variable at correction steps
		real(8) :: ndf ! nomber of degrees of freesom
		


		real(8):: dt
		integer:: md_i, correc_step_i, iat
		!Assign masses to each atom (for MD)

		!start MD parameters
		kb_HaK=((8.617343d-5)/27.21138386)
		r_t=matmul(latvec_in, xred_in)
		v_t= vel_in
		s_t= s_in
		s_t_dot= s_in_dot
		ndf= 3.d0*nat
		dt=dtion_md !define time step
		call get_af_t(fcart_in, amu, nat, af_t)
		call get_asv(v_t, amu, nat, asvt) 
        do md_i= 1, md_steps
			ar_t= af_t/(s_t**2.0) - 2.0*(s_t_dot/s_t**2.0)*v_t
			as_t= (asvt/s_t - (ndf+1)*kb_HaK*temp/s_t)/Qmass
			r_t_dt= r_t + dt*(v_t/s_t) + 0.5*dt*dt*ar_t
			s_t_dt= s_t + dt*s_t_dot + 0.5*dt*dt*as_t
			
			!prediction correction 
			!prediction, calculate input for correction data
			v_c= s_t_dt*v_t/s_t + dt*s_t_dt*ar_t
			s_c_dot= s_t_dot + dt*as_t
			call get_asv(v_c, amu, nat, asvc)
			
			!correction steps
			do correc_step_i=1, correc_steps
				ar_c= af_t/(s_t_dt**2.0) - 2.0*(s_c_dot/s_t_dt**2.0)*v_c
				as_c= (asvc/s_t_dt - (ndf+1)*kb_HaK*temp/s_t_dt)/Qmass    
				v_c= v_t + 0.5*dt*(ar_c + ar_t)
				s_c_dot= s_t_dot + 0.5*dt*(as_c + as_t)
				call get_asv(v_c, amu, nat, asvc) 
			end do !correction
			call atomic_constrains_atom_fix(r_t_dt,dt, bond_valu, angl_valu,atom_fix_valu,atom_fix_cord,&
						&bond_const,angl_const,atom_fix_const,&
						&bool_bond_cons,bool_angl_cons,bool_atom_fix_cons,&
						&nat,numb_angl_cons, numb_bond_cons,numb_atom_fix_cons,&
						&r_t)
			!r_t= r_t_dt
			s_t= s_t_dt
			v_t= v_c
			s_t_dot= s_c_dot
		end do !md steps
		s_out= s_t
		s_out_dot= s_t_dot
		call get_x_t(latvec_in, r_t, nat, x_out)
		!r_out= r_t
		v_out= v_t
		return
	end subroutine md_nvt_constrains 
	
!FUNCTIONS FOR TESTING
	subroutine testing_ax(afx_t, v_t,h_t, g_t_tot, s_t, s_t_dot, nat, ax)
		implicit none
		integer, intent(in) :: nat
		real(8),intent(in) :: afx_t(3,nat)
		real(8),intent(in) :: v_t(3,nat)
		real(8),intent(in) :: g_t_tot(3,3)
		real(8),intent(in) :: h_t(3,3)
		real(8),intent(in) :: s_t
		real(8),intent(in) :: s_t_dot
		real(8),intent(out) :: ax(3,nat)
		real(8) :: x_t_dot(3,nat)
		real(8) :: h_t_inv(3,3)
		call invertmat(h_t,h_t_inv,3)
		x_t_dot=matmul(h_t_inv,v_t)
		ax= afx_t/(s_t**2.0) - 2.0*(s_t_dot/s_t)*x_t_dot &
			& - matmul(g_t_tot,x_t_dot)
		return
	end subroutine testing_ax

	subroutine testing_ah(vol_t, bmass, target_pressure_habohr,&
               &  sigma, v_t_mat, eta_t, nat, ah_t)
		implicit none
		integer,intent(in) :: nat
		real(8),intent(in) :: vol_t
		real(8),intent(in) :: bmass
		real(8),intent(in) :: target_pressure_habohr
		real(8),intent(in) :: sigma(3,3)
		real(8),intent(in) :: v_t_mat(3,3)
		real(8),intent(in) :: eta_t(3,3)
		real(8),intent(out) :: ah_t(3,3)
		real(8), parameter :: kb_HaK=8.617343d-5/27.21138386d0
		real(8):: pressure_md(3,3)
		real(8):: unitmat(3,3)
		integer :: i
		unitmat=0.d0
		do i=1,3
			unitmat(i,i)=1.d0
		enddo
		pressure_md(:,:)= target_pressure_habohr*unitmat(:,:)
		ah_t= (v_t_mat)-sigma - pressure_md
		ah_t= matmul(ah_t,eta_t)/bmass
		return
	end subroutine testing_ah

	subroutine testing_as(Qmass, temp, s_t, asvt, nat, as_t)
		implicit none
		integer,intent(in) :: nat
		real(8),intent(in) :: Qmass 
		real(8),intent(in) :: temp
		real(8),intent(in) :: s_t
		real(8),intent(in) :: asvt
		real(8),intent(out) :: as_t
		real(8):: ndf
		real(8), parameter :: kb_HaK=8.617343d-5/27.21138386d0
		ndf= 3.d0*nat
		as_t= (asvt/s_t - (ndf+1)*kb_HaK*temp/s_t)/Qmass
		return
	end subroutine testing_as
end module npt_MD_suite
