program test_oopeos
	use comp_list_mod !which does use pureprop_mod
	use pengrob_mod
	use cub_eos_mod
	use ice_mod
	use MR_vdw1f_symmetric_mod
	use hydrates_mod
	use hydrates_cond_mod

	use Gex_based_phase_model_mod
	use Gex_model_mod
	use NRTL_mod
	use pur_liq_sat_model_mod
	
	implicit none
	
	call sub1
	call sub2
	print*, 'end program'
	contains
	
	subroutine sub1
	implicit none
		character(100) :: string, phase_file_name, paused
		integer :: main_option, phase_file_number
		class(phase_model_c), pointer :: obj, obj_ice, obj_pur_mr
		class(pur_cub_eos_par_c), pointer :: pur_obj_par
		class(hydrate_model_c), pointer :: obj_hyd
		class(hydrate_cond_c), pointer :: obj_hyd_cond
		integer :: i, fnum, pur_index
		real(8) :: real1, real2 !acho que no escopo dessa sub2 esse r sobrescreve o R global da complist
		real(8), allocatable :: v(:)
		string = 'input/oopeostest/comp_list.dat'
		call load_comp_list(string)
		print*, 'compilou'
		allocate(pur_obj_par,source=pengrob_i())
		print*, 'alocou puro'
		!vamos testar::
		select type(pur_obj_par)
		type is (pengrob_c)
			call pur_obj_par%debug_pur_ceos_par
			!ou call debug_PR(pur_obj_par) !mesma coisa
			real1 = 273.15 !K
			call pur_obj_par%calc_T_dep(r)
			call pur_obj_par%debug_pur_ceos_par
		end select
		phase_file_name = 'input/oopeostest/Lw.dat'
		open(newunit=phase_file_number,file=trim(phase_file_name),status='old',action='read')
		allocate(obj,source=MR_vdw1f_symmetric_i(phase_file_number,pengrob_i()))
		close(phase_file_number)
		print*, 'alocou mistura'
		!vamos testar::
		select type(obj)
		type is (MR_vdw1f_symmetric)
			call obj%pur_cub_eos_par_obj%debug_pur_ceos_par
			real1 = 273.15 !K
			call obj%pur_cub_eos_par_obj%calc_T_dep(r)
			call obj%pur_cub_eos_par_obj%debug_pur_ceos_par
			call obj%debug_phase_model
			!inicialização ok
			!calcular phi (e pobjos intermediarios)
			real1 = 273.15 !K
			real2 = 50.d5 !Pa
			call obj%calc_MixingParam()
			call obj%debug_phase_model
			call obj%calc_volume(real1, real2)
			call obj%debug_phase_model
			call obj%calc_phi(real1, real2)
			call obj%debug_phase_model
			open(unit=fnum,file='output/plotPvsV',status='replace',action='write')
				call obj%plot_PvsV(r,fnum)
			close(fnum)
		end select
!ceos_pure_c
		pur_index = 2
		allocate(obj_pur_mr,source=ceos_pure_i(pengrob_i()))

		select type(obj_pur_mr)
		type is (ceos_pure_c)
			call obj_pur_mr%set_ceos_pure_index(pur_index)
		
			call obj_pur_mr%pur_cub_eos_par_obj%debug_pur_ceos_par
			real1 = 273.15 !K
			call obj_pur_mr%pur_cub_eos_par_obj%calc_T_dep(r)
			call obj_pur_mr%pur_cub_eos_par_obj%debug_pur_ceos_par
print*, '...'
			call obj_pur_mr%debug_phase_model
!print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused
			!inicialização ok
			!calcular phi (e pobjos intermediarios)
			real1 = 273.15 !K
			real2 = 50.d5 !Pa
			obj_pur_mr%Theta_m = obj_pur_mr%mix_Theta() !(z%x,z%thetai)
			obj_pur_mr%b_m = obj_pur_mr%mix_b() !(z%x,z%bi)
!print*, ''
			call obj_pur_mr%debug_phase_model
!print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused
			call obj_pur_mr%calc_volume(real1, real2)
print*, '...'
			call obj_pur_mr%debug_phase_model
!print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused
			call obj_pur_mr%calc_phi(real1, real2)
print*, '...'
			call obj_pur_mr%debug_phase_model
!print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused
			open(unit=fnum,file='output/plotPvsV',status='replace',action='write')
			call obj_pur_mr%plot_PvsV(r,fnum)
		end select
		
		!compara
print*,'...'; print*, 'comparação do puro'
		obj%x = (/0.d0,1.d0,0.d0/)
		call obj%calc_phi(real1, real2)
		call obj%debug_phase_model
!print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused

!ICE
		phase_file_name = 'input/oopeostest/Ice.dat'
		open(newunit=phase_file_number,file=trim(phase_file_name),status='old',action='read')
		allocate(obj_ice,source=ice_sublim_i(phase_file_number))
		close(phase_file_number)
		call obj_ice%debug_phase_model
		!print*, 'alocou gelo'
		select type(obj_ice)
		type is (ice_sublim_c)
			real1 = 273.15 !K
			real2 = 50.d5 !Pa
			call obj_ice%calc_phi(real1, real2)
			call obj_ice%debug_phase_model
		end select
		
!HYD
		phase_file_name = 'input/oopeostest/Hyd.dat'
		open(newunit=phase_file_number,file=trim(phase_file_name),status='old',action='read')
		allocate(obj_hyd,source=hydrate_model_i(phase_file_number))
		close(phase_file_number)
		call obj_hyd%debug_phase_model
		!print*, 'alocou hyd_basico'
		select type(obj_hyd)
		type is (hydrate_model_c)
			real1 = 273.15 !K
			real2 = 50.d5 !Pa
!			call obj_hyd%calc_phi(real1, real2) !DEFEreal2ED
			call obj_hyd%debug_phase_model
		end select
		
!HYDcond
		phase_file_name = 'input/oopeostest/Hyd.dat'
		open(newunit=phase_file_number,file=trim(phase_file_name),status='old',action='read')
		allocate(obj_hyd_cond,source=hydrate_cond_i(phase_file_number))
		close(phase_file_number)
		call obj_hyd_cond%debug_phase_model
		!print*, 'alocou hyd_basico'
		select type(obj_hyd_cond)
		type is (hydrate_cond_c)
			real1 = 273.15 !K
			real2 = 50.d5 !Pa
			call obj_hyd_cond%calc_phi(real1, real2)
			call obj_hyd_cond%debug_phase_model
		end select
		
		
		
		!print*, 'alocou hyd_condensed'
		!print*, 'alocou hyd_sublimation'
	end subroutine sub1

	subroutine sub2
		
		character(100) :: string, phase_file_name, paused
		integer :: main_option, phase_file_number
		
		class(phase_model_c), pointer :: Gex_based_phase_model_obj
		class(Gex_model_c), pointer :: Gex_model_obj
		class(NRTL_c), pointer :: NRTL_obj
		class(pur_liq_sat_model_c), pointer :: pur_liq_sat_model_obj
		
		integer :: i, fnum, pur_index

		real(8), allocatable :: x_in
		real(8) :: T_in, P_in

		string = 'input/oopeostest/comp_list.dat'
		call load_comp_list(string)

		allocate(Gex_based_phase_model_obj,source=Gex_based_phase_model_i())
		call Gex_based_phase_model_obj%debug_phase_model

		allocate(Gex_model_obj,source=Gex_model_i())
		call Gex_model_obj%debug_Gex_model

		allocate(NRTL_obj,source=NRTL_i())
		call NRTL_obj%debug_Gex_model

		allocate(pur_liq_sat_model_obj,source=pur_liq_sat_model_i())
		call pur_liq_sat_model_obj%debug_phase_model

		nullify(Gex_based_phase_model_obj) !let it leak

		open(newunit=fnum,file=trim('input/oopeostest/phases/phase.dat'),status='old',action='read')
		allocate(Gex_based_phase_model_obj,source=Gex_based_phase_model_i(fnum,NRTL_i()))
		close(fnum)
		call Gex_based_phase_model_obj%debug_phase_model

!***todos os testes ok
	!faltam calculos: parametros de psat(T), parametros de volume(T,P)
	
	T_in = 273.15d0
	P_in = 1.0d5
	
	select type (Gex_based_phase_model_obj)
	type is (Gex_based_phase_model_c)
		call Gex_based_phase_model_obj%calc_Volume_L(T_in,P_in)
		print*,;	print*,;	print*,;	print*,;	print*,;	print*,;	print*,;	print*,;
		call Gex_based_phase_model_obj%debug_phase_model
!		print*,;print*, 'paused at ', __FILE__, __LINE__
!		read(*,*) paused
	end select
	
	call Gex_based_phase_model_obj%calc_phi(T_in,P_in)
	
	print*,;	print*,;	print*,;	print*,;	print*,;	print*,;	print*,;	print*,;
	call Gex_based_phase_model_obj%debug_phase_model
	
	end subroutine sub2
	
end program

