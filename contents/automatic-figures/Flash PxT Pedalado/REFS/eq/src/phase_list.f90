module phase_list_mod
	use comp_list_mod !which does use pureprop_mod
	use pengrob_mod
	use cub_eos_mod
	use ice_mod
	use MR_vdw1f_symmetric_mod
	use hydrates_mod
	use hydrates_cond_mod
!	use hydrates_sublim_mod
	use Gex_based_phase_model_mod
	use Gex_model_mod
	use NRTL_mod
	use pur_liq_sat_model_mod

	implicit none
	!
	integer :: nphas0
	character(2) :: nphas_str
	character(100) :: sim_case! declarado aqui, alimentado no "main"
	!
	type phase_list_c
		class(phase_model_c), pointer :: phase !pointer para iniciar com comando allocate
	end type phase_list_c
	!
	type(phase_list_c), pointer :: phase_list(:) ! a lista está criada aqui porque só se espera ter uma lista mesmo
	!														!pointer para ser uma lista
	contains
	!
	subroutine load_phase_list(phase_list_file)
		character(100) :: phase_list_file
		character(20), allocatable :: phase_model_name(:), phase_id(:)
		character(100), allocatable :: phase_files(:)
		integer :: i, load_phase_list_fnum, phase_model_fnum
		!implementation
		open(newunit=load_phase_list_fnum,file=trim(phase_list_file),status='old',action='read')
		read(load_phase_list_fnum,*); read(load_phase_list_fnum,*) nphas0
		write(nphas_str,'(I2.1)') nphas0
		allocate(phase_id(nphas0))
		allocate(phase_list(nphas0))
		allocate(phase_files(nphas0))
		allocate(phase_model_name(nphas0))
		read(load_phase_list_fnum,*) !cabeçalho
		do i=1,nphas0;
			read(load_phase_list_fnum,*) phase_id(i), phase_model_name(i), phase_files(i)
			phase_files(i) = 'input/'//trim(sim_case)//'/phases/'//trim(phase_files(i))
		end do
		close(load_phase_list_fnum)
		!montar fases
		do i = 1, nphas0
			open(newunit=phase_model_fnum,file=trim(phase_files(i)),status='old',action='read')
			!estrutura de arquivo de cada fase
			select case(trim(phase_model_name(i)))
			case('vdw1f_sym_pr')
				print*, i, ' phase is ',phase_model_name(i)
				allocate(phase_list(i)%phase, source=MR_vdw1f_symmetric_i(phase_model_fnum,pengrob_i()))
			case('ice')
				print*, i, ' phase is ',phase_model_name(i)
				allocate(phase_list(i)%phase, source=ice_sublim_i(phase_model_fnum))
			case('hyd_cond')
				print*, i, ' phase is ',phase_model_name(i)
				allocate(phase_list(i)%phase, source=hydrate_cond_i(phase_model_fnum))
!			case('Hyd_sublim')
!				print*, i, ' phase is ',phase_model_name(i)
!				allocate(phase_list(i)%phase, source=hydrate_sublim_i(phase_file_number))
			case('NRTL')
				print*, i, ' phase is ',phase_model_name(i)
				allocate(phase_list(i)%phase, source=Gex_based_phase_model_i(phase_model_fnum,NRTL_i()))
			case default
				print*, 'unexpected case argument in ', __FILE__, __LINE__; stop
			end select
			close(phase_model_fnum)
		end do
!		print*, 'end subroutine load_phase_list'; !pause
	end subroutine load_phase_list
!
end module phase_list_mod
