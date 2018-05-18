! Author: iuri soter viana segtovich
!
!Created on October 11, 2013, 10:05 PM
!
program controle
!
	use multiflash_mod
!	
	implicit none
!	
	character(100) :: string
	character(20) :: intstr, strtemp
	integer :: main_option
	class(phase_model_c), pointer :: ass
	integer :: i, j, fnum_in, i1,i2,i3
	real(8) :: r8temp
	
	!INICIALIZAÇÃO DO PROGRAMA DE CONTROLE
	open(unit=10,file='input/select_sim_case',status='old',action='read')
		read(10,*) sim_case !nome da pasta com intruções para simulação
	close(10)
	
	open(unit=10,file='input/'//trim(sim_case)//'/acompanhamento.dat',status='old',action='read')
		read(10,*) acompanhamentoscreen, acompanhamentolog, acompanhamentolog2
	close(10)
	
	if (acompanhamentolog == 1) then
		open(unit=20,file='output/acompanhamento.dat',status='replace',action='write')
	end if
	
	if (acompanhamentolog2 == 1) then
		open(unit=21,file='output/it_log.dat',status='replace',action='write')
		open(unit=22,file='output/echo_log.dat',status='replace',action='write')
	end if

	string = 'input/'//trim(sim_case)//'/main.dat'
	open(unit=10,file=trim(string),status='old',action='read')
		read(10,*) main_option
	close(10)
	
	select case(main_option)
	
	case(0) !rascunho
			print*, 'rascunho'
			call rascunho
		
	case(1)
		print*, 'multiflash'
		
		string = 'input/'//trim(sim_case)//'/comp_list.dat'
		call load_comp_list(string)
		
		string = 'input/'//trim(sim_case)//'/phase_list.dat'
		call load_phase_list(string)
		
		string = 'input/'//trim(sim_case)//'/flash.dat'
		call load_flash(string)
		
		print*, 'INITIALIZATION COMPLETE'
		

		call multiphaseflash !calc
		call multiphaseflash_report !result
!		print*, Tflash, Pflash !/1.d5 !com todas as casas decimais
		do !infinito
	!		stop
	!		pause
			print*, '(1) RESS COLLAPSED PHASES AND RUN'
			print*, '(2) CHANGE SPEC RUN'
			print*, '(3) CHANGE CONDITION RUN'
			print*, '(4) Plot PvsV at solution'
			print*, '(5) Switch phases order'
			print*, '(0) Stop'
			read(*,*) paused
			if( trim(paused) == '1' ) then
				print*, 'ress'
				call ressurrect_phases
				print*, 'calc'
				call multiphaseflash !calc
				print*, 'report'
				call multiphaseflash_report !result
!				print*, Tflash, Pflash !/1.d5 !com todas as casas decimais
			elseif( trim(paused) == '11' ) then
				print*, 'calc'
				call multiphaseflash !calc
				print*, 'report'
				call multiphaseflash_report !result
			elseif( trim(paused) == '2' ) then !CHANGE SPEC
				print*, '(1) P_.T_'
				print*, '(2) P_.F1'
				print*, '(3) T_.F1'
				print*, '(4) F1.F2'
				print*, '(5) call diagram'
				print*, '(6) call binary seq hyd'
				read(*,*) paused
				if( trim(paused) == '1' ) then
					FLASH_SPEC(1) = 'P_'
					FLASH_SPEC(2) = 'T_'
!TRY1, TO AVOID SINGULAR JAC WHEN SETTING A FLASHPT AFTER A FLASH_B OR FLASH_BB, SEEMS FINE
!REMOVES INCIPIENT PHASE AND POSSIBLY OTHERS, LEAVING A MAX OF PRESENT PHASE EQUAL TO NCOMP, WHICH IS THE MAXIMUM FOR GL .GE. 2
					stab_var(1:nphas-ncomp) = .1d0 !arbitrary positive value ***try too_small later
					phase_frac(1:nphas-ncomp) = 0.d0 !absent phase
					phase_flick(1:nphas-ncomp) = 0
!DEBUG END
				elseif( trim(paused) == '2' .or. trim(paused) == '3' ) then
!TRY1, TO AVOID SINGULAR JAC WHEN SETTING A FLASH_B AFTER A FLASH_BB, SEEMS FINE.
!REMOVES 2ND INCIPIENT PHASE AND POSSIBLY OTHERS, LEAVING A MAX OF PRESENT PHASE EQUAL TO NCOMP+1, WHICH IS THE MAXIMUM FOR GL .GE. 1
					stab_var(2:nphas-ncomp) = .1d0 !arbitrary positive value ***try too_small later
					phase_frac(2:nphas-ncomp) = 0.d0 !absent phase
					phase_flick(2:nphas-ncomp) = 0 !absent phase
					if( trim(paused) == '2' ) then
						FLASH_SPEC(1) = 'P_'
						FLASH_SPEC(2) = 'F1'
					elseif( trim(paused) == '3' ) then
						FLASH_SPEC(1) = 'T_'
						FLASH_SPEC(2) = 'F1'
					endif
				elseif( trim(paused) == '4' ) then
					FLASH_SPEC(1) = 'F1'
					FLASH_SPEC(2) = 'F2'
				elseif( trim(paused) == '5' ) then
					if (FLASH_SPEC(2) .ne. 'F1') then
						print*, 'for a diagram calculation, flash spec must be T_F1 or P_F1'
					elseif (FLASH_SPEC(1) == 'T_') then
						print*, 'Tlim[K]'
						read(*,*) Tlim
						call diagram_loop_extendedNR
								call multiphaseflash_report !result
					elseif(FLASH_SPEC(1) == 'P_') then
						print*, 'Plim[Pa]'
						read(*,*) Plim
						call diagram_loop_extendedNR
								call multiphaseflash_report !result
					endif
!					stop
				elseif( trim(paused) == '6' ) then
					call binary_sequential_hyd
					stop
				endif
!				call ressurrect_phases
!				call multiphaseflash !calc
!				call multiphaseflash_report !result
!			elseif( trim(paused) == '3' ) then !DISTURB CONDITION and retry
			elseif( trim(paused) == '3' ) then !DISTURB CONDITION and retry
				print*, '("1."(?)) :. "P=Px(1+"\1".d-2)"'
				print*, '("2."(?)) :. "P=P/(1+"\1".d-2)"'
				print*, '("3."(?)) :. "T=T+"\1".d-1"'
				print*, '("4."(?)) :. "T=T-"\1".d-1"'
				read(*,*) paused
				
				if (len(trim(paused)) == 1) then
					paused = trim(paused)//'.1'
				end if
				strtemp=paused(3:len(trim(paused)))
!				print*, 'pause= ', paused, 'strtemp= ', strtemp, __FILE__, __LINE__; stop
				read(strtemp,*) r8temp; print*, r8temp
!				print*, 'r8temp', r8temp, __FILE__, __LINE__; stop
				if( paused(1:1) == '1' ) then
!					Pflash = Pflash*1.01d0

					Pflash = Pflash*(1+r8temp*1.d-2)
				elseif( paused(1:1) == '2' ) then
					Pflash = Pflash/(1+r8temp*1.d-2)
				elseif( paused(1:1) == '3' ) then
					Tflash = Tflash + r8temp*1.d-1
				elseif( paused(1:1) == '4' ) then
					Tflash = Tflash - r8temp*1.d-1
				end if
!				call ressurrect_phases
!				call multiphaseflash !calc
!				call multiphaseflash_report !result
!				print*, Tflash, Pflash !/1.d5 !com todas as casas decimais
			elseif( trim(paused) == '4' ) then
			print*, 'plotPvsV'
			do j = 1, nphas
				ass => phase_list(j)%phase
				select type(ass)
				class is (cub_eos_c)
					write(intstr,'(I1.1)') j
					open(newunit=fnum_in,file='output/'//'plot_PvsV'//trim(intstr)//'.dat',status='replace',action='write')
		!escrever cabeçalho identificando a fase
					write(fnum_in,'(A,I1.1,A)') 'cubic_eos isoplete for solution composition for phase ',j,' from input phase list'
					write(fnum_in,*) ' '
		!LITERAL
					do i = 1, ncomp
						write(fnum_in,'(A,I1.1,A)', advance = 'NO') '    z_comp(',i,')'
					end do
					write(fnum_in,'(A)', advance = 'NO') ' loop_x'
					write(fnum_in,'(2(A))', advance = 'NO') '       Tflash', '       Pflash'

						!essa é pra conferir:
						write(fnum_in,'(A,I1.1,A)', advance = 'NO') '      cond(',j,')'
						write(fnum_in,'(2(A,I1.1,A))', advance = 'NO') '   stabvar(',j,')', '  phasfrac(',j,')'
						do i = 1, ncomp
							write(fnum_in,'(A,I2.2,A,I1.1,A)', advance = 'NO') '      x(',i,';',j,')'
						end do !i = 1, ncomp

					write(fnum_in,'(2(A))', advance = 'NO') '      conv_FE', '       conv_x'
					write(fnum_in,'(5(A))', advance = 'NO') '  phfr<0', '   swref', '  stvr<0', '   toosm', '  collps'
					write(fnum_in,'(3(A))', advance = 'NO') '         x_rf', '        NR_rf', '       MI_aux'
					write(fnum_in,*)
!NUMÉRICO
					do i = 1, ncomp
						write(fnum_in,'(ES13.5)', advance = 'NO') z_comp(i) !se eu for dar open com nome de arquivo, o open deve ser no main, e o filenumber deve chegar aqui, pode ser variavel do modulo multiflash
					end do
					write(fnum_in,'(I7.1)', advance = 'NO') loop_x
					write(fnum_in,'(2(ES13.5))', advance = 'NO') Tflash, Pflash

						if( ph_index_inv(j) /= 0 ) then
							!essa é pra conferir:
							write(fnum_in,'(2(A))', advance = 'NO') '        ', adjustr(phase_list(ph_index(ph_index_inv(j)))%phase%phase_cond)
							write(fnum_in,'(2(ES13.5))', advance = 'NO') stab_var(ph_index_inv(j)), phase_frac(ph_index_inv(j))
							do i = 1, ncomp
								write(fnum_in,'(ES13.5)', advance = 'NO') phase_list(ph_index(ph_index_inv(j)))%phase%x(i)
							end do !i = 1, ncomp
						else
							!essa é pra conferir:
							write(fnum_in,'(A)', advance = 'NO') adjustr('            X')
							write(fnum_in,'(2(ES13.5))', advance = 'NO') 0.d0, 0.d0
							do i = 1, ncomp
								write(fnum_in,'(ES13.5)', advance = 'NO') 0.d0
							end do !i = 1, ncomp
						end if
					write(fnum_in,*) ' '
					
					write(fnum_in,*) ' '
!CALCULAR ISOPLETA
					call ass%Plot_PvsV(Tflash,fnum_in)
					
					close(fnum_in)
				end select
			end do !i = 1, nphas
			print*, 'plotPvsV concluded', __FILE__, __LINE__
			
			elseif( trim(paused) == '5' ) then
				print*, 'i_1st, \n i_2nd, \n i_ref'
				read(*,*) i1
				read(*,*) i2
				read(*,*) i3
				call switch_phases(i1,i2,i3)
			elseif( trim(paused) == '0' ) then
				print*, 'MULTIFLASH COMPLETE'
				stop
			end if
		end do
		
	case(2)
		print*, 'G_chart'
		
		string = 'input/'//trim(sim_case)//'/comp_list.dat'
		call load_comp_list(string)
		
		string = 'input/'//trim(sim_case)//'/phase_list.dat'
		call load_phase_list(string)
		
		print*, 'INITIALIZATION COMPLETE'
		
		call mixture_G_chart !tabela de phi de cada componente para binário
		
		print*, 'mixture_G_chart COMPLETE'
		
	case(3)
		print*, 'incipient_flash_old_deactivated'
		
!		string = 'input/'//trim(sim_case)//'/comp_list.dat'
!		call load_comp_list(string)
!		
!		string = 'input/'//trim(sim_case)//'/phase_list.dat'
!		call load_phase_list(string)
!		
!		string = 'input/'//trim(sim_case)//'/flash.dat'
!		call load_flash(string)
!		
!		print*, 'INITIALIZATION COMPLETE'
!		
!		string = 'input/'//trim(sim_case)//'/incipient.dat'
!		call load_incipient(string)
!		
!		call solve_for_incipient_phase
!		
!		print*, 'incipient complete'
		
	case(4)
		print*, 'diagram_loop_old_deactivated'
		
!		string = 'input/'//trim(sim_case)//'/comp_list.dat'
!		call load_comp_list(string)
!		
!		string = 'input/'//trim(sim_case)//'/phase_list.dat'
!		call load_phase_list(string)
!		
!		string = 'input/'//trim(sim_case)//'/flash.dat'
!		call load_flash(string)
!		
!		print*, 'INITIALIZATION COMPLETE'
!		
!		string = 'input/'//trim(sim_case)//'/incipient.dat'
!		call load_incipient(string)
!		
!		call ddiagram_loop
!		
!		print*, 'diagram complete, check fort.92'
!		
	case(5)
		print*, 'PR_test'
		
		string = 'input/'//trim(sim_case)//'/comp_list.dat'
		call load_comp_list(string)
		
		string = 'input/'//trim(sim_case)//'/phase_list.dat'
		call load_phase_list(string)
		
		call teste1 !teste do peng robinson
		
		stop
		
	case(6)
		string = 'input/'//trim(sim_case)//'/comp_list.dat'
		call load_comp_list(string)
		
		print*, 'nothing to do here ', __FILE__, __LINE__; stop
		
	case(7)
		print*, 'psatplot'
		call psatplot
		
	case(8)
		print*, 'diagram_loop extendedNR'
		
		string = 'input/'//trim(sim_case)//'/comp_list.dat'
		call load_comp_list(string)
		
		string = 'input/'//trim(sim_case)//'/phase_list.dat'
		call load_phase_list(string)
		
		string = 'input/'//trim(sim_case)//'/flash.dat'
		call load_flash(string)
		
		print*, 'INITIALIZATION COMPLETE'
		
		if( flash_spec(1) == 'P_' .and. flash_spec(2) == 'T_') then
			print*, 'choose a P_F1 or T_F1 flash spec for diagram'; stop
		elseif ( flash_spec(2) == 'F1' ) then
			if( flash_spec(1) == 'P_' .or. flash_spec(1) == 'T_' ) then
				!ok
			end if
		elseif ( flash_spec(1) == 'F1' .and. flash_spec(2) == 'F2' ) then
			print*, 'choose a P_F1 or T_F1 flash spec for diagram'; stop
		else
			print*, 'unexpected flash_spec at', __FILE__, __LINE__; stop
		end if
		
		call diagram_loop_extendedNR
		
		print*, 'diagram complete, check fort.92'
		
	case(9)
		print*, 'binary sequential calculation'
		
		string = 'input/'//trim(sim_case)//'/comp_list.dat'
		call load_comp_list(string)
	
		string = 'input/'//trim(sim_case)//'/phase_list.dat'
		call load_phase_list(string)
	
		string = 'input/'//trim(sim_case)//'/flash.dat'
		call load_flash(string)
	
		print*, 'INITIALIZATION COMPLETE'
		
		call binary_sequential_calculation
		
	case(10) !sequential simulation baseada em dadosexp.dat
		
	case(11) !sequential calculation para hidrato com dois formadores, variando a proporção entre os formadores em excesso de água, e dado T
		print*, 'binary sequential hyd'
		
		string = 'input/'//trim(sim_case)//'/comp_list.dat'
		call load_comp_list(string)
	
		string = 'input/'//trim(sim_case)//'/phase_list.dat'
		call load_phase_list(string)
	
		string = 'input/'//trim(sim_case)//'/flash.dat'
		call load_flash(string)
	
		print*, 'INITIALIZATION COMPLETE'
		
		call binary_sequential_hyd
		
	case default
		print*, 'select calculation function in input/'//trim(sim_case)//'//main.dat'; stop
	end select
	
	if (acompanhamentoscreen == 1) then
		print*, (time_2-time_1), 'time'
	end if
	if (acompanhamentolog == 1) then
		write(20,*) (time_2-time_1), 'time'
		close(20)
	end if
	if (acompanhamentolog2 == 1) then
		close(21)
		close(22)
	end if
	
	write(*,*) 
	write(*,*)  'RUN COMPLETE'

	contains

subroutine psatplot
	real(8) :: zerofunc(2), dzdp
	integer :: i, j
	Tflash = 390.15d0
	Pflash = 1.d5
	
	string = 'input/'//trim(sim_case)//'/comp_list.dat'
	call load_comp_list(string)
	
	ncomp = 1
	
	string = 'input/'//trim(sim_case)//'/phase_list_psat.dat'
	call load_phase_list(string)		!que as duas fases iniciais do arquivo sejam PR
	
	nphas = 2
	phase_list(1)%phase%x(1) = 1.d0
	phase_list(1)%phase%phase_cond = 'L'
	
	phase_list(2)%phase%x(1) = 1.d0
	phase_list(2)%phase%phase_cond = 'V'
	
	
	print*, 'INITIALIZATION COMPLETE'
	
	open(unit=10,file='output/psat.dat',status='replace',action='write')
	
	do i = 1, 120
		do j = 1, 1000
			!ponto central
			

			call phase_list(1)%phase%calc_phi(Tflash,Pflash)
			

			call phase_list(2)%phase%calc_phi(Tflash,Pflash)
			
			print*, Pflash, Tflash
			print*, dexp(phase_list(1)%phase%phi(1)), phase_list(2)%phase%phi(1)
			print*, phase_list(1)%phase%V, phase_list(2)%phase%V
!			read(*,*) Pflash
			zerofunc(1) = dabs(phase_list(1)%phase%phi(1)-phase_list(2)%phase%phi(1))
			print*, zerofunc(1)
			
!			pause
			
			if (zerofunc(1) < 1d-4) then
				write(10,*) Tflash, Pflash
				print*, 'convergence' !pause
				exit
			end if
			
			Pflash = Pflash + 1.d0 ![Pa]
			

			call phase_list(1)%phase%calc_phi(Tflash,Pflash)
			

			call phase_list(2)%phase%calc_phi(Tflash,Pflash)
			
			print*, Pflash, Tflash
			print*, phase_list(1)%phase%phi(1), phase_list(2)%phase%phi(1)
			print*, phase_list(1)%phase%V, phase_list(2)%phase%V
!			read(*,*) Pflash
			zerofunc(2) = dabs(phase_list(1)%phase%phi(1)-phase_list(2)%phase%phi(1))
			print*, zerofunc(2)
			
			Pflash = Pflash - 1.d0
			
			dzdp = (zerofunc(2)-zerofunc(1))/2.d0 ![Pa]
			
			Pflash = Pflash - .1d0*zerofunc(1)/dzdp
			if (j==999) stop
		end do !j = 1, 1000
		Tflash = Tflash - 1.d0 ![K]
	end do !i = 1, 100
	close(10)
end subroutine psatplot

subroutine teste1 !calcular phi arbitrariamente
real(8) :: T, P
T = 273.2d0
P = 2577007.6d0
print*, phase_list(1)%phase%phase_cond
!call phases_list(1)%phase%debug_phase_model
call phase_list(1)%phase%calc_phi(T,P)
print*, phase_list(1)%phase%V, 'V'
print*, phase_list(1)%phase%phi(:), 'phi 1' !acesso direto a variável
print*, P*phase_list(1)%phase%phi(:), 'fug'
!demais fases

end subroutine teste1

!!!!subroutine ddiagram_loop !PT

!!!!	integer :: i, npontos
!!!!	real(8) Tinicial, Tfinal, passoT
!!!!	real(8) Pinicial, Pfinal, passoP

!!!!!	!Loop para diagrama

!!!!!	!!Se acontecer um colapso de fase (do tipo raiz fase instável deixando de existir completamente), a fase colapsada tem que ser de índice superior à fase incipiente buscada

!!!!	npontos = 100
!!!!	do i = 1, npontos
!!!!		
!!!!		if (fv == 1) then		!!P GRID
!!!!			Pinicial = Pflash !a última pressão usada
!!!!			Pfinal = 10.d0*Pinicial
!!!!!			Pfinal = Pinicial/10.d0
!!!!			passoP = (Pfinal-Pinicial)/dfloat(npontos)
!!!!			call solve_for_incipient_phase
!!!!			write(92,*) Tflash, Pflash
!!!!		
!!!!		elseif( fv == 2 ) then	!	!!T GRID
!!!!			Tinicial = Tflash !a temperatura do flash de inicialização
!!!!			Tfinal = 2.d0*Tinicial
!!!!!			Tfinal = Tinicial/2.d0
!!!!			passoT = (Tfinal-Tinicial)/dfloat(npontos)
!!!!			call solve_for_incipient_phase
!!!!			write(92,*) Tflash, Pflash
!!!!		end if
!!!!		!pause
!!!!		!faz analise de sensibilidade pra ver se continua chutando T ou passa pra P
!!!!			!#########
!!!!			!
!!!!			!#########
!!!!			!
!!!!			!dá o passo
!!!!		if (fv == 1) then
!!!!			Pflash = Pflash + passoP
!!!!			Tflash= Tflash - 1.d0		!pra tentar garantir de jogar a estimativa inicial de T para dentro do envelope caso esteja ficando retrogrado
!!!!		elseif (fv == 2) then
!!!!			Tflash = Tflash + passoT
!!!!		end if
!!!!		
!!!!	end do !i = 1, npontos

!!!!end subroutine ddiagram_loop

subroutine diagram_loop_extendedNR

!	!Loop para diagrama

!	!!Se acontecer um colapso de fase (do tipo raiz fase instável deixando de existir completamente), a fase colapsada tem que ser de índice superior à fase incipiente buscada

!	VERSÃO NOVA 140624
	integer :: i, ii, j, npontos
	real(8) :: sign_
	real(8) :: passoT = .1d0 !2.d0 !2.d0 !0.1d0 !0.01d0
	real(8) :: passoP
	real(8) :: passoPrel = .01d0 !.1d0 !.005d0 !0.001d0
	npontos = 10000 !máximo

		if ( flash_spec(1) == 'P_' ) then
			if( Pflash < Plim ) then
				sign_ = 1.d0
			else
				sign_ = -1.d0
			end if
		elseif ( flash_spec(1) == 'T_' ) then
			if( Tflash < Tlim ) then
				sign_ = 1.d0
			else
				sign_ = -1.d0
			end if
		end if

!		write(92,'(2(A))', advance = 'NO') '       Tflash', '       Pflash'
!		do j = 1, nphas
!			write(92,'(A)', advance = 'NO') '   phase_cond'
!			write(92,'(A)', advance = 'NO') '      stab_var'
!			write(92,'(A)', advance = 'NO') '   phase_frac'
!			do ii = 1, ncomp
!				write(92,'(A)', advance = 'NO') '            x'
!			end do !i = 1, ncomp
!				write(92,'(A)', advance = 'NO') '            V'
!		end do !j = 1, nphas0
!		write(92,*)

	do i = 1, npontos
!		call ressurrect_phases
		call multiphaseflash
		if(loop_x .GE. 100) then !*** <100> equivale ao loop_x_max da sub multiflash
			print*, 'I guess it is ', '"Too many steps in loop_x"'
			return
		end if
		
!		!REPORT DIAGRAM POINT
!		write(92,'(2(ES13.5))', advance = 'NO') Tflash, Pflash
!		do j = 1, nphas
!			write(92,'(2(A))', advance = 'NO') '        ', adjustr(phase_list(ph_index(j))%phase%phase_cond)
!			write(92,'(A, ES12.5)', advance = 'NO') '  ', stab_var(j)
!			write(92,'(ES13.5)', advance = 'NO') phase_frac(j)
!			do ii = 1, ncomp
!				write(92,'(ES13.5)', advance = 'NO') phase_list(ph_index(j))%phase%x(ii)
!			end do !i = 1, ncomp
!				write(92,'(ES13.5)', advance = 'NO') phase_list(ph_index(j))%phase%V
!		end do !j = 1, nphas0
!		write(92,*)
		
		if ( flash_spec(1) == 'P_' ) then
			passoP = passoPrel*Pflash*sign_
			Pflash = Pflash + passoP
			if( (Pflash-Plim)*sign_ > 0.d0 ) then
				exit
			end if
!			Tflash = Tflash + !mantém valor do último cálculo
		elseif ( flash_spec(1) == 'T_' ) then
			Tflash = Tflash + passoT*sign_
			if( (Tflash-Tlim)*sign_ > 0.d0 ) then
				exit
			end if
		end if
		
	end do !i = 1, npontos

end subroutine diagram_loop_extendedNR

subroutine mixture_G_chart

!grid vars
real(8) :: fracz(2) !seria o w
real(8) :: phiL(2), phiV(2)
real(8) :: miL(2), miV(2)
real(8) :: GL, GV, G, GmixL, GmixV

!bulk for TPD
real(8) :: ztest1(2), ztest2(2), ztest3(2)
real(8) :: miLbulk1(2), miVbulk1(2), miLbulk2(2), miVbulk2(2), miLbulk3(2), miVbulk3(2)
real(8) :: TPD1L, TPD2L, TPD3L, TPD1V, TPD2V, TPD3V, dTPD1, SdTPD1 !para fase incipiente L ou então V

integer :: j
real(8) :: T, P

!pure components
real(8) :: G1L, G2L, G1V, G2V, G1, G2

T = 300.d0
P = 1.d5

SdTPD1 = 0.d0
SdTPD1 = 2.063747d1 + 6.249122d-1

write(101,*) 'fracz(1), phiL(1), phiL(2), phiV(1), phiV(2), miL(1), miL(2), miV(1), miV(2), GL, GV, GmixL, GmixV'
write(102,'(A)') 'fracz(1), TPD1L, TPD2L, TPD3L, TPD1V, TPD2V, TPD3V'

!puros
fracz = 1.d0
fracz(2) = 1.d0 - fracz(1)

phase_list(1)%phase%x(:) = fracz(:)

phase_list(1)%phase%phase_cond = 'L'

call phase_list(1)%phase%calc_phi(T,P)
phiL(:) = phase_list(1)%phase%phi(:)

phase_list(1)%phase%phase_cond = 'V'

call phase_list(1)%phase%calc_phi(T,P)
phiV(:) = phase_list(1)%phase%phi(:)

miL(1) = R*T*dlog(fracz(1)*phiL(1))

miV(1) = R*T*dlog(fracz(1)*phiV(1))

GL = miL(1)
GV = miV(1)

G1 = min(GL,GV)

!write(101,'(12(ES15.6))') fracz(1), phiL(1), phiL(2), phiV(1), phiV(2), miL(1), miL(2), miV(1), miV(2), GL, GV, G1

fracz = 0.d0
fracz(2) = 1.d0 - fracz(1)

phase_list(1)%phase%x(:) = fracz(:)

phase_list(1)%phase%phase_cond = 'L'

call phase_list(1)%phase%calc_phi(T,P)
phiL(:) = phase_list(1)%phase%phi(:)

phase_list(1)%phase%phase_cond = 'V'

call phase_list(1)%phase%calc_phi(T,P)
phiV(:) = phase_list(1)%phase%phi(:)

miL(2) = R*T*dlog(fracz(2)*phiL(2))

miV(2) = R*T*dlog(fracz(2)*phiV(2))

GL = miL(2)
GV = miV(2)

G2 = min(GL,GV)

!write(101,'(12(ES15.6))') fracz(1), phiL(1), phiL(2), phiV(1), phiV(2), miL(1), miL(2), miV(1), miV(2), GL, GV, G2

print*, 'G1 e G2 ok'

!bulk ztests for TPD

!sejam ztest os z bulk e fracz os w incipiente

!bulk1
ztest1(1) = .01d0
ztest1(2) = 1.d0 - ztest1(1)

!tem que ver se o menor G é GL ou GV pra decidir se pega o miL ou miV
!ou sempre que o GL for o menor é porque o miL e miV também são menores
!e aí seria só pegar o menor mi
!e ae?

phase_list(1)%phase%x(:) = ztest1(:)

phase_list(1)%phase%phase_cond = 'L'

call phase_list(1)%phase%calc_phi(T,P)
phiL(:) = phase_list(1)%phase%phi(:)

phase_list(1)%phase%phase_cond = 'V'

call phase_list(1)%phase%calc_phi(T,P)
phiV(:) = phase_list(1)%phase%phi(:)

miLbulk1(:) = R*T*dlog(ztest1(:)*phiL(:))

miVbulk1(:) = R*T*dlog(ztest1(:)*phiV(:))

!bulk 2
ztest2(1) = .25d0
ztest2(2) = 1.d0 - ztest2(1)

phase_list(1)%phase%x(:) = ztest2(:)

phase_list(1)%phase%phase_cond = 'L'

call phase_list(1)%phase%calc_phi(T,P)
phiL(:) = phase_list(1)%phase%phi(:)

phase_list(1)%phase%phase_cond = 'V'

call phase_list(1)%phase%calc_phi(T,P)
phiV(:) = phase_list(1)%phase%phi(:)

miLbulk2(:) = R*T*dlog(ztest2(:)*phiL(:))

miVbulk2(:) = R*T*dlog(ztest2(:)*phiV(:))

!bulk3
ztest3(1) = .75d0
ztest3(2) = 1.d0 - ztest3(1)

phase_list(1)%phase%x(:) = ztest3(:)

phase_list(1)%phase%phase_cond = 'L'

call phase_list(1)%phase%calc_phi(T,P)
phiL(:) = phase_list(1)%phase%phi(:)

phase_list(1)%phase%phase_cond = 'V'

call phase_list(1)%phase%calc_phi(T,P)
phiV(:) = phase_list(1)%phase%phi(:)

miLbulk3(:) = R*T*dlog(ztest3(:)*phiL(:))

miVbulk3(:) = R*T*dlog(ztest3(:)*phiV(:))

!grid

fracz = 0.d0
do j = 1, 999
fracz(1) = fracz(1) + 0.001d0 !0.001d0
fracz(2) = 1.d0 - fracz(1)

phase_list(1)%phase%x(:) = fracz(:)

phase_list(1)%phase%phase_cond = 'L'

call phase_list(1)%phase%calc_phi(T,P)
phiL(:) = phase_list(1)%phase%phi(:)

phase_list(1)%phase%phase_cond = 'V'

call phase_list(1)%phase%calc_phi(T,P)
phiV(:) = phase_list(1)%phase%phi(:)

miL(:) = R*T*dlog(fracz(:)*phiL(:))

miV(:) = R*T*dlog(fracz(:)*phiV(:))

GL = fracz(1)*miL(1)+fracz(2)*miL(2)
GV = fracz(1)*miV(1)+fracz(2)*miV(2)

GmixL = GL - fracz(1)*G1 - fracz(2)*G2

GmixV = GV - fracz(1)*G1 - fracz(2)*G2

!TPD (analytical)

!um TPD para dada fase bulk, para cada fase incipiente, supondo cada modelo possível para modelar ela
!nesse caso: TPD pra fase z bulk Liq e fase w incipiente Liq tb :: nenhuma fase estacionaria
!e se fase incipiente vapor:
TPD1L = fracz(1)*(miL(1)-miVbulk1(1)) + fracz(2)*(miL(2)-miVbulk1(2))		!SUPONDO BULK VAPOR
TPD1V = fracz(1)*(miV(1)-miVbulk1(1)) + fracz(2)*(miV(2)-miVbulk1(2))		!SUPONDO BULK VAPOR


TPD2L = fracz(1)*(miL(1)-miLbulk2(1)) + fracz(2)*(miL(2)-miLbulk2(2))		!SUPONDO BULK LIQ
TPD2V = fracz(1)*(miV(1)-miLbulk2(1)) + fracz(2)*(miV(2)-miLbulk2(2))		!SUPONDO BULK LIQ

TPD3L = fracz(1)*(miL(1)-miLbulk3(1)) + fracz(2)*(miL(2)-miLbulk3(2))		!SUPONDO BULK LIQ
TPD3V = fracz(1)*(miV(1)-miLbulk3(1)) + fracz(2)*(miV(2)-miLbulk3(2))		!SUPONDO BULK LIQ

write(101,'(16(ES15.6))') fracz(1), phiL(1), phiL(2), phiV(1), phiV(2), miL(1), miL(2), miV(1), miV(2), GL, GV, GmixL, GmixV

write(102,'(7(ES15.6))') fracz(1), TPD1L, TPD2L, TPD3L, TPD1V, TPD2V, TPD3V


end do

end subroutine mixture_G_chart

subroutine binary_sequential_calculation
	integer :: i
	class(phase_model_c), pointer :: ass

	ass => phase_list(1)%phase
	select type(ass)
	class is (MR_vdw1f)
	
	ass%kij(1,2) = -0.18
	ass%kij(2,1) = ass%kij(1,2)
	
	write(92,*) trim(comp_list(1)%comp%name_), '(1) ', trim(comp_list(2)%comp%name_), '(2) ', 'kij= ', ass%kij(1,2)
	write(92,*) 'Tflash, Pflash, x1, x2, y1, y2, z1, z2'
	end select
	
	z_comp = (/1.d0, 0.d0/)
	
	do i = 1, 101

	select type(ass)
	class is (MR_vdw1f)
	ass%kij(1,2) = -0.05d0-25.9d0+0.03459d0*Tflash+4801.4d0/Tflash
	ass%kij(2,1) = ass%kij(1,2)
	end select

		call multiphaseflash !calc
		write(92,*) Tflash, Pflash, phase_list(1)%phase%x(:), phase_list(2)%phase%x(:), z_comp(:)
		call multiphaseflash_report !result
!			pause
		z_comp(1) = z_comp(1) - 1.d-2
		z_comp(2) = 1.d0 - z_comp(1)
	end do !i = 1, 100
end subroutine binary_sequential_calculation

subroutine binary_sequential_hyd
	integer :: i
	
	do i = 1, 2001
		call multiphaseflash !calc
		!hyd
		!composição do vapor, pressão, fase hidrato, composição dela, pressão de novo(para a tabela do excel)
!		write(92,*) phase_list(5)%phase%x(:), Pflash, phase_list(1)%phase%phase_cond, phase_list(1)%phase%x(:), Pflash
		!fluid
		!composição do vapor, pressão, e fração de cada fase hidrato presente.
!		write(92,*) phase_list(1)%phase%x(:), Pflash, phase_list(2)%phase%phase_cond, phase_frac(2),phase_list(3)%phase%phase_cond, phase_frac(3)
		call multiphaseflash_report !result
!			print*, 'paused at', __FILE__, __LINE__; read(*,*) paused
		z_comp(2) = z_comp(2) + 0.00001d0
		z_comp(3) = z_comp(3) - 0.00001d0
	end do
end subroutine binary_sequential_hyd

subroutine rascunho
		!VARIÁVEIS RASCUNHO
		integer :: i
		Real(8) :: Tvec(12)

		string = 'input/'//trim(sim_case)//'/comp_list.dat'
		call load_comp_list(string)
		
		string = 'input/'//trim(sim_case)//'/phase_list.dat'
		call load_phase_list(string)
		
		string = 'input/'//trim(sim_case)//'/flash.dat'
		call load_flash(string)
		
		print*, 'INITIALIZATION COMPLETE'
		
!		Tvec = (/273.15000, 278.15000, 283.15000,293.15000,303.15000,313.15000,323.15000,333.15000,343.15000,353.15000,363.15000,373.15000/)
!		
!	do i=1,12
!			Tflash = Tvec(i)
!			call multiphaseflash !calc
!			print*, phase_list(1)%phase%V*1.d3
!	end do
		Tflash = 273.d0
	do i=1,100
		Tflash = Tflash-1.d0
		call multiphaseflash !calc
		print*, Tflash, Pflash, phase_list(1)%phase%phi(1)
	end do
end subroutine rascunho

subroutine rascunho2
		print*, 'multiflash'
		
		string = 'input/'//trim(sim_case)//'/comp_list.dat'
		call load_comp_list(string)
		
		string = 'input/'//trim(sim_case)//'/phase_list.dat'
		call load_phase_list(string)
		
		string = 'input/'//trim(sim_case)//'/flash.dat'
		call load_flash(string)
		
		print*, 'INITIALIZATION COMPLETE'
		
!		do
		call multiphaseflash !calc
		call multiphaseflash_report !result
		print*, Tflash, Pflash !/1.d5 !com todas as casas decimais
!		pause
		print*, 'ress collapsed phases and try again? (y/n)'; read(*,*) paused
		if( trim(paused) == 'y' ) then
			print*, 'ok'
			call ressurrect_phases
			call multiphaseflash !calc
			call multiphaseflash_report !result
			print*, Tflash, Pflash !/1.d5 !com todas as casas decimais
		end if
!		end do
		
		do 
!		print*, 'newT'
!		read(*,*) Tflash
!		print*, 'newP'
!		read(*,*) Pflash
		print*, 'newZ3'
		read(*,*) z_comp(3)
		z_comp(2) = 1.d0 - z_comp(1) - z_comp(3)
		call multiphaseflash !calc
		call multiphaseflash_report !result
		end do
		
		print*, 'MULTIFLASH COMPLETE'
end subroutine rascunho2

end program controle

