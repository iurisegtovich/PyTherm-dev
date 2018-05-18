module multiflash_mod
!versão pós defesa M.Sc.
!edição iniciada em 15.01.02
!version 15.01.07, tentativa de algoritmo bi-gen, funcional mas com problemas
!15.01.08 - refinamentos do bi-gen, parece que ficou muito bom, ver "anotaçoes.15.01.08.txt++"
!15.01.09 - desenvolvimento do NR-parcial
	use comp_list_mod !which does use pureprop_mod
	use phase_list_mod !which does use eos_mod
	use acompanhamento_mod
	!
	implicit none
	!
	integer :: nphas
	integer :: serial_report_fnum
	real(8) :: Tflash, Pflash !temperature [K], pressure [Pa]
	real(8) :: Tflash_spec, Pflash_spec, Tlim, Plim 												!MULTI-SPEC
	real(8), allocatable :: z_comp(:) !cada componente; composição global !global component molar composition in stream
	real(8), allocatable :: phase_frac(:), stab_var(:) !cada fase, fração molar da fase na corrente, e variavel de estabilidade
	integer, allocatable :: phase_flick(:) !0 means its equal to zero, therefore stab_var should be calculated, 1 means the phase is present therefore stab_var is zero, and phase_frac itslf should be calculated
	
	real(8), allocatable :: dist_coef(:,:) !cada componente; cada fase
	integer, allocatable :: ph_index(:), ph_index_inv(:)
	character(2) :: flash_spec(2) !elem. 1: P_ ou T_ ou F1; elem. 2: T_ ou F1 ou F2 	!MULTI-SPEC
								!para fazer o F1-F2 ambos T e Pdevem estar livres, logo tem que haver ainda mais uma linha e coluna no jacobiano ainda não programadas
									!no jac, em vez de ultimo elemento entram os 2x2 elementos do canto inferior diretio da matriz: dTdP;dTdT;dPdT; dPdP (0;1;0;1?)
									
!ITERATION COUNTERS
	integer :: loop_x
	integer :: loop_FE
	integer :: loop_FE_total = 0
!variaveis de convergencia do flash
	real(8) :: conv_x
	real(8) :: conv_FE
!NUMERICAL METHOD CONSTANTS/reguladores de convergência/parâmetros de tolerância numérica
	integer :: loop_x_max = 1000
	integer :: loop_FE_max = 1000
	real(8) :: conv_x_TOL = 1.d-8
	real(8) :: conv_FE_TOL = 1.d-8
	real(8) :: Too_Small = 1.d-9		
	real(8) :: GTzero_TOL = 1.d-9		!usado para 'maior que 0' ou 'menor ou igual a 0'
	real(8) :: LTzero_TOL = -1.d-9
	real(8) :: collapse_TOL = 1.d-4
	contains
	!
	subroutine load_flash(flash_file)
		!arguments
		character(100) :: flash_file
		!local
		integer :: i, j, jj
		!implementation
		nphas=nphas0
		open(unit=10,file=trim(flash_file),status='old',action='read')
!		read(10,*); read(10,*) Tflash 													!SINGLE-SPEC
!		read(10,*); read(10,*) Pflash 													!SINGLE-SPEC
		read(10,*); read(10,*) Tflash_spec, Tlim												!MULTI-SPEC
		read(10,*); read(10,*) Pflash_spec, Plim												!MULTI-SPEC
		Tflash = Tflash_spec 															!MULTI-SPEC
		Pflash = Pflash_spec 															!MULTI-SPEC
		allocate(z_comp(ncomp))
		allocate(phase_frac(nphas))
		allocate(stab_var(nphas))
		allocate(phase_flick(nphas))
		allocate(dist_coef(ncomp,nphas))
		allocate(ph_index(nphas))
		allocate(ph_index_inv(nphas))
! 		read(10,*); do i=1,ncomp; read(10,*) z_comp(i); enddo !EM COLUNA
		read(10,*); read(10,*) (z_comp(i),i=1,ncomp); ! EM LINHA
		read(10,*); read(10,*) (phase_frac(i), i=1,nphas)

		if( phase_frac(nphas) < 1.d-9 ) then !não permitir inicialização com referencia ausente
			phase_frac(nphas) = 1.d0
		endif
		
		phase_frac(:) = phase_frac(:)/sum(phase_frac)
		
		!inicialização das variaveis de estabilidade e de interruptor do NR-parcial.
		
		do j = 1, nphas
			if( phase_frac(j) < 1.d-9 ) then
				phase_frac(j) = 0.d0
				phase_flick(j) = 0
				stab_var(j) = too_small !1d-1 !valor positivo arbitrário, testar 1d-9 depois que o multiflash3 estiver pronto***
			else
!				phase_frac(j) igual ao que foi lido
				phase_flick(j) = 1
				stab_var(j) = 0.d0
			endif
		enddo !j = 1, nphas
		
		read(10,*); read(10,*) flash_spec(1); read(10,*) flash_spec(2);					!MULTI-SPEC
		close(10)
		if( flash_spec(1) == 'P_' .and. flash_spec(2) == 'T_') then
			!Flash-P-T
		elseif ( flash_spec(2) == 'F1' ) then
			if( flash_spec(1) == 'P_' .or. flash_spec(1) == 'T_' ) then
				if( nphas < 2 ) then
					print*, 'for a P_F1 or T_F1 flash, nphas must be .GE. 2'; stop
				endif !nphas < 2
			else
				print*, 'unexpected flash_spec at', __FILE__, __LINE__; stop
			endif
		elseif ( flash_spec(1) == 'F1' .and. flash_spec(2) == 'F2' ) then
			if( nphas < 3 ) then
				print*, 'for a P_F1 or T_F1 flash, nphas must be .GE. 3'; stop
			endif !nphas < 3
		else
			print*, 'unexpected flash_spec at', __FILE__, __LINE__; stop
		endif

		
		
!		!inicialização padrão das variáveis do flash
!		stab_var(:) = 0.d0 !todas as fases consideradas estão presentes
!		phase_frac(:) = 1.d0/nphas !fases igualmente distribuídas
!		!a fase de referência será sempre a última fase dos vetores
		!
		!inicializar o flash com as fases na ordem de inicialização
		do j = 1, nphas
			ph_index(j) = j
		enddo !j = 1, nphas
		
		open(newunit=serial_report_fnum, file='output/serial_report.dat', status='replace', action='write')
		
		do i = 1, ncomp
			write(serial_report_fnum,'(A,I1.1,A)', advance = 'NO') '    z_comp(',i,')'
		enddo
		write(serial_report_fnum,'(A)', advance = 'NO') ' loop_x'
		write(serial_report_fnum,'(2(A))', advance = 'NO') '       Tflash', '       Pflash'
		do j = 1, nphas0
			!essa é pra conferir:
			write(serial_report_fnum,'(A,I1.1,A)', advance = 'NO') '      cond(',j,')'
			write(serial_report_fnum,'(2(A,I1.1,A))', advance = 'NO') '   stabvar(',j,')', '  phasfrac(',j,')'
			do i = 1, ncomp
				write(serial_report_fnum,'(A,I2.2,A,I1.1,A)', advance = 'NO') '      x(',i,';',j,')'
			enddo !i = 1, ncomp
		enddo !j = 1, nphas0
		write(serial_report_fnum,'(2(A))', advance = 'NO') '      conv_FE', '       conv_x'
		write(serial_report_fnum,'(5(A))', advance = 'NO') '  phfr<0', '   swref', '  stvr<0', '   toosm', '  collps'
		write(serial_report_fnum,'(3(A))', advance = 'NO') '         x_rf', '        NR_rf', '       MI_aux'
		write(serial_report_fnum,*) ' '
		
	end subroutine load_flash
	!
	subroutine unload_flash
		!declaration
		
		!implementation
		close(serial_report_fnum)
		
	end subroutine unload_flash
	!
	subroutine multiphaseflash
		!DECLARATIONS
		!NO ARGUMENTS
		!
		!LOCAL
		real(8), allocatable :: x(:,:), x_ol(:,:) !pra poder jogar em forma de matriz na função update x, tendo que depois atualizar os x em cada objeto.
		real(8), allocatable :: phase_frac_ol(:), stab_var_ol(:) !cada fase
		real(8), allocatable :: phase_flick_ol(:) !cada fase
		!
		!variáveis livres para testes/depuração
		real(8) :: dummy(40), dummyummy(20,20)

		real(8) :: sum1

		integer :: i, j, jj, ref_index, re_set_theta_index
		!
		integer :: checkpoint = 10
		!
		character(100) :: format_
		!
		integer :: collapse_tag(2) !indice da fase que crescerá e índice da fase que será eliminada do sistema de equações
		!
		integer :: nphas_present
		!
		time_1 = 0.d0
		time_2 = 0.d0
		checkpoint = 10
		call cpu_time(time_1)
		!
		!implementation
		allocate(x(ncomp,nphas),x_ol(ncomp,nphas))
		allocate(phase_frac_ol(nphas))
		allocate(stab_var_ol(nphas))
		allocate(phase_flick_ol(nphas))

		!
!		inicializa variável local x a partir dos objetos fase
		do j = 1, nphas
			x(:,j) = phase_list(ph_index(j))%phase%x(:)
		enddo !j = 1, nphas
		!
		!ECHO
		if (ACOMPANHAMENTOLOG == 1) then
			write(20,*) TFlash, 'Tflash'
			write(20,*) PFlash, 'Pflash'
		endif
!		if (ACOMPANHAMENTOLOG2 == 1) then
!		NENHUM ECHO VAI FICAR MELHOR QUE OS PRÓPRIOS ARQUIVOS DE INPUT
!			E ELES JÁ ESTÃO ORGANIZADOS EM PASTAS, ENTÃO É MELHOR COPIAR A PASTA
!				USANDO CALL SYSTEM(cp blablabla), ou manualmente
!			
!			!tudo sem formatação, várias casas decimais
!			!cabeçalho (algoritimicos, termodinamicos sistema, termodinamicos por comp)
!			write(22,*) 'Tflash', 'Pflash'
!			!dados
!			write(22,*) flash_spec(1), flash_spec(2), ncomp, nphas
!			write(22,*) (trim(comp_list(i)%comp%name_),i=1,ncomp)
!			write(22,*) (z_comp,i=1,ncomp)
!			write(22,*) TFlash, PFlash
!			write(22,*) (phase_list(j)%phase%phase_cond,j=1,nphas)
!		endif
		!
		
		if( ACOMPANHAMENTOLOG2 == 1 ) then
				write(21,'(A)', advance = 'NO') ' loop_x'
				write(21,'(2(A))', advance = 'NO') '       Tflash', '       Pflash'
				do j = 1, nphas0
						!essa é pra conferir:
						write(21,'(A,I1.1,A)', advance = 'NO') '      cond(',j,')'
						write(21,'(2(A,I1.1,A))', advance = 'NO') '   stabvar(',j,')', '  phasfrac(',j,')'
						do i = 1, ncomp
							write(21,'(A,I2.2,A,I1.1,A)', advance = 'NO') '      x(',i,';',j,')'
						enddo !i = 1, ncomp
				enddo !j = 1, nphas0
				write(21,'(2(A))', advance = 'NO') '      conv_FE', '       conv_x'
				write(21,'(5(A))', advance = 'NO') '  phfr<0', '   swref', '  stvr<0', '   toosm', '  collps'
				write(21,'(3(A))', advance = 'NO') '         x_rf', '        NR_rf', '       MI_aux'
				write(21,*)
		endif
		
		!INÌCIO DO ALGORITMO
		do loop_x = 1, loop_x_max
			conv_x = 0.d0
			x_ol(:,:) = x(:,:)
			if (ACOMPANHAMENTOLOG == 1) then
				write(20,*) __FILE__, __LINE__
!				do j = 1, nphas
!					write(20,*) x(:,j), 'x(:,', j
!!					write(20,*) ph_index(j)
!				enddo !j = 1, nphas
			endif
			call update_K(dist_coef,Tflash,Pflash,x)
			
!			if (ACOMPANHAMENTOLOG == 1) then !BOTEI O PRINT LÁ DENTRO
!				do j = 1, nphas
!					write(20,*) dist_coef(:,j), 'K(:,', j
!				enddo !j = 1, nphas
!			endif
			do loop_FE = 1, loop_FE_max
				conv_FE = 0.d0
				!salva as variáveis antes do newton_step para permitir rebobinar um passo
				phase_frac_ol(:) = phase_frac(:)
				stab_var_ol(:) = stab_var(:)
				phase_flick_ol(:) = phase_flick(:)

				!ACOMPANHAMENTOS
				if (ACOMPANHAMENTOLOG == 1) then; write(20,*) loop_FE, 'loop_FE count'; write(20,*) stab_var(:), 'stab var'; write(20,*) phase_frac(:), 'phase frac'; write(20,*) phase_flick(:), 'phase_flick' ;write(20,*) 'step:'; endif
				if (ACOMPANHAMENTOSCREEN == 1) then; print*, stab_var(:), 'stab var'; print*, phase_frac(:), 'phase frac'; print*, phase_flick(:), 'phase_flick'; print*, 'step:'; endif
				!NEWTON STEP:
				if( flash_spec(1) == 'P_' .and. flash_spec(2) == 'T_') then
					call newton_step_PT(z_comp,dist_coef,stab_var,phase_frac,phase_flick)		!SINGLE-SPEC
				elseif ( flash_spec(2) == 'F1' ) then
					phase_frac(1) = 0.d0
					stab_var(1) = 0.d0
					phase_flick(1) = 1 !fisicamente, nesse cálculo, a fase está em igualdade de fugacidades, então phase_flick = 1, mas ele não importa pois essas variáveis vão para o NR como constantes e nenhuma das duas será calculada
					if( flash_spec(1) == 'P_' .or. flash_spec(1) == 'T_' ) then
!!DEBUG
!	format_ = '(' // trim(nphas_str) // '(" ", A, "           "),A,)'; !print*, format_ !debug
!	write(*,format_) (phase_list(ph_index(j))%phase%phase_cond,j=1,nphas), ' tipo_de_fase'
!	format_ = '(' // trim(nphas_str) // '("       ",I1.1, "         ")," ",I1.1,A,)'; !print*, format_ !debug
!	write(*,format_) phase_flick(:), sum(phase_flick), '_fases_presentes'
!!DEBUGEND
						call newton_step_B(z_comp,dist_coef,stab_var,phase_frac,phase_flick)	!MULTI-SPEC
					else
						print*, 'unexpected flash_spec at', __FILE__, __LINE__; stop
					endif
				elseif ( flash_spec(1) == 'F1' .and. flash_spec(2) == 'F2' ) then
					phase_frac(1) = 0.d0
					stab_var(1) = 0.d0
					phase_flick(1) = 1 !fisicamente, nesse cálculo, a fase está em igualdade de fugacidades, então phase_flick = 1, mas ele não importa
!					:pois essas variáveis vão para o NR como constantes e nenhuma das duas será calculada
					phase_frac(2) = 0.d0
					stab_var(2) = 0.d0
					phase_flick(2) = 1 !fisicamente, nesse cálculo, a fase está em igualdade de fugacidades, então phase_flick = 1, mas ele não importa
!					:pois essas variáveis vão para o NR como constantes e nenhuma das duas será calculada
					call newton_step_BB(z_comp,dist_coef,stab_var,phase_frac,phase_flick)	!MULTI-SPEC
				else
					print*, 'unexpected flash_spec at', __FILE__, __LINE__; stop
				endif
				!:NEWTON STEPPED
				!ACOMPANHAMENTOS:
				if (ACOMPANHAMENTOLOG == 1) then; write(20,*) Tflash, 'Tflash', Pflash, 'Pflash'; write(20,*) stab_var(:), 'stab var'; write(20,*) phase_frac(:), 'phase frac'; write(20,*) phase_flick(:), 'phase_flick'; endif
				if (ACOMPANHAMENTOSCREEN == 1) then; print*, Tflash, 'Tflash', Pflash, 'Pflash' ;print*, stab_var(:), 'stab var'; print*, phase_frac(:), 'phase frac'; print*, phase_flick(:), 'phase_flick'; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
				!
				!se a T ou P muda, o K tem que mudar, mesmo sem mudar o x, pra ser consistente	!MULTI-SPEC
!				if (ACOMPANHAMENTOLOG == 1) then
!					write(20,*) __FILE__, __LINE__
!					do j = 1, nphas
!						write(20,*) x(:,j), 'x(:,', j
!					enddo !j = 1, nphas
!				endif

				if( flash_spec(1) == 'P_' .and. flash_spec(2) == 'T_') then
					!não precisa calcular K aqui para o flash PT
				else !flash_B ou flash_BB
					call update_K(dist_coef,Tflash,Pflash,x)										!MULTI-SPEC
					x_ol(:,:) = x(:,:)
				endif
				!obs: o update_K chama as rotinas que calculam phi, e para a fase hidrato, essas rotinas vão alterar o phase%x,
				!e, desde a versão 15.01.07 SSbi-generic, o x do algoritmo também

!				if (ACOMPANHAMENTOLOG == 1) then
!					do j = 1, nphas
!						write(20,*) dist_coef(:,j), 'K(:,', j
!					enddo !j = 1, nphas
!				endif
!
				!verificar consistência dos novos valores de solução
				!i.e., alpha e theta não negativos
				!
				!avaliar tolerância na definição de positivo ou negativo:
				!x > 0 pode ser avaliado por x > 1.d-7 e x < 0 pode ser x < -1.d-7 com entre -1d-7 e +1d-7 sendo uma zona de tolerância para igual a zero
				!
!PRIMEIRA VERIFICAÇÃO DE CONSISTÊNCIA DA CONVERGÊNCIA
				!desaparecimento de fase qualquer
				do j = 1, nphas-1
					if( phase_frac(j) < LTzero_TOL ) then
						!ACOMPANHAMENTOS:
						if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'alpha(',j,')<0'; endif
						if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'alpha(',j,')<0'; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
						event_flag(1) = ph_index(j)
						!eliminar fase
						phase_frac(j) = 0.d0
						phase_flick(j) = 0.d0
						!normaliza demais fases
						sum1 = 0.d0
						do jj = 1, nphas
							sum1 = sum1 + phase_frac(jj)
						enddo !j = 1, nphas
						phase_frac(:) = phase_frac(:)/sum1
						!re_set_theta
						!OPçãO1a
!						stab_var(j) = 1.d-1
						!OPçãO1b
						stab_var(j) = too_small
!						print*, dist_coef(1,:), __FILE__, __LINE__; pause !OFF
						!OPçãO2
!						stab_var(j) = re_set_theta(z_comp,dist_coef,stab_var,phase_frac,j)
!						if( stab_var(j) < GTzero_TOL ) then
!							!OPçãO1a
!							stab_var(j) = 1.d-1 !*** testar 1.d-9 qndo o multiflash3 estiver pronto
!							!OPçãO1b
!!						stab_var(j) = too_small !atrapalhava a convergência no multiflash antigo
!!						stab_var(j) = 0.d0				!os dois iguais a zero dava NaNaNaNaN no multiflash antigo
!						end if

						!exige nova iteração
						conv_FE = 1.d0
						conv_x = 1.d0 !em x também
					endif
				enddo !j = 1, nphas-1
				!
!SEGUNDA VERIFICAÇÃO DE CONSISTÊNCIA DA CONVERGÊNCIA
!!!				if( conv_FE > conv_FE_TOL) then !se foi feita exigência de nova iteração,
!!!					cycle !loopFE !fazer imediantamente #CYCLEMUSTDIE
!!!				endif
				if( conv_FE < conv_FE_TOL) then
					!desaparecimento da fase de referência
					if( phase_frac(nphas) < LTzero_TOL ) then
						!ACOMPANHAMENTOS:
						if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'alpha_ref<0'; endif
						if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'alpha_ref<0'; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
						!rebobinar todas as fases
						phase_frac(:) = phase_frac_ol(:)
						stab_var(:) = stab_var_ol(:)
						!eliminar fase (ref)
						phase_frac(nphas) = 0.d0
						phase_flick(nphas) = 0.d0
						!re_set_theta será feito APÓS REORDER POIS O RESULTADO DEPENDE DA FASE DE REFERÊNCIA
						
						!normaliza demais fases
						sum1 = 0.d0
						do jj = 1, nphas
							sum1 = sum1 + phase_frac(jj)
						enddo !j = 1, nphas
						phase_frac(:) = phase_frac(:)/sum1
						!
						!reorder (nova ref é a fase com maior alpha, ou seguindo alguma lista de prioridades)
						!seleciona a nova fase
						if( flash_spec(1) == 'P_' .and. flash_spec(2) == 'T_') then
							ref_index = 1 !#reveza #rodízio de fases de referência: 123 -> 231 -> 312 -> 123
						elseif ( flash_spec(2) == 'F1' ) then
							if( flash_spec(1) == 'P_' .or. flash_spec(1) == 'T_' ) then
								ref_index = 2 !não mexe na fase 1, para não atrapalhar algoritmo de busca de fase incipiente #reveza #rodízio de fases de referência: 123 -> 231 -> 312 -> 123
							else
								print*, 'unexpected flash_spec at', __FILE__, __LINE__; stop
							endif
						elseif ( flash_spec(1) == 'F1' .and. flash_spec(2) == 'F2' ) then
							ref_index = 3
						else
							print*, 'unexpected flash_spec at', __FILE__, __LINE__; stop
						endif
!						print*, 'teste de proibição de referencia @ ', __FILE__, __LINE__; ref_index = 4;
						!reorder
						!ACOMPANHAMENTOS:
						if (ACOMPANHAMENTOLOG == 1) then; write(20,*) stab_var(:), 'stab var'; write(20,*) phase_frac(:), 'phase frac'; write(20,*) 'switch ref:'; endif
						if (ACOMPANHAMENTOSCREEN == 1) then; print*, stab_var(:), 'stab var'; print*, phase_frac(:), 'phase frac'; print*, 'switch ref:'; endif
						
						event_flag(2) = ph_index(ref_index)
						
						call switch_reference(x,x_ol,dist_coef,stab_var,phase_frac,phase_flick,ph_index,ref_index)
						!a stab_var_ol e phase_frac_ol não foram reordenadas, mas tudo bem porque vamos dar o cycl imediatamente
						
						!ACOMPANHAMENTOS:
						if (ACOMPANHAMENTOLOG == 1) then; write(20,*) stab_var(:), 'stab var'; write(20,*) phase_frac(:), 'phase frac'; endif
						if (ACOMPANHAMENTOSCREEN == 1) then; print*, stab_var(:), 'stab var'; print*, phase_frac(:), 'phase frac'; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
						!
						!recalcula os coeficientes de distribuição com base na nova referência
!						if (ACOMPANHAMENTOLOG == 1) then
!							write(20,*) __FILE__, __LINE__
!							do j = 1, nphas
!								write(20,*) x(:,j), 'x(:,', j
!							enddo !j = 1, nphas
!							do j = 1, nphas
!								write(20,*) dist_coef(:,j), 'K(:,', j
!							enddo !j = 1, nphas
!							write(20,*) Tflash, Pflash
!						endif

						call update_K(dist_coef,Tflash,Pflash,x)
						x_ol(:,:) = x(:,:)
						
!						if (ACOMPANHAMENTOLOG == 1) then
!							do j = 1, nphas
!								write(20,*) dist_coef(:,j), 'K(:,', j
!							enddo !j = 1, nphas
!						endif
						!
						!re_set_theta (APÓS REORDER POIS O RESULTADO DEPENDE DA FASE DE REFERÊNCIA)
						re_set_theta_index = nphas-1 !índice da fase que precisar ser recalculada após o reordenamento (o novo índice da fase que era a fase de referência)
						
						!OPçãO1a
!						stab_var(re_set_theta_index) = 1.d-1 !*** testar 1.d-9 qndo o multiflash3 estiver pronto
						!OPçãO1b
						stab_var(re_set_theta_index) = too_small
!						print*, dist_coef(1,:), __FILE__, __LINE__; pause !OFF
						!OPçãO2
!						stab_var(re_set_theta_index) = re_set_theta(z_comp,dist_coef,stab_var,phase_frac,re_set_theta_index)
!						if( stab_var(re_set_theta_index) < GTzero_TOL ) then
!							!OPçãO1a
!							stab_var(re_set_theta_index) = 1.d-1 !*** testar 1.d-9 qndo o multiflash3 estiver pronto
!							!OPçãO1b
!!						stab_var(j) = too_small !atrapalhava a convergência no multiflash antigo
!!						stab_var(j) = 0.d0				!os dois iguais a zero dava NaNaNaNaN no multiflash antigo
!						end if
						
						stab_var(re_set_theta_index) = re_set_theta(z_comp,dist_coef,stab_var,phase_frac,re_set_theta_index) !não é mais esse j após ter reordenado, é o j-1 se estiver fazendo o rodízio (se index_ref usado foi = 1); se j era 1, j-1 é zero, e deverá ser usado j=nphas
						!exige nova iteração
						conv_FE = 1.d0
						conv_x = 1.d0 !em x também
					endif
				endif !conv_FE < conv_FE_TOL
				!
!TERCEIRA VERIFICAÇÃO DE CONSISTÊNCIA DA CONVERGÊNCIA
!!!				if( conv_FE > conv_FE_TOL) then !se foi feita exigência de nova iteração,
!!!				cycle !fazer imediantamente #CYCLEMUSTDIE
!!!				endif
				if( conv_FE < conv_FE_TOL) then
					!identificação de fase teoricamente mais estável que a fase de referência
						!ESSA ATITUDE ESTÁ DRÁSTICA DEMAIS, IMPEDINDO CONVERGÊNCIA PRÓXIMO A PONTOS DE ORVALHO
						!TALVEZ ELA DEVA SER TOMADA, SE PERTINENTE, APENAS APÓS A CONVERGÊNCIA OU ESTOURO DE PASSOS DO NEWTON
						!E AQUI NO MEIO DOS LOOPS ISSO SERIA PERMITIDO
							!MESMO SEM ESSAS VERIFICAÇÕES DEU PRA CONVERGIR A 300, 350, 400, 500 E 200K
							!A VERIFICAÇÃO DE ALFA(REF) DEVE ESTAR SENDO SUFICIENTE.
					do j = 1, nphas-1
						if( stab_var(j) < LTzero_TOL ) then
							!ACOMPANHAMENTOS:
							if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'theta(',j,')<0'; endif
							if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'theta(',j,')<0'; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
							!
							event_flag(3) = ph_index(j)
							!
	!						!rebobinar
	!						phase_frac(:) = phase_frac_ol(:)
	!						stab_var(:) = stab_var_ol(:)
	!						!eliminar fase (ref)
	!						phase_frac(nphas) = 0.d0
	!						!reorder (nova ref é a fase com theta negativo ou a próxima de uma lista de prioridades)
	!						!seleciona a nova fase
	!						ref_index = 1 !rodízio de fases de referência: 123 -> 231 -> 312 -> 123
	!						!reorder
	!						print*,; print*, stab_var(:); print*, phase_frac(:)
	!						call switch_reference(x,x_ol,dist_coef,stab_var,phase_frac,ph_index,ref_index)
	!						x_ol(:,:) = x(:,:) !precisa reordenar o xol tb
	!						print*,; print*, stab_var(:); print*, phase_frac(:); !if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
	!						!re_set_theta (após reorder pois o resultado depende da fase de referência)
	!						if( j == 1 ) then !índice da fase que precisar ser recalculada após o reordenamento
	!							re_set_theta_index = nphas
	!						else !(j > 1)
	!							re_set_theta_index = j-1
	!						endif
	!						call update_K(dist_coef,Tflash,Pflash,x)
	!						stab_var(re_set_theta_index) = re_set_theta(z_comp,dist_coef,stab_var,phase_frac,re_set_theta_index) !não é mais esse j após ter reordenado, é o j-1 se estiver fazendo o rodízio (se index_ref usado foi = 1); se j era 1, j-1 é zero, e deverá ser usado j=nphas
	!						!recalcula os coeficientes de distribuição com base na nova referência
	!						!exige nova iteração
	!						conv_FE = 1.d0
		!					conv_x = 1.d0 !em x também
	!						exit !parar de analisar as fases
						
!PRESENT PHASE COUNTING
								nphas_present = 0
								do jj = 1, nphas
!OPçãO1
!									if( dabs(phase_frac(jj)) > too_small/2.d0 ) then
!									!talvez eu devesse usar a seguinte verificação em vez dessa:
!!!!!!!								if( (stab_var(jj)) < too_small/2.d0 ) then
!									!preciso testar
!										nphas_present = nphas_present + 1
									!endif
!OPçãO2
									if( phase_flick(jj) == 1 ) then
										nphas_present = nphas_present + 1
									endif
								enddo !jj = 1, nphas

					!verifica regra das fases de Gibbs para não gerar Jacobiano singular no caso do flash T_P
							if( flash_spec(1) == 'P_' .and. flash_spec(2) == 'T_') then

								if( nphas_present < (ncomp+1) ) then !número máximo de fases simultaneamente presentes nas iterações do flash PT deve ser NC+1
									!NC para um sistema com 2 graus de liberdade, NC mais um para poder haver troca de fases, se não o algoritmo fica preso
									!e converge para uma tabela com fases com theta negativo que não foram permitidas de aparecer
									stab_var(j) = 0.d0
									phase_flick(j) = 1
									phase_frac(j) = too_small
								else
									!não permitir o aparecimento da fase
								endif
								
							elseif ( flash_spec(2) == 'F1' ) then
								if( flash_spec(1) == 'P_' .or. flash_spec(1) == 'T_' ) then
									if( nphas_present < (ncomp+2) ) then
										stab_var(j) = 0.d0
										phase_flick(j) = 1
										phase_frac(j) = too_small
									endif !nphas_present < (ncomp+1)
								else
									print*, 'unexpected flash_spec at', __FILE__, __LINE__; stop
								endif
							elseif ( flash_spec(1) == 'F1' .and. flash_spec(2) == 'F2' ) then
								if( nphas_present < (ncomp+3) ) then
									stab_var(j) = 0.d0
									phase_flick(j) = 1
									phase_frac(j) = too_small
								endif !nphas_present < (ncomp+2)
							else
								print*, 'unexpected flash_spec at', __FILE__, __LINE__; stop
							endif

	!						stab_var(j) = too_small		!os dois iguais a too_small volta para cá infinitamente
						
	!						print*, phase_frac(:)
	!						print*, stab_var(:), 'stab_var'
	!						print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused
						
	!						phase_frac(j) = 0				!os dois iguais a zero dá NaNaNaNaN
						endif !stab_var(j) < LTzero_TOL
					enddo !j = 1, nphas-1
				endif !conv_FE < conv_FE_TOL
				
!QUARTA VERIFICAÇÃO DE CONSISTÊNCIA DA CONVERGÊNCIA, DESATIVADA
!!!				if( conv_FE > conv_FE_TOL) then !se foi feita exigência de nova iteração,
!!!				cycle !fazer imediantamente #CYCLEMUSTDIE
!!!				endif
!				if( conv_FE < conv_FE_TOL) then
					!é uma recomendação do artigo que não está sendo útil para mim
				
					!se os valores de estabilidade e quantidade são numericamente muito pequenos simultaneamente, perturbá-los (e exigir nova iteração?).
					!OBS: isso vai ser ativado a cada iteração se houver sido realizada um colapso de fases iguais
						!pois vão ter duas fases iguais, logo estáveis, theta=0, mas uma delas vai ter sido esvaziada, alpha = 0.
						!e essa solução é válida, não é um realmente um problema.
							!solução: só fazer isso se além de alfa e theta serem zero, não tiver outra fase com xizes iguais?
						
					!OBS2: agora o meu colapso remove a fase do sistema numérico
						!posso voltar a testar a utilidade desse bloco depois
						
	!				do j = 1, nphas-1
	!					if( dabs(stab_var(j)) < Too_Small .and. dabs(phase_frac(j)) < Too_Small ) then
	!					event_flag(4) = ph_index(j)
	!						print*,; print*, j,'too small'; !if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif !ACOMPANHAMENTO
	!						stab_var(j) = Too_Small
	!						phase_frac(j) = Too_Small
	!!!!!						print*,; print*, 'x'; print*, x(1,:); print*, x(2,:)
	!!						conv_FE = 1.d0 !exige nova iteração
	!						conv_x = 1.d0 !em x também
	!					endif
	!				enddo !j = 1, nphas-1
!				endif !conv_FE < conv_FE_TOL

!CONCLUSÃO DO LOOP FE
!!!				if( conv_FE > conv_FE_TOL) then !se foi feita exigência de nova iteração,
!!!				cycle !fazer imediantamente #CYCLEMUSTDIE
!!!				endif
				if( conv_FE < conv_FE_TOL) then
				!avaliar a convergência apenas se já não se tiver exigido uma nova iteração
				!i.e. sendo os valores fisicamente consistentes
				!
!REGISTRAR VARIÁVEIS SE DESEJÁEL
					!ACOMPANHAMENTOS:
					!se aqui der x NaN ou infinity para fase hidrato, se deve a problema no cálculo desse x dentro do módulo termodinâmico, e não no update_x
					if (ACOMPANHAMENTOLOG == 1) then; write(20,*) __FILE__, __LINE__; do j = 1, nphas; write(20,*) phase_list(ph_index(j))%phase%x(:), 'ph%x(:,', j; enddo; do j = 1, nphas; write(20,*) phase_list(ph_index(j))%phase%phi(:)*Pflash*phase_list(ph_index(j))%phase%x(:), 'fug(:,', j; enddo; endif
					if (ACOMPANHAMENTOSCREEN == 1) then; do j = 1, nphas; print*, phase_list(ph_index(j))%phase%x(:), 'ph%x(:,', j; enddo; do j = 1, nphas; print*, phase_list(ph_index(j))%phase%phi(:)*Pflash*phase_list(ph_index(j))%phase%x(:), 'fug(:,', j; enddo; endif
					!
!calcular NLEs e verificar convergencia do loop interno
					do j = 1, nphas-1
						conv_FE = conv_FE + dabs(NLE_E_j(z_comp,dist_coef,stab_var,phase_frac,j))
					enddo !j = 1, nphas-1
					!
				endif !conv_FE < conv_FE_TOL
				!
!REGISTRAR VARIÁVEIS SE DESEJÁEL
				!ACOMPANHAMENTOS:
				if (ACOMPANHAMENTOLOG == 1) then; write(20,*) conv_FE, 'conv_FE'; endif
				if (ACOMPANHAMENTOSCREEN == 1) then; print*, conv_FE, 'conv_FE'; endif
!				print*, Pflash/1.d5, __FILE__, __LINE__
!CONVERGIU?
				if (conv_FE < conv_FE_TOL) then
					!*** ANOTAR LOOP_FE EM ALGUM LUGAR? IMPIRMIR NO REPORT, JUNTO COM O LOOP_X, SUM(LOOP_FE)/LOOP_X OU ENTÃO MAX(LOOP_FE)?
					Loop_FE_total = Loop_FE_total + loop_FE
					exit !loop_FE !convergiu !eu uso exit porque minha variavel de contador é um valor máximo de iterações para não permitir loop_infinito
				endif
				!
!NÃO CONVERGIU
				if (loop_FE == loop_FE_max) then
						!ACOMPANHAMENTOS: 
						if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'Too many steps in loop_E'; endif
						if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'Too many steps in loop_E'; endif
!					print*, 'Too many steps in loop_E'; stop	!FLAG UPDATE K IN EVERY STEP
					Loop_FE_total = Loop_FE_total + loop_FE
					conv_x = 1.d0
					exit !loop_FE !(nesse bloco, o exit na verdade aconteceria independentemente de eu escrever exit pois é o maximo do contador) para mudar os x e dar outra chance, mas já condenando o loopx a fazer outra iteração
				endif
!				call update_x(x,z_comp,dist_coef,stab_var,phase_frac)		!isso não pode ficar ligado pois aqui não tem o controle para colapasar as fases e etc
			enddo !loop_FE = 1, 10000 !############################################################################################################################
			!
!LOOP EXTERNO
			!ACOMPANHAMENTOS:
			if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'loop_FE_convergiu', loop_FE; write(20,*) stab_var(:), 'stab var'; write(20,*) phase_frac(:), 'phase frac'; endif
			if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'loop_FE_convergiu', loop_FE; print*, stab_var(:), 'stab var'; print*, phase_frac(:), 'phase frac'; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
			!
			call update_x(x,z_comp,dist_coef,stab_var,phase_frac)
			!:atualiza o x das fases do grupo 1 (o grupo das fluidas), e atualiza a fug do grupo 2
			!
			!ACOMPANHAMENTOS:
			if (ACOMPANHAMENTOSCREEN == 1) then; do j = 1, nphas; print*, x(:,j), 'x(:,', j; enddo; endif
			if (ACOMPANHAMENTOLOG == 1) then; write(20,*) loop_x, 'loop_x count', __FILE__, __LINE__; do j = 1, nphas; write(20,*) x(:,j), 'x(:,', j; write(20,*) ph_index(j); enddo; endif
			!
!COLAPSAR FASES IGUAIS ***proibimos que colapse a fase ref, e vamos tentar evitar que colapse a fase 1 e fase 2
				!vamos supor que em dada passada do algoritmo por aqui
					!só vai aparecer no máximo um par de fases para colapsar,
						!e então colapsar ela ( a uma fase detectada após a verificação de trivialidade dentre todas as fases)
			collapse_tag(:) = 0
			do j = 1, nphas-1			!naturalmente só roda se nphas > 1
				do jj = j+1, nphas		!naturalmente só roda se nphas > 1
					!
!					if( phase_frac(j) > GTzero_TOL .AND. phase_frac(jj) > GTzero_TOL ) then !aplicar ação de colapsar apenas caso as duas fases de fato estejam presentes
						!condicional desativada, fases não presentes também podem deixar o jacobiano singular se convergirem para a mesma solução
!					endif
					sum1 = 0.d0
					do i = 1, ncomp
						sum1 = sum1 + dabs(x(i,j)-x(i,jj))
					enddo !i = 1, ncomp
						sum1 = sum1 + 1.d3*dabs(phase_list(ph_index(j))%phase%V-phase_list(ph_index(jj))%phase%V)					!não colapsar fases com volumes diferentes (Equilibrio Liq-Vapor de substancias puras e azeótropos)
					if( sum1 < collapse_TOL ) then						!avaliar a tolerância para colapsar e a para convergir,
!					print*, j, jj, sum1, dabs(x(1,j)-x(1,jj)), dabs(phase_list(ph_index(j))%phase%V-phase_list(ph_index(jj))%phase%V); pause
																				!para que o algoritmo não julgue convergencia sem ter colapsado algo que deva ser colapsado
																				!nem colapse à toa parando o algoritmo em cálculos próximo ao ponto crítico,
																					!o colapso deve ser feito antes que o Jac se torne singular
																						!a tolerancia de convergencia verifica todos os xij entre uma iteração e outra
																							!essa aqui verifica todos os xi de uma fase j e outra jj
						!só tentar colapsar fases fluidas (Liq ou Vap, que vem da equação de estado)
						if( phase_list(ph_index(j))%phase%phase_cond == 'V' .or. phase_list(ph_index(j))%phase%phase_cond == 'L' ) then
							if( phase_list(ph_index(jj))%phase%phase_cond == 'V' .or. phase_list(ph_index(jj))%phase%phase_cond == 'L' ) then
								!COLAPSAR essas fases e eliminar j ou jj !a collapse_tag(2) é eliminada
								if( jj == nphas ) then !não eliminar a última
									collapse_tag(1) = jj; collapse_tag(2) = j
								elseif( j == 1 ) then !evitar eliminar a 1a
									collapse_tag(1) = j; collapse_tag(2) = jj
								elseif( j == 2 ) then !evitar eliminar a 2a
									collapse_tag(1) = j; collapse_tag(2) = jj
								elseif( phase_frac(j) > phase_frac(jj) ) then !eliminar a menor
									collapse_tag(1) = j; collapse_tag(2) = jj
								else !elseif( phase_frac(jj) .GE. phase_frac(j) ) then
									collapse_tag(1) = jj; collapse_tag(2) = j
								endif
								!forçar exigência de mais um loop
		!						conv_x = 1.d0												!desativada
							endif !phase_cond jj
						endif !phase_cond j
					endif
				enddo !jj = j+1, nphas
			enddo !j = 1, nphas-1
			!
			if( collapse_tag(1) .NE. 0 ) then
				event_flag(5) = ph_index(collapse_tag(2)) !fase removida
				phase_frac(collapse_tag(1)) = phase_frac(collapse_tag(1)) + phase_frac(collapse_tag(2))

				call collapse_phase(x,x_ol,dist_coef,stab_var,phase_frac,phase_flick,ph_index,collapse_tag(2))					!<---REDUÇÂO DA ORDEM DO NEWTON
				!
				!ACOMPANHAMENTOS:
				if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'collapsed, novo nphas = ', nphas; endif
				if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'collapsed, novo nphas = ', nphas; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
				!
				!redimensionar os backups de alfa e theta, não precisa mexer nos valores pois eles já serão atualizados no próximo loop
				deallocate(phase_frac_ol)
				deallocate(stab_var_ol)
				deallocate(phase_flick_ol)															!<---REDUÇÂO DA ORDEM DO NEWTON
				allocate(phase_frac_ol(nphas))
				allocate(stab_var_ol(nphas))
				allocate(phase_flick_ol(nphas))												!<---REDUÇÂO DA ORDEM DO NEWTON
			endif
			
			!avaliação da convergência em x
			if( conv_x < conv_x_TOL ) then !avaliar a convergência apenas se já não se tiver exigido uma nova iteração
				do i = 1, ncomp
					do j = 1, nphas
						conv_x = conv_x + dabs(x(i,j)-x_ol(i,j))
					enddo !j = 1, nphas
				enddo !i = 1, ncomp
			endif
			
			!ACOMPANHAMENTOS:
			!CHECKPOINTS
			if( ACOMPANHAMENTOSCREEN /= 2 ) then; if( loop_x == checkpoint ) then; write(*,*), 'loop_x is at ', loop_x; checkpoint = checkpoint*10; endif; endif

			!it_log.dat e serial report
!			if( ACOMPANHAMENTOLOG2 == 1 .or. serial_report_FLAG==1 ) then
				!se ph_index_inv(j) não aparecer de 1 até nphas, deixa valor zero.
				ph_index_inv(:) = 0
				j = 1
				jj = 1
				do
					if (j > nphas0) exit
					if( j == ph_index(jj) ) then
						ph_index_inv(j) = jj
						jj = 1
						j = j + 1
					else
						jj = jj+1
						if (jj > nphas) then
							jj = 1
							j=j+1
						endif
					endif
				enddo
!			endif !( ACOMPANHAMENTOLOG2 == 1 .or. serial_report_FLAG==1 ) then

			!serial_report
			if( conv_x < conv_x_TOL ) then
				do i = 1, ncomp
					write(serial_report_fnum,'(ES13.5)', advance = 'NO') z_comp(i) !se eu for dar open com nome de arquivo, o open deve ser no main, e o filenumber deve chegar aqui, pode ser variavel do modulo multiflash
				enddo
				write(serial_report_fnum,'(I7.1)', advance = 'NO') loop_x
				write(serial_report_fnum,'(2(ES13.5))', advance = 'NO') Tflash, Pflash
				do j = 1, nphas0
					if( ph_index_inv(j) /= 0 ) then
						!essa é pra conferir:
						write(serial_report_fnum,'(2(A))', advance = 'NO') '        ', adjustr(phase_list(ph_index(ph_index_inv(j)))%phase%phase_cond)
						write(serial_report_fnum,'(2(ES13.5))', advance = 'NO') stab_var(ph_index_inv(j)), phase_frac(ph_index_inv(j))
						do i = 1, ncomp
							write(serial_report_fnum,'(ES13.5)', advance = 'NO') x(i,ph_index_inv(j))
						enddo !i = 1, ncomp
					else
						!essa é pra conferir:
						write(serial_report_fnum,'(A)', advance = 'NO') adjustr('            X')
						write(serial_report_fnum,'(2(ES13.5))', advance = 'NO') 0.d0, 0.d0
						do i = 1, ncomp
							write(serial_report_fnum,'(ES13.5)', advance = 'NO') 0.d0
						enddo !i = 1, ncomp
					endif
				enddo !j = 1, nphas0
				write(serial_report_fnum,'(2(ES13.5))', advance = 'NO') conv_FE, conv_x
				write(serial_report_fnum,'(5(I8.1))', advance = 'NO') (event_flag(j),j=1,5)
				write(serial_report_fnum,'(3(ES13.5))', advance = 'NO') (met_num_par_val(j),j=1,3)
				write(serial_report_fnum,*) ' '
			endif !( conv_x < conv_x_TOL ) then !serial_report

			!it_log.dat:
			if( ACOMPANHAMENTOLOG2 == 1 ) then
				write(21,'(I7.1)', advance = 'NO') loop_x
				write(21,'(2(ES13.5))', advance = 'NO') Tflash, Pflash
				do j = 1, nphas0
					if( ph_index_inv(j) /= 0 ) then
						!essa é pra conferir:
						write(21,'(2(A))', advance = 'NO') '        ', adjustr(phase_list(ph_index(ph_index_inv(j)))%phase%phase_cond)
						write(21,'(2(ES13.5))', advance = 'NO') stab_var(ph_index_inv(j)), phase_frac(ph_index_inv(j))
						do i = 1, ncomp
							write(21,'(ES13.5)', advance = 'NO') x(i,ph_index_inv(j))
						enddo !i = 1, ncomp
					else
						!essa é pra conferir:
						write(21,'(A)', advance = 'NO') adjustr('            X')
						write(21,'(2(ES13.5))', advance = 'NO') 0.d0, 0.d0
						do i = 1, ncomp
							write(21,'(ES13.5)', advance = 'NO') 0.d0
						enddo !i = 1, ncomp
					endif
				enddo !j = 1, nphas0
				write(21,'(2(ES13.5))', advance = 'NO') conv_FE, conv_x
				write(21,'(5(I8.1))', advance = 'NO') (event_flag(j),j=1,5)
				!resetar a flag até que outros eventos ocorram
				event_flag(:) = 0
				write(21,'(3(ES13.5))', advance = 'NO') (met_num_par_val(j),j=1,3)
				write(21,*)
			endif !ACOMPANHAMENTOLOG2 > it_log.dat
			
			!***NÃO BOTEI O FLICK NO SERIAL REPORT e it_log.dat

!CONVERGIU?
			if (conv_x < conv_x_TOL) then
				exit !convergiu
			endif
!NÃO CONVERGIU
			if (loop_x == loop_x_max) then
				print*, 'Too many steps in loop_x' !msg importante independente do acompanhamento, pois ela fecha o programa
				!fechar programa
!				stop
!				voltar ao menu em vez de fechar
				!recuperar variaveis pré divergir ?
					!seriam as da ultima simulação convergida, presentes no serial_report, se essa não for a primeira simulação do programa
				
				!retornar
				call cpu_time(time_2)
				return
			endif
			
			!ACOMPANHAMENTO:
			if (ACOMPANHAMENTOLOG == 1) then; write(20,*) loop_x, 'loop_x'; endif
		enddo !loop_x = 1, loop_x_max
		
!		print*, conv_x, conv_FE; pause
		
!		!informa o x novo para os objeto fase, não precisa, o update_x e o update_K já estão cuidando disso.
!		do j = 1, nphas
!			phase_list(ph_index(j))%phase%x(:) = x(:,j)
!		enddo !j = 1, nphas
		
!		RESULTS
!		print*,;print*,;print*,;print*,;print*,;print*,;print*,;print*,;print*,;print*,
!		print*, 'Flash is complete'
!		print*, 'tipo de fase'
!		format_ = '(' // trim(nphas_str) // '(" ", A, "          "))'; !print*, format_ !debug
!		write(*,format_) (phase_list(ph_index(j))%phase%phase_cond,j=1,nphas)
!		format_ = '(' // trim(nphas_str) // '(" ", F7.4, "          "))'; !print*, format_ !debug
!		write(*,format_) (phase_list(ph_index(j))%phase%V*1.d3,j=1,nphas)
!		print*, 'estabilidade'
!		format_ = '(' // trim(nphas_str) // '(F7.4, "         "))'; !print*, format_ !debug
!		write(*,format_) stab_var(:)
!		print*, 'fração de fase'
!		write(*,format_) phase_frac(:)
!		print*, 'composição'
!		format_ = '(' // trim(nphas_str) // '(F7.4, "         "), A)'; !print*, format_ !debug
!		do i = 1, ncomp
!			write(*,format_) x(i,:), comp_list(i)%comp%name_
!!			write(*,format_) (phase_list(ph_index(j))%phase%x(i),j=1,nphas), comp_list(i)%comp%name_ !CONFERE
!		enddo !i = 1, ncomp
!		
!		print*, Tflash, Pflash
		
!		print*, phase_list(2)%phase%phi(1), x(1,2)
!		print*, phase_list(1)%phase%phi(1), x(1,1)
!		
!		print*, phase_list(1)%phase%phi(1)*x(1,1) - phase_list(2)%phase%phi(1)*x(1,2)
!		print*, phase_list(1)%phase%phi(2)*x(2,1) - phase_list(2)%phase%phi(2)*x(2,2)
		
!!!!		!PERDEDOR DE TEMPO
!!!!		do i = 1, 100000
!!!!			do j = 1, 1000
!!!!				jj=i+j
!!!!			enddo !1 = 1, 100000
!!!!		enddo !1 = 1, 100000

!		print*, ph_index, __FILE__, __LINE__
!		print*, ph_index_inv, __FILE__, __LINE__
!		read(*,*) paused

		call cpu_time(time_2)
		print*, Tflash, Pflash, __FILE__, __LINE__
!		print*, time_2-time_1, 'seconds @ ', __FILE__, __LINE__
	end subroutine multiphaseflash

!!!!!SUBROUTINES DO FLASH
	subroutine multiphaseflash_report
		character(100) :: format_
		integer :: i, j
		call system('clear') !'limpa' tela do terminal antes de imprimir os resultados
!		write(*,*)
		write(*,*) '################ ################ ################ ################ ################ ################'
!		write(*,*) 'Tflash_[K]       ', 'Pflash_[Pa]      ', ('z[', adjustr(comp_list(i)%comp%name_(1:7)), '] ',i=1,ncomp)
		write(*,*) 'Tflash_[K]       ', 'Pflash_[Pa]      ', ('z_', adjustl(comp_list(i)%comp%name_(1:15)),i=1,ncomp)
		format_ = '(' // 'F8.2, "          ", ES9.2, "       ",' // trim(ncomp_str) // '(F8.4,"         "))'
		write(*,format_), Tflash, Pflash, z_comp
		write(*,*) '################ ################ ################ ################ ################ ################'
		format_ = '(' // trim(nphas_str) // '(" ", A, "           "),A,)'; !print*, format_ !debug
		write(*,format_) (phase_list(ph_index(j))%phase%phase_cond,j=1,nphas), ' tipo_de_fase'
		format_ = '(' // trim(nphas_str) // '(ES13.5, "    "),A,)'; !print*, format_ !debug
		write(*,format_) (phase_list(ph_index(j))%phase%V*1.d3,j=1,nphas), ' volume_[L/mol]'
		format_ = '(' // trim(nphas_str) // '(F8.4, "         "),A,)'; !print*, format_ !debug
		write(*,format_) stab_var(:), ' estabilidade'
		write(*,format_) phase_frac(:), ' fração_de_fase'
		format_ = '(' // trim(nphas_str) // '("       ",I1.1, "         ")," ",I1.1,A,)'; !print*, format_ !debug
		write(*,format_) phase_flick(:), sum(phase_flick), '_fases_presentes'
		write(*,*) 'composição______ ________________ ________________ ________________ ________________ ________________'
		format_ = '(' // trim(nphas_str) // '(F8.4, "         ")," ",A)'; !print*, format_ !debug
		do i = 1, ncomp
			write(*,format_) (phase_list(ph_index(j))%phase%x(i),j=1,nphas), trim(comp_list(i)%comp%name_) !CONFERE
		enddo !i = 1, ncomp
		write(*,*) 'fugacidade_[Pa]_ ________________ ________________ ________________ ________________ ________________'
		format_ = '(' // trim(nphas_str) // '(ES13.5, "    ")," ",A)'; !print*, format_ !debug
		do i = 1, ncomp
			write(*,format_) (Pflash*phase_list(ph_index(j))%phase%phi(i)*phase_list(ph_index(j))%phase%x(i),j=1,nphas), trim(comp_list(i)%comp%name_) !CONFERE
		enddo !i = 1, ncomp
!		print*, phase_list(ph_index(nphas))%phase%phi(1) !*phase_list(ph_index(nphas))%phase%x(1)*Pflash
!!		print*, phase_list(ph_index(nphas-1))%phase%phi(1) !*phase_list(ph_index(nphas-1))%phase%x(1)*Pflash
		
!		phase_list(ph_index(4))%phase%x(:) = 0.d0
!		phase_list(ph_index(4))%phase%x(2) = 1.d0
!		call phase_list(ph_index(4))%phase%calc_phi(Tflash,Pflash)
!		print*, phase_list(ph_index(4))%phase%phi(2)
!		
!		phase_list(ph_index(4))%phase%x(:) = .0002d0
!		phase_list(ph_index(4))%phase%x(2) = .9998d0
!		call phase_list(ph_index(4))%phase%calc_phi(Tflash,Pflash)
!		print*, phase_list(ph_index(4))%phase%phi(2)
		write(*,*) '________________ ________________ ________________ ________________ ________________ ________________'
		write(*,'(I7.6,A,F9.2,A)') loop_x, ' loop_x iterations in (', time_2-time_1, ') seconds'
		write(*,*) '________________ ________________ ________________ ________________ ________________ ________________'
		write(*,*) Tflash, Pflash, 'brought to you by ', __FILE__
		write(*,*) loop_FE_total, 'loop_FE_total'
		write(*,*) conv_x, 'conv_x'
	end subroutine multiphaseflash_report

	subroutine update_K(K_,T_,P_,x_)
		real(8), intent(out) :: K_(ncomp,nphas)
		real(8), intent(in) :: T_, P_
		real(8), intent(inout) :: x_(ncomp,nphas) !x_ é argumento daqui pois no update_K pra a fase hidrato, o x é alterado
!		integer, intent(in) :: cond_(nphas)
		real(8) :: phi_(ncomp,nphas)
		integer :: i, j
		integer :: class_flag
		class(phase_model_c), pointer :: associatename
		
		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'update_K:'; write(20,*) 'T', T_, 'P', P_; do j = 1, nphas; write(20,*) 'x', x_(:,j), phase_list(ph_index(j))%phase%phase_cond; enddo; endif

		!MODO COM SELECT TYPE para algoritmo Bi-Generic
		do j = 1, nphas
			associatename => phase_list(ph_index(j))%phase
			select type (associatename)
			class is (cub_eos_c)
				class_flag=1
			class is (hydrate_model_c)
				class_flag=2
			class is (ice_sublim_c)
				class_flag=1
				!vou deixar como antes, analogo a fase fluida, sendo que phi será um valor muito grande e logo x naturalmente será aproximadamente zero para todos os componentes não presentes por definição.
			class is (Gex_based_phase_model_c)
				class_flag=1
			class default
				print*, 'unexpected class/type argument in (', __FILE__, __LINE__, ') for phase ', j; stop
			end select
			
			if(class_flag==1)then
				call associatename%calc_phi(T_,P_)
				phi_(:,j) = associatename%phi(:) !esse phi é calculado a partir do phase(j)%x atual, que pode ser modificado pela sub update_x
			elseif (class_flag==2)then
				call associatename%calc_phi(T_,P_)
				phi_(:,j) = associatename%phi(:) !esse phi é calculado a partir do phase(j)%fug_aux atual, que pode ser modificado pela sub update_x
				x_(:,j) = associatename%x(:) !estou atualizando o x do algoritmo com o x da fase
			endif
			
		enddo !j = 1, nphas
		!END MODO COM SELECT TYPE
!
!CALCULAR K-values
		do i = 1, ncomp
			do j = 1, nphas
				K_(i,j) = phi_(i,nphas)/phi_(i,j)
!				print*, K_(i,j), __FILE__, __LINE__
			enddo !j = 1, nphas
		enddo !i = 1, ncomp
!		print*, K_(:,2)
!		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; do j = 1, nphas; write(20,*) 'K', K_(:,j); enddo; endif
	end subroutine update_K
	
	subroutine update_x(x_,z_,K_,th_,al_) !essa rotina vai fazer o trabalho de atualizar as fug para a fase hidrato sem emxer no x dela.
		real(8), intent(inout) :: x_(ncomp,nphas)
		real(8), intent(in) :: z_(ncomp),K_(ncomp,nphas),th_(nphas),al_(nphas)

		real(8) :: x_nu(ncomp,nphas)
		real(8) :: x_ol(ncomp,nphas)

		real(8) :: step(ncomp,nphas), max_step(ncomp,nphas), rel_step(ncomp,nphas), reducing_factor
		integer :: i, j
		integer :: class_flag
		real(8) :: sum1
		class(phase_model_c), pointer :: associatename
		
		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'update_x'; write(20,*) 'th', th_(:); write(20,*) 'al', al_(:); do j = 1, nphas; write(20,*) 'K', K_(:,j); enddo; do j = 1, nphas; write(20,*) 'x', x_(:,j), phase_list(ph_index(j))%phase%phase_cond; enddo; endif
		
!		print*, x_, 'update_x_in

!MODO COM SELECT TYPE para algoritmo Bi-Generic

!Hidrato ref: fug novo = "x_novo" da equação X, vezes phi velho, disponível no objeto fase, vezes Pressão
!						que é equivalente a dizer fug_novo = fug_velho*(x_novo/x_velho)
!						o x real da fase não é alterado por essa rotina para ficar igual ao x_novo, ele continua como o x_velho calculado pelo update_K
!Fluid ref: x_novo = x_novo da equação X
!Hidrato j: fug_novo = "x_novo" da equação X, vezes phi velho, disponível no objeto fase, vezes Pressão
!						como x_novo para fase j é igual a x_ref_novo vezes K vezes exp(theta),
!							então é equivalente a fug_novo = fug_ref_j(x_novo,phi_velho,P)*exp(theta_j)
!Fluid j: x_novo = x_novo da equação X
!!						que é igual a x_ref_novo vezes K vezes exp(theta),

!backup
	x_ol(:,:) = x_(:,:)

!FASE DE REFERENCIA
		!MODO COM SELECT TYPE para algoritmo Bi-Generic
		associatename => phase_list(ph_index(nphas))%phase
		select type (associatename)
		class is (cub_eos_c)
			class_flag = 1
		class is (hydrate_model_c)
			class_flag = 2
		class is (ice_sublim_c) !futuramente pode ser o caso para sólidos puros em geral
			class_flag = 1
			!vou deixar como antes, analogo a fase fluida, sendo que phi será um valor muito grande e logo x naturalmente será aproximadamente zero para todos os componentes não presentes por definição.
		class is (Gex_based_phase_model_c)
			class_flag=1
		class default
			print*, 'unexpected class/type argument in (', __FILE__, __LINE__, ') for phase ', j; stop
		end select
		if(class_flag==1)then

			do i = 1, ncomp
				sum1 = 0.d0
				do j = 1, nphas-1
					sum1 = sum1 + al_(j)*((K_(i,j)*dexp(th_(j)))-1.d0)
				enddo !j = 1, nphas
				x_nu(i,nphas) = z_(i)/(1.d0+sum1) 
				!DEBUG !impedir x negativo
				if(x_nu(i,nphas)<0)then
					!***usar alguma flag de acompanhamento aqui				
					print*, 'x_nu ', i, ' nphas', ' <0', ' @', __FILE__, __LINE__
					x_nu(i,nphas)=0
				endif
			enddo
			
			x_(:,nphas)=x_nu(:,nphas)

		elseif(class_flag==2)then

!ANTES
!calcula um hidrato na fase de referencia baseado em uma fase 1
!			do i = 1, ncomp
!				associatename%fug_aux(i)  = x_nu(i,1)*phase_list(ph_index(1))%phase%phi(i)*Pflash*dexp(-1.d0*stab_var(1)) !usando a primeira fase para definir a fugacidade dos hospedes, já que o próprio hyd é fase de ref
!				!usei Pflash do módulo aqui, não tenho ela como argumento
!			enddo !i = 1, ncomp

!AGORA
!usa a equação X para respeitar o balanço de massa
!"x_novo"
			do i = 1, ncomp
				sum1 = 0.d0
				do j = 1, nphas-1
					sum1 = sum1 + al_(j)*((K_(i,j)*dexp(th_(j)))-1.d0)
				enddo !j = 1, nphas
				x_nu(i,nphas) = z_(i)/(1.d0+sum1)
				!DEBUG !impedir x negativo
				if(x_nu(i,nphas)<0)then
					!***usar alguma flag de acompanhamento aqui				
					print*, 'x_nu ', i, ' nphas', ' <0', ' @', __FILE__, __LINE__
					x_nu(i,nphas)=0
				endif
			enddo

			do i = 1, ncomp
				associatename%fug_aux(i)  = x_nu(i,nphas)*phase_list(ph_index(nphas))%phase%phi(i)*Pflash
				!x pela equação nova, phi velho e P velho
			enddo !i = 1, ncomp

!			x_ continua sendo igual ao x_ol

		endif

!REMAINING PHASES
		!MODO COM SELECT TYPE para algoritmo Bi-Generic
		do j = 1, nphas-1
			associatename => phase_list(ph_index(j))%phase
			select type (associatename)
			class is (cub_eos_c)
				class_flag=1
			class is (hydrate_model_c)
				class_flag=2
			class is (ice_sublim_c) !futuramente pode ser o caso para sólidos puros em geral
				class_flag=1
				!vou deixar como antes, analogo a fase fluida, sendo que phi será um valor muito grande e logo x naturalmente será aproximadamente zero para todos os componentes não presentes por definição.
			class is (Gex_based_phase_model_c)
				class_flag=1
			class default
				print*, 'unexpected class/type argument in (', __FILE__, __LINE__, ') for phase ', j; stop
			end select
			if(class_flag==1)then

				do i = 1, ncomp
					x_nu(i,j) = K_(i,j)*x_nu(i,nphas)*dexp(th_(j))
	!				print*, K_(i,j), x_nu(i,nphas), dexp(th_(j)), __FILE__, __LINE__
				enddo !i = 1, ncomp

				x_(:,j)=x_nu(:,j)

			elseif(class_flag==2)then

!ANTES
!				do i = 1, ncomp
!					associatename%fug_aux(i)  = x_nu(i,nphas)*phase_list(ph_index(nphas))%phase%phi(i)*Pflash*dexp(stab_var(j)) !usando a fase de referencia para definir a fugacidade dos hospedes
!				!usei Pflash do módulo aqui, não tenho ela como argumento
!				enddo !i = 1, ncomp

!AGORA
!"x_novo"
!				do i = 1, ncomp
!					x_nu(i,j) = K_(i,j)*x_nu(i,nphas)*dexp(th_(j))
!				enddo !i = 1, ncomp

!				do i = 1, ncomp
!					associatename%fug_aux(i)  = x_nu(i,j)*phase_list(ph_index(j))%phase%phi(i)*Pflash
!					!x pela equação nova, phi velho e P velho
!				enddo !i = 1, ncomp

!EQUIVALENTE A
				do i = 1, ncomp
						associatename%fug_aux(i)  = x_nu(i,nphas)*phase_list(ph_index(nphas))%phase%phi(i)*Pflash*dexp(stab_var(j)) !usando a fase de referencia para definir a fugacidade dos hospedes
	!				!usei Pflash do módulo aqui, não tenho ela como argumento
				enddo !i = 1, ncomp

!				x_ continua sendo igual ao x_ol

			endif
		enddo !j = 1, nphas-

!RESIDUO DE QUANDO EU FAZIA QUESTÃO DE FAZER AS FASES DE HIDRATO ANTES DAS FLUIDAS:
!!!!!!FASE DE REFERENCIA
!!!!!		associatename => phase_list(ph_index(j))%phase
!!!!!		select type (associatename)
!!!!!		class is (cub_eos_c)
!!!!!			class_flag = 1
!!!!!		class is (hydrate_model_c)
!!!!!			class_flag = 2
!!!!!		class is (ice_sublim_c) !futuramente pode ser o caso para sólidos puros em geral
!!!!!			class_flag = 1
!!!!!			!vou deixar como antes, analogo a fase fluida, sendo que phi será um valor muito grande e logo x naturalmente será aproximadamente zero para todos os componentes não presentes por definição.
!!!!!		class default
!!!!!			print*, 'unexpected class/type argument in (', __FILE__, __LINE__, ') for phase ', j; stop
!!!!!		end select
!!!!!		if(class_flag==1)then
!!!!!			do i = 1, ncomp
!!!!!				sum1 = 0.d0
!!!!!				do j = 1, nphas-1
!!!!!					sum1 = sum1 + al_(j)*((K_(i,j)*dexp(th_(j)))-1.d0)
!!!!!				enddo !j = 1, nphas
!!!!!				x_nu(i,nphas) = z_(i)/(1.d0+sum1) 
!!!!!				!DEBUG !impedir x negativo
!!!!!				if(x_nu(i,nphas)<0)then
!!!!!					print*, 'x_nu ', i, ' nphas', ' <0', ' @', __FILE__, __LINE__
!!!!!					x_nu(i,nphas)=0
!!!!!				endif
!!!!!			enddo
!!!!!		elseif(class_flag==2)then
!!!!!			!FEITO
!!!!!		endif
!!!!!		!END MODO COM SELECT TYPE

!!!!!!REMAINING PHASES
!!!!!		!MODO COM SELECT TYPE para algoritmo Bi-Generic
!!!!!		do j = 1, nphas-1
!!!!!			associatename => phase_list(ph_index(j))%phase
!!!!!			select type (associatename)
!!!!!			class is (cub_eos_c)
!!!!!				class_flag=1
!!!!!			class is (hydrate_model_c)
!!!!!				class_flag=2
!!!!!			class is (ice_sublim_c) !futuramente pode ser o caso para sólidos puros em geral
!!!!!				class_flag=1
!!!!!				!vou deixar como antes, analogo a fase fluida, sendo que phi será um valor muito grande e logo x naturalmente será aproximadamente zero para todos os componentes não presentes por definição.
!!!!!			class default
!!!!!				print*, 'unexpected class/type argument in (', __FILE__, __LINE__, ') for phase ', j; stop
!!!!!			end select
!!!!!			if(class_flag==1)then
!!!!!				do i = 1, ncomp
!!!!!					x_nu(i,j) = K_(i,j)*x_nu(i,nphas)*dexp(th_(j))
!!!!!	!				print*, K_(i,j), x_nu(i,nphas), dexp(th_(j)), __FILE__, __LINE__
!!!!!				enddo !i = 1, ncomp
!!!!!			elseif(class_flag==2)then
!!!!!			!FEITO
!!!!!			endif
!!!!!			
!!!!!		enddo !j = 1, nphas-


		!END MODO COM SELECT TYPE

!		print*, x_, 'update_x_raw'

		!!!!RESTRIÇÂO DE TAMANHO DE PASSO
!		passo é x_ (em que atualizei para x_nu apenas as fases do grupo 1, e mantive x_ol para as demais)
		step(:,:) = x_(:,:) - x_ol(:,:)			!baseada no valor de x que entra como argumento (inout)
		max_step(:,:) = dabs(.01d0*x_ol(:,:))+1.d-4 !termo independente define o mínimo tamanho máximo de passo
		rel_step(:,:) = dabs(step(:,:))/max_step(:,:)
		reducing_factor = maxval(rel_step(:,:))
		met_num_par_val(1) = reducing_factor
		if( reducing_factor > 1.d0 ) then
!			print*, x_ol(:,:), maxval(step(:,:)), maxval(step(:,:))/reducing_factor
!			print*, x_ol(1,1), step(1,1), step(1,1)/reducing_factor
			step(:,:) = step(:,:)/reducing_factor
			x_(:,:) = x_ol(:,:) + step(:,:)
		else
!			print*, maxval(step(:,:)) !, step(1,1)
!			print*, x_ol(1,1), step(1,1), step(1,1)
		endif
!		!redução de passo por exaustão além do reducing_factor
!		if (loop_x > 1000) then
!			x_(:,:) = x_ol(:,:) + step(:,:)/(dfloat(loop_x)/1.d3)
!		endif
		
		!!!Conclusão da rotina
!		print*, x_, 'update_x_restricted'

		!normalizar se necessário
		do j = 1, nphas
			sum1 = 0.d0
			do i = 1, ncomp
!			print*, i,j,x_(i,j), __FILE__, __LINE__
				sum1 = sum1 + x_(i,j)
			enddo !i = 1, ncomp
			!
			if((sum1-1)>1d-5)then
				!***usar alguma flag de acompanhamento aqui
!				print*, 'sum(x_) ', ' nphas', j, ' =', sum1, ' @', __FILE__, __LINE__, phase_list(ph_index(j))%phase%phase_cond
				do i = 1, ncomp
					x_(i,j) = x_(i,j)/sum1
	!				print*, 'n',i,j,x_(i,j), __FILE__, __LINE__
				enddo !i = 1, ncomp
			endif
		enddo !j = 1, nphas
		
!atualizar o x da fase com o x do algoritmo (que só mudou para as do grupo 1 (fluidas e gelo))
		do j = 1, nphas
			phase_list(ph_index(j))%phase%x(:) = x_(:,j)
		enddo
		
!		print*, x_, 'update_x_normalized'
		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'update_x'; do j = 1, nphas; write(20,*) 'x', x_(:,j), phase_list(ph_index(j))%phase%phase_cond; enddo; endif
		
	end subroutine update_x
	
	function NLE_E_j(z_,K_,th_,al_,j_)
		real(8) :: NLE_E_j
		real(8) :: z_(ncomp),K_(ncomp,nphas),th_(nphas),al_(nphas)
		integer :: j_
		integer :: i, j
		real(8) :: sum1, sum2
		sum1 = 0.d0
		do i = 1, ncomp
			sum2 = 0.d0
			do j = 1, nphas-1
				sum2 = sum2 + al_(j)*((K_(i,j)*dexp(th_(j)))-1.d0)
			enddo !j = 1, nphas-1
			sum1 = sum1 + (z_(i)*((K_(i,j_)*dexp(th_(j_)))-1.d0))/(1.d0+sum2)
		enddo !i = 1, ncomp
		NLE_E_j = sum1
		return
	end function NLE_E_j
	
	function dFEjdtheta(z_,K_,th_,al_,j_,jj_) !derivada da função F(j) em relação ao parâmetro theta(jj)
		real(8) :: dFEjdtheta
		real(8) :: z_(ncomp),K_(ncomp,nphas),th_(nphas),al_(nphas)
		integer :: j_, jj_
		integer :: i, j
		real(8) :: sum1, sum2, term1, term2
		!termo1
		sum1 = 0.d0
		do i = 1, ncomp
			sum2 = 0.d0
			do j = 1, nphas-1
				sum2 = sum2 + al_(j)*((K_(i,j)*dexp(th_(j)))-1.d0)
			enddo !j = 1, nphas-1
			sum1 = sum1 + (z_(i)*((K_(i,j_)*dexp(th_(j_)))-1.d0)*al_(jj_)*(K_(i,jj_)*dexp(th_(jj_))))/((1.d0+sum2)**2)
		enddo !i = 1, ncomp
		term1 = -1.d0*sum1
		if( j_.NE.jj_ ) then
			dFEjdtheta = term1
		else !(j == jj_)
			!termo2
			sum1 = 0.d0
			do i = 1, ncomp
				sum2 = 0.d0
				do j = 1, nphas-1
					sum2 = sum2 + al_(j)*((K_(i,j)*dexp(th_(j)))-1.d0)
				enddo !j = 1, nphas-1
				sum1 = sum1 + (z_(i)*K_(i,jj_)*dexp(th_(jj_)))/(1.d0+sum2)
			enddo !i = 1, ncomp
			term2 = sum1
			dFEjdtheta = term1 + term2
		endif
		return
	end function dFEjdtheta
	
	function dFEjdalpha(z_,K_,th_,al_,j_,jj_)
		real(8) :: dFEjdalpha
		real(8) :: z_(ncomp),K_(ncomp,nphas),th_(nphas),al_(nphas)
		integer :: j_, jj_
		integer :: i, j
		real(8) :: sum1, sum2
		sum1 = 0.d0
		do i = 1, ncomp
			sum2 = 0.d0
			do j = 1, nphas-1
				sum2 = sum2 + al_(j)*((K_(i,j)*dexp(th_(j)))-1.d0)
			enddo !j = 1, nphas-1
			sum1 = sum1 + (z_(i)*((K_(i,j_)*dexp(th_(j_)))-1.d0)*((K_(i,jj_)*dexp(th_(jj_)))-1.d0))/((1.d0+sum2)**2)
		enddo !i = 1, ncomp
		dFEjdalpha = -1.d0*sum1
		return
	end function dFEjdalpha
	
!	function NLE_stability(th_,al_) !para um par de valores th_ e al_ de uma fase específica
!		real(8) :: NLE_stability
!		real(8) :: th_, al_
!		NLE_stability = (th_*al_)/(th_+al_)
!		return
!	end function NLE_stability
	
!	function dFsdtheta(al_,th_,j_,jj_) !derivada da função indice j_ em relação a variável indice jj_, para j_ == jj_ ( a derivada é zero se j_ /= jj_)
!		real(8) :: dFsdtheta
!		real(8) :: al_, th_
!		integer :: j_, jj_
!		
!		if( j_ .NE. jj_ ) then
!			dFsdtheta = 0.d0
!			dFsdtheta = 0.d0
!		else !(j == jj)
!!			dFsdtheta = (al_)/(th_+al_) !não está incluso o segundo termo da derivada, o qual só é diferente de zero de tanto al_ quanto th_ forem diferentes de zero
!										!se al_ == th_ == número pequeno, então essa fórmula dá 1/2, enquanto a fórmula completa daria 1/4.
!			dFsdtheta = (al_)/(th_+al_) - (al_*th_)/((th_+al_)**2) !equação comleta, se al_ ou th_ são zero, ela dá 0 ou 1,
!																	!se al_ = th_ = too small, ela dá 1/4 em vez de 1/2, mas fazendo isso parece estragar a convergência
!		endif
!		return
!	end function dFsdtheta
	
!	function dFsdalpha(al_,th_,j_,jj_)
!		real(8) :: dFsdalpha
!		real(8) :: al_, th_
!		integer :: j_, jj_
!		
!		if( j_ .NE. jj_ ) then
!			dFsdalpha = 0.d0
!			dFsdalpha = 0.d0
!		else !(j == jj)
!!			dFsdalpha = (th_)/(th_+al_)
!			dFsdalpha = (th_)/(th_+al_) - (al_*th_)/((th_+al_)**2) !equação comleta, se al_ ou th_ são zero, ela dá 0 ou 1,
!!			write(93,*) dFsdalpha 									!se al_ = th_ = too small, ela dá 1/4 em vez de 1/2, mas fazendo isso parece estragar a convergência
!		endif
!		return
!	end function dFsdalpha
	
	subroutine Jacobiano_flash_PT(Jac_,z_,K_,th_,al_,al_flick_)		!SINGLE-SPEC
		!ARGUMENTS
		real(8), allocatable :: Jac_(:,:)
		real(8), allocatable :: z_(:)
		real(8), allocatable :: K_(:,:)
		real(8), allocatable :: th_(:)
		real(8), allocatable :: al_(:)
		integer, allocatable :: al_flick_(:)
		!LOCAL
		integer :: j, jj
		integer :: NLESsize

		NLESsize = nphas-1

		do j = 1, nphas-1
			do jj = 1, nphas-1
				if(al_flick_(jj)==1) then
					Jac_(j,jj) = dFEjdalpha(z_,K_,th_,al_,j,jj)
				else !elseif (al_flick_(jj)==0) then
					Jac_(j,jj) = dFEjdtheta(z_,K_,th_,al_,j,jj)
				endif
			enddo !j = 1, nphas-1
		enddo !jj = 1, nphas-1
	end subroutine Jacobiano_flash_PT
	
	subroutine Jacobiano_flash_B(Jac_,z_,K_,th_,al_,al_flick_)		!MULTI-SPEC
		!ARGUMENTS
		real(8), allocatable :: Jac_(:,:)
		real(8), allocatable :: z_(:)
		real(8), allocatable :: K_(:,:)
		real(8), allocatable :: th_(:)
		real(8), allocatable :: al_(:)
		integer, allocatable :: al_flick_(:)
		!LOCAL
		integer :: i, j, jj
		integer :: NLESsize
		real(8) :: dFEjdT(nphas-1), dKijdT(ncomp,nphas), dFEjdP(nphas-1), dKijdP(ncomp,nphas)
		real(8) :: K_mais_T(ncomp,nphas), K_menos_T(ncomp,nphas), T_mais, T_menos
		real(8) :: K_mais_P(ncomp,nphas), K_menos_P(ncomp,nphas), P_mais, P_menos
		real(8) :: passinho_T, passinho_P
		real(8) :: sum1, sum2
		real(8) :: x_l(ncomp,nphas)
		!
		!ALLOCATION
		!Jac_ is already allocated as in/out argument
		NLESsize = nphas-1
!		allocate()
		!IMPLEMENTATION
		!FIRST COLUMN, ALL LINES
		select case(flash_spec(1))
		case('P_') !pressão fixa, temperatura livre
			!derivadas das NLE_ E e S em relação à temperatura, usando derivada numérica para o termo dKdT
			!calcula dKijdT
			passinho_T = 1.d-3
			T_mais = Tflash + passinho_T; !print*, 'T_mais', T_mais
			T_menos = Tflash - passinho_T
	!		print*, 'derivada numérica'
	!		print*, 'K_', K_, 'K_'; !pause
			call update_K(K_mais_T,T_mais,Pflash,x_l) !esse x_l irá conter a atualização da composição da fase hidrato, o que não é interessante de guardar pois esse passo não representa uma mudança da condição do sistema, mas apenas um cálculo auxiliar para a derivada numérica
			call update_K(K_menos_T,T_menos,Pflash,x_l)
	!		print*, 'K_mais', K_mais, 'K_mais'; !pause
	!		print*, 'K_menos', K_menos, 'K_menos'; !pause
			dKijdT(:,:) = (K_mais_T(:,:)-K_menos_T(:,:))/(2.d0*passinho_T)
	!		print*, dkijdt(:,:), 'dkijdt(:,:)'; pause
			!
			!calcula dFEjdT
			do j = 1, nphas-1
			dFEjdT(j) = 0.d0
				do i = 1, ncomp
					sum1 = 0.d0
					sum2 = 0.d0
					do jj = 1, nphas-1
						sum1 = sum1 + al_(jj)*dKijdT(i,jj)*dexp(th_(jj))
						sum2 = sum2 + al_(jj)*(K_(i,jj)*dexp(th_(jj))-1.d0)
					enddo !jj = 1, nphas-1
					!obs: essa equação tem a fórmula explicita do x embutida, para não usar a variavel x desatualizada
					!está como deduzida diretamente da NLE_E_j (como no meu texto), sem substituir por x para simplificar a escrita dela (como no texto do ballard)
					dFEjdT(j) = dFEjdT(j) + (z_(i)/(1.d0+sum2))*((dKijdT(i,j)*dexp(th_(j)))-(((K_(i,j)*dexp(th_(j))-1.d0)*sum1)/(1.d0+sum2)))
				enddo !i = 1, ncomp
			enddo !j = 1, nphas-1
		case('T_')
			!derivadas das NLE_ E e S em relação à pressão, usando derivada numérica para o termo dKdP
			!calcula dKijdP
			passinho_P = 1.d-5*Pflash
			P_mais = Pflash + passinho_P; !print*, 'P_mais', P_mais
			P_menos = Pflash - passinho_P
	!		print*, 'derivada numérica'
	!		print*, 'K_', K_, 'K_'; !pause
			call update_K(K_mais_P,Tflash,P_mais,x_l) !esse x_l irá conter a atualização da composição da fase hidrato, o que não é interessante de guardar pois esse passo não representa uma mudança da condição do sistema, mas apenas um cálculo auxiliar para a derivada numérica
			call update_K(K_menos_P,Tflash,P_menos,x_l)
!*** esses update_K com T_mais/menos, ou P_mais/menos irão mudar o phi da fase hidrato, e esse phi velho é utilizado no loop de SS para atualizar o Fug
!	como a perturbação é muito pequena não deve influenciar na funcionalidade do SS, mas é válido notar essa imperfeição na implementação
!	para corrigir eu poderia
!	1) faze backup das variáveis das fases e devolver para elas após a pertrubação (o que não seria muito trabalhoso, mas também, provavelmente desnecessário)
!	2) fazer essas derivadas parciais usando uma classe slave (o que seria bastante trabalhoso para implementar as declarações e inicializações, e aumentaria o consumo de memória do código)

	!		print*, 'K_mais', K_mais, 'K_mais'; !pause
	!		print*, 'K_menos', K_menos, 'K_menos'; !pause
			dKijdP(:,:) = (K_mais_P(:,:)-K_menos_P(:,:))/(2.d0*passinho_P)
	!		print*, dkijdt(:,:), 'dkijdt(:,:)'; pause
			!
			!calcula dFEjdP
			do j = 1, nphas-1
			dFEjdP(j) = 0.d0
				do i = 1, ncomp
					sum1 = 0.d0
					sum2 = 0.d0
					do jj = 1, nphas-1
						sum1 = sum1 + al_(jj)*dKijdP(i,jj)*dexp(th_(jj))
						sum2 = sum2 + al_(jj)*(K_(i,jj)*dexp(th_(jj))-1.d0)
					enddo !jj = 1, nphas-1
					!obs: essa equação tem a fórmula explicita do x embutida, para não usar a variavel x desatualizada
					!está como deduzida diretamente da NLE_E_j (como no meu texto), sem substituir por x para simplif
					dFEjdP(j) = dFEjdP(j) + (z_(i)/(1.d0+sum2))*((dKijdP(i,j)*dexp(th_(j)))-(((K_(i,j)*dexp(th_(j))-1.d0)*sum1)/(1.d0+sum2)))
				enddo !i = 1, ncomp
			enddo !j = 1, nphas-1
		case default
			print*, 'unexpected flash_spec(1) case argument in ', __FILE__, __LINE__; stop
		end select
		!
!IMPLICIT LOOP
		select case(flash_spec(1))
		case('P_') !pressão fixa, temperatura livre
			Jac_(1:nphas-1,1) = (/(dFEjdT(j),j=1,nphas-1)/) !dFEdT
		case('T_') !temperatura fixa, pressão livre.
			Jac_(1:nphas-1,1) = (/(dFEjdP(j),j=1,nphas-1)/) !dFEdP
		case default
			print*, 'unexpected flash_spec(1) case argument in ', __FILE__, __LINE__; stop
		end select

!SECOND TO LAST COLUMN, ALL LINES

		do j = 1, nphas-1
			do jj = 2, nphas-1
				if(al_flick_(jj)==1) then
					Jac_(j,jj) = dFEjdalpha(z_,K_,th_,al_,j,jj)
				else !elseif (al_flick_(jj)==0) then
					Jac_(j,jj) = dFEjdtheta(z_,K_,th_,al_,j,jj)
				endif
			enddo !j = 1, nphas-1
		enddo !jj = 2, nphas-1

	end subroutine Jacobiano_flash_B
	
	subroutine Jacobiano_flash_BB(Jac_,z_,K_,th_,al_,al_flick_)		!MULTI-SPEC
		!ARGUMENTS
		real(8), allocatable :: Jac_(:,:)
		real(8), allocatable :: z_(:)
		real(8), allocatable :: K_(:,:)
		real(8), allocatable :: th_(:)
		real(8), allocatable :: al_(:)
		integer, allocatable :: al_flick_(:)
		!LOCAL
		integer :: i, j, jj
		integer :: NLESsize
		real(8) :: dFEjdT(nphas-1), dKijdT(ncomp,nphas), dFEjdP(nphas-1), dKijdP(ncomp,nphas)
		real(8) :: K_mais_T(ncomp,nphas), K_menos_T(ncomp,nphas), T_mais, T_menos
		real(8) :: K_mais_P(ncomp,nphas), K_menos_P(ncomp,nphas), P_mais, P_menos
		real(8) :: passinho_T, passinho_P
		real(8) :: sum1, sum2
		real(8) :: x_l(ncomp,nphas)
		!
		!ALLOCATION
		!Jac_ is already allocated as in/out argument
		NLESsize = nphas-1
!		allocate()
		!IMPLEMENTATION
		!FIRST COLUMN, ALL LINES
		!derivadas das NLE_ E e S em relação à temperatura, usando derivada numérica para o termo dKdT
		!calcula dKijdT
		passinho_T = 1.d-3
		T_mais = Tflash + passinho_T; !print*, 'T_mais', T_mais
		T_menos = Tflash - passinho_T
	!	print*, 'derivada numérica'
	!	print*, 'K_', K_, 'K_'; !pause
		call update_K(K_mais_T,T_mais,Pflash,x_l)
		call update_K(K_menos_T,T_menos,Pflash,x_l)
	!	print*, 'K_mais', K_mais, 'K_mais'; !pause
	!	print*, 'K_menos', K_menos, 'K_menos'; !pause
		dKijdT(:,:) = (K_mais_T(:,:)-K_menos_T(:,:))/(2.d0*passinho_T)
!		print*, dkijdt(:,:), 'dkijdt(:,:)'; pause
		!
		!calcula dFEjdT
		do j = 1, nphas-1
		dFEjdT(j) = 0.d0
			do i = 1, ncomp
				sum1 = 0.d0
				sum2 = 0.d0
				do jj = 1, nphas-1
					sum1 = sum1 + al_(jj)*dKijdT(i,jj)*dexp(th_(jj))
					sum2 = sum2 + al_(jj)*(K_(i,jj)*dexp(th_(jj))-1.d0)
				enddo !jj = 1, nphas-1
				!obs: essa equação tem a fórmula explicita do x embutida, para não usar a variavel x desatualizada
				dFEjdT(j) = dFEjdT(j) + (z_(i)/(1.d0+sum2))*((dKijdT(i,j)*dexp(th_(j)))-(((K_(i,j)*dexp(th_(j))-1.d0)*sum1)/(1.d0+sum2)))
			enddo !i = 1, ncomp
		enddo !j = 1, nphas-1

		!derivadas das NLE_ E e S em relação à pressão, usando derivada numérica para o termo dKdP
		!calcula dKijdP
		passinho_P = 1.d-5*Pflash
		P_mais = Pflash + passinho_P; !print*, 'P_mais', P_mais
		P_menos = Pflash - passinho_P
!		print*, 'derivada numérica'
!		print*, 'K_', K_, 'K_'; !pause
		call update_K(K_mais_P,Tflash,P_mais,x_l)
		call update_K(K_menos_P,Tflash,P_menos,x_l)
!		print*, 'K_mais', K_mais, 'K_mais'; !pause
!		print*, 'K_menos', K_menos, 'K_menos'; !pause
		dKijdP(:,:) = (K_mais_P(:,:)-K_menos_P(:,:))/(2.d0*passinho_P)
!		print*, dkijdt(:,:), 'dkijdt(:,:)'; pause
		!
		!calcula dFEjdP
		do j = 1, nphas-1
		dFEjdP(j) = 0.d0
			do i = 1, ncomp
				sum1 = 0.d0
				sum2 = 0.d0
				do jj = 1, nphas-1
					sum1 = sum1 + al_(jj)*dKijdP(i,jj)*dexp(th_(jj))
					sum2 = sum2 + al_(jj)*(K_(i,jj)*dexp(th_(jj))-1.d0)
				enddo !jj = 1, nphas-1
				!obs: essa equação tem a fórmula explicita do x embutida, para não usar a variavel x desatualizada
				dFEjdP(j) = dFEjdP(j) + (z_(i)/(1.d0+sum2))*((dKijdP(i,j)*dexp(th_(j)))-(((K_(i,j)*dexp(th_(j))-1.d0)*sum1)/(1.d0+sum2)))
			enddo !i = 1, ncomp
		enddo !j = 1, nphas-1
		!
		!PREENCHE JAC:
		!derivada para T
		Jac_(1:nphas-1,1) = (/(dFEjdT(j),j=1,nphas-1)/) !dFEdT
		Jac_(nphas:NLESsize,1) = 0.d0 !dFSdT
		!derivada para P
		Jac_(1:nphas-1,2) = (/(dFEjdP(j),j=1,nphas-1)/) !dFEdP
		Jac_(nphas:NLESsize,2) = 0.d0 !dFSdP

!IMPLICIT LOOP

		do j = 1, nphas-1
			do jj = 3, nphas-1
				if(al_flick_(jj)==1) then
					Jac_(j,jj) = dFEjdalpha(z_,K_,th_,al_,j,jj)
				else !elseif (al_flick_(jj)==0) then
					Jac_(j,jj) = dFEjdtheta(z_,K_,th_,al_,j,jj)
				endif
			enddo !j = 1, nphas-1
		enddo !jj = 2, nphas-1
!
	end subroutine Jacobiano_flash_BB
	
	subroutine newton_step_PT(z_,K_,th_,al_,al_flick_)		!SINGLE-SPEC
											!atualiza as variaveis de estabilidade e frações de fase com 1 passo de newton
		!ARGUMENTS
		real(8), allocatable :: z_(:)
		real(8), allocatable :: K_(:,:)
		real(8), allocatable :: th_(:)
		real(8), allocatable :: al_(:)
		integer, allocatable :: al_flick_(:)
		!LOCAL
		real(8), allocatable :: IV_(:)
		real(8), allocatable ::  NLE_(:)
		real(8), allocatable ::  Jac_(:,:)
		real(8), allocatable ::  Jac_inv(:,:)
		integer :: j, jj
		integer :: NLESsize
		real(8) :: sum1, reducing_factor, aux
		real(8), allocatable :: step(:)
		real(8), allocatable :: max_step(:)
		real(8), allocatable :: rel_step(:)
		
		!ALLOCATE
		NLESsize = nphas-1
		allocate(IV_(NLESsize))
		allocate(NLE_(NLESsize))
		allocate(Jac_(NLESsize,NLESsize))
		allocate(Jac_inv(NLESsize,NLESsize))
		allocate(step(NLESsize))
		allocate(max_step(NLESsize))
		allocate(rel_step(NLESsize))

		!IMPLEMENTATION
		do j = 1, nphas-1
			if(al_flick_(j) == 1) then
					IV_(j) = al_(j)
				else !elseif(al_flick_(j) == 0)
					IV_(j) = th_(j)
				endif
		enddo !j = 1, nphas-1
		
		do j = 1, nphas-1
			NLE_(j) = NLE_E_j(z_,K_,th_,al_,j)
		enddo !j = 1, nphas
		
		call Jacobiano_flash_PT(Jac_,z_,K_,th_,al_,al_flick_)

		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'JAC'; do j = 1, NLESsize; write(20,*) Jac_(j,:); enddo; endif
		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'JAC'; do j = 1, NLESsize; print*, Jac_(j,:); enddo; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
		!
		call matrix_inversion(NLESsize,Jac_,Jac_inv,aux)
		
		met_num_par_val(3) = aux
		!
		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'JACinv'; do j = 1, NLESsize; write(20,*) Jac_inv(j,:); enddo; endif
		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'JACinv'; do j = 1, 2*nphas-2; print*, Jac_inv(j,:); enddo; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
		!
		step = matmul(Jac_inv,NLE_)
		!
		!!!!RESTRIÇÂO DE TAMANHO DE PASSO
		max_step(:) = dabs(.1d0*IV_(:))+1.d-2 !termo independente define o mínimo tamanho máximo de passo
		rel_step(:) = dabs(step(:))/max_step(:)
		reducing_factor = maxval(rel_step(:))
!		print*, step
!		print*, max_step
		met_num_par_val(2) = reducing_factor
		if( reducing_factor > 1.d0 ) then
			step(:) = step(:)/reducing_factor
		endif
		!!!Conclusão da rotina
		IV_(:) = IV_(:) - step(:)
		do j = 1, nphas-1
			if(al_flick_(j)==1) then
				al_(j) = IV_(j)
			else !elseif (al_flick_(j)==0) then
				th_(j) = IV_(j)
			endif
		enddo !j = 1, nphas-1
		
		th_(nphas) = 0.d0 !UNNECESSARY LINE AS TH_(NPHAS) IS NEVER RECALCULATED
		
		sum1 = 0.d0
		do j = 1, nphas-1
			sum1 = sum1 + al_(j)
		enddo !j = 1, nphas-1
		al_(nphas) = 1.d0 - sum1
	end subroutine newton_step_PT
	
	subroutine newton_step_B(z_,K_,th_,al_,al_flick_) !atualiza as variaveis de estabilidade e frações de fase com 1 passo de newton
		!ARGUMENTS
		real(8), allocatable :: z_(:)
		real(8), allocatable :: K_(:,:)
		real(8), allocatable :: th_(:)
		real(8), allocatable :: al_(:)
		integer, allocatable :: al_flick_(:)
		
		!LOCAL
		integer :: NLESsize
		real(8), allocatable :: IV_(:)
		real(8), allocatable :: NLE_(:)
		real(8), allocatable :: Jac_(:,:)
		real(8), allocatable :: Jac_inv(:,:)
		integer :: j, jj
		real(8), allocatable :: step(:)
		real(8), allocatable :: max_step(:)
		real(8), allocatable :: rel_step(:)
		real(8) :: sum1
		real(8) :: reducing_factor
		real(8) :: TPstep
		real(8) :: aux
		!ALLOCATION
		NLESsize = nphas-1
		allocate (IV_(NLESsize))
		allocate (NLE_(NLESsize))
		allocate (Jac_(NLESsize,NLESsize))
		allocate (Jac_inv(NLESsize,NLESsize))
		allocate (step(NLESsize))
		allocate (max_step(NLESsize))
		allocate (rel_step(NLESsize))
		
		!implementation
		select case(flash_spec(1))
		case('P_') !pressão fixa, temperatura livre
			IV_(1) = Tflash
		case('T_') !temperatura fixa, pressão livre
			IV_(1) = Pflash
		case default
			print*, 'unexpected flash_spec(1) case argument in ', __FILE__, __LINE__; stop
		end select
!!DEBUG
!print*, 'DEBUG', __FILE__, __LINE__
!th_(1:nphas) = (/(10*j,j=1,nphas)/)
!al_(1:nphas) = (/(j,j=1,nphas)/)
!IV_ = 0.d0
!!DEBUG_END

!EQUIVALENT IMPLICIT LOOP
		do j = 2, nphas-1
			if(al_flick_(j) == 1) then
					IV_(j) = al_(j)
				else !elseif(al_flick_(j) == 0)
					IV_(j) = th_(j)
				endif
		enddo !j = 1, nphas-1

!!DEBUG
!print*, 'DEBUG', __FILE__, __LINE__
!print*, IV_(:), __FILE__, __LINE__; stop
!!DEBUG END

!!DEBUG
!print*, 'DEBUG', __FILE__, __LINE__
!print*, NLE_(:), __FILE__, __LINE__
!!DEBUG END

!EQUIVALENT IMPLLICIT LOOP
		NLE_(1:nphas-1)=(/(NLE_E_j(z_,K_,th_,al_,j),j=1,nphas-1)/)

!!DEBUG
!print*, 'DEBUG', __FILE__, __LINE__
!print*, NLE_(:), __FILE__, __LINE__; stop
!!DEBUG END

		!
		call Jacobiano_flash_B(Jac_,z_,K_,th_,al_,al_flick_)
		
		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'JAC'; do j = 1, NLESsize; write(20,*) Jac_(j,:); enddo; endif
		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'JAC'; do j = 1, NLESsize; print*, Jac_(j,:); enddo; endif
!		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
		!
		call matrix_inversion(NLESsize,Jac_,Jac_inv,aux)
		
		met_num_par_val(3) = aux
		!
		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'JACinv'; do j = 1, NLESsize; write(20,*) Jac_inv(j,:); enddo; endif
		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'JACinv'; do j = 1, NLESsize; print*, Jac_inv(j,:); enddo; endif
!		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif

		step(:) = matmul(Jac_inv,NLE_)
		!!!!RESTRIÇÂO DE TAMANHO DE PASSO
		max_step(:) = dabs(.1d0*IV_(:))+1.d-2 !esse termo independente define o mínimo tamanho máximo de passo
		rel_step(:) = dabs(step(:))/max_step(:)
		reducing_factor = maxval(rel_step(:))
		met_num_par_val(2) = reducing_factor
		if( reducing_factor > 1.d0 ) then
			step(:) = step(:)/reducing_factor
		endif
!		print*, step
!		print*, max_step
!		print*, reducing_factor, 'NR'
		!redução de passo por exaustão
!		if (loop_x > 1000) then
!			step(:) = step(:)/(dfloat(loop_x)/1.d4)
!		endif

		!!!Conclusão da rotina
		IV_(:) = IV_(:) - step(:)

!OUTPUT TRANSLATION
!FOR TH and AL
		do j = 2, nphas-1
			if(al_flick_(j)==1) then
				al_(j) = IV_(j)
			else !elseif (al_flick_(j)==0) then
				th_(j) = IV_(j)
			endif
		enddo !j = 1, nphas-1

!FOR T OR P
		!restrição de passo adicional em T e P
		select case(flash_spec(1))
		case('P_') !pressão fixa, temperatura livre
			!restrição do passo em T
			TPstep = Tflash-IV_(1)
			if( dabs(TPstep) > 1.d0 ) then
				TPstep = 1.d0*(TPstep/dabs(TPstep)) !step máximo vezes direção
			endif
			Tflash = Tflash - TPstep
		case('T_') !temperatura fixa, pressão livre
			!restrição do passo em P
			TPstep = Pflash-IV_(1)
			if( dabs(TPstep) > .05d0*Pflash ) then
				TPstep = .05d0*Pflash*(TPstep/dabs(TPstep)) !step máximo vezes direção
			endif
			Pflash = Pflash - TPstep
		case default
			print*, 'unexpected flash_spec(1) case argument in ', __FILE__, __LINE__; stop
		end select
		!
		th_(nphas) = 0.d0 !UNNECESSARY LINE, TH_(NPHAS) IS NEVER MODIFIED, AND ALWAYS ZERO
		sum1 = 0.d0
		do j = 1, nphas-1
			sum1 = sum1 + al_(j)
		enddo !j = 1, nphas-1
		al_(nphas) = 1.d0 - sum1 !NECESSARY LINE FOR UPDATING REF's PHASE FRACTION
	end subroutine newton_step_B
	
	subroutine newton_step_BB(z_,K_,th_,al_,al_flick_) !atualiza as variaveis de estabilidade e frações de fase com 1 passo de newton
		!ARGUMENTS
		!ARGUMENTS
		real(8), allocatable :: z_(:)
		real(8), allocatable :: K_(:,:)
		real(8), allocatable :: th_(:)
		real(8), allocatable :: al_(:)
		integer, allocatable :: al_flick_(:)
		!LOCAL
		integer :: NLESsize
		real(8), allocatable :: IV_(:), NLE_(:), Jac_(:,:), Jac_inv(:,:)
		integer :: j, jj
		real(8), allocatable :: step(:), max_step(:), rel_step(:)
		real(8) :: sum1, reducing_factor, TPstep, aux
		!allocation
		NLESsize = nphas-1
		allocate (IV_(NLESsize))
		allocate (NLE_(NLESsize))
		allocate (Jac_(NLESsize,NLESsize))
		allocate (Jac_inv(NLESsize,NLESsize))
		allocate (step(NLESsize))
		allocate (max_step(NLESsize))
		allocate (rel_step(NLESsize))
		!IMPLEMENTATION
		!T
		IV_(1) = Tflash
		!P
		IV_(2) = Pflash
		!Remaining
		do j = 3, nphas-1
			if(al_flick_(j) == 1) then
					IV_(j) = al_(j)
				else !elseif(al_flick_(j) == 0)
					IV_(j) = th_(j)
				endif
		enddo !j = 1, nphas-1

!IMPLICIT LOOP for NLEs
		NLE_(1:nphas-1)=(/(NLE_E_j(z_,K_,th_,al_,j),j=1,nphas-1)/)
!
		call Jacobiano_flash_BB(Jac_,z_,K_,th_,al_,al_flick_)
!
		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'JAC'; do j = 1, NLESsize; write(20,*) Jac_(j,:); enddo; endif
		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'JAC'; do j = 1, NLESsize; print*, Jac_(j,:); enddo; endif
!		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif

		call matrix_inversion(NLESsize,Jac_,Jac_inv,aux)
		
		met_num_par_val(3) = aux
!
		!ACOMPANHAMENTOS:
		if (ACOMPANHAMENTOLOG == 1) then; write(20,*) 'JACinv'; do j = 1, NLESsize; write(20,*) Jac_inv(j,:); enddo; endif
		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'JACinv'; do j = 1, NLESsize; print*, Jac_inv(j,:); enddo; endif
!		if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused endif
!
		step(:) = matmul(Jac_inv,NLE_)
		!!!!RESTRIÇÂO DE TAMANHO DE PASSO
		max_step(:) = dabs(.1d0*IV_(:))+1.d-2 !esse termo independente define o mínimo tamanho máximo de passo
		rel_step(:) = dabs(step(:))/max_step(:)
		reducing_factor = maxval(rel_step(:))
		met_num_par_val(2) = reducing_factor
		if( reducing_factor > 1.d0 ) then
			step(:) = step(:)/reducing_factor
		endif
!		print*, step
!		print*, max_step
!
		!!!Conclusão da rotina
		IV_(:) = IV_(:) - step(:)
!OUTPUT TRANSLATION
!FOR TH and AL
		do j = 3, nphas-1
			if(al_flick_(j)==1) then
				al_(j) = IV_(j)
			else !elseif (al_flick_(j)==0) then
				th_(j) = IV_(j)
			endif
		enddo !j = 1, nphas-1
!FOR T AND P
!
		!restrição de passo adicional em T e P
		!step em T
		TPstep = Tflash-IV_(1) !obs: Tflash contem o valor da variável antes do step e Sol contém o valor após o step com a restrição de passo de 25%
		if( dabs(TPstep) > 1.d0 ) then
			TPstep = 1.d0*(TPstep/dabs(TPstep)) !step máximo vezes direção
		endif
		Tflash = Tflash - TPstep
		!
		!step em P limitado a no máximo 5%
		TPstep = Pflash-IV_(2)  !obs: Pflash contem o valor da variável antes do step e Sol contém o valor após o step com a restrição de passo de 25%
		if( dabs(TPstep) > .05d0*Pflash ) then
			TPstep = .05d0*Pflash*(TPstep/dabs(TPstep)) !step máximo vezes direção
		endif
		Pflash = Pflash - TPstep
!
!UPDATE LAST PHASE FRACTION (DEPENDENT) USING ITS OWN DEFINIION EQUATION:
		th_(nphas) = 0.d0 !UNNECESSARY LINE, TH_(NPHAS) IS NEVER MODIFIED, AND ALWAYS ZERO
		sum1 = 0.d0
		do j = 1, nphas-1
			sum1 = sum1 + al_(j)
		enddo !j = 1, nphas-1
		al_(nphas) = 1.d0 - sum1 !NECESSARY LINE FOR UPDATING REF's PHASE FRACTION
!
	end subroutine newton_step_BB
	
	subroutine matrix_inversion(ND,A,B,C)
		IMPLICIT NONE
		INTEGER i, j						! contadores
		INTEGER, INTENT(IN)  :: ND			! dimensão das matrizes
		REAL(8), INTENT(IN)  :: A(ND,ND)	! matriz original
		REAL(8), INTENT(OUT) :: B(ND,ND)	! matriz invertida
		REAL(8), INTENT(OUT) :: C			! min(aux)
		REAL(8) M(ND,2*ND), V(2*ND), aux	! variáveis locais auxiliares
		!Inicializando a matriz M <-- [A I]
		M = 0.d0
		M(1:ND,1:ND) = A
		FORALL (i=1:ND) M(i,i+ND) = 1.d0
		!inicia o procedimento de inversão
		DO i = 1,ND
			!pivotação
			aux = M(i,i)
			DO j = i+1,ND
				IF (DABS(M(j,i)) <= DABS(aux)) CYCLE
				V(:) = M(j,:)
				M(j,:) = M(i,:)
				M(i,:) = V(:)
				aux = M(i,i)
			enddo
			!verifica a possibilidade da matriz ser não inversível
			if( i==1 ) then
				C=dabs(aux)
			else
				C=min(C,dabs(aux))
			endif
			IF (DABS(aux) < 1.d-10) then
!				WRITE(*,*) 'Possibilidade de matriz NAO inversivel'
!				print*, dabs(aux), Tflash, Pflash
!				stop
			endif
			!utiliza o método de Gauss para a inversão
			M(i,:) = M(i,:)/aux
			DO j = 1,ND
				IF (i == j) CYCLE
				aux = M(j,i)
				M(j,:) = M(j,:) - aux*M(i,:)
			enddo
		enddo
		!Fornecendo a matriz inversa: M --> [I B]
		B = M(:,ND+1:2*ND)
		RETURN
	end subroutine matrix_inversion
	
	function re_set_theta(z_,K_,th_,al_,j_)	!algumas estimativas de composição da fase de referência geram thetas negativos
									!para as fases instáveis removidas
									!nesse caso seria melhor calcular o theta por outro meio:
									!aqui testamos o calculo do theta novo a partir da equação de definição, usando cada componente
									!foi observado no computador chinês mathcad que usando-se menor valor positivo de theta calculado
									!(em 1 exemplo, ainda a verificar-se em mais casos)
									!gerava-se posteriormente uma menor função objetivo NLE_EK para a fase manipulada
		real(8) :: z_(ncomp),K_(ncomp,nphas),th_(nphas),al_(nphas)
		real(8) :: re_set_theta
!		real(8) :: theta_k_25,  !baseado na eq 25 do gupta 1991
!		real(8) :: theta_k_def(ncomp) !baseado na definição por componente
		real(8) :: x_l(ncomp,nphas) !composição teórica atualizada localmente no algoritmo só para uso nessa função
		integer :: i, j, j_
		real(8) :: sum1, term1(ncomp)
		

		
		if (ACOMPANHAMENTOLOG == 1) then
			write(20,*) __FILE__, __LINE__
		endif
		
			do i = 1, ncomp
				sum1 = 0.d0
				do j = 1, nphas-1
					sum1 = sum1 + al_(j)*((K_(i,j)*dexp(th_(j)))-1.d0)
				enddo !j = 1, nphas
				x_l(i,nphas) = z_(i)/(1.d0+sum1) 
				!DEBUG !impedir x negativo
				if(x_l(i,nphas)<0)then
					!***usar alguma flag de acompanhamento aqui				
					print*, 'x_nu ', i, ' nphas', ' <0', ' @', __FILE__, __LINE__
					x_l(i,nphas)=0
				endif
			enddo
				do i = 1, ncomp
					x_l(i,j_) = K_(i,j_)*x_l(i,nphas)*dexp(th_(j_))
				enddo !i = 1, ncomp


!		print*, x_l_; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; 
!		print*, x_l_; if (ACOMPANHAMENTOSCREEN == 1) then; print*, 'paused at ', __FILE__, __LINE__; read(*,*) paused; endif
		
		!aqui eu atualizei o x usando todos os theta,
		!mas se o alfa da fase sendo atualizada já é zero o theta dela própria não faz efeito na conta.
		!então é coerente
		!
		
		!estrofe para versão da function com xizes desatualizados, como eu fazia antes
!		x_l_(:,:) = x(:,:)
		!continua
		
		!equação 25
		
		do i = 1, ncomp
!		print*, K_(i,j_), x_l_(i,nphas), __FILE__, __LINE__; pause !OFF
			term1(i) = K_(i,j_)*x_l(i,nphas); !print*, term1(i), 'resettheta term1'
		enddo !i = 1, ncomp
		sum1=0
		do i = 1, ncomp
			sum1 = sum1 + term1(i)
		enddo !i = 1, ncomp
!		print*, sum1; pause !OFF
		re_set_theta = dlog(1.d0/sum1)
!		
!		theta_k_25 = dlog(1.d0/sum1)
!		!equação de definição, por componente
!		do i = 1, ncomp
!			theta_k_def(i) = dlog(x_l_(i,j_)/term1(i))
!		enddo !i = 1, ncomp
!		print*, 're_set_theta', theta_k_25, theta_k_def
!		write(95,*) j_ !re_set_theta_debug
!		write(95,*) 're_set_theta', theta_k_25, theta_k_def !re_set_theta_debug
!		re_set_theta = theta_k_25
!		!resultado como mínimo valor positivo dentre os resultados das equações alternativas
!		do i = 1, ncomp
!			if( theta_k_def(i) > GTzero_TOL ) then !ao obter um valor positivo
!				if( re_set_theta < GTzero_TOL .or. re_set_theta > theta_k_def(i)) then !caso o anterior seja negativo (muito pequeno), ou maior que o novo (muito grande)
!					re_set_theta = theta_k_def(i) !aceito esse valor
!				else
!				!re_set_theta preservado como theta_k_25
!				endif
!			else
!			!re_set_theta preservado como theta_k_25
!			endif
!		enddo !i = 1, ncomp
!		print*, 're_set_theta', re_set_theta

	end function re_set_theta
	
	subroutine switch_reference(x_,xol_,K_,th_,al_,al_flick_,phindex_,ref_index)
	!apenas alteramos a ordem dos valores nos arrays dos problemas,
	!pois nesse programa todas as fases são calculadas analogamente (mesmo modelo com argumentos diferentes)
		real(8), intent(inout) :: x_(ncomp,nphas),xol_(ncomp,nphas),K_(ncomp,nphas),th_(nphas),al_(nphas)
		integer, intent(inout) :: al_flick_(nphas)
		integer, intent(inout) :: phindex_(nphas)
		!
		integer, intent(in) :: ref_index
		!
		real(8) :: x_n(ncomp,nphas),xol_n(ncomp,nphas),K_n(ncomp,nphas),th_n(nphas),al_n(nphas),al_flick_n(nphas)
		integer :: phindex_n(nphas)
!		integer :: cond_(nphas), cond_n(nphas)
		integer :: j, jj

		x_n(:,nphas) = x_(:,ref_index)
		xol_n(:,nphas) = xol_(:,ref_index)
		K_n(:,nphas) = K_(:,ref_index)
		th_n(nphas) = th_(ref_index)
		al_n(nphas) = al_(ref_index)
		al_flick_n(nphas) = al_flick_(ref_index)
!		cond_n(nphas) = cond_(ref_index)
			!ao invés de reordenar o cond, reordenar a lista de fases com a cond e a metodologia adequada para resolução
!		phase_list_n(nphas) = phase_list(ph_index(ref_index))
!			ao invés de reordenar a lista de fases renomear o ph_index
		phindex_n(nphas) = phindex_(ref_index)
!														SOON
!
		jj=1
		do j = 1, nphas
			if( j.NE.ref_index ) then
				x_n(:,jj) = x_(:,j)
				xol_n(:,jj) = xol_(:,j)
				K_n(:,jj) = K_(:,j)
				th_n(jj) = th_(j)
				al_n(jj) = al_(j)
				al_flick_n(jj) = al_flick_(j)
!				cond_n(jj) = cond_(j)
!				phase_list_n(jj) = phase_list(ph_index(j))
				phindex_n(jj) = phindex_(j)
				jj=jj+1
			endif
		enddo !j = 1, nphas
		x_(:,:) = x_n(:,:)
		xol_(:,:) = xol_n(:,:)
		K_(:,:) = K_n(:,:)
		th_(:) = th_n(:)
		al_(:) = al_n(:)
		al_flick_(:) = al_flick_n(:)
		
!		cond_(:) = cond_n(:)
!		phase_list(:) = phase_list_n(:)
		phindex_(:) = phindex_n(:)
		
		!obs: se eu peguei uma fase não presente pra freferencia:
!		al_(nphas) = too_small
		th_(nphas) = 0.d0
		al_flick_(nphas) = 1
		!:***será que vale a pena fazer a verificação da regra das fases aqui?
			!eu preciso que essa fase seja presente, então a verificação da regra das fases seria feito um pouco antes para proibir uma fase não presente de virar referencia caso já existissem fases demais presentes
		
	end subroutine switch_reference
	
	subroutine collapse_phase(x_,xol_,K_,th_,al_,al_flick_,phindex_,j_) !add j to jj, then remove j.
		!só deve ser chamada pela sub multiflash
		real(8), allocatable, intent(inout) :: x_(:,:),xol_(:,:),K_(:,:),th_(:),al_(:) !allocatable pra eu poder redimensionar
		integer, allocatable, intent(inout) :: phindex_(:)
		integer, allocatable, intent(inout) :: al_flick_(:)

		!colapsa as variáveis do module, não recebe cada uma como argumento
		!recebe apenas o indice que identifica as fases
		
		integer, intent(in) :: j_ !índices

		real(8) :: x_n(ncomp,nphas), xol_n(ncomp,nphas), K_n(ncomp,nphas), th_n(nphas), al_n(nphas), al_flick_n(nphas) !variáveis dummy com valores das antigas no tamanho inicial
		integer :: phindex_n(nphas)
		integer :: i, j, jj !contador

		!backup dos valores
		x_n(:,:) = x_(:,:)
		xol_n(:,:) = xol_(:,:)
		K_n(:,:) = K_(:,:) 
		th_n(:) = th_(:)
		al_n(:) = al_(:)
		al_flick_n(:) = al_flick_(:)
		phindex_n(:) = phindex_(:)
		
		!Resize variáveis do módulo (que entraram via argumento)
		deallocate(x_,xol_,K_,th_,al_,al_flick_,phindex_)
		nphas = nphas - 1
		write(nphas_str,'(I2.1)') nphas
		allocate(x_(ncomp,nphas),xol_(ncomp,nphas),K_(ncomp,nphas),th_(nphas),al_(nphas),al_flick_(nphas),phindex_(nphas))
		
		!recuperação dos valores das fases pertinentes
		!loop em j//jj para j diferente de j_
		jj=1
		do j = 1, nphas + 1 !o loop tem que ser do tamanho das variaveis antigas
			if( j.NE.j_ ) then
				x_(:,jj) = x_n(:,j)
				xol_(:,jj) = xol_n(:,j)
				K_(:,jj) = K_n(:,j)
				th_(jj) = th_n(j)
				al_(jj) = al_n(j)
				al_flick_(jj) = al_flick_n(j)
				phindex_(jj) = phindex_n(j)
				jj=jj+1
			endif
		enddo !j = 1, nphas
!!		!
!		!CONFERINDO
!!		print*, 'conferindo'
!!do i = 1, ncomp
!!		print*, i
!!		print*, 		x_n(i,:) 
!!		print*, 		x_(i,:)
!!enddo
!!		
!!		!print*, 		xol_n(:,:) 
!!		!print*, 		 xol_(:,:)
!!		!print*, 		K_n(:,:) 
!!		!print*, 		 K_(:,:) 
!!		!print*, 		th_n(:) 
!!		!print*, 		 th_(:)
!!		print*, 'conferindo'
!!		print*, 		al_n(:) 
!!		print*, 'conferindo'
!!		print*, 		 al_(:)
!!		
!!		print*, 		phindex_n(:) 
!!		print*, 		 phindex_(:)
!!		
	end subroutine collapse_phase

	subroutine switch_phases(index_1st,index_2nd,index_ref) !tem que ser baseada em nphas0
		!rodando fora da sub multiflash só precisa se importar com as var do módulo, nao as var da sub
		!apenas alteramos a ordem dos valores nos arrays dos problemas,
	!	!pois nesse programa todas as fases são calculadas analogamente (mesmo modelo com argumentos diferentes)
		!DECLARATIONS
		!ARGUMENTS
		integer :: index_1st
		integer :: index_2nd
		integer :: index_ref
		integer :: j
		integer :: jj
		!LOCAL
		real(8), allocatable :: phase_frac_ol(:)
		real(8), allocatable :: phase_flick_ol(:)
		real(8), allocatable :: stab_var_ol(:)
		real(8), allocatable :: dist_coef_ol(:,:)
		integer, allocatable :: ph_index_ol(:) !relaciona a ordem atual dentre nphas com a ordem original dentre nphas0

		if (index_1st > nphas .or. index_2nd > nphas .or. index_ref > nphas) then !***novo no multiflash3, não tem no 1 e 2
			print*, 'invalid switch_phases index, greater than nphas'
			return
		endif

		!ALLOCATE
		allocate(phase_frac_ol(nphas))
		allocate(phase_flick_ol(nphas))
		allocate(stab_var_ol(nphas))
		allocate(dist_coef_ol(ncomp,nphas))
		allocate(ph_index_ol(nphas))

		!BACKING-UP
		phase_frac_ol(:) = phase_frac(:)
		phase_flick_ol(:) = phase_flick(:)
		stab_var_ol(:) = stab_var(:)
		dist_coef_ol(:,:) = dist_coef(:,:)
		ph_index_ol(:) = ph_index(:)

!!***
!!DEBUG
!print*, phase_frac_ol(:), __FILE__, __LINE__
!print*, stab_var_ol(:), __FILE__, __LINE__
!print*, dist_coef_ol(:,:), __FILE__, __LINE__
!print*, ph_index_ol(:), __FILE__, __LINE__
!!DEBUG END
!!***
		!RE-FILLING SPECIAL SLOTS
		phase_frac(1) = phase_frac_ol(index_1st)
		phase_flick(1) = phase_flick_ol(index_1st)
		stab_var(1) = stab_var_ol(index_1st)
		dist_coef(:,1) = dist_coef_ol(:,index_1st)
		ph_index(1) = ph_index_ol(index_1st)
		
		phase_frac(2) = phase_frac_ol(index_2nd)
		phase_flick(2) = phase_flick_ol(index_2nd)
		stab_var(2) = stab_var_ol(index_2nd)
		dist_coef(:,2) = dist_coef_ol(:,index_2nd)
		ph_index(2) = ph_index_ol(index_2nd)
		
		phase_frac(nphas) = phase_frac_ol(index_ref)
		phase_flick(nphas) = phase_flick_ol(index_ref)
		stab_var(nphas) = stab_var_ol(index_ref)
		dist_coef(:,nphas) = dist_coef_ol(:,index_ref)
		ph_index(nphas) = ph_index_ol(index_ref)
		
		!RE_FILLING REMAINING
		jj=3
		do j = 1, nphas
			if( j .NE. index_1st .and. j .NE. index_2nd .and. j .NE. index_ref ) then
				phase_frac(jj) = phase_frac_ol(j)
				phase_flick(jj) = phase_flick_ol(j)
				stab_var(jj) = stab_var_ol(j)
				dist_coef(:,jj) = dist_coef_ol(:,j)
				ph_index(jj) = ph_index_ol(j)
				jj=jj+1 !j = 1 até nphas dif 1st 2nd ref, logo nphas-3 elementnos, logo jj vai de 3 até (3) + (nphas-3) (-1) = nphas-1
			endif
		enddo !j = 1, nph as
		
!!***
!!DEBUG
!print*, phase_frac_ol(:), __FILE__, __LINE__
!print*, stab_var_ol(:), __FILE__, __LINE__
!print*, dist_coef_ol(:,:), __FILE__, __LINE__
!print*, ph_index_ol(:), __FILE__, __LINE__
!!DEBUG END
!!***
	end subroutine switch_phases

	subroutine ressurrect_phases
		!essa rotina reintegra ao conjunto de fases do flash as fases que foram eliminadas pelo collapse_phase
		!obs: se a convergência foi para uma fase tipo vapor (vide volume molar), mas ocorreu na fase com phase_cond 'L',
			!então essa rotina vai reintegrar ao código a fase com phase_cond 'V' e reinicializar ela, que vai obviamente colapsar de novo
			!mas o interessante seria passar os valores dessa fase tipo vapor para a fase 'V' e reinicializar a fase 'L'
		integer :: i, j, jj
		integer :: Ress_fnum
		real(8), allocatable :: phase_frac_ol(:)
		real(8), allocatable :: phase_flick_ol(:)
		real(8), allocatable :: stab_var_ol(:)
		real(8), allocatable :: dist_coef_ol(:,:)
		integer, allocatable :: ph_index_ol(:)
		!#1 select phases
	!	print*, ph_index_inv(:)
		do j = 1, nphas0
			allocate(phase_frac_ol(nphas))
			allocate(phase_flick_ol(nphas))
			allocate(stab_var_ol(nphas))
			allocate(dist_coef_ol(ncomp,nphas))
			allocate(ph_index_ol(nphas))
			if( ph_index_inv(j) == 0) then
				!backup
				phase_frac_ol(:) = phase_frac(:)
				phase_flick_ol(:) = phase_flick(:)
				stab_var_ol(:) = stab_var(:)
				dist_coef_ol(:,:) = dist_coef(:,:)
				ph_index_ol(:) = ph_index(:)
				!resize
				deallocate(phase_frac,phase_flick,stab_var,dist_coef,ph_index)
				nphas = nphas + 1
				write(nphas_str,'(I2.1)') nphas
				allocate(phase_frac(nphas))
				allocate(phase_flick(nphas))
				allocate(stab_var(nphas))
				allocate(dist_coef(ncomp,nphas))
				allocate(ph_index(nphas))
				!restore
				!nphas
				phase_frac(nphas) = phase_frac_ol(nphas-1)
				phase_flick(nphas) = phase_flick_ol(nphas-1)
				stab_var(nphas) = stab_var_ol(nphas-1)
				dist_coef(:,nphas) = dist_coef_ol(:,nphas-1)
				ph_index(nphas) = ph_index_ol(nphas-1)
				!1,nphas-2
				do jj = 1, nphas-2
					phase_frac(jj) = phase_frac_ol(jj)
					phase_flick(jj) = phase_flick_ol(jj)
					stab_var(jj) = stab_var_ol(jj)
					dist_coef(:,jj) = dist_coef_ol(:,jj)
					ph_index(jj) = ph_index_ol(jj)
				enddo !j = 1, nphas-2
				!re-initialize
				!nphas-1
				phase_frac(nphas-1) = 0.d0 !1.d-9
				phase_flick(nphas-1) = 0
				stab_var(nphas-1) = 0.1d0 !***ou toosmall ou re_set_theta
				
	!			dist_coef(:,nphas-1) = dist_coef_ol(:,nphas) !não precisa dar valor aqui
				ph_index(nphas-1) = j
!				print*, 'ressurrected phase has default composition', __FILE__, __LINE__ !***
!				phase_list(j)%phase%x(:) = (/.1d0,.1d0,.8d0/) !só falta botar esse valor pra vir do mesmo arquivo que inicializou a fase ou de outro
				!pronto
				print*, 'ressurrected phase initialized according to file ', 'input/'//trim(sim_case)//'/phases/'//'Ress.dat', __FILE__, __LINE__ !***
				open(newunit=Ress_fnum, file='input/'//trim(sim_case)//'/phases/'//'Ress.dat', status='old', action='read')
				read(Ress_fnum,*); read(Ress_fnum,*) (phase_list(ph_index(nphas-1))%phase%x(i),i=1,ncomp); !EM LINHA
				read(Ress_fnum,*); read(Ress_fnum,*) phase_list(ph_index(nphas-1))%phase%phase_cond
				close(Ress_fnum)
			endif
			deallocate(phase_frac_ol)
			deallocate(phase_flick_ol)
			deallocate(stab_var_ol)
			deallocate(dist_coef_ol)
			deallocate(ph_index_ol)
		enddo !i = 1, nphas0
	
		!temporário composition
	!	phase_list(nphas-1)%phase%x(:) = (/.8807d0,.111d0,.008d0/) !só falta botar esse valor pra vir do mesmo arquivo que inicializou a fase ou de outro
	!	phase_list(nphas-2)%phase%x(:) = (/.10d0,.09d0,.910d0/) !só falta botar esse valor pra vir do mesmo arquivo que inicializou a fase ou de outro
	!	
	end subroutine ressurrect_phases

!!!!!ENDSUBROUTINES DO FLASH
!
end module multiflash_mod
!eos agora não pode ter T nem P, tem que ter seu x(:), apenas receber T e P para os procedimentos como argumento)

