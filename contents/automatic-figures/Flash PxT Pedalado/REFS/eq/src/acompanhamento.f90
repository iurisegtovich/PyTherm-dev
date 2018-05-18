module acompanhamento_mod
	implicit none
	character(100) :: paused
	integer :: acompanhamentoscreen, acompanhamentolog, acompanhamentolog2 !1 enabe pauses 0 just run !1 log em arquivo, 0 just run
	real(8) :: time_1, time_2
	integer :: event_flag(5)
	real(8) :: met_num_par_val(3)
end module acompanhamento_mod
