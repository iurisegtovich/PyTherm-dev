program editcase
	implicit none
	character(100) :: sim_case
	open(unit=10,file='input/select_sim_case',status='old',action='read')
		read(10,*) sim_case !nome da pasta com intruções para simulação
	close(10)
	
	call system('gedit input/'//trim(sim_case)//'/*.dat')
	call system('gedit input/'//trim(sim_case)//'/*/*.dat')
end program editcase
