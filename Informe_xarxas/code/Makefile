MKDIR_P := mkdir -p
OUT_DIR := output

.PHONY: directories all clean

all: ${OUT_DIR} dynamic_process subroutine.so ass4.so
subroutine.so : subroutine.f90 
	f2py3 -m -c network subroutine.f90

ass4.so : assignment4.f90 
	f2py3 -m -c network subroutine.f90

dynamic_process : module_dp.o dynamic_process.o
	gfortran -o dynamic_process.x dynamic_process.o module_dp.o 



module_dp.o : module_dp.f90 
	gfortran -c module_dp.f90

dynamic_process.o: module_dp.f90 dynamic_process.f90
	gfortran -c dynamic_process.f90


directories: $(OUT_DIR)
${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}
clean:
	rm -f *.mod
	rm -f *.o 
	rm -f *.x
	rm -f *.so

run_dynam:
	./dynamic_process.x

run_ass1:
	python3 assignment1.py

run_ass4:
	python3 assignment4.py

run_ass5:
	python3 assignment5.py
	
run_plots_dyn:
	python3 plots_dynamics.py