dg.exe: main.f90 module_gau_elimi.o  module_gauss_tri.o module_regular_mesh.o module_globals.o module_prob.o module_polynomials_tri.o module_output_tri.o module_regular_mesh_data.o 
	ifort main.f90 module_gau_elimi.o  module_gauss_tri.o module_regular_mesh.o module_globals.o module_prob.o module_polynomials_tri.o module_output_tri.o module_regular_mesh_data.o -r8 -O3 -o dg.exe
module_regular_mesh.mod: module_regular_mesh.o module_regular_mesh.f90  
	ifort module_regular_mesh.f90 -c -r8
module_regular_mesh.o: module_regular_mesh.f90 module_prob.mod module_regular_mesh_data.mod module_gauss_tri.mod module_gau_elimi.mod 
	ifort module_regular_mesh.f90 -c -r8 
module_output_tri.o: module_output_tri.f90 module_regular_mesh_data.mod 
	ifort -c -r8 module_output_tri.f90 
module_gauss_tri.mod: module_gauss_tri.o module_gauss_tri.f90
	ifort -c -r8 module_gauss_tri.f90
module_gauss_tri.o: module_gauss_tri.f90 module_prob.mod
	ifort -c -r8 module_gauss_tri.f90 
module_prob.mod: module_prob.o module_prob.f90
	ifort -c -r8 module_prob.f90
module_prob.o: module_prob.f90
	ifort -c -r8 module_prob.f90
module_globals.mod: module_globals.o module_globals.f90
	ifort -c -r8 module_globals.f90
module_globals.o: module_globals.f90
	ifort -c -r8 module_globals.f90
module_regular_mesh_data.mod: module_regular_mesh_data.o module_regular_mesh_data.f90
	ifort -c -r8 module_regular_mesh_data.f90
module_regular_mesh_data.o: module_regular_mesh_data.f90
	ifort -c -r8 module_regular_mesh_data.f90
module_polynomials_tri.mod: module_polynomials_tri.f90 module_polynomials_tri.o
	ifort -c -r8 module_polynomials_tri.f90
module_polynomials_tri.o: module_polynomials_tri.f90
	ifort -c -r8 module_polynomials_tri.f90
module_gau_elimi.mod: module_gau_elimi.f90 module_gau_elimi.o
	ifort -c -r8 module_gau_elimi.f90
module_gau_elimi.o: module_gau_elimi.f90
	ifort -c -r8 module_gau_elimi.f90
clean:
	rm *.o *.mod
