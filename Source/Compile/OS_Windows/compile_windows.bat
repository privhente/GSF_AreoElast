:: Batch script to compile DeSiO source code with intel fortran compiler in windows 
:: Author: Christian Hente, ISD
:: Date: 26.10.2022

:: remove old files from directory
del *.obj 
del *.mod 
del *.exe 
del *.optrpt

:: set ifort, mkl enviroment variables
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" --config config.txt
set INCLUDE=%INCLUDE%;%cd%

:: set path directory of DeSiO-Source
set path_desio=C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Source

:: set path for DeSiO solvers
set path_structure="%path_desio%\Structure"
set path_aero="%path_desio%\Aero"
set path_fsi="%path_desio%\FSI"
set path_shared="%path_desio%\shared_routines"
set path_controller="%path_desio%\Controller_routines"
set path_hydro="%path_desio%\Hydro"

:: compile common files
ifort /c    %path_shared%\mkl_dss.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_shared%\sparse_function.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_shared%\SingPrec.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_shared%\my_FileIO_functions.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16

:: compile DeSiO-Structure
ifort /c    %path_structure%\class_load_6.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_sparse_matrix.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_load_12.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\my_constants_structure.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_nullspace.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_point_mass.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\my_math_structure.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_boundary.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\my_materials.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_node_structure.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_beam_element_dissipation.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_constraint.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_beam_element.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_matrix12.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\solver_variables.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_node_6.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_amplitude.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_constraint_6_to_12.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_node_12.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_condensed12.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_spring12.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_shell_element_dissipation.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_constraint_6.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_single_layer_material.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_constraint_12.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_rigid_body.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_beam.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_damping12.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_multi_layer_material.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_shell_element.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_shell.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\class_model_structure.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\solver_functions.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\invariants.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\main_structure.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\static_solver_arc.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\modal_solver.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\buckling_solver.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\static_solver.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_structure%\dynamic_solver.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16

:: compile DeSiO-Aero
ifort /c    %path_aero%\my_math_aero.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\my_constants_aero.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_inflow_aero.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_section.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_node_aero.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_node3_aero.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\my_vlm_functions.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_vortex_segment.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_vortex_ring.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_vortex_sheet.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_bounded_vortex_sheet.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_unbounded_vortex_sheet.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\class_model_aero.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\num_diff_q.f90      /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\main_aero.f90       /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\num_diff_v.f90      /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_aero%\num_diff_DG.f90     /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16

:: compile DeSiO-Controller_routines
ifort /c    %path_controller%\IWinAPI_m.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_controller%\DLLHelper_m.f90 /O2 /fp:strict /Qopenmp /Qmkl:parallel /heap-arrays16
ifort /c    %path_controller%\DISCON_m.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_controller%\class_controller.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16

:: compile DeSIO-Hydro routines
ifort /c    %path_hydro%\my_constants_hydro.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_hydro%\my_math_hydro.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_hydro%\class_node_hydro.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_hydro%\class_node3_hydro.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_hydro%\class_morison_segment.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_hydro%\class_model_hydro.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_hydro%\dynamic_solver_hydro.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_hydro%\main_hydro.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
:: ifort /c    %path_hydro%\main_program.f90     /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16

:: compile DeSiO-FSI
ifort /c    %path_fsi%\my_fsi_functions.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_fsi%\class_model_fsi.f90  /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_fsi%\fsi_dynamic_solver.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_fsi%\fsi_static_solver.f90    /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_fsi%\fsi_kinematic_solver.f90 /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_fsi%\main_fsi_nonlinear.f90   /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16
ifort /c    %path_fsi%\main_program.f90     /O2 /fp:strict  /Qopenmp    /Qmkl:parallel  /heap-arrays16

:: merge object files to .exe file
ifort /exe:DeSiO.exe *.obj

pause
