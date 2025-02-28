#!/bin/bash
# Batch script to compile DeSiO source code with Intel ifx compiler"
# Author: Christian Hente,
# Date: 08.02.2025
clear
#S
# Import Intel Fortran and MKL
# source ./opt/intel/oneapi/setvars.sh
#
# Set DeSiO source path
path_desio="/home/christian/DeSiO/Source"
path_mod=$path_desio/mod
#
path_structure=$path_desio/Structure
path_structure_obj=$path_structure/obj
#
path_aero=$path_desio/Aero
path_aero_obj=$path_aero/obj
#
path_fsi=$path_desio/FSI
path_fsi_obj=$path_fsi/obj
#
path_shared=$path_desio/shared_routines
path_shared_obj=$path_shared/obj
#
path_controller=$path_desio/Controller_routines
path_controller_obj=$path_controller/obj
#
rm $path_structure_obj $path_aero_obj $path_fsi_obj $path_shared_obj $path_controller_obj $path_mod -rf
mkdir $path_mod
mkdir $path_structure_obj
mkdir $path_aero_obj
mkdir $path_fsi_obj
mkdir $path_shared_obj
mkdir $path_controller_obj
#
# Compiler settings:
FC=ifx
FFLAGS="-O2 -fp-model strict -fpp -fiopenmp -qmkl=parallel -heap-arrays -module $path_mod" #-DController_ADDON
#
LIBS="-qmkl=parallel"
#
# compile shared routines
$FC $FFLAGS $LIBS -c -o  "$path_shared_obj/mkl_dss.o"                          $path_shared/mkl_dss.f90
$FC $FFLAGS $LIBS -c -o  "$path_shared_obj/sparse_function.o"                  $path_shared/sparse_function.f90
$FC $FFLAGS $LIBS -c -o  "$path_shared_obj/SingPrec.o"                         $path_shared/SingPrec.f90
$FC $FFLAGS $LIBS -c -o  "$path_shared_obj/my_FileIO_functions.o"              $path_shared/my_FileIO_functions.f90
#
# DeSiO-Structure                                
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_load_6.o"                                  $path_structure/class_load_6.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_sparse_matrix.o"                           $path_structure/class_sparse_matrix.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_load_12.o"                                 $path_structure/class_load_12.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/my_constants_structure.o"                        $path_structure/my_constants_structure.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_nullspace.o"                               $path_structure/class_nullspace.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_point_mass.o"                             $path_structure/class_point_mass.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/my_math_structure.o"                              $path_structure/my_math_structure.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_boundary.o"                                $path_structure/class_boundary.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/my_materials.o"                                   $path_structure/my_materials.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_node_structure.o"                           $path_structure/class_node_structure.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_beam_element_dissipation.o"                $path_structure/class_beam_element_dissipation.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_constraint.o"                           $path_structure/class_constraint.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_beam_element.o"                          $path_structure/class_beam_element.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_matrix12.o"                                $path_structure/class_matrix12.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/solver_variables.o"                             $path_structure/solver_variables.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_node_6.o"                                 $path_structure/class_node_6.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_amplitude.o"                             $path_structure/class_amplitude.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_constraint_6_to_12.o"                    $path_structure/class_constraint_6_to_12.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_node_12.o"                                  $path_structure/class_node_12.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_condensed12.o"                            $path_structure/class_condensed12.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_spring12.o"                                $path_structure/class_spring12.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_shell_element_dissipation.o"               $path_structure/class_shell_element_dissipation.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_constraint_6.o"                            $path_structure/class_constraint_6.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_single_layer_material.o"                  $path_structure/class_single_layer_material.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_constraint_12.o"                           $path_structure/class_constraint_12.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_rigid_body.o"                              $path_structure/class_rigid_body.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_beam.o"                                    $path_structure/class_beam.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_damping12.o"                               $path_structure/class_damping12.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_multi_layer_material.o"                    $path_structure/class_multi_layer_material.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_shell_element.o"                           $path_structure/class_shell_element.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_shell.o"                                   $path_structure/class_shell.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/class_model_structure.o"                         $path_structure/class_model_structure.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/solver_functions.o"                              $path_structure/solver_functions.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/invariants.o"                                    $path_structure/invariants.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/main_structure.o"                                $path_structure/main_structure.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/static_solver_arc.o"                             $path_structure/static_solver_arc.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/modal_solver.o"                                  $path_structure/modal_solver.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/buckling_solver.o"                               $path_structure/buckling_solver.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/static_solver.o"                                 $path_structure/static_solver.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/dynamic_solver.o"                                $path_structure/dynamic_solver.f90
$FC $FFLAGS $LIBS -c -o "$path_structure_obj/main_program.o"                                   $path_structure/main_program.f90
#
# DeSiO-Aero                                    
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/my_math_aero.o"                                 $path_aero/my_math_aero.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/my_types_aero.o"                                $path_aero/my_types_aero.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/my_constants_aero.o"                            $path_aero/my_constants_aero.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_inflow_aero.o"                            $path_aero/class_inflow_aero.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_section.o"                                $path_aero/class_section.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_node_aero.o"                              $path_aero/class_node_aero.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_node3_aero.o"                             $path_aero/class_node3_aero.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/my_vlm_functions.o"                             $path_aero/my_vlm_functions.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_vortex_segment.o"                         $path_aero/class_vortex_segment.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_vortex_ring.o"                            $path_aero/class_vortex_ring.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_vortex_sheet.o"                           $path_aero/class_vortex_sheet.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_bounded_vortex_sheet.o"                   $path_aero/class_bounded_vortex_sheet.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_unbounded_vortex_sheet.o"                 $path_aero/class_unbounded_vortex_sheet.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/class_model_aero.o"                             $path_aero/class_model_aero.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/num_diff_q.o"                                   $path_aero/num_diff_q.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/main_aero.o"                                    $path_aero/main_aero.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/num_diff_v.o"                                   $path_aero/num_diff_v.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/num_diff_DG.o"                                  $path_aero/num_diff_DG.f90
$FC $FFLAGS $LIBS -c -o "$path_aero_obj/main_program.o"                                 $path_aero/main_program.f90
# 
# compile DeSiO-Controller_routines
#$FC $FFLAGS $LIBS -c -o IWinAPI_m.o                         $path_controller/IWinAPI_m.f90
#$FC $FFLAGS $LIBS -c -o DLLHelper_m.o                       $path_controller/DLLHelper_m.f90
#$FC $FFLAGS $LIBS -c -o DISCON_m.o                          $path_controller/DISCON_m.f90
#$FC $FFLAGS $LIBS -c -o DISCON_m_linux.o                    $path_controller/DISCON_m_linux.f90
$FC $FFLAGS $LIBS -c -o "$path_controller_obj/class_controller.o"                   $path_controller/class_controller.f90
#
# DeSiO-FSI                                      
$FC $FFLAGS $LIBS -c -o "$path_fsi_obj/my_fsi_functions.o"                            $path_fsi/my_fsi_functions.f90
$FC $FFLAGS $LIBS -c -o "$path_fsi_obj/class_model_fsi.o"                             $path_fsi/class_model_fsi.f90
$FC $FFLAGS $LIBS -c -o "$path_fsi_obj/fsi_dynamic_solver.o"                          $path_fsi/fsi_dynamic_solver.f90
$FC $FFLAGS $LIBS -c -o "$path_fsi_obj/fsi_static_solver.o"                           $path_fsi/fsi_static_solver.f90
$FC $FFLAGS $LIBS -c -o "$path_fsi_obj/fsi_kinematic_solver.o"                        $path_fsi/fsi_kinematic_solver.f90
$FC $FFLAGS $LIBS -c -o "$path_fsi_obj/main_fsi_nonlinear.o"                          $path_fsi/main_fsi_nonlinear.f90
$FC $FFLAGS $LIBS -c -o "$path_fsi_obj/main_program.o"                                $path_fsi/main_program.f90
#
# Merge object files into executable
OBJ_FILES_AERO=$(find "$path_shared_obj" "$path_aero_obj" -name "*.o")
OBJ_FILES_STRUC=$(find "$path_shared_obj" "$path_structure_obj" -name "*.o")
OBJ_FILES_FSI=$(find "$path_shared_obj" "$path_structure_obj" "$path_aero_obj" "$path_controller_obj" "$path_fsi_obj" -name "*.o" | grep -v -E "$path_structure_obj/main_program.o|$path_aero_obj/main_program.o")

$FC $FFLAGS $LIBS $OBJ_FILES_AERO -o DeSiO_Aero
$FC $FFLAGS $LIBS $OBJ_FILES_STRUC -o DeSiO_Structure
$FC $FFLAGS $LIBS $OBJ_FILES_FSI -o DeSiO
#
