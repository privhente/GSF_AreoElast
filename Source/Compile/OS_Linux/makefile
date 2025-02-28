# Intel Fortran Compiler Makefile for DeSiO
# Christian Hente - 24.02.2025
#
# initiate variables
CONTROLLER_ADDON ?= 0

# Paths
path_desio = /home/christian/DeSiO/Source
path_shared = $(path_desio)/shared_routines
path_aero = $(path_desio)/Aero
path_structure = $(path_desio)/Structure
path_controller = $(path_desio)/Controller
path_fsi = $(path_desio)/FSI
path_mod = $(path_desio)/mod

# Compiler settings
FC = ifx
FFLAGS = -O2 -fp-model strict -fpp -fiopenmp -qmkl=parallel -heap-arrays -module $(path_mod)
ifeq ($(CONTROLLER_ADDON),1)
  FFLAGS += -DController_ADDON
endif
LIBS = -qmkl=parallel

# Debug-Flag: Setze BUILD=DEBUG, um Debug-Einstellungen zu aktivieren.
ifeq ($(BUILD),DEBUG)
  FFLAGS = -O0 -g -traceback -check all -fpp -fiopenmp -qmkl=parallel -heap-arrays -module $(path_mod)
endif

# Automatische Erfassung von .f90-Dateien
SRC_SHARED = $(wildcard $(path_shared)/*.f90)
SRC_AERO = $(wildcard $(path_aero)/*.f90)
SRC_AERO_NO_MAIN = $(filter-out $(path_aero)/main_program.f90, $(SRC_AERO))
SRC_STRUC = $(wildcard $(path_structure)/*.f90)
SRC_STRUC_NO_MAIN = $(filter-out $(path_structure)/main_program.f90, $(SRC_STRUC))
SRC_CONTROLLER = $(wildcard $(path_controller)/class_controller.f90)
SRC_FSI = $(wildcard $(path_fsi)/*.f90)
SRC_FSI_NO_MAIN = $(filter-out $(path_fsi)/main_program.f90, $(SRC_FSI))

# Objektdateien definieren
OBJ_SHARED = $(patsubst $(path_shared)/%.f90, $(path_shared)/obj/%.o, $(SRC_SHARED))
OBJ_AERO = $(patsubst $(path_aero)/%.f90, $(path_aero)/obj/%.o, $(SRC_AERO_NO_MAIN))
OBJ_AERO_NO_MAIN = $(patsubst $(path_aero)/%.f90, $(path_aero)/obj/%.o, $(SRC_AERO_NO_MAIN))
OBJ_STRUC = $(patsubst $(path_structure)/%.f90, $(path_structure)/obj/%.o, $(SRC_STRUC_NO_MAIN))
OBJ_STRUC_NO_MAIN = $(patsubst $(path_structure)/%.f90, $(path_structure)/obj/%.o, $(SRC_STRUC_NO_MAIN))
OBJ_CONTROLLER = $(patsubst $(path_controller)/%.f90, $(path_controller)/obj/%.o, $(SRC_CONTROLLER))
OBJ_FSI = $(patsubst $(path_fsi)/%.f90, $(path_fsi)/obj/%.o, $(SRC_FSI_NO_MAIN))

# main_program.o als letztes hinzufügen
OBJ_AERO += $(path_aero)/obj/main_program.o
TARGET_AERO = DeSiO_Aero

OBJ_STRUC += $(path_structure)/obj/main_program.o
TARGET_STRUC = DeSiO_Structure

OBJ_FSI += $(path_fsi)/obj/main_program.o
TARGET_FSI = DeSiO

# Default target
all: dirs $(TARGET_FSI) $(TARGET_STRUC) $(TARGET_AERO)

# Kompilieren von Objektdateien mit automatischer Dependency-Erzeugung
$(path_shared)/obj/%.o: $(path_shared)/%.f90
	$(FC) $(FFLAGS) $(LIBS) -c -o $@ $< -gen-dep=make

$(path_structure)/obj/%.o: $(path_structure)/%.f90
	$(FC) $(FFLAGS) $(LIBS) -c -o $@ $< -gen-dep=make

$(path_aero)/obj/%.o: $(path_aero)/%.f90
	$(FC) $(FFLAGS) $(LIBS) -c -o $@ $< -gen-dep=make

$(path_controller)/obj/%.o: $(path_controller)/%.f90
	$(FC) $(FFLAGS) $(LIBS) -c -o $@ $< -gen-dep=make

$(path_fsi)/obj/%.o: $(path_fsi)/%.f90
	$(FC) $(FFLAGS) $(LIBS) -c -o $@ $< -gen-dep=make

# Linken der Executables
$(TARGET_FSI): $(OBJ_SHARED) $(OBJ_STRUC_NO_MAIN) $(OBJ_AERO_NO_MAIN) $(OBJ_CONTROLLER) $(OBJ_FSI)
	$(FC) $(FFLAGS) $(LIBS) $^ -o $@

$(TARGET_STRUC): $(OBJ_SHARED) $(OBJ_STRUC)
	$(FC) $(FFLAGS) $(LIBS) $^ -o $@

$(TARGET_AERO): $(OBJ_SHARED) $(OBJ_AERO)
	$(FC) $(FFLAGS) $(LIBS) $^ -o $@

# Bereinigen
clean:
	rm -f $(path_shared)/obj/*.o $(path_shared)/obj/*.d \
          $(path_structure)/obj/*.o $(path_structure)/obj/*.d \
          $(path_aero)/obj/*.o $(path_aero)/obj/*.d \
          $(path_controller)/obj/*.o $(path_controller)/obj/*.d \
          $(path_fsi)/obj/*.o $(path_fsi)/obj/*.d \ $(path_mod)/*.mod \
          $(TARGET_STRUC) $(TARGET_AERO) $(TARGET_FSI)

# Erstkompilierung mit Skript
first_build:
	bash compile_ifx.sh

.PHONY: dirs
dirs:
	mkdir -p $(path_shared)/obj $(path_aero)/obj $(path_structure)/obj $(path_controller)/obj $(path_fsi)/obj

# Einbinden der automatisch generierten Dependency-Dateien
-include $(wildcard $(path_shared)/obj/*.d) \
         $(wildcard $(path_aero)/obj/*.d) \
         $(wildcard $(path_structure)/obj/*.d) \
         $(wildcard $(path_controller)/obj/*.d) \
         $(wildcard $(path_fsi)/obj/*.d)
