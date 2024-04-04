# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/maipe/aspect/crustal-thermal

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/maipe/aspect/crustal-thermal

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/maipe/aspect/crustal-thermal/CMakeFiles /home/maipe/aspect/crustal-thermal//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/maipe/aspect/crustal-thermal/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named melt_petrol

# Build rule for target.
melt_petrol: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 melt_petrol
.PHONY : melt_petrol

# fast build rule for target.
melt_petrol/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/melt_petrol.dir/build.make CMakeFiles/melt_petrol.dir/build
.PHONY : melt_petrol/fast

#=============================================================================
# Target rules for targets named c_heating

# Build rule for target.
c_heating: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 c_heating
.PHONY : c_heating

# fast build rule for target.
c_heating/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/c_heating.dir/build.make CMakeFiles/c_heating.dir/build
.PHONY : c_heating/fast

c_heating.o: c_heating.cc.o
.PHONY : c_heating.o

# target to build an object file
c_heating.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/c_heating.dir/build.make CMakeFiles/c_heating.dir/c_heating.cc.o
.PHONY : c_heating.cc.o

c_heating.i: c_heating.cc.i
.PHONY : c_heating.i

# target to preprocess a source file
c_heating.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/c_heating.dir/build.make CMakeFiles/c_heating.dir/c_heating.cc.i
.PHONY : c_heating.cc.i

c_heating.s: c_heating.cc.s
.PHONY : c_heating.s

# target to generate assembly for a file
c_heating.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/c_heating.dir/build.make CMakeFiles/c_heating.dir/c_heating.cc.s
.PHONY : c_heating.cc.s

melt_petrol.o: melt_petrol.cc.o
.PHONY : melt_petrol.o

# target to build an object file
melt_petrol.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/melt_petrol.dir/build.make CMakeFiles/melt_petrol.dir/melt_petrol.cc.o
.PHONY : melt_petrol.cc.o

melt_petrol.i: melt_petrol.cc.i
.PHONY : melt_petrol.i

# target to preprocess a source file
melt_petrol.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/melt_petrol.dir/build.make CMakeFiles/melt_petrol.dir/melt_petrol.cc.i
.PHONY : melt_petrol.cc.i

melt_petrol.s: melt_petrol.cc.s
.PHONY : melt_petrol.s

# target to generate assembly for a file
melt_petrol.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/melt_petrol.dir/build.make CMakeFiles/melt_petrol.dir/melt_petrol.cc.s
.PHONY : melt_petrol.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... c_heating"
	@echo "... melt_petrol"
	@echo "... c_heating.o"
	@echo "... c_heating.i"
	@echo "... c_heating.s"
	@echo "... melt_petrol.o"
	@echo "... melt_petrol.i"
	@echo "... melt_petrol.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

