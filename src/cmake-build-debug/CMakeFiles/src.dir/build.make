# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /snap/clion/58/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/58/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/src.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/src.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/src.dir/flags.make

CMakeFiles/src.dir/kmer_categorizer.cpp.o: CMakeFiles/src.dir/flags.make
CMakeFiles/src.dir/kmer_categorizer.cpp.o: ../kmer_categorizer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/src.dir/kmer_categorizer.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/src.dir/kmer_categorizer.cpp.o -c "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/kmer_categorizer.cpp"

CMakeFiles/src.dir/kmer_categorizer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/src.dir/kmer_categorizer.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/kmer_categorizer.cpp" > CMakeFiles/src.dir/kmer_categorizer.cpp.i

CMakeFiles/src.dir/kmer_categorizer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/src.dir/kmer_categorizer.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/kmer_categorizer.cpp" -o CMakeFiles/src.dir/kmer_categorizer.cpp.s

# Object files for target src
src_OBJECTS = \
"CMakeFiles/src.dir/kmer_categorizer.cpp.o"

# External object files for target src
src_EXTERNAL_OBJECTS =

src: CMakeFiles/src.dir/kmer_categorizer.cpp.o
src: CMakeFiles/src.dir/build.make
src: /usr/lib/x86_64-linux-gnu/libboost_program_options.so
src: CMakeFiles/src.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable src"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/src.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/src.dir/build: src

.PHONY : CMakeFiles/src.dir/build

CMakeFiles/src.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/src.dir/cmake_clean.cmake
.PHONY : CMakeFiles/src.dir/clean

CMakeFiles/src.dir/depend:
	cd "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src" "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src" "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/cmake-build-debug" "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/cmake-build-debug" "/home/whiskas/Documents/MEGAsync/FMFI/2. ročník Mgr/Diplomovka/Hybrid-Genome-Assembler/src/cmake-build-debug/CMakeFiles/src.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/src.dir/depend

