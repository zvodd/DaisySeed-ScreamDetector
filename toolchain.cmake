# This file is a CMake toolchain file for the ARM GCC toolchain.
# It should be placed in the project's root directory.

# 1. Set the system name to "Generic" for cross-compilation
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

# 2. Define the path to your ARM GCC toolchain directory
#    Replace this with the actual path to your toolchain.
if(DEFINED ENV{DAISY_TOOLCHAIN_DIR})
    # If it is, set the CMake variable to its value
    set(DAISY_TOOLCHAIN_DIR "$ENV{DAISY_TOOLCHAIN_DIR}")
    message(STATUS "DAISY_TOOLCHAIN_DIR is set to: ${DAISY_TOOLCHAIN_DIR}")
else()
    message(FATAL_ERROR "DAISY_TOOLCHAIN_DIR environment variable is not set. Please set it before running CMake.")
endif()

# 3. Specify the compilers and other tools
set(CMAKE_C_COMPILER_WORKS 1)
set(CMAKE_CXX_COMPILER_WORKS 1)

set(CMAKE_C_COMPILER   "${DAISY_TOOLCHAIN_DIR}/bin/arm-none-eabi-gcc.exe")
set(CMAKE_CXX_COMPILER "${DAISY_TOOLCHAIN_DIR}/bin/arm-none-eabi-g++.exe")
set(CMAKE_ASM_COMPILER "${DAISY_TOOLCHAIN_DIR}/bin/arm-none-eabi-as.exe")

# 4. Define the linker, archiver, and other utilities
set(CMAKE_OBJCOPY   "${DAISY_TOOLCHAIN_DIR}/bin/arm-none-eabi-objcopy" CACHE INTERNAL "Objcopy")
set(CMAKE_SIZE      "${DAISY_TOOLCHAIN_DIR}/bin/arm-none-eabi-size"    CACHE INTERNAL "Size")
set(CMAKE_OBJDUMP   "${DAISY_TOOLCHAIN_DIR}/bin/arm-none-eabi-objdump" CACHE INTERNAL "Objdump")
set(CMAKE_AR        "${DAISY_TOOLCHAIN_DIR}/bin/arm-none-eabi-ar"      CACHE INTERNAL "Archiver")
set(CMAKE_STRIP     "${DAISY_TOOLCHAIN_DIR}/bin/arm-none-eabi-strip"   CACHE INTERNAL "Strip")

# 5. Tell CMake not to look for the C and C++ standard libraries in the host system
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)