# Tested to run under MSYS2 bash with cmake in Windows 10
# This points to the location of the DaisyToolchain (Windows) install
set DAISY_TOOLCHAIN_DIR="/d/Code/LIBS/DaisyToolchain/"
cmake -DCMAKE_TOOLCHAIN_FILE=toolchain.cmake -B build -S .