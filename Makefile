# Project Name
TARGET = DAISY_SEED_CLAPPER_RMS_FIRMWARE

# Sources
CPP_SOURCES = main_clapper_rms.cpp

# Library Locations
LIBDAISY_DIR = ../DaisyExamples/libDaisy/
DAISYSP_DIR = ../DaisyExamples/DaisySP/

# Core location, and generic Makefile.
SYSTEM_FILES_DIR = $(LIBDAISY_DIR)/core
include $(SYSTEM_FILES_DIR)/Makefile