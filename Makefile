# Project Name
TARGET = DAISY_SEED_SCREAM_STREAM

# Sources
CPP_SOURCES = ./src/steaming_mfcc.cpp

# Library Locations
LIBDAISY_DIR = ../DaisyExamples/libDaisy/
DAISYSP_DIR = ../DaisyExamples/DaisySP/

# C_SOURCES = $(wildcard ../DaisyExamples/libDaisy/Drivers/CMSIS/DSP/Source/*/*.c)
# C_INCLUDES += \
# -I$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Include
# # ASM_SOURCES = \
# # $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_bitreversal2.s
# C_DEFS += -DARM_TABLE_TWIDDLECOEF_F32_2048 -DARM_TABLE_BITREVIDX_FLT_2048 -DARM_TABLE_TWIDDLECOEF_RFFT_F32_4096

# ASM_SOURCES = \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_bitreversal2.S
# C_SOURCES = \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_rfft_fast_f32.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_rfft_fast_init_f32.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/CommonTables/arm_common_tables.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/CommonTables/arm_const_structs.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_cfft_f32.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_bitreversal.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_cfft_radix8_f32.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_cfft_radix4_f32.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/ComplexMathFunctions/arm_cmplx_mag_f32.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/ComplexMathFunctions/arm_cmplx_mag_squared_f32.c \
# $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/StatisticsFunctions/arm_max_f32.c
# C_INCLUDES += \
# -I$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Include

# C_DEFS += -DARM_DSP_CONFIG_TABLES -DARM_TABLE_TWIDDLECOEF_F32_1024 -DARM_TABLE_BITREVIDX_FLT_1024 -DARM_TABLE_TWIDDLECOEF_RFFT_F32_2048


# Try this simplified list
C_SOURCES = \
    $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_rfft_fast_init_f32.c \
    $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_rfft_fast_f32.c \
    $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/ComplexMathFunctions/arm_cmplx_mag_squared_f32.c \
    $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/CommonTables/arm_common_tables.c \
    $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/CommonTables/arm_const_structs.c
CPPFLAGS += -DARM_DSP_CONFIG_TABLES -DARM_TABLE_REAL_FFT_F32_2048



# Core location, and generic Makefile.
SYSTEM_FILES_DIR = $(LIBDAISY_DIR)/core
include $(SYSTEM_FILES_DIR)/Makefile