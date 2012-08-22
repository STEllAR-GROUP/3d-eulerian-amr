################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/src/assert.cpp \
../src/src/indexer_2d.cpp \
../src/src/indexer_2d_by2.cpp \
../src/src/main.cpp \
../src/src/program.cpp \
../src/src/reconstruct.cpp 

OBJS += \
./src/src/assert.o \
./src/src/indexer_2d.o \
./src/src/indexer_2d_by2.o \
./src/src/main.o \
./src/src/program.o \
./src/src/reconstruct.o 

CPP_DEPS += \
./src/src/assert.d \
./src/src/indexer_2d.d \
./src/src/indexer_2d_by2.d \
./src/src/main.d \
./src/src/program.d \
./src/src/reconstruct.d 


# Each subdirectory must supply rules for building sources it contributes
src/src/%.o: ../src/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	icpc -g -O0 -I"/home/dmarce1/workspace/amr/src" -openmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


