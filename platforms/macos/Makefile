ROOT := ../..
BUILD_DIR := $(ROOT)/out/macos
CMAKE_COMMON := -S "$(ROOT)" -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=ON

all: serial

serial:
	cmake $(CMAKE_COMMON) -DODFR4_BUILD_SERIAL=ON -DODFR4_BUILD_MPI=OFF
	cmake --build "$(BUILD_DIR)" --parallel --target ODFR4_3_serial

mpi:
	cmake $(CMAKE_COMMON) -DODFR4_BUILD_SERIAL=ON -DODFR4_BUILD_MPI=ON
	cmake --build "$(BUILD_DIR)" --parallel --target ODFR4_3_mpi

test-serial: serial
	ctest --test-dir "$(BUILD_DIR)" -L serial --output-on-failure

test-mpi: mpi
	ctest --test-dir "$(BUILD_DIR)" -L mpi --output-on-failure

clean:
	cmake -E remove_directory "$(BUILD_DIR)"

.PHONY: all serial mpi test-serial test-mpi clean
