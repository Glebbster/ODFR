cmake_minimum_required(VERSION 3.20)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" ODFR4_VERSION
  LIMIT_COUNT 1
)
string(STRIP "${ODFR4_VERSION}" ODFR4_VERSION)
if(NOT ODFR4_VERSION MATCHES
   "^[0-9]+\\.[0-9]+\\.[0-9]+(-[0-9A-Za-z.-]+)?$")
  message(FATAL_ERROR
    "VERSION must contain a version such as 4.3.0 or 4.3.1-dev"
  )
endif()
string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+"
  ODFR4_PROJECT_VERSION "${ODFR4_VERSION}"
)

project(ODFR VERSION "${ODFR4_PROJECT_VERSION}" LANGUAGES Fortran)

option(ODFR4_BUILD_SERIAL "Build the single-process executable without an MPI installation" ON)
option(ODFR4_BUILD_MPI "Build the real MPI executable" OFF)

if(NOT ODFR4_BUILD_SERIAL AND NOT ODFR4_BUILD_MPI)
  message(FATAL_ERROR "Enable ODFR4_BUILD_SERIAL, ODFR4_BUILD_MPI, or both")
endif()

set(ODFR4_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(ODFR4_GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${ODFR4_GENERATED_DIR}")
configure_file(
  "${ODFR4_SOURCE_DIR}/version.inc.in"
  "${ODFR4_GENERATED_DIR}/version.inc"
  @ONLY
)

set(ODFR4_SOURCES
  "${ODFR4_SOURCE_DIR}/optimization_steps.for"
  "${ODFR4_SOURCE_DIR}/anxlk.f"
  "${ODFR4_SOURCE_DIR}/bessel.f"
  "${ODFR4_SOURCE_DIR}/ccrints.f"
  "${ODFR4_SOURCE_DIR}/cd2km.f"
  "${ODFR4_SOURCE_DIR}/cfs.f"
  "${ODFR4_SOURCE_DIR}/cgltrimy.f"
  "${ODFR4_SOURCE_DIR}/choiceExp.for"
  "${ODFR4_SOURCE_DIR}/covar.f"
  "${ODFR4_SOURCE_DIR}/covrpt.f"
  "${ODFR4_SOURCE_DIR}/cscg.f"
  "${ODFR4_SOURCE_DIR}/enorm.f"
  "${ODFR4_SOURCE_DIR}/eprls.f"
  "${ODFR4_SOURCE_DIR}/fitl.f"
  "${ODFR4_SOURCE_DIR}/fitNL2.for"
  "${ODFR4_SOURCE_DIR}/fz.f"
  "${ODFR4_SOURCE_DIR}/gconvl.f"
  "${ODFR4_SOURCE_DIR}/inkor.for"
  "${ODFR4_SOURCE_DIR}/INPNamDR.for"
  "${ODFR4_SOURCE_DIR}/ipar.f"
  "${ODFR4_SOURCE_DIR}/lfunR.f"
  "${ODFR4_SOURCE_DIR}/librat.for"
  "${ODFR4_SOURCE_DIR}/lmnls.f"
  "${ODFR4_SOURCE_DIR}/lmpar.f"
  "${ODFR4_SOURCE_DIR}/MADR.f"
  "${ODFR4_SOURCE_DIR}/matrll.f"
  "${ODFR4_SOURCE_DIR}/momdls.f"
  "${ODFR4_SOURCE_DIR}/Nl2def1.for"
  "${ODFR4_SOURCE_DIR}/nl2sol_2015.f90"
  "${ODFR4_SOURCE_DIR}/NLSINIT.f"
  "${ODFR4_SOURCE_DIR}/ODFint.for"
  "${ODFR4_SOURCE_DIR}/ODFR4.f"
  "${ODFR4_SOURCE_DIR}/ODFread.for"
  "${ODFR4_SOURCE_DIR}/ODFwrite.for"
  "${ODFR4_SOURCE_DIR}/PARreadN.for"
  "${ODFR4_SOURCE_DIR}/PARwrN.for"
  "${ODFR4_SOURCE_DIR}/plgndr.f"
  "${ODFR4_SOURCE_DIR}/qrfac.f"
  "${ODFR4_SOURCE_DIR}/qrsolv.f"
  "${ODFR4_SOURCE_DIR}/REALFT.f"
  "${ODFR4_SOURCE_DIR}/scmvm.f"
  "${ODFR4_SOURCE_DIR}/SpecOut.for"
  "${ODFR4_SOURCE_DIR}/specread.for"
  "${ODFR4_SOURCE_DIR}/stats.f"
  "${ODFR4_SOURCE_DIR}/stvect.f"
  "${ODFR4_SOURCE_DIR}/tchange.f"
  "${ODFR4_SOURCE_DIR}/w3j.f"
  "${ODFR4_SOURCE_DIR}/WignerL.for"
  "${ODFR4_SOURCE_DIR}/Wignermy.f"
  "${ODFR4_SOURCE_DIR}/xpack.f"
  "${ODFR4_SOURCE_DIR}/zaxpy.f"
  "${ODFR4_SOURCE_DIR}/zaypx.f"
  "${ODFR4_SOURCE_DIR}/zdotu.f"
)

function(odfr4_configure_target target_name)
  target_include_directories("${target_name}" PRIVATE
    "${ODFR4_GENERATED_DIR}"
    "${ODFR4_SOURCE_DIR}"
  )
  set_target_properties("${target_name}" PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/${target_name}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )

  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options("${target_name}" PRIVATE
      -fbacktrace
      -ffixed-line-length-132
      -std=legacy
      -fallow-invalid-boz
      -fallow-argument-mismatch
      -fdec
      $<$<CONFIG:Debug>:-fcheck=all>
      $<$<CONFIG:Debug>:-ffpe-trap=invalid,zero,overflow>
    )
  endif()
endfunction()

if(ODFR4_BUILD_SERIAL)
  add_executable(ODFR4_3_serial
    ${ODFR4_SOURCES}
    "${CMAKE_CURRENT_SOURCE_DIR}/compat/serial/mpi_stubs.f90"
  )
  target_include_directories(ODFR4_3_serial BEFORE PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/compat/serial"
  )
  set_target_properties(ODFR4_3_serial PROPERTIES OUTPUT_NAME ODFR5_serial)
  odfr4_configure_target(ODFR4_3_serial)
endif()

if(ODFR4_BUILD_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  add_executable(ODFR4_3_mpi ${ODFR4_SOURCES})
  target_link_libraries(ODFR4_3_mpi PRIVATE MPI::MPI_Fortran)
  set_target_properties(ODFR4_3_mpi PROPERTIES OUTPUT_NAME ODFR5_mpi)
  odfr4_configure_target(ODFR4_3_mpi)
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(optimization_steps_probe
    "${ODFR4_SOURCE_DIR}/optimization_steps.for"
    "${ODFR4_SOURCE_DIR}/nl2sol_2015.f90"
    "${ODFR4_SOURCE_DIR}/Nl2def1.for"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/optimization_steps/probe.for"
  )
  odfr4_configure_target(optimization_steps_probe)
  find_package(Python3 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME optimization_steps_observer
      COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/optimization_steps/check_probe.py"
        "$<TARGET_FILE:optimization_steps_probe>"
        "${CMAKE_CURRENT_BINARY_DIR}/tests/observer"
    )
    set_tests_properties(optimization_steps_observer PROPERTIES TIMEOUT 60)
  endif()

  if(ODFR4_BUILD_SERIAL)
    add_test(NAME smoke_serial
      COMMAND "${CMAKE_COMMAND}"
        "-DODFR4_PROGRAM=$<TARGET_FILE:ODFR4_3_serial>"
        "-DODFR4_INPUT_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/input"
        "-DODFR4_RUN_ROOT=${CMAKE_CURRENT_BINARY_DIR}/tests"
        "-DODFR4_TEST_NAME=serial"
        "-DODFR4_EXPECTED_VERSION=${ODFR4_VERSION}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/run_smoke.cmake"
    )
    set_tests_properties(smoke_serial PROPERTIES LABELS "smoke;serial" TIMEOUT 900)

    add_test(NAME hgo_h2_rt_serial
      COMMAND "${CMAKE_COMMAND}"
        "-DODFR4_PROGRAM=$<TARGET_FILE:ODFR4_3_serial>"
        "-DODFR4_INPUT_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests/hgo_h2_rt/input"
        "-DODFR4_RUN_ROOT=${CMAKE_CURRENT_BINARY_DIR}/tests"
        "-DODFR4_TEST_NAME=hgo_h2_rt_serial"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/hgo_h2_rt/run_hgo_h2_rt.cmake"
    )
    set_tests_properties(hgo_h2_rt_serial PROPERTIES
      LABELS "regression;hgo_h2_rt;serial"
      TIMEOUT 900
    )
  endif()

  if(ODFR4_BUILD_MPI)
    add_test(NAME smoke_mpi_2
      COMMAND "${CMAKE_COMMAND}"
        "-DODFR4_PROGRAM=$<TARGET_FILE:ODFR4_3_mpi>"
        "-DODFR4_INPUT_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/input"
        "-DODFR4_RUN_ROOT=${CMAKE_CURRENT_BINARY_DIR}/tests"
        "-DODFR4_TEST_NAME=mpi_2"
        "-DODFR4_EXPECTED_VERSION=${ODFR4_VERSION}"
        "-DODFR4_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DODFR4_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DODFR4_NUMPROCS=2"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/run_smoke.cmake"
    )
    set_tests_properties(smoke_mpi_2 PROPERTIES LABELS "smoke;mpi" TIMEOUT 900)

    add_test(NAME hgo_h2_rt_mpi_2
      COMMAND "${CMAKE_COMMAND}"
        "-DODFR4_PROGRAM=$<TARGET_FILE:ODFR4_3_mpi>"
        "-DODFR4_INPUT_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests/hgo_h2_rt/input"
        "-DODFR4_RUN_ROOT=${CMAKE_CURRENT_BINARY_DIR}/tests"
        "-DODFR4_TEST_NAME=hgo_h2_rt_mpi_2"
        "-DODFR4_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DODFR4_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DODFR4_NUMPROCS=2"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/hgo_h2_rt/run_hgo_h2_rt.cmake"
    )
    set_tests_properties(hgo_h2_rt_mpi_2 PROPERTIES
      LABELS "regression;hgo_h2_rt;mpi"
      TIMEOUT 900
    )
  endif()
endif()
