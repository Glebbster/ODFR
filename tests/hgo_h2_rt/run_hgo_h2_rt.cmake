cmake_minimum_required(VERSION 3.20)

foreach(required_var ODFR4_PROGRAM ODFR4_INPUT_DIR ODFR4_RUN_ROOT ODFR4_TEST_NAME)
  if(NOT DEFINED ${required_var} OR "${${required_var}}" STREQUAL "")
    message(FATAL_ERROR "Missing required variable: ${required_var}")
  endif()
endforeach()

if(NOT EXISTS "${ODFR4_PROGRAM}")
  message(FATAL_ERROR "Executable does not exist: ${ODFR4_PROGRAM}")
endif()

set(run_dir "${ODFR4_RUN_ROOT}/${ODFR4_TEST_NAME}")
file(REMOVE_RECURSE "${run_dir}")
file(MAKE_DIRECTORY "${run_dir}")
file(COPY "${ODFR4_INPUT_DIR}/" DESTINATION "${run_dir}")

set(odfr4_command "${ODFR4_PROGRAM}")
if(DEFINED ODFR4_MPIEXEC AND NOT "${ODFR4_MPIEXEC}" STREQUAL "")
  if(NOT DEFINED ODFR4_NUMPROC_FLAG OR NOT DEFINED ODFR4_NUMPROCS)
    message(FATAL_ERROR "MPI regression test requires ODFR4_NUMPROC_FLAG and ODFR4_NUMPROCS")
  endif()
  set(odfr4_command
    "${ODFR4_MPIEXEC}"
    "${ODFR4_NUMPROC_FLAG}"
    "${ODFR4_NUMPROCS}"
    "${ODFR4_PROGRAM}"
  )
endif()

execute_process(
  COMMAND ${odfr4_command}
  WORKING_DIRECTORY "${run_dir}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr
  TIMEOUT 840
)

file(WRITE "${run_dir}/stdout.log" "${stdout}")
file(WRITE "${run_dir}/stderr.log" "${stderr}")

if(NOT result EQUAL 0)
  message(FATAL_ERROR
    "HGO H2 RT regression test failed with result '${result}'.\n"
    "stdout:\n${stdout}\n"
    "stderr:\n${stderr}"
  )
endif()

foreach(expected_file history ODFres.par final1.par final1.dat)
  if(NOT EXISTS "${run_dir}/${expected_file}")
    message(FATAL_ERROR "Expected output was not created: ${expected_file}")
  endif()
endforeach()

file(READ "${run_dir}/history" history_text)
file(READ "${run_dir}/final1.dat" spectrum_text LIMIT 4096)
set(report_text "${stdout}\n${history_text}")

if(NOT report_text MATCHES "Final loss reporting:")
  message(FATAL_ERROR "Final loss report was not found in stdout/history")
endif()
if(NOT report_text MATCHES "SSE_W= *1\\.[0-9]+E-02")
  message(FATAL_ERROR "HGO weighted SSE left the expected 1.xE-02 range")
endif()
if(NOT report_text MATCHES "SSE_UNW=" OR
   NOT report_text MATCHES "RED_W=" OR
   NOT report_text MATCHES "RED_UNW=")
  message(FATAL_ERROR "Weighted/unweighted HGO loss fields were not found")
endif()
if(NOT report_text MATCHES "N_TOTAL= *1024")
  message(FATAL_ERROR "HGO report does not contain the expected 1024 points")
endif()
if(NOT report_text MATCHES "Spectrum contributions:")
  message(FATAL_ERROR "Per-spectrum HGO loss report was not found")
endif()
if(NOT spectrum_text MATCHES "s1\\.esr" OR
   NOT spectrum_text MATCHES "uglerodnik\\.par" OR
   NOT spectrum_text MATCHES "NPTS= *1024")
  message(FATAL_ERROR "final1.dat does not describe the expected HGO data set")
endif()

file(SIZE "${run_dir}/final1.dat" spectrum_size)
if(spectrum_size LESS 100000)
  message(FATAL_ERROR "final1.dat is unexpectedly short: ${spectrum_size} bytes")
endif()

message(STATUS
  "HGO H2 RT regression test '${ODFR4_TEST_NAME}' passed; outputs: ${run_dir}"
)
