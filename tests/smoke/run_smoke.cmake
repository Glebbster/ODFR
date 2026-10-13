cmake_minimum_required(VERSION 3.20)

foreach(required_var ODFR4_PROGRAM ODFR4_INPUT_DIR ODFR4_RUN_ROOT
                     ODFR4_TEST_NAME ODFR4_EXPECTED_VERSION)
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
string(REPEAT "5\n" 512 smoke_weights)
file(WRITE "${run_dir}/weight" "${smoke_weights}")

set(odfr4_command "${ODFR4_PROGRAM}")
if(DEFINED ODFR4_MPIEXEC AND NOT "${ODFR4_MPIEXEC}" STREQUAL "")
  if(NOT DEFINED ODFR4_NUMPROC_FLAG OR NOT DEFINED ODFR4_NUMPROCS)
    message(FATAL_ERROR "MPI smoke test requires ODFR4_NUMPROC_FLAG and ODFR4_NUMPROCS")
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
    "ODFR4 smoke test failed with result '${result}'.\n"
    "stdout:\n${stdout}\n"
    "stderr:\n${stderr}"
  )
endif()

foreach(expected_file history ODFR4_out.par ODFR4F_out.par smoke_out.dat)
  if(NOT EXISTS "${run_dir}/${expected_file}")
    message(FATAL_ERROR "Expected output was not created: ${expected_file}")
  endif()
endforeach()

file(READ "${run_dir}/history" history_text)
file(READ "${run_dir}/smoke_out.dat" spectrum_text LIMIT 4096)
set(report_text "${stdout}\n${history_text}")

set(expected_version_line "ODFR5 version ${ODFR4_EXPECTED_VERSION}")
string(FIND "${history_text}" "${expected_version_line}" version_position)
if(version_position EQUAL -1)
  message(FATAL_ERROR
    "Expected ODFR4 version was not written to history: ${ODFR4_EXPECTED_VERSION}"
  )
endif()
if(NOT report_text MATCHES "Weights OK")
  message(FATAL_ERROR "The generated 512-point weight file was not accepted")
endif()
if(NOT report_text MATCHES "Final loss reporting:")
  message(FATAL_ERROR "Final loss report was not found in stdout/history")
endif()
if(NOT report_text MATCHES "SSE_W=" OR NOT report_text MATCHES "SSE_UNW=")
  message(FATAL_ERROR "Weighted/unweighted SSE fields were not found")
endif()
if(NOT report_text MATCHES "Spectrum contributions:")
  message(FATAL_ERROR "Per-spectrum loss report was not found")
endif()
if(NOT spectrum_text MATCHES "SSE_W=" OR NOT spectrum_text MATCHES "SSE_UNW=")
  message(FATAL_ERROR "Loss fields were not written to smoke_out.dat")
endif()

message(STATUS "ODFR4 smoke test '${ODFR4_TEST_NAME}' passed; outputs: ${run_dir}")
