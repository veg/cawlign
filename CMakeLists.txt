cmake_minimum_required(VERSION 3.10)

project(cawlign)

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)

set(CMAKE_CONFIGURATION_TYPES Release)


# SET VERSION FROM FILE
file (STRINGS "VERSION.txt" VERSION_NUMBER)
add_definitions (-DVERSION_NUMBER=\"${VERSION_NUMBER}\")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION.txt")

include_directories(
    src/
)

add_executable(
    cawlign
    src/alignment.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse.cpp
    src/cawlign.cpp
    src/configparser.cpp
    src/scoring.cpp
)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64|ARM64")
    target_compile_options (cawlign PRIVATE -fsigned-char -O3 -std=c++14 -mcpu=native)
    target_link_options (cawlign PRIVATE -O3 -mcpu=native -std=c++14 -fsigned-char)
else()
    target_compile_options (cawlign PRIVATE -fsigned-char -O3 -std=c++14 -march=native)
    target_link_options (cawlign PRIVATE -O3 -march=native -std=c++14 -fsigned-char)
endif()


add_executable(
    cawlign_debug EXCLUDE_FROM_ALL
    src/alignment.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse.cpp
    src/cawlign.cpp
    src/configparser.cpp
    src/scoring.cpp
    
)

target_compile_options (cawlign_debug PRIVATE -O0 -std=c++14  -g -fsanitize=address)
target_link_options (cawlign_debug PRIVATE -O0 -std=c++14  -g -fsanitize=address)


add_executable(
    cawlign_prof EXCLUDE_FROM_ALL
    src/alignment.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse.cpp
    src/cawlign.cpp
    src/configparser.cpp
    src/scoring.cpp
)

target_compile_options (cawlign_prof PRIVATE -g -O3 -std=c++14  -funroll-loops)
target_link_options (cawlign_prof PRIVATE -g -O3 -march=native -std=c++14 -funroll-loops)


set_property(
    TARGET cawlign
    APPEND PROPERTY COMPILE_DEFINITIONS LIBRARY_PATH="${INSTALL_PREFIX}/share/cawlign"
)





#add_definitions(${DEFAULT_COMPILE_OPTIONS})

find_package(OpenMP)

if(${OPENMP_FOUND})
   add_definitions(${OpenMP_CXX_FLAGS})
   target_link_libraries(cawlign PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(cawlign_debug PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(cawlign_prof PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
endif(${OPENMP_FOUND})


#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------

option(ENABLE_COVERAGE "Enable code coverage" OFF)

if(ENABLE_COVERAGE)
    target_compile_options(cawlign PRIVATE --coverage -O0)
    target_link_options(cawlign PRIVATE --coverage)
    
    find_program(LCOV_PATH lcov)
    find_program(GENHTML_PATH genhtml)
    
    if(LCOV_PATH AND GENHTML_PATH)
        add_custom_target(coverage
            COMMAND ${LCOV_PATH} --capture --directory . --output-file coverage.info --ignore-errors mismatch --ignore-errors unused
            COMMAND ${LCOV_PATH} --remove coverage.info '/usr/*' '*/_deps/*' '*/test/*' --output-file coverage.filtered.info
            COMMAND ${GENHTML_PATH} coverage.filtered.info --output-directory coverage_report
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating coverage report"
        )
    endif()
endif()

enable_testing()

set(TEST_PYTHON_EXE python3)
set(TEST_RUNNER ${CMAKE_SOURCE_DIR}/test/test_runner.py)
set(BASELINES_DIR ${CMAKE_SOURCE_DIR}/test/baselines)
set(CASES_DIR ${CMAKE_SOURCE_DIR}/test/cases)
set(RES_DIR ${CMAKE_SOURCE_DIR}/res)

# Define tests
add_test(NAME nuc_default COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_default.fa $<TARGET_FILE:cawlign> -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404.fa)
add_test(NAME nuc_linear COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_linear.fa $<TARGET_FILE:cawlign> -S linear -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404.fa)
add_test(NAME nuc_local COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_local.fa $<TARGET_FILE:cawlign> -l local -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404.fa)
add_test(NAME nuc_custom_scoring COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_custom_scoring.fa $<TARGET_FILE:cawlign> -s ${RES_DIR}/scoring/NUC4.4 -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404.fa)

add_test(NAME prot_default COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/prot_default.fa $<TARGET_FILE:cawlign> -t protein -r ${RES_DIR}/references/HXB2_pr -s ${RES_DIR}/scoring/BLOSUM62_prot -q ${CASES_DIR}/OQ948666-pr.fa)

add_test(NAME codon_default COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_default.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/universal -r ${RES_DIR}/references/HXB2_pr -s ${RES_DIR}/scoring/BLOSUM62 -q ${CASES_DIR}/OQ948666-pr.fa)
add_test(NAME codon_mtdna COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_mtdna.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/vertebrate-mtdna -r ${RES_DIR}/references/COXI-mouse -s ${RES_DIR}/scoring/BLOSUM62 -q ${CASES_DIR}/R_norvegicus_COXI.fa)

add_test(NAME nuc_rc COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_rc.fa $<TARGET_FILE:cawlign> -R annotated -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404-rc.fa)

add_test(NAME nuc_refalign COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_refalign.fa $<TARGET_FILE:cawlign> -f refalign -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404.fa)
add_test(NAME nuc_pairwise COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_pairwise.fa $<TARGET_FILE:cawlign> -f pairwise -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404.fa)
add_test(NAME nuc_refmap_with_ref COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_refmap_with_ref.fa $<TARGET_FILE:cawlign> -I -f refmap -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404.fa)

#-------------------------------------------------------------------------------
# Negative & Edge Case Tests
#-------------------------------------------------------------------------------

# Invalid Datatype
add_test(NAME negative_invalid_datatype COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} --expect-fail --expect-error "invalid data type: fake_type" dummy.fa $<TARGET_FILE:cawlign> -t fake_type -q ${CASES_DIR}/KX988404.fa)

# Invalid Genetic Code
add_test(NAME negative_invalid_genetic_code COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} --expect-fail --expect-error "unknown NCBI genetic code number 999" dummy.fa $<TARGET_FILE:cawlign> -t codon -c 999 -q ${CASES_DIR}/OQ948666-pr.fa)
set_tests_properties(negative_invalid_genetic_code PROPERTIES ENVIRONMENT "CAWLIGN_LIBRARY_PATH=${RES_DIR}")

# Incompatible Options (Protein + RC)
add_test(NAME negative_protein_rc COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} --expect-fail --expect-error "Reverse complement options other than 'none' are not compatible with the protein data type" dummy.fa $<TARGET_FILE:cawlign> -t protein -R annotated -q ${CASES_DIR}/OQ948666-pr.fa)

# Codon Mode - Reference length not divisible by 3
add_test(NAME negative_codon_div3 COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} --expect-fail --expect-error "The reference sequence must have length divisible by 3" dummy.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${CMAKE_SOURCE_DIR}/test/edge_cases/div3_fail.fa -q ${CASES_DIR}/OQ948666-pr.fa)

# Protein Mode - codon scoring file (no [ALPHABET] section)
add_test(NAME negative_protein_codon_scoring COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} --expect-fail --expect-error "This scoring file has no [ALPHABET] section" dummy.fa $<TARGET_FILE:cawlign> -t protein -s ${RES_DIR}/scoring/BLOSUM62 -r ${RES_DIR}/references/HXB2_pr -q ${CASES_DIR}/OQ948666-pr.fa)

# Codon Mode - Stop codon in reference
add_test(NAME negative_codon_stop_ref COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} --expect-fail --expect-error "The reference sequence must not have stop codons in it" dummy.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${CMAKE_SOURCE_DIR}/test/edge_cases/stop_ref.fa -q ${CASES_DIR}/OQ948666-pr.fa)

# Edge Case - Empty input file (should handle gracefully)
add_test(NAME edge_empty_input COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} --allow-empty ${CMAKE_SOURCE_DIR}/test/edge_cases/empty.fa $<TARGET_FILE:cawlign> -q ${CMAKE_SOURCE_DIR}/test/edge_cases/empty.fa)

# Edge Case - Very short sequence
add_test(NAME edge_short_seq COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/short_seq.fa $<TARGET_FILE:cawlign> -r ${RES_DIR}/references/HXB2_pol -q ${CMAKE_SOURCE_DIR}/test/edge_cases/short_qry.fa)

#-------------------------------------------------------------------------------
# Advanced Coverage Tests
#-------------------------------------------------------------------------------

# Codon Slippage (24-move state machine stress test)
add_test(NAME codon_slippage COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_slippage.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/int-slippage.fa)

# Non-Affine Scoring (-a)
add_test(NAME nuc_non_affine COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_non_affine.fa $<TARGET_FILE:cawlign> -a -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404.fa)

# Nucleotide Ambiguity Handling
add_test(NAME nuc_ambig COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/nuc_ambig.fa $<TARGET_FILE:cawlign> -r ${CMAKE_SOURCE_DIR}/test/edge_cases/ambig_nuc.fa -q ${CMAKE_SOURCE_DIR}/test/edge_cases/ambig_nuc.fa)

# Codon Quiet Mode stderr Suppression Check (verify that -q produces no stderr output)
add_test(NAME nuc_quiet COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} --check-empty-stderr ${BASELINES_DIR}/nuc_default.fa $<TARGET_FILE:cawlign> -q -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/KX988404.fa)

#-------------------------------------------------------------------------------
# Complex Codon Frameshift Tests (State Machine Coverage)
#-------------------------------------------------------------------------------

# +1 Insertion in query
add_test(NAME codon_fs_plus_1 COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_fs_plus_1.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${CMAKE_SOURCE_DIR}/test/edge_cases/codon_fs_ref.fa -q ${CMAKE_SOURCE_DIR}/test/edge_cases/codon_fs_plus_1.fa)

# -1 Deletion in query
add_test(NAME codon_fs_minus_1 COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_fs_minus_1.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${CMAKE_SOURCE_DIR}/test/edge_cases/codon_fs_ref.fa -q ${CMAKE_SOURCE_DIR}/test/edge_cases/codon_fs_minus_1.fa)

# +4 Insertion in query (Large transition)
add_test(NAME codon_fs_plus_4 COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_fs_plus_4.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${CMAKE_SOURCE_DIR}/test/edge_cases/codon_fs_ref.fa -q ${CMAKE_SOURCE_DIR}/test/edge_cases/codon_fs_plus_4.fa)

# Pathological State Machine Stress Tests (109 edge cases against HXB2_pr and HXB2_pol)
add_test(NAME codon_pathological_pr COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/pathological_cases.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${RES_DIR}/references/HXB2_pr -q ${CMAKE_SOURCE_DIR}/test/pathological/pathological_cases.fa)
add_test(NAME codon_pathological_pol COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/pathological_cases_pol.fa $<TARGET_FILE:cawlign> -t codon -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${RES_DIR}/references/HXB2_pol -q ${CMAKE_SOURCE_DIR}/test/pathological/pathological_cases.fa)

# Regression tests for sparse/divergent corridor fallback
add_test(NAME codon_kx988404_pol COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_kx988404_pol.fa $<TARGET_FILE:cawlign> -q -t codon -s ${RES_DIR}/scoring/BLOSUM62 -c ${RES_DIR}/genetic_codes/universal -r ${RES_DIR}/references/HXB2_pol ${CASES_DIR}/KX988404_codon.fa)
add_test(NAME codon_kf576595_rt COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_kf576595_rt.fa $<TARGET_FILE:cawlign> -q -t codon -s ${RES_DIR}/scoring/BLOSUM62 -c ${RES_DIR}/genetic_codes/universal -r ${RES_DIR}/references/HXB2_rt ${CASES_DIR}/KF576595_rt.fa)
add_test(NAME codon_ins123_pol COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_ins123_pol.fa $<TARGET_FILE:cawlign> -q -t codon -s ${RES_DIR}/scoring/BLOSUM62 -c ${RES_DIR}/genetic_codes/universal -r ${RES_DIR}/references/HXB2_pol ${CASES_DIR}/ins123_pol.fa)
add_test(NAME codon_margin400_pol COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_margin400_pol.fa $<TARGET_FILE:cawlign> -q -t codon -s ${RES_DIR}/scoring/BLOSUM62 -c ${RES_DIR}/genetic_codes/universal -r ${RES_DIR}/references/HXB2_pol ${CASES_DIR}/margin400_pol.fa)

# Regression tests for corridor soundness: unanchored leading segments and internal inversions (baselines from main)
add_test(NAME codon_band_lead_inv_pol COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_band_lead_inv_pol.fa $<TARGET_FILE:cawlign> -q -t codon -s ${RES_DIR}/scoring/BLOSUM62 -c ${RES_DIR}/genetic_codes/universal -r ${RES_DIR}/references/HXB2_pol ${CASES_DIR}/band_lead_inv_pol.fa)
add_test(NAME codon_band_lead_inv_gag COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_band_lead_inv_gag.fa $<TARGET_FILE:cawlign> -q -t codon -s ${RES_DIR}/scoring/BLOSUM62 -c ${RES_DIR}/genetic_codes/universal -r ${RES_DIR}/references/HXB2_gag ${CASES_DIR}/band_lead_inv_gag.fa)

# Codon global/local output must not depend on thread count, record order, or malloc fill pattern
add_test(NAME codon_affine_determinism COMMAND ${CMAKE_SOURCE_DIR}/test/codon_affine_determinism.sh $<TARGET_FILE:cawlign> ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/test-tmp)

# Regression tests for codon global/local affine col-0 insertion init (baselines from main)
add_test(NAME codon_global_affine COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_global_affine.fa $<TARGET_FILE:cawlign> -t codon -l global -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/lanl-pol-10.fa)
add_test(NAME codon_local_affine COMMAND ${TEST_PYTHON_EXE} ${TEST_RUNNER} ${BASELINES_DIR}/codon_local_affine.fa $<TARGET_FILE:cawlign> -t codon -l local -c ${RES_DIR}/genetic_codes/universal -s ${RES_DIR}/scoring/BLOSUM62 -r ${RES_DIR}/references/HXB2_pol -q ${CASES_DIR}/lanl-pol-10.fa)


#-------------------------------------------------------------------------------
# Installation
#-------------------------------------------------------------------------------

install(
    TARGETS cawlign
	RUNTIME DESTINATION bin
	OPTIONAL
)

install(DIRECTORY res/ DESTINATION share/cawlign)

