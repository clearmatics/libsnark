set(SNARK_EXTRALIBS)
if(${CURVE} STREQUAL "BN128")
  set(
    SNARK_EXTRALIBS

    ${SNARK_EXTRALIBS}
    ${PROCPS_LIBRARIES}
    zm
  )
endif()

# TODO(madars): Revisit the decision to only offer a static
# library. In theory, it might be nice to offer a shared object, but
# it may be impractical for users, since many of the build-time
# options are mutually exclusive. At least as a static library with a
# simple, reproducible build, users of this library can pick the
# appropriate options themselves for embedding.
add_library(
  snark
  STATIC

  common/data_structures/set_commitment.cpp
  gadgetlib1/constraint_profiling.cpp
)
if("${USE_LINKED_LIBRARIES}")
  target_link_libraries(
    snark

    ff
    ${GMP_LIBRARIES}
    ${GMPXX_LIBRARIES}
    ${CRYPTO_LIBRARIES}
    ${PROCPS_LIBRARIES}
    ${LIBFF_LIBRARIES}
    ${SNARK_EXTRALIBS}
  )
else()
  target_link_libraries(
    snark

    ff
    ${GMP_LIBRARIES}
    ${GMPXX_LIBRARIES}
    ${CRYPTO_LIBRARIES}
    ${PROCPS_LIBRARIES}
    ${SNARK_EXTRALIBS}
  )
  target_include_directories(
    snark

    PUBLIC
    ${DEPENDS_DIR}/libff
    ${DEPENDS_DIR}/libfqfft
  )
endif()

install(
  DIRECTORY "" DESTINATION "include/libsnark"
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.tcc"
  PATTERN "tests" EXCLUDE
  PATTERN "examples" EXCLUDE
)

install(
  TARGETS
  snark

  DESTINATION lib
)

# Profiling
add_executable(
  profile_r1cs_ppzksnark
  EXCLUDE_FROM_ALL

  zk_proof_systems/ppzksnark/r1cs_ppzksnark/profiling/profile_r1cs_ppzksnark.cpp
)
target_link_libraries(
  profile_r1cs_ppzksnark

  snark
)

add_executable(
  profile_r1cs_se_ppzksnark
  EXCLUDE_FROM_ALL

  zk_proof_systems/ppzksnark/r1cs_se_ppzksnark/profiling/profile_r1cs_se_ppzksnark.cpp
)
target_link_libraries(
  profile_r1cs_se_ppzksnark

  snark
)

add_executable(
  profile_r1cs_gg_ppzksnark
  EXCLUDE_FROM_ALL

  zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/profiling/profile_r1cs_gg_ppzksnark.cpp
)
target_link_libraries(
  profile_r1cs_gg_ppzksnark

  snark
)

# Tests
if ("${IS_LIBSNARK_PARENT}")
  add_executable(
    gadgetlib1_simple_test
    EXCLUDE_FROM_ALL

    gadgetlib1/tests/gadgetlib1_test.cpp
  )
  target_link_libraries(
    gadgetlib1_simple_test

    snark
    gtest_main
  )

  add_executable(
    gadgetlib1_r1cs_ppzksnark_verifier_gadget_test
    EXCLUDE_FROM_ALL

    gadgetlib1/gadgets/verifiers/tests/test_r1cs_ppzksnark_verifier_gadget.cpp
  )
  target_link_libraries(
    gadgetlib1_r1cs_ppzksnark_verifier_gadget_test

    snark
  )

  add_executable(
    relations_qap_test
    EXCLUDE_FROM_ALL

    relations/arithmetic_programs/qap/tests/test_qap.cpp
  )
  target_link_libraries(
    relations_qap_test

    snark
  )

  add_executable(
    zk_proof_systems_r1cs_ppzksnark_test
    EXCLUDE_FROM_ALL

    zk_proof_systems/ppzksnark/r1cs_ppzksnark/tests/test_r1cs_ppzksnark.cpp
  )
  target_link_libraries(
    zk_proof_systems_r1cs_ppzksnark_test

    snark
  )

  add_executable(
    zk_proof_systems_r1cs_se_ppzksnark_test
    EXCLUDE_FROM_ALL

    zk_proof_systems/ppzksnark/r1cs_se_ppzksnark/tests/test_r1cs_se_ppzksnark.cpp
  )
  target_link_libraries(
    zk_proof_systems_r1cs_se_ppzksnark_test

    snark
  )

  add_executable(
    zk_proof_systems_r1cs_gg_ppzksnark_test
    EXCLUDE_FROM_ALL

    zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/tests/test_r1cs_gg_ppzksnark.cpp
  )
  target_link_libraries(
    zk_proof_systems_r1cs_gg_ppzksnark_test

    snark
  )

  add_executable(
    test_knapsack_gadget
    EXCLUDE_FROM_ALL

    gadgetlib1/gadgets/hashes/knapsack/tests/test_knapsack_gadget.cpp
  )
  target_link_libraries(
    test_knapsack_gadget

    snark
  )

  add_executable(
    test_merkle_tree_gadgets
    EXCLUDE_FROM_ALL

    gadgetlib1/gadgets/merkle_tree/tests/test_merkle_tree_gadgets.cpp
  )
  target_link_libraries(
    test_merkle_tree_gadgets

    snark
  )

  add_executable(
    test_set_commitment_gadget
    EXCLUDE_FROM_ALL

    gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp
  )
  target_link_libraries(
    test_set_commitment_gadget

    snark
  )

  add_executable(
    test_sha256_gadget
    EXCLUDE_FROM_ALL

    gadgetlib1/gadgets/hashes/sha256/tests/test_sha256_gadget.cpp
  )
  target_link_libraries(
    test_sha256_gadget

    snark
  )

  include(CTest)
  add_test(
    NAME gadgetlib1_simple_test
    COMMAND gadgetlib1_simple_test
  )
  add_test(
    NAME gadgetlib1_r1cs_ppzksnark_verifier_gadget_test
    COMMAND gadgetlib1_r1cs_ppzksnark_verifier_gadget_test
  )
  add_test(
    NAME relations_qap_test
    COMMAND relations_qap_test
  )
  add_test(
    NAME zk_proof_systems_r1cs_ppzksnark_test
    COMMAND zk_proof_systems_r1cs_ppzksnark_test
  )
  add_test(
    NAME zk_proof_systems_r1cs_se_ppzksnark_test
    COMMAND zk_proof_systems_r1cs_se_ppzksnark_test
  )
  add_test(
    NAME zk_proof_systems_r1cs_gg_ppzksnark_test
    COMMAND zk_proof_systems_r1cs_gg_ppzksnark_test
  )
  add_test(
    NAME test_knapsack_gadget
    COMMAND test_knapsack_gadget
  )
  add_test(
    NAME test_merkle_tree_gadgets
    COMMAND test_merkle_tree_gadgets
  )
  add_test(
    NAME test_set_commitment_gadget
    COMMAND test_set_commitment_gadget
  )
  add_test(
    NAME test_sha256_gadget
    COMMAND test_sha256_gadget
  )

  add_dependencies(check gadgetlib1_simple_test)
  add_dependencies(check gadgetlib1_r1cs_ppzksnark_verifier_gadget_test)
  add_dependencies(check relations_qap_test)
  add_dependencies(check zk_proof_systems_r1cs_ppzksnark_test)
  add_dependencies(check zk_proof_systems_r1cs_se_ppzksnark_test)
  add_dependencies(check zk_proof_systems_r1cs_gg_ppzksnark_test)
  add_dependencies(check test_knapsack_gadget)
  add_dependencies(check test_merkle_tree_gadgets)
  add_dependencies(check test_set_commitment_gadget)
  add_dependencies(check test_sha256_gadget)
endif()

# TODO (howardwu): Resolve runtime on targets:
# gadgetlib1_fooram_test, zk_proof_systems_r1cs_mp_ppzkpcd_test, zk_proof_systems_r1cs_sp_ppzkpcd_test, zk_proof_systems_ram_zksnark_test

# add_test(
#   NAME gadgetlib1_fooram_test
#   COMMAND gadgetlib1_fooram_test
# )
# add_test(
#   NAME zk_proof_systems_r1cs_mp_ppzkpcd_test
#   COMMAND zk_proof_systems_r1cs_mp_ppzkpcd_test
# )
# add_test(
#   NAME zk_proof_systems_r1cs_sp_ppzkpcd_test
#   COMMAND zk_proof_systems_r1cs_sp_ppzkpcd_test
# )
# add_test(
#   NAME zk_proof_systems_ram_zksnark_test
#   COMMAND zk_proof_systems_ram_zksnark_test
# )
# add_dependencies(check gadgetlib1_fooram_test)
# add_dependencies(check zk_proof_systems_r1cs_mp_ppzkpcd_test)
# add_dependencies(check zk_proof_systems_r1cs_sp_ppzkpcd_test)
# add_dependencies(check zk_proof_systems_ram_zksnark_test)
