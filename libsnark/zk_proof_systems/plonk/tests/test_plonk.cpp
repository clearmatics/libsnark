/** @file
 *****************************************************************************
 Test program that exercises the Plonk protocol (first setup, then
 prover, then verifier) on a synthetic R1CS instance.

 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <gtest/gtest.h>

#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>

#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>
#include <libfqfft/polynomial_arithmetic/naive_evaluate.hpp>

#include <libsnark/zk_proof_systems/plonk/plonk.hpp>
#include <libsnark/zk_proof_systems/plonk/common_input.hpp>
#include <libsnark/zk_proof_systems/plonk/srs.hpp>
#include <libsnark/zk_proof_systems/plonk/prover.hpp>
#include <libsnark/zk_proof_systems/plonk/verifier.hpp>

// maximum polynomial degree (resp. maximum number of gates)
const size_t MAX_DEGREE = 254;

namespace libsnark
{

  template<typename ppT> void test_plonk()
  {
#ifndef DEBUG
    printf("[%s:%d] DEBUG info disabled. Please, enable. Terminate...\n", __FILE__, __LINE__);
    assert(0);    
#endif // #ifndef DEBUG
    
    // Execute all tests for the given curve.
    ppT::init_public_params();

    using Field = libff::Fr<ppT>;

    // --- TEST VECTORS ---    
    // load test vector values from example circuit
    plonk_example<ppT> example;
    
    // --- SETUP ---
    printf("[%s:%d] Common preprocessed input...\n", __FILE__, __LINE__);
    // common preprocessed input (CPI)
    common_preprocessed_input<ppT> common_input;
    // setup the CPI using an example circuit
    common_input.setup_from_example(example);

    // random hidden element secret (toxic waste). we fix it to a
    // constant in order to match against the test vectors
    Field secret = example.secret;
#ifdef DEBUG
    printf("[%s:%d] secret ", __FILE__, __LINE__);
    secret.print();
#endif // #ifdef DEBUG

    // --- SRS ---
    printf("[%s:%d] SRS...\n", __FILE__, __LINE__);
    // create SRS object
    srs<ppT> srs;
    // compute SRS = powers of secret times G1: 1*G1, secret^1*G1,
    // secret^2*G1, ... and secret times G2: 1*G2, secret^1*G2
    srs.derive_from_secret(secret, common_input.num_gates);    
    // sompare SRS against reference test values
#ifdef DEBUG
    for (int i = 0; i < (int)common_input.num_gates + 3; ++i) {
      printf("secret_power_G1[%2d] ", i);
      srs.secret_powers_g1[i].print();
      // test from generator
      libff::G1<ppT> srs_secret_powers_g1_i(srs.secret_powers_g1[i]);
      srs_secret_powers_g1_i.to_affine_coordinates();
      assert(srs_secret_powers_g1_i.X == example.secret_powers_g1[i][0]);
      assert(srs_secret_powers_g1_i.Y == example.secret_powers_g1[i][1]);
    }
    for (int i = 0; i < 2; ++i) {
      printf("secret_power_G2[%2d] ", i);
      srs.secret_powers_g2[i].print();
    }
#endif // #ifdef DEBUG
    
    // --- PROVER ---
    printf("[%s:%d] Prover...\n", __FILE__, __LINE__);        
    // initialize prover
    plonk_prover<ppT> prover(common_input);
    // compute proof
    plonk_proof<ppT> proof = prover.compute_proof(srs, common_input);
    // compare proof against test vector values (debug)
#ifdef DEBUG
    assert(proof.a_zeta == example.a_zeta);
    assert(proof.b_zeta == example.b_zeta);
    assert(proof.c_zeta == example.c_zeta);
    assert(proof.S_0_zeta == example.S_0_zeta);
    assert(proof.S_1_zeta == example.S_1_zeta);
    assert(proof.z_poly_xomega_zeta == example.z_poly_xomega_zeta);
    assert(proof.r_zeta == example.r_zeta);
    //    assert(proof.W_polys_blinded_at_secret_g1 == example.W_polys_blinded_at_secret_g1);
    for (int i = 0; i < (int)NUM_HGEN; ++i) {
      printf("W_polys_at_secret_g1[%d]\n", i);
      proof.W_polys_blinded_at_secret_g1[i].print();
      libff::G1<ppT> W_polys_blinded_at_secret_g1_i(proof.W_polys_blinded_at_secret_g1[i]);
      W_polys_blinded_at_secret_g1_i.to_affine_coordinates();
      assert(W_polys_blinded_at_secret_g1_i.X == example.W_polys_blinded_at_secret_g1[i][0]);
      assert(W_polys_blinded_at_secret_g1_i.Y == example.W_polys_blinded_at_secret_g1[i][1]);
    }
    //    assert(proof.z_poly_at_secret_g1 == example.z_poly_at_secret_g1);
    proof.z_poly_at_secret_g1.print();
    libff::G1<ppT> z_poly_at_secret_g1_aff(proof.z_poly_at_secret_g1);
    z_poly_at_secret_g1_aff.to_affine_coordinates();
    assert(z_poly_at_secret_g1_aff.X == example.z_poly_at_secret_g1[0]);
    assert(z_poly_at_secret_g1_aff.Y == example.z_poly_at_secret_g1[1]);
    //    assert(proof.t_poly_at_secret_g1 == example.t_poly_at_secret_g1);
    for (int i = 0; i < (int)NUM_HGEN; ++i) {
      printf("[%s:%d] t_poly_at_secret_g1[%d]\n", __FILE__, __LINE__, i);
      proof.t_poly_at_secret_g1[i].print();
      libff::G1<ppT> t_poly_at_secret_g1_i(proof.t_poly_at_secret_g1[i]);
      t_poly_at_secret_g1_i.to_affine_coordinates();
      assert(t_poly_at_secret_g1_i.X == example.t_poly_at_secret_g1[i][0]);
      assert(t_poly_at_secret_g1_i.Y == example.t_poly_at_secret_g1[i][1]);
    }
    //    assert(proof.W_zeta_at_secret == example.W_zeta_at_secret);
    proof.W_zeta_at_secret.print();
    libff::G1<ppT> W_zeta_at_secret_aff(proof.W_zeta_at_secret);
    W_zeta_at_secret_aff.to_affine_coordinates();
    assert(W_zeta_at_secret_aff.X == example.W_zeta_at_secret[0]);
    assert(W_zeta_at_secret_aff.Y == example.W_zeta_at_secret[1]);    
    //    assert(proof.W_zeta_omega_at_secret == example.W_zeta_omega_at_secret);
    proof.W_zeta_omega_at_secret.print();
    libff::G1<ppT> W_zeta_omega_at_secret_aff(proof.W_zeta_omega_at_secret);
    W_zeta_omega_at_secret_aff.to_affine_coordinates();
    assert(W_zeta_omega_at_secret_aff.X == example.W_zeta_omega_at_secret[0]);
    assert(W_zeta_omega_at_secret_aff.Y == example.W_zeta_omega_at_secret[1]);
#endif // #ifdef DEBUG
    
    // --- VERIFIER ---
    printf("[%s:%d] Verifier...\n", __FILE__, __LINE__);
    // initialize verifier
    plonk_verifier<ppT> verifier;
    // verify proof
    bool b_valid_proof = verifier.verify_proof(proof, srs, common_input);
    assert(b_valid_proof);
    
    // end 
    printf("[%s:%d] Test OK\n", __FILE__, __LINE__);
  }

  TEST(TestPlonk, BLS12_381)
  {
    test_plonk<libff::bls12_381_pp>();
  }

} // namespace libsnark
