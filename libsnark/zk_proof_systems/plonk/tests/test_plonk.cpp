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
#include <libsnark/zk_proof_systems/plonk/common_input.hpp>
#include <libsnark/zk_proof_systems/plonk/prover.hpp>
#include <libsnark/zk_proof_systems/plonk/verifier.hpp>
//#include <libsnark/zk_proof_systems/plonk/plonk.hpp>
//#include <libsnark/zk_proof_systems/plonk/srs.hpp>

namespace libsnark
{
// Manipulate elements of a valid proof to assert that proof
// verification fails
//
// Plonk proof Pi is composed of the following elements:
//
// Pi ([a]_1, [b]_1, [c]_1, [z]_1,
//     [t_lo]_1, [t_mi]_1, [t_hi]_1,
//     \bar{a}, \bar{b}, \bar{c},
//     \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
//     [W_zeta]_1, [W_{zeta omega_roots}]_1
//     r_zeta (*))
//
template<typename ppT>
void test_verify_invalid_proof(
    const plonk_proof<ppT> valid_proof,
    const srs<ppT> srs,
    const common_preprocessed_input<ppT> common_input)
{
    // initialize verifier
    plonk_verifier<ppT> verifier;
    bool b_accept = true;

    // random element on the curve initialized to zero
    libff::G1<ppT> G1_noise = libff::G1<ppT>::zero();
    // random element in the scalar field initialized to zero
    libff::Fr<ppT> Fr_noise = libff::Fr<ppT>::zero();
    // initialize the manipulated proof to the valid one
    plonk_proof<ppT> proof = valid_proof;

    // manipulate [a]_1, [b]_1, [c]_1
    for (size_t i = 0; i < valid_proof.W_polys_blinded_at_secret_g1.size();
         ++i) {
        // re-initialize the manipulated proof
        proof = valid_proof;
        G1_noise = libff::G1<ppT>::random_element();
        proof.W_polys_blinded_at_secret_g1[i] =
            proof.W_polys_blinded_at_secret_g1[i] + G1_noise;
        b_accept = verifier.verify_proof(proof, srs, common_input);
        ASSERT_FALSE(b_accept);
    }
    // manipulate [z]_1
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.z_poly_at_secret_g1 = proof.z_poly_at_secret_g1 + G1_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
    // manipulate [t_lo]_1, [t_mi]_1, [t_hi]_1
    for (size_t i = 0; i < valid_proof.t_poly_at_secret_g1.size(); ++i) {
        // re-initialize the manipulated proof
        proof = valid_proof;
        G1_noise = libff::G1<ppT>::random_element();
        proof.t_poly_at_secret_g1[i] = proof.t_poly_at_secret_g1[i] + G1_noise;
        b_accept = verifier.verify_proof(proof, srs, common_input);
        ASSERT_FALSE(b_accept);
    }
    // manipulate \bar{a}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.a_zeta = proof.a_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{b}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.b_zeta = proof.b_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{c}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.c_zeta = proof.c_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{S_sigma1}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.S_0_zeta = proof.S_0_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{S_sigma2}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.S_1_zeta = proof.S_1_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{z_w}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.z_poly_xomega_zeta = proof.z_poly_xomega_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
    // manipulate [W_zeta]_1
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.W_zeta_at_secret = proof.W_zeta_at_secret + G1_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
    // manipulate [W_{zeta omega_roots}]_1
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.W_zeta_omega_at_secret = proof.W_zeta_omega_at_secret + G1_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
    // manipulate r_zeta
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.r_zeta = proof.r_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, common_input);
    ASSERT_FALSE(b_accept);
}

template<typename ppT> void test_plonk()
{
#ifndef DEBUG
    printf("[%s:%d] WARNING! DEBUG info disabled.\n", __FILE__, __LINE__);
    //    assert(0);
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

    printf("[%s:%d] SRS...\n", __FILE__, __LINE__);
    // --- USRS ---
    // create USRS object
    usrs<ppT> usrs;
    // compute SRS = powers of secret times G1: 1*G1, secret^1*G1,
    // secret^2*G1, ... and secret times G2: 1*G2, secret^1*G2
    usrs.derive_from_secret(secret);
    // --- SRS ---
    // circuit description. consists only of number of gates for now.
    circuit_t<ppT> circuit;
    circuit.num_gates = common_input.num_gates;
    // create SRS object
    srs<ppT> srs(circuit);
    // derive SRS from USRS and the circuit description
    srs.derive(usrs);
    // sompare SRS against reference test values
    for (int i = 0; i < (int)common_input.num_gates + 3; ++i) {
        printf("secret_power_G1[%2d] ", i);
        srs.secret_powers_g1[i].print();
        // test from generator
        libff::G1<ppT> srs_secret_powers_g1_i(srs.secret_powers_g1[i]);
        srs_secret_powers_g1_i.to_affine_coordinates();
        ASSERT_EQ(srs_secret_powers_g1_i.X, example.secret_powers_g1[i][0]);
        ASSERT_EQ(srs_secret_powers_g1_i.Y, example.secret_powers_g1[i][1]);
    }
    for (int i = 0; i < 2; ++i) {
        printf("secret_power_G2[%2d] ", i);
        srs.secret_powers_g2[i].print();
    }

    // --- PROVER ---
    printf("[%s:%d] Prover...\n", __FILE__, __LINE__);

    // initialize prover
    plonk_prover<ppT> prover;
    // compute proof
    plonk_proof<ppT> proof = prover.compute_proof(srs, common_input);
    // compare proof against test vector values (debug)
    ASSERT_EQ(proof.a_zeta, example.a_zeta);
    ASSERT_EQ(proof.b_zeta, example.b_zeta);
    ASSERT_EQ(proof.c_zeta, example.c_zeta);
    ASSERT_EQ(proof.S_0_zeta, example.S_0_zeta);
    ASSERT_EQ(proof.S_1_zeta, example.S_1_zeta);
    ASSERT_EQ(proof.z_poly_xomega_zeta, example.z_poly_xomega_zeta);
    ASSERT_EQ(proof.r_zeta, example.r_zeta);
    for (int i = 0; i < (int)NUM_HGEN; ++i) {
        printf("W_polys_at_secret_g1[%d]\n", i);
        proof.W_polys_blinded_at_secret_g1[i].print();
        libff::G1<ppT> W_polys_blinded_at_secret_g1_i(
            proof.W_polys_blinded_at_secret_g1[i]);
        W_polys_blinded_at_secret_g1_i.to_affine_coordinates();
        ASSERT_EQ(
            W_polys_blinded_at_secret_g1_i.X,
            example.W_polys_blinded_at_secret_g1[i][0]);
        ASSERT_EQ(
            W_polys_blinded_at_secret_g1_i.Y,
            example.W_polys_blinded_at_secret_g1[i][1]);
    }
    proof.z_poly_at_secret_g1.print();
    libff::G1<ppT> z_poly_at_secret_g1_aff(proof.z_poly_at_secret_g1);
    z_poly_at_secret_g1_aff.to_affine_coordinates();
    ASSERT_EQ(z_poly_at_secret_g1_aff.X, example.z_poly_at_secret_g1[0]);
    ASSERT_EQ(z_poly_at_secret_g1_aff.Y, example.z_poly_at_secret_g1[1]);
    for (int i = 0; i < (int)NUM_HGEN; ++i) {
        printf("[%s:%d] t_poly_at_secret_g1[%d]\n", __FILE__, __LINE__, i);
        proof.t_poly_at_secret_g1[i].print();
        libff::G1<ppT> t_poly_at_secret_g1_i(proof.t_poly_at_secret_g1[i]);
        t_poly_at_secret_g1_i.to_affine_coordinates();
        ASSERT_EQ(t_poly_at_secret_g1_i.X, example.t_poly_at_secret_g1[i][0]);
        ASSERT_EQ(t_poly_at_secret_g1_i.Y, example.t_poly_at_secret_g1[i][1]);
    }
    proof.W_zeta_at_secret.print();
    libff::G1<ppT> W_zeta_at_secret_aff(proof.W_zeta_at_secret);
    W_zeta_at_secret_aff.to_affine_coordinates();
    ASSERT_EQ(W_zeta_at_secret_aff.X, example.W_zeta_at_secret[0]);
    ASSERT_EQ(W_zeta_at_secret_aff.Y, example.W_zeta_at_secret[1]);
    proof.W_zeta_omega_at_secret.print();
    libff::G1<ppT> W_zeta_omega_at_secret_aff(proof.W_zeta_omega_at_secret);
    W_zeta_omega_at_secret_aff.to_affine_coordinates();
    ASSERT_EQ(W_zeta_omega_at_secret_aff.X, example.W_zeta_omega_at_secret[0]);
    ASSERT_EQ(W_zeta_omega_at_secret_aff.Y, example.W_zeta_omega_at_secret[1]);

    // --- VERIFIER ---
    printf("[%s:%d] Verifier...\n", __FILE__, __LINE__);
    // initialize verifier
    plonk_verifier<ppT> verifier;
    // verify proof
    bool b_valid_proof = verifier.verify_proof(proof, srs, common_input);
    ASSERT_TRUE(b_valid_proof);
    // assert that proof verification fails when the proof is
    // manipulated. must be executed when DEBUG is not defined to
    // disable certain assert-s that may fail before the verify
    // invalid proof test
#ifndef DEBUG
    test_verify_invalid_proof(proof, srs, common_input);
#endif // #ifndef DEBUG
#ifndef DEBUG
    printf("[%s:%d] WARNING! DEBUG info was disabled.\n", __FILE__, __LINE__);
#endif // #ifndef DEBUG
    printf("[%s:%d] Test OK\n", __FILE__, __LINE__);
    // end
}

TEST(TestPlonk, BLS12_381) { test_plonk<libff::bls12_381_pp>(); }

} // namespace libsnark
