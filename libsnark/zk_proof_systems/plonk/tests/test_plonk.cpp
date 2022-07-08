/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libsnark/zk_proof_systems/plonk/circuit.hpp"
#include "libsnark/zk_proof_systems/plonk/prover.hpp"
#include "libsnark/zk_proof_systems/plonk/verifier.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <gtest/gtest.h>

/// Test program that exercises the Plonk protocol (first setup, then
/// prover, then verifier) on a synthetic R1CS instance.

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
template<typename ppT>
void test_verify_invalid_proof(
    const plonk_proof<ppT> &valid_proof, const srs<ppT> &srs)
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
        b_accept = verifier.verify_proof(proof, srs);
        ASSERT_FALSE(b_accept);
    }
    // manipulate [z]_1
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.z_poly_at_secret_g1 = proof.z_poly_at_secret_g1 + G1_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
    // manipulate [t_lo]_1, [t_mi]_1, [t_hi]_1
    for (size_t i = 0; i < valid_proof.t_poly_at_secret_g1.size(); ++i) {
        // re-initialize the manipulated proof
        proof = valid_proof;
        G1_noise = libff::G1<ppT>::random_element();
        proof.t_poly_at_secret_g1[i] = proof.t_poly_at_secret_g1[i] + G1_noise;
        b_accept = verifier.verify_proof(proof, srs);
        ASSERT_FALSE(b_accept);
    }
    // manipulate \bar{a}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.a_zeta = proof.a_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{b}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.b_zeta = proof.b_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{c}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.c_zeta = proof.c_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{S_sigma1}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.S_0_zeta = proof.S_0_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{S_sigma2}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.S_1_zeta = proof.S_1_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{z_w}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.z_poly_xomega_zeta = proof.z_poly_xomega_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
    // manipulate [W_zeta]_1
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.W_zeta_at_secret = proof.W_zeta_at_secret + G1_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
    // manipulate [W_{zeta omega_roots}]_1
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.W_zeta_omega_at_secret = proof.W_zeta_omega_at_secret + G1_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
    // manipulate r_zeta
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.r_zeta = proof.r_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs);
    ASSERT_FALSE(b_accept);
}

/// Compute or fill-in ciruit specific data from example.
template<typename ppT>
circuit_t<ppT> plonk_circuit_description_from_example(
    const plonk_example example)
{
    // the example class is defined specifically for the BLS12-381
    // curve, so make sure we are using this curve TODO: remove when
    // the implementation is stable and tested
    try {
        plonk_exception_assert_curve_bls12_381<ppT>();
    } catch (const std::domain_error &e) {
        std::cout << "Error: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }

    using Field = libff::Fr<ppT>;

    // public input (PI)
    Field PI_value = example.public_input;
    // index of the row of the PI in the non-transposed gates_matrix
    int PI_index = example.public_input_index;
    // Transposed gates matrix: each row is a q-vector
    std::vector<std::vector<Field>> gates_matrix_transpose =
        example.gates_matrix_transpose;
    // wire permutation
    std::vector<size_t> wire_permutation = example.wire_permutation;
    // Generate domains on which to evaluate the witness
    // polynomials. k1,k2 can be random, but we fix them for debug to
    // match against the test vector values
    libff::Fr<ppT> k1 = example.k1;
    libff::Fr<ppT> k2 = example.k2;
#ifdef DEBUG_PLONK
    printf("[%s:%d] k1 ", __FILE__, __LINE__);
    k1.print();
    printf("[%s:%d] k2 ", __FILE__, __LINE__);
    k2.print();
#endif // #ifdef DEBUG_PLONK

    size_t num_gates = example.num_gates;
    // TODO: throw exception
#ifdef DEBUG_PLONK
    // ensure that num_gates is not 0
    assert(num_gates);
    // ensure num_gates is power of 2
    assert((num_gates & (num_gates - 1)) == 0);
#endif // #ifdef DEBUG_PLONK

    size_t num_qpolys = example.num_qpolys;

    // We represent the constraints q_L, q_R, q_O, q_M, q_C and the
    // witness w_L, w_R, w_O as polynomials in the roots of unity
    // e.g. f_{q_L}(omega_i) = q_L[i], 0\le{i}<8
    // compute Lagrange basis
    std::vector<polynomial<Field>> L_basis;
    L_basis.resize(num_gates, polynomial<Field>(num_gates));
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(num_gates);
    plonk_compute_lagrange_basis<Field>(num_gates, L_basis);

    // compute public input (PI) polynomial
    polynomial<Field> PI_poly;
    std::vector<Field> PI_points(num_gates, Field(0));
    PI_points[PI_index] = Field(-PI_value);
    plonk_compute_public_input_polynomial(PI_points, PI_poly);

    // compute the selector polynomials (q-polynomials) from the
    // transposed gates matrix over the Lagrange basis q_poly = \sum_i
    // q[i] * L[i] where q[i] is a coefficient (a scalar Field
    // element) and L[i] is a polynomial with Field coefficients
    std::vector<polynomial<Field>> Q_polys;
    Q_polys.resize(num_qpolys, polynomial<Field>(num_gates));
    plonk_compute_selector_polynomials<Field>(gates_matrix_transpose, Q_polys);

    // number of generators for H, Hk1, Hk2
    int num_hgen = NUM_HSETS;

    // omega[0] are the n roots of unity, omega[1] are omega[0]*k1,
    // omega[2] are omega[0]*k2
    std::vector<std::vector<Field>> omega_roots;
    plonk_compute_roots_of_unity_omega(num_gates, k1, k2, omega_roots);
    // H_gen contains the generators of H, k1 H and k2 H in one place
    // ie. circuit.omega_roots, circuit.omega_roots_k1 and
    // circuit.omega_roots_k2
    std::vector<Field> H_gen;
    plonk_roots_of_unity_omega_to_subgroup_H(omega_roots, H_gen);

    // TODO: write unit test for plonk_roots_of_unity_omega_to_subgroup_H
#ifdef DEBUG_PLONK
    printf("[%s:%d] H_gen\n", __FILE__, __LINE__);
    print_vector(H_gen);
    for (int i = 0; i < (int)H_gen.size(); ++i) {
        assert(H_gen[i] == example.H_gen[i]);
    }
#endif // #ifdef DEBUG_PLONK

    // permute circuit.H_gen according to the wire permutation
    std::vector<Field> H_gen_permute;
    H_gen_permute.resize(num_hgen * num_gates, Field(0));
    plonk_permute_subgroup_H<Field>(H_gen, wire_permutation, H_gen_permute);
    // TODO: write unit test for plonk_permute_subgroup_H
#ifdef DEBUG_PLONK
    printf("[%s:%d] H_gen_permute\n", __FILE__, __LINE__);
    print_vector(H_gen_permute);
    for (size_t i = 0; i < H_gen_permute.size(); ++i) {
        assert(H_gen_permute[i] == example.H_gen_permute[i]);
    }
#endif // #ifdef DEBUG_PLONK

    // compute the permutation polynomials S_sigma_1, S_sigma_2,
    // S_sigma_3 (see [GWC19], Sect. 8.1) (our indexing starts from 0)
    std::vector<polynomial<Field>> S_polys;
    S_polys.resize(num_hgen, polynomial<Field>(num_gates));
    plonk_compute_permutation_polynomials<Field>(
        H_gen_permute, num_gates, S_polys);

    circuit_t<ppT> circuit(
        std::move(num_gates),
        std::move(num_qpolys),
        std::move(L_basis),
        std::move(PI_poly),
        std::move(Q_polys),
        std::move(S_polys),
        std::move(omega_roots),
        std::move(H_gen),
        std::move(H_gen_permute),
        std::move(k1),
        std::move(k2));
    return circuit;
}

template<typename ppT> void test_plonk()
{
#ifndef DEBUG_PLONK
    printf("[%s:%d] WARNING! DEBUG_PLONK info disabled.\n", __FILE__, __LINE__);
    //    assert(0);
#endif // #ifndef DEBUG_PLONK

    try {
        plonk_exception_assert_curve_bls12_381<ppT>();
    } catch (const std::domain_error &e) {
        std::cout << "Error: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }

    // Execute all tests for the given curve.
    ppT::init_public_params();

    using Field = libff::Fr<ppT>;

    // --- TEST VECTORS ---
    // the example class is defined specifically for the BLS12-381
    // curve, so make sure we are using this curve
    try {
        plonk_exception_assert_curve_bls12_381<ppT>();
    } catch (const std::domain_error &e) {
        std::cout << "Error: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    // load test vector values from example circuit
    plonk_example example;

    // --- SETUP ---
    printf("[%s:%d] Setup...\n", __FILE__, __LINE__);
    // random hidden element secret (toxic waste). we fix it to a
    // constant in order to match against the test vectors
    Field secret = example.secret;
#ifdef DEBUG_PLONK
    printf("[%s:%d] secret ", __FILE__, __LINE__);
    secret.print();
#endif // #ifdef DEBUG_PLONK

    printf("[%s:%d] SRS...\n", __FILE__, __LINE__);
    // --- USRS ---
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = 254;
    // compute SRS = powers of secret times G1: 1*G1, secret^1*G1,
    // secret^2*G1, ... and secret times G2: 1*G2, secret^1*G2
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    // --- circuit ---
    circuit_t<ppT> circuit =
        plonk_circuit_description_from_example<ppT>(example);
    // --- SRS ---
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(usrs, circuit);
    // sompare SRS against reference test values
    for (int i = 0; i < (int)srs.num_gates + 3; ++i) {
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
    plonk_proof<ppT> proof = prover.compute_proof(srs);
    // compare proof against test vector values (debug)
    ASSERT_EQ(proof.a_zeta, example.a_zeta);
    ASSERT_EQ(proof.b_zeta, example.b_zeta);
    ASSERT_EQ(proof.c_zeta, example.c_zeta);
    ASSERT_EQ(proof.S_0_zeta, example.S_0_zeta);
    ASSERT_EQ(proof.S_1_zeta, example.S_1_zeta);
    ASSERT_EQ(proof.z_poly_xomega_zeta, example.z_poly_xomega_zeta);
    ASSERT_EQ(proof.r_zeta, example.r_zeta);
    for (int i = 0; i < (int)NUM_HSETS; ++i) {
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
    for (int i = 0; i < (int)NUM_HSETS; ++i) {
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
    bool b_valid_proof = verifier.verify_proof(proof, srs);
    ASSERT_TRUE(b_valid_proof);
    // assert that proof verification fails when the proof is
    // manipulated. must be executed when DEBUG_PLONK is not defined to
    // disable certain assert-s that may fail before the verify
    // invalid proof test
#ifndef DEBUG_PLONK
    test_verify_invalid_proof(proof, srs);
#endif // #ifndef DEBUG_PLONK
#ifndef DEBUG_PLONK
    printf(
        "[%s:%d] WARNING! DEBUG_PLONK info was disabled.\n",
        __FILE__,
        __LINE__);
#endif // #ifndef DEBUG_PLONK
    printf("[%s:%d] Test OK\n", __FILE__, __LINE__);
    // end
}

TEST(TestPlonk, BLS12_381) { test_plonk<libff::bls12_381_pp>(); }

} // namespace libsnark
