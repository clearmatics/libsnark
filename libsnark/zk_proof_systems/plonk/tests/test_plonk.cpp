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

#define PLONK_MAX_DEGREE 245

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
    const plonk_proof<ppT> &valid_proof,
    const srs<ppT> &srs,
    transcript_hasher<ppT> &hasher)
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
        b_accept = verifier.verify_proof(proof, srs, hasher);
        ASSERT_FALSE(b_accept);
    }
    // manipulate [z]_1
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.z_poly_at_secret_g1 = proof.z_poly_at_secret_g1 + G1_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate [t_lo]_1, [t_mi]_1, [t_hi]_1
    for (size_t i = 0; i < valid_proof.t_poly_at_secret_g1.size(); ++i) {
        // re-initialize the manipulated proof
        proof = valid_proof;
        G1_noise = libff::G1<ppT>::random_element();
        proof.t_poly_at_secret_g1[i] = proof.t_poly_at_secret_g1[i] + G1_noise;
        b_accept = verifier.verify_proof(proof, srs, hasher);
        ASSERT_FALSE(b_accept);
    }
    // manipulate \bar{a}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.a_zeta = proof.a_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{b}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.b_zeta = proof.b_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{c}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.c_zeta = proof.c_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{S_sigma1}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.S_0_zeta = proof.S_0_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{S_sigma2}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.S_1_zeta = proof.S_1_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{z_w}
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.z_poly_xomega_zeta = proof.z_poly_xomega_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate [W_zeta]_1
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.W_zeta_at_secret = proof.W_zeta_at_secret + G1_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate [W_{zeta omega_roots}]_1
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.W_zeta_omega_at_secret = proof.W_zeta_omega_at_secret + G1_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate r_zeta
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.r_zeta = proof.r_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, hasher);
    ASSERT_FALSE(b_accept);
}

/// Compute or fill-in ciruit specific data from example.
/// \attention the example class is defined specifically for the BLS12-381
/// curve, so make sure we are using this curve
template<typename ppT>
circuit_t<ppT> plonk_circuit_description_from_example(
    const plonk_example example)
{
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
    // ensure that num_gates is not 0
    assert(num_gates > 0);
    // ensure num_gates is power of 2
    assert((num_gates & (num_gates - 1)) == 0);

    size_t num_qpolys = example.num_qpolys;

    // We represent the constraints q_L, q_R, q_O, q_M, q_C and the
    // witness w_L, w_R, w_O as polynomials in the roots of unity
    // e.g. f_{q_L}(omega_i) = q_L[i], 0\le{i}<8
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(num_gates);

    // compute public input (PI) polynomial
    polynomial<Field> PI_poly;
    std::vector<Field> PI_points(num_gates, Field(0));
    PI_points[PI_index] = Field(-PI_value);
    plonk_compute_public_input_polynomial(PI_points, PI_poly);

    // compute the selector polynomials (q-polynomials) from the
    // transposed gates matrix over the Lagrange basis q_poly = \sum_i
    // q[i] * L[i] where q[i] is a coefficient (a scalar Field
    // element) and L[i] is a polynomial with Field coefficients
    std::vector<polynomial<Field>> Q_polys =
        plonk_compute_selector_polynomials<Field>(
            num_gates, num_qpolys, gates_matrix_transpose);

    // omega[0] are the n roots of unity, omega[1] are omega[0]*k1,
    // omega[2] are omega[0]*k2
    std::vector<std::vector<Field>> omega_roots;
    plonk_compute_roots_of_unity_omega(num_gates, k1, k2, omega_roots);

    // H_gen contains the generators of H, k1 H and k2 H in one place
    // ie. circuit.omega_roots, circuit.omega_roots_k1 and
    // circuit.omega_roots_k2
    std::vector<Field> H_gen;
    plonk_compute_cosets_H_k1H_k2H(num_gates, k1, k2, H_gen);

    // TODO: write unit test for plonk_roots_of_unity_omega_to_subgroup_H
#ifdef DEBUG_PLONK
    printf("[%s:%d] H_gen\n", __FILE__, __LINE__);
    libff::print_vector(H_gen);
    for (int i = 0; i < (int)H_gen.size(); ++i) {
        assert(H_gen[i] == example.H_gen[i]);
    }
#endif // #ifdef DEBUG_PLONK

    // permute circuit.H_gen according to the wire permutation
    std::vector<Field> H_gen_permute =
        plonk_permute_subgroup_H<Field>(H_gen, wire_permutation, num_gates);
    // TODO: write unit test for plonk_permute_subgroup_H
#ifdef DEBUG_PLONK
    printf("[%s:%d] H_gen_permute\n", __FILE__, __LINE__);
    libff::print_vector(H_gen_permute);
    for (size_t i = 0; i < H_gen_permute.size(); ++i) {
        assert(H_gen_permute[i] == example.H_gen_permute[i]);
    }
#endif // #ifdef DEBUG_PLONK

    // compute the permutation polynomials S_sigma_1, S_sigma_2,
    // S_sigma_3 (see [GWC19], Sect. 8.1) (our indexing starts from 0)
    std::vector<polynomial<Field>> S_polys =
        plonk_compute_permutation_polynomials<Field>(H_gen_permute, num_gates);

    circuit_t<ppT> circuit(
        std::move(num_gates),
        std::move(num_qpolys),
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

template<typename ppT>
void test_plonk_compute_accumulator(
    const plonk_example &example,
    const libff::Fr<ppT> &beta,
    const libff::Fr<ppT> &gamma,
    const std::vector<libff::Fr<ppT>> &witness,
    const srs<ppT> &srs)
{
    using Field = libff::Fr<ppT>;
    // A[0] = 1; ... A[i] = computed from (i-1)
    std::vector<Field> A_vector = plonk_compute_accumulator(
        srs.num_gates, beta, gamma, witness, srs.H_gen, srs.H_gen_permute);
    polynomial<Field> A_poly(srs.num_gates);
    plonk_interpolate_polynomial_from_points<Field>(A_vector, A_poly);

    // initialize hard-coded values from example circuit
    printf("[%s:%d] A_poly\n", __FILE__, __LINE__);
    libff::print_vector(A_poly);
    ASSERT_EQ(A_poly, example.A_poly);
}

template<typename ppT>
void test_plonk_prover_round_one(
    const plonk_example &example,
    const round_zero_out_t<ppT> &round_zero_out,
    const std::vector<libff::Fr<ppT>> &witness,
    const srs<ppT> &srs,
    transcript_hasher<ppT> &hasher)
{
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    round_one_out_t<ppT> round_one_out = plonk_prover<ppT>::round_one(
        round_zero_out, blind_scalars, witness, srs, hasher);
    for (int i = 0; i < (int)NUM_HSETS; ++i) {
        printf("[%s:%d] this->W_polys[%d]\n", __FILE__, __LINE__, (int)i);
        libff::print_vector(round_one_out.W_polys[i]);
        ASSERT_EQ(round_one_out.W_polys[i], example.W_polys[i]);
    }
    for (int i = 0; i < (int)NUM_HSETS; ++i) {
        printf("[%s:%d] W_polys_blinded[%d]\n", __FILE__, __LINE__, i);
        libff::print_vector(round_one_out.W_polys_blinded[i]);
        ASSERT_EQ(round_one_out.W_polys_blinded[i], example.W_polys_blinded[i]);
    }
    printf("[%s:%d] Output from Round 1\n", __FILE__, __LINE__);
    for (int i = 0; i < (int)NUM_HSETS; ++i) {
        printf("W_polys_at_secret_g1[%d]\n", i);
        round_one_out.W_polys_blinded_at_secret_g1[i].print();
        libff::G1<ppT> W_polys_blinded_at_secret_g1_i(
            round_one_out.W_polys_blinded_at_secret_g1[i]);
        W_polys_blinded_at_secret_g1_i.to_affine_coordinates();
        ASSERT_EQ(
            W_polys_blinded_at_secret_g1_i.X,
            example.W_polys_blinded_at_secret_g1[i][0]);
        ASSERT_EQ(
            W_polys_blinded_at_secret_g1_i.Y,
            example.W_polys_blinded_at_secret_g1[i][1]);
    }
}

template<typename ppT>
void test_plonk_prover_round_two(
    const plonk_example &example,
    const libff::Fr<ppT> &beta,
    const libff::Fr<ppT> &gamma,
    const round_zero_out_t<ppT> &round_zero_out,
    const std::vector<libff::Fr<ppT>> &blind_scalars,
    const std::vector<libff::Fr<ppT>> &witness,
    const srs<ppT> &srs,
    transcript_hasher<ppT> &hasher)
{
    round_two_out_t<ppT> round_two_out = plonk_prover<ppT>::round_two(
        beta, gamma, round_zero_out, blind_scalars, witness, srs, hasher);
    printf("[%s:%d] z_poly\n", __FILE__, __LINE__);
    libff::print_vector(round_two_out.z_poly);
    ASSERT_EQ(round_two_out.z_poly, example.z_poly);
    printf("[%s:%d] Output from Round 2\n", __FILE__, __LINE__);
    printf("[%s:%d] z_poly_at_secret_g1\n", __FILE__, __LINE__);
    round_two_out.z_poly_at_secret_g1.print();
    libff::G1<ppT> z_poly_at_secret_g1_aff(round_two_out.z_poly_at_secret_g1);
    z_poly_at_secret_g1_aff.to_affine_coordinates();
    ASSERT_EQ(z_poly_at_secret_g1_aff.X, example.z_poly_at_secret_g1[0]);
    ASSERT_EQ(z_poly_at_secret_g1_aff.Y, example.z_poly_at_secret_g1[1]);
}

template<typename ppT>
void test_plonk_prover_round_three(
    const plonk_example &example,
    const libff::Fr<ppT> &alpha,
    const libff::Fr<ppT> &beta,
    const libff::Fr<ppT> &gamma,
    const round_zero_out_t<ppT> &round_zero_out,
    const round_one_out_t<ppT> &round_one_out,
    const round_two_out_t<ppT> &round_two_out,
    const srs<ppT> &srs,
    transcript_hasher<ppT> &hasher)
{
    round_three_out_t<ppT> round_three_out = plonk_prover<ppT>::round_three(
        alpha,
        beta,
        gamma,
        round_zero_out,
        round_one_out,
        round_two_out,
        srs,
        hasher);
    printf("[%s:%d] Output from Round 3\n", __FILE__, __LINE__);
    printf("[%s:%d] t_poly_long\n", __FILE__, __LINE__);
    libff::print_vector(round_three_out.t_poly_long);
    ASSERT_EQ(round_three_out.t_poly_long, example.t_poly_long);
    for (int i = 0; i < (int)NUM_HSETS; ++i) {
        printf("[%s:%d] t_poly[%d]\n", __FILE__, __LINE__, i);
        libff::print_vector(round_three_out.t_poly[i]);
        ASSERT_EQ(round_three_out.t_poly[i], example.t_poly[i]);
    }
    for (int i = 0; i < (int)NUM_HSETS; ++i) {
        printf("[%s:%d] t_poly_at_secret_g1[%d]\n", __FILE__, __LINE__, i);
        round_three_out.t_poly_at_secret_g1[i].print();
        libff::G1<ppT> t_poly_at_secret_g1_i(
            round_three_out.t_poly_at_secret_g1[i]);
        t_poly_at_secret_g1_i.to_affine_coordinates();
        ASSERT_EQ(t_poly_at_secret_g1_i.X, example.t_poly_at_secret_g1[i][0]);
        ASSERT_EQ(t_poly_at_secret_g1_i.Y, example.t_poly_at_secret_g1[i][1]);
    }
}

template<typename ppT>
void test_plonk_prover_round_four(
    const plonk_example &example,
    const libff::Fr<ppT> &zeta,
    const round_one_out_t<ppT> &round_one_out,
    const round_three_out_t<ppT> &round_three_out,
    const srs<ppT> &srs,
    transcript_hasher<ppT> &hasher)
{
    round_four_out_t<ppT> round_four_out = plonk_prover<ppT>::round_four(
        zeta, round_one_out, round_three_out, srs, hasher);
    // Prover Round 4 output check against test vectors
    printf("[%s:%d] Output from Round 4\n", __FILE__, __LINE__);
    printf("a_zeta ");
    round_four_out.a_zeta.print();
    ASSERT_EQ(round_four_out.a_zeta, example.a_zeta);
    printf("b_zeta ");
    round_four_out.b_zeta.print();
    ASSERT_EQ(round_four_out.b_zeta, example.b_zeta);
    printf("c_zeta ");
    round_four_out.c_zeta.print();
    ASSERT_EQ(round_four_out.c_zeta, example.c_zeta);
    printf("S_0_zeta ");
    round_four_out.S_0_zeta.print();
    ASSERT_EQ(round_four_out.S_0_zeta, example.S_0_zeta);
    printf("S_1_zeta ");
    round_four_out.S_1_zeta.print();
    ASSERT_EQ(round_four_out.S_1_zeta, example.S_1_zeta);
    printf("t_zeta ");
    round_four_out.t_zeta.print();
    ASSERT_EQ(round_four_out.t_zeta, example.t_zeta);
    printf("z_poly_xomega_zeta ");
    round_four_out.z_poly_xomega_zeta.print();
    ASSERT_EQ(round_four_out.z_poly_xomega_zeta, example.z_poly_xomega_zeta);
}

template<typename ppT>
void test_plonk_prover_round_five(
    const plonk_example &example,
    const libff::Fr<ppT> &alpha,
    const libff::Fr<ppT> &beta,
    const libff::Fr<ppT> &gamma,
    const libff::Fr<ppT> &zeta,
    const libff::Fr<ppT> &nu,
    const round_zero_out_t<ppT> &round_zero_out,
    const round_one_out_t<ppT> &round_one_out,
    const round_two_out_t<ppT> &round_two_out,
    const round_three_out_t<ppT> &round_three_out,
    const round_four_out_t<ppT> &round_four_out,
    const srs<ppT> &srs,
    transcript_hasher<ppT> &hasher)
{
    round_five_out_t<ppT> round_five_out = plonk_prover<ppT>::round_five(
        alpha,
        beta,
        gamma,
        zeta,
        nu,
        round_zero_out,
        round_one_out,
        round_two_out,
        round_three_out,
        round_four_out,
        srs,
        hasher);

    printf("[%s:%d] Outputs from Prover round 5\n", __FILE__, __LINE__);
    printf("r_zeta ");
    round_five_out.r_zeta.print();
    ASSERT_EQ(round_five_out.r_zeta, example.r_zeta);
    printf("[%s:%d] W_zeta_at_secret \n", __FILE__, __LINE__);
    round_five_out.W_zeta_at_secret.print();
    libff::G1<ppT> W_zeta_at_secret_aff(round_five_out.W_zeta_at_secret);
    W_zeta_at_secret_aff.to_affine_coordinates();
    ASSERT_EQ(W_zeta_at_secret_aff.X, example.W_zeta_at_secret[0]);
    ASSERT_EQ(W_zeta_at_secret_aff.Y, example.W_zeta_at_secret[1]);
    printf("[%s:%d] W_zeta_omega_at_secret \n", __FILE__, __LINE__);
    round_five_out.W_zeta_omega_at_secret.print();
    libff::G1<ppT> W_zeta_omega_at_secret_aff(
        round_five_out.W_zeta_omega_at_secret);
    W_zeta_omega_at_secret_aff.to_affine_coordinates();
    ASSERT_EQ(W_zeta_omega_at_secret_aff.X, example.W_zeta_omega_at_secret[0]);
    ASSERT_EQ(W_zeta_omega_at_secret_aff.Y, example.W_zeta_omega_at_secret[1]);
}

/// \attention the example class is defined specifically for the BLS12-381
/// curve, so make sure we are using this curve
template<typename ppT> void test_plonk_prover_rounds()
{
    using Field = libff::Fr<ppT>;

    ppT::init_public_params();

    // load test vector values from example circuit
    plonk_example example;
    // random hidden element secret (toxic waste)
    Field secret = example.secret;
    // example witness
    std::vector<Field> witness = example.witness;
    // example circuit
    circuit_t<ppT> circuit =
        plonk_circuit_description_from_example<ppT>(example);
    // hard-coded values for the "random" blinding constants from
    // example circuit
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    // prepare srs
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(usrs, circuit);

    // initialize hasher
    std::vector<uint8_t> buffer;
    transcript_hasher<ppT> hasher(buffer);

    // Prover Round 0 (initialization)
    round_zero_out_t<ppT> round_zero_out = plonk_prover<ppT>::round_zero(srs);

    // --- Unit test Prover Round 1 ---
    // reset buffer at the start of the round (needed for testing only)
    printf("[%s:%d] Unit test Prover Round 1...\n", __FILE__, __LINE__);
    test_plonk_prover_round_one<ppT>(
        example, round_zero_out, witness, srs, hasher);

    // --- Unit test Prover Round 2 ---
    // reset buffer at the start of the round (needed for testing only)
    printf("[%s:%d] Unit test Prover Round 2...\n", __FILE__, __LINE__);
    round_one_out_t<ppT> round_one_out = plonk_prover<ppT>::round_one(
        round_zero_out, blind_scalars, witness, srs, hasher);
    // clear hash buffer
    hasher.buffer_clear();
    // add outputs from Round 1 to the hash buffer
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[a]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[b]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[c]);
    const libff::Fr<ppT> beta = hasher.get_hash();
    hasher.add_element(libff::Fr<ppT>::one());
    const libff::Fr<ppT> gamma = hasher.get_hash();
    test_plonk_prover_round_two<ppT>(
        example,
        beta,
        gamma,
        round_zero_out,
        blind_scalars,
        witness,
        srs,
        hasher);

    // Unit test plonk_compute_accumulator
    test_plonk_compute_accumulator<ppT>(example, beta, gamma, witness, srs);

    // --- Unit test Prover Round 3 ---
    // reset buffer at the start of the round (needed for testing only)
    printf("[%s:%d] Prover Round 3...\n", __FILE__, __LINE__);
    round_two_out_t<ppT> round_two_out = plonk_prover<ppT>::round_two(
        beta, gamma, round_zero_out, blind_scalars, witness, srs, hasher);
    // clear hash buffer
    hasher.buffer_clear();
    // add outputs from Round 1 to the hash buffer
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[a]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[b]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[c]);
    hasher.add_element(libff::Fr<ppT>::one());
    // add outputs from Round 2 to the hash buffer
    hasher.add_element(round_two_out.z_poly_at_secret_g1);
    const libff::Fr<ppT> alpha = hasher.get_hash();
    test_plonk_prover_round_three<ppT>(
        example,
        alpha,
        beta,
        gamma,
        round_zero_out,
        round_one_out,
        round_two_out,
        srs,
        hasher);

    // --- Unit test Prover Round 4 ---
    printf("[%s:%d] Prover Round 4...\n", __FILE__, __LINE__);
    round_three_out_t<ppT> round_three_out = plonk_prover<ppT>::round_three(
        alpha,
        beta,
        gamma,
        round_zero_out,
        round_one_out,
        round_two_out,
        srs,
        hasher);
    // clear hash buffer
    hasher.buffer_clear();
    // add outputs from Round 1 to the hash buffer
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[a]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[b]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[c]);
    hasher.add_element(libff::Fr<ppT>::one());
    // add outputs from Round 2 to the hash buffer
    hasher.add_element(round_two_out.z_poly_at_secret_g1);
    // add outputs from Round 3 to the hash buffer
    hasher.add_element(round_three_out.t_poly_at_secret_g1[lo]);
    hasher.add_element(round_three_out.t_poly_at_secret_g1[mid]);
    hasher.add_element(round_three_out.t_poly_at_secret_g1[hi]);
    const libff::Fr<ppT> zeta = hasher.get_hash();
    test_plonk_prover_round_four<ppT>(
        example, zeta, round_one_out, round_three_out, srs, hasher);

    // --- Unit test Prover Round 5 ---
    printf("[%s:%d] Unit test Prover Round 5...\n", __FILE__, __LINE__);
    round_four_out_t<ppT> round_four_out = plonk_prover<ppT>::round_four(
        zeta, round_one_out, round_three_out, srs, hasher);
    // clear hash buffer
    hasher.buffer_clear();
    // add outputs from Round 1 to the hash buffer
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[a]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[b]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[c]);
    hasher.add_element(libff::Fr<ppT>::one());
    // add outputs from Round 2 to the hash buffer
    hasher.add_element(round_two_out.z_poly_at_secret_g1);
    // add outputs from Round 3 to the hash buffer
    hasher.add_element(round_three_out.t_poly_at_secret_g1[lo]);
    hasher.add_element(round_three_out.t_poly_at_secret_g1[mid]);
    hasher.add_element(round_three_out.t_poly_at_secret_g1[hi]);
    // add outputs from Round 4 to the hash buffer
    hasher.add_element(round_four_out.a_zeta);
    hasher.add_element(round_four_out.b_zeta);
    hasher.add_element(round_four_out.c_zeta);
    hasher.add_element(round_four_out.S_0_zeta);
    hasher.add_element(round_four_out.S_1_zeta);
    hasher.add_element(round_four_out.z_poly_xomega_zeta);
    const libff::Fr<ppT> nu = hasher.get_hash();
    test_plonk_prover_round_five<ppT>(
        example,
        alpha,
        beta,
        gamma,
        zeta,
        nu,
        round_zero_out,
        round_one_out,
        round_two_out,
        round_three_out,
        round_four_out,
        srs,
        hasher);
}

/// \attention the example class is defined specifically for the BLS12-381
/// curve, so make sure we are using this curve
template<typename ppT> void test_plonk_srs()
{
    using Field = libff::Fr<ppT>;

    ppT::init_public_params();
    // load test vector values from example circuit
    plonk_example example;
    // random hidden element secret (toxic waste)
    Field secret = example.secret;
    // example circuit
    circuit_t<ppT> circuit =
        plonk_circuit_description_from_example<ppT>(example);
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    // --- USRS ---
    // compute SRS = powers of secret times G1: 1*G1, secret^1*G1,
    // secret^2*G1, ... and secret times G2: 1*G2, secret^1*G2
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    // --- SRS ---
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(usrs, circuit);
    // compare SRS against reference test values
    printf("[%s:%d] secret ", __FILE__, __LINE__);
    secret.print();
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
}

/// \attention the example class is defined specifically for the BLS12-381
/// curve, so make sure we are using this curve
template<typename ppT> void test_plonk_prover()
{
    using Field = libff::Fr<ppT>;

    ppT::init_public_params();
    // load test vector values from example circuit
    plonk_example example;
    // random hidden element secret (toxic waste)
    Field secret = example.secret;
    // example witness
    std::vector<Field> witness = example.witness;
    // example circuit
    circuit_t<ppT> circuit =
        plonk_circuit_description_from_example<ppT>(example);
    // hard-coded values for the "random" blinding constants from
    // example circuit
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    // prepare srs
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(usrs, circuit);

    // initialize hasher
    std::vector<uint8_t> buffer;
    transcript_hasher<ppT> hasher(buffer);

    // initialize prover
    plonk_prover<ppT> prover;
    // compute proof
    plonk_proof<ppT> proof =
        prover.compute_proof(srs, witness, blind_scalars, hasher);
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
}

/// \attention the example class is defined specifically for the BLS12-381
/// curve, so make sure we are using this curve
template<typename ppT>
void test_plonk_verifier_preprocessed_input(
    const plonk_example &example, const srs<ppT> &srs)
{
    // compute verifier preprocessed input
    const verifier_preprocessed_input_t<ppT> preprocessed_input =
        plonk_verifier<ppT>::preprocessed_input(srs);

    for (int i = 0; i < (int)srs.Q_polys.size(); ++i) {
        printf("srs.Q_polys_at_secret_G1[%d] \n", i);
        preprocessed_input.Q_polys_at_secret_g1[i].print();
        libff::G1<ppT> Q_poly_at_secret_g1_i(
            preprocessed_input.Q_polys_at_secret_g1[i]);
        Q_poly_at_secret_g1_i.to_affine_coordinates();
        ASSERT_EQ(Q_poly_at_secret_g1_i.X, example.Q_polys_at_secret_g1[i][0]);
        ASSERT_EQ(Q_poly_at_secret_g1_i.Y, example.Q_polys_at_secret_g1[i][1]);
    }
    for (int i = 0; i < (int)srs.S_polys.size(); ++i) {
        printf("S_polys_at_secret_G1[%d] \n", i);
        preprocessed_input.S_polys_at_secret_g1[i].print();
        libff::G1<ppT> S_poly_at_secret_g1_i(
            preprocessed_input.S_polys_at_secret_g1[i]);
        S_poly_at_secret_g1_i.to_affine_coordinates();
        ASSERT_EQ(S_poly_at_secret_g1_i.X, example.S_polys_at_secret_g1[i][0]);
        ASSERT_EQ(S_poly_at_secret_g1_i.Y, example.S_polys_at_secret_g1[i][1]);
    }
}

template<typename ppT>
void test_plonk_verifier_step_five(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const srs<ppT> &srs)
{
    const step_five_out_t<ppT> step_five_out =
        plonk_verifier<ppT>::step_five(step_four_out, srs);
    printf("[%s:%d] zh_zeta ", __FILE__, __LINE__);
    step_five_out.zh_zeta.print();
    ASSERT_EQ(step_five_out.zh_zeta, example.zh_zeta);
}

template<typename ppT>
void test_plonk_verifier_step_six(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const srs<ppT> &srs)
{
    const step_six_out_t<ppT> step_six_out =
        plonk_verifier<ppT>::step_six(step_four_out, srs);
    printf("L_0_zeta ");
    step_six_out.L_0_zeta.print();
    ASSERT_EQ(step_six_out.L_0_zeta, example.L_0_zeta);
}

template<typename ppT>
void test_plonk_verifier_step_seven(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const srs<ppT> &srs)
{
    const step_seven_out_t<ppT> step_seven_out =
        plonk_verifier<ppT>::step_seven(step_four_out, srs);
    printf("PI_zeta ");
    step_seven_out.PI_zeta.print();
    ASSERT_EQ(step_seven_out.PI_zeta, example.PI_zeta);
}

template<typename ppT>
void test_plonk_verifier_step_eight(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const step_five_out_t<ppT> &step_five_out,
    const step_six_out_t<ppT> &step_six_out,
    const step_seven_out_t<ppT> &step_seven_out,
    const plonk_proof<ppT> &proof)
{
    const step_eight_out_t<ppT> step_eight_out =
        plonk_verifier<ppT>::step_eight(
            step_four_out, step_five_out, step_six_out, step_seven_out, proof);
    ASSERT_EQ(step_eight_out.r_prime_zeta, example.r_prime_zeta);
}

template<typename ppT>
void test_plonk_verifier_step_nine(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const step_six_out_t<ppT> &step_six_out,
    const plonk_proof<ppT> &proof,
    const verifier_preprocessed_input_t<ppT> &preprocessed_input,
    const srs<ppT> &srs)
{
    step_nine_out_t<ppT> step_nine_out = plonk_verifier<ppT>::step_nine(
        step_four_out, step_six_out, proof, preprocessed_input, srs);
    step_nine_out.D1.print();
    libff::G1<ppT> D1_aff(step_nine_out.D1);
    D1_aff.to_affine_coordinates();
    ASSERT_EQ(D1_aff.X, example.D1[0]);
    ASSERT_EQ(D1_aff.Y, example.D1[1]);
}

template<typename ppT>
void test_plonk_verifier_step_ten(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const step_nine_out_t<ppT> &step_nine_out,
    const plonk_proof<ppT> &proof,
    const verifier_preprocessed_input_t<ppT> &preprocessed_input,
    const srs<ppT> &srs)
{
    step_ten_out_t<ppT> step_ten_out = plonk_verifier<ppT>::step_ten(
        step_four_out, step_nine_out, proof, preprocessed_input, srs);
    printf("[%s:%d] F1\n", __FILE__, __LINE__);
    step_ten_out.F1.print();
    libff::G1<ppT> F1_aff(step_ten_out.F1);
    F1_aff.to_affine_coordinates();
    ASSERT_EQ(F1_aff.X, example.F1[0]);
    ASSERT_EQ(F1_aff.Y, example.F1[1]);
}

template<typename ppT>
void test_plonk_verifier_step_eleven(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const step_eight_out_t<ppT> &step_eight_out,
    const plonk_proof<ppT> &proof)
{
    const step_eleven_out_t<ppT> step_eleven_out =
        plonk_verifier<ppT>::step_eleven(step_four_out, step_eight_out, proof);
    printf("[%s:%d] E1\n", __FILE__, __LINE__);
    step_eleven_out.E1.print();
    libff::G1<ppT> E1_aff(step_eleven_out.E1);
    E1_aff.to_affine_coordinates();
    ASSERT_EQ(E1_aff.X, example.E1[0]);
    ASSERT_EQ(E1_aff.Y, example.E1[1]);
}

/// test verifier step 12
/// \attention the example class is defined specifically for the BLS12-381
/// curve, so make sure we are using this curve
template<typename ppT>
void test_plonk_verifier_pairing(
    const plonk_example &example,
    const libff::Fr<ppT> &zeta, // step 4
    const libff::Fr<ppT> &u,    // step 4
    const libff::G1<ppT> &F1,   // step 10
    const libff::G1<ppT> &E1,   // step 11
    const plonk_proof<ppT> &proof,
    const srs<ppT> &srs)
{
    using Field = libff::Fr<ppT>;
    std::vector<libff::G1<ppT>> curve_points_lhs{
        proof.W_zeta_at_secret, proof.W_zeta_omega_at_secret};
    std::vector<libff::Fr<ppT>> scalar_elements_lhs{Field(1), u};
    libff::G1<ppT> pairing_first_lhs =
        plonk_multi_exp_G1<ppT>(curve_points_lhs, scalar_elements_lhs);
    libff::G2<ppT> pairing_second_lhs = srs.secret_powers_g2[1];
    std::vector<libff::G1<ppT>> curve_points_rhs{
        proof.W_zeta_at_secret, proof.W_zeta_omega_at_secret, F1, E1};
    std::vector<libff::Fr<ppT>> scalar_elements_rhs{
        // Warning! raise to the power of -1 to check e() * e()^-1 = 1
        Field(-1) * zeta,
        Field(-1) * u * zeta * srs.omega_roots[base][1],
        Field(-1) * Field(1),
        Field(-1) * Field(-1)};

    libff::G1<ppT> pairing_first_rhs =
        plonk_multi_exp_G1<ppT>(curve_points_rhs, scalar_elements_rhs);
    libff::G2<ppT> pairing_second_rhs = srs.secret_powers_g2[0];

    const libff::G1_precomp<ppT> _A = ppT::precompute_G1(pairing_first_lhs);
    const libff::G2_precomp<ppT> _B = ppT::precompute_G2(pairing_second_lhs);
    const libff::G1_precomp<ppT> _C = ppT::precompute_G1(pairing_first_rhs);
    const libff::G2_precomp<ppT> _D = ppT::precompute_G2(pairing_second_rhs);
    const libff::Fqk<ppT> miller_result =
        ppT::double_miller_loop(_A, _B, _C, _D);
    const libff::GT<ppT> result = ppT::final_exponentiation(miller_result);
    bool b_accept = (result == libff::GT<ppT>::one());

    printf("[%s:%d] pairing_first_lhs\n", __FILE__, __LINE__);
    pairing_first_lhs.print();
    libff::G1<ppT> pairing_first_lhs_aff(pairing_first_lhs);
    pairing_first_lhs_aff.to_affine_coordinates();
    ASSERT_EQ(pairing_first_lhs_aff.X, example.pairing_first_lhs[0]);
    ASSERT_EQ(pairing_first_lhs_aff.Y, example.pairing_first_lhs[1]);

    printf("[%s:%d] pairing_first_rhs\n", __FILE__, __LINE__);
    pairing_first_rhs.print();
    libff::G1<ppT> pairing_first_rhs_aff(pairing_first_rhs);
    pairing_first_rhs_aff.to_affine_coordinates();
    ASSERT_EQ(pairing_first_rhs_aff.X, example.pairing_first_rhs[0]);
    ASSERT_EQ(pairing_first_rhs_aff.Y, example.pairing_first_rhs[1]);

    ASSERT_TRUE(b_accept);
}

/// \attention the example class is defined specifically for the BLS12-381
/// curve, so make sure we are using this curve
template<typename ppT> void test_plonk_verifier_steps()
{
    using Field = libff::Fr<ppT>;

    ppT::init_public_params();
    // load test vector values from example circuit
    plonk_example example;
    // random hidden element secret (toxic waste)
    Field secret = example.secret;
    // example witness
    std::vector<Field> witness = example.witness;
    // example circuit
    circuit_t<ppT> circuit =
        plonk_circuit_description_from_example<ppT>(example);
    // hard-coded values for the "random" blinding constants from
    // example circuit
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    // prepare srs
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(usrs, circuit);

    // initialize hasher
    std::vector<uint8_t> buffer;
    transcript_hasher<ppT> hasher(buffer);

    // initialize prover
    plonk_prover<ppT> prover;
    // compute proof
    plonk_proof<ppT> proof =
        prover.compute_proof(srs, witness, blind_scalars, hasher);

    // Unit test verifier preprocessed input
    test_plonk_verifier_preprocessed_input(example, srs);

    // unit test verifier step 5
    const step_four_out_t<ppT> step_four_out =
        plonk_verifier<ppT>::step_four(proof, hasher);
    test_plonk_verifier_step_five(example, step_four_out, srs);

    // unit test verifier step 6
    test_plonk_verifier_step_six(example, step_four_out, srs);

    // unit test verifier step 7
    test_plonk_verifier_step_seven(example, step_four_out, srs);

    // unit test verifier step 8
    const step_five_out_t<ppT> step_five_out =
        plonk_verifier<ppT>::step_five(step_four_out, srs);
    const step_six_out_t<ppT> step_six_out =
        plonk_verifier<ppT>::step_six(step_four_out, srs);
    const step_seven_out_t<ppT> step_seven_out =
        plonk_verifier<ppT>::step_seven(step_four_out, srs);
    test_plonk_verifier_step_eight(
        example,
        step_four_out,
        step_five_out,
        step_six_out,
        step_seven_out,
        proof);

    // unit test verifier step 9
    const verifier_preprocessed_input_t<ppT> preprocessed_input =
        plonk_verifier<ppT>::preprocessed_input(srs);
    test_plonk_verifier_step_nine(
        example, step_four_out, step_six_out, proof, preprocessed_input, srs);

    // unit test verifier step 10
    step_nine_out_t<ppT> step_nine_out = plonk_verifier<ppT>::step_nine(
        step_four_out, step_six_out, proof, preprocessed_input, srs);
    test_plonk_verifier_step_ten(
        example, step_four_out, step_nine_out, proof, preprocessed_input, srs);

    // unit test verifier step 11
    const step_eight_out_t<ppT> step_eight_out =
        plonk_verifier<ppT>::step_eight(
            step_four_out, step_five_out, step_six_out, step_seven_out, proof);
    test_plonk_verifier_step_eleven(
        example, step_four_out, step_eight_out, proof);

    // unit test verifier pairing (step 12)
    step_ten_out_t<ppT> step_ten_out = plonk_verifier<ppT>::step_ten(
        step_four_out, step_nine_out, proof, preprocessed_input, srs);
    const step_eleven_out_t<ppT> step_eleven_out =
        plonk_verifier<ppT>::step_eleven(step_four_out, step_eight_out, proof);
    test_plonk_verifier_pairing(
        example,
        step_four_out.zeta,
        step_four_out.u,
        step_ten_out.F1,
        step_eleven_out.E1,
        proof,
        srs);
}

/// \attention the example class is defined specifically for the BLS12-381
/// curve, so make sure we are using this curve
template<typename ppT> void test_plonk_verifier()
{
    using Field = libff::Fr<ppT>;

    ppT::init_public_params();
    // load test vector values from example circuit
    plonk_example example;
    // random hidden element secret (toxic waste)
    Field secret = example.secret;
    // example witness
    std::vector<Field> witness = example.witness;
    // example circuit
    circuit_t<ppT> circuit =
        plonk_circuit_description_from_example<ppT>(example);
    // hard-coded values for the "random" blinding constants from
    // example circuit
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    // prepare srs
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(usrs, circuit);

    // initialize hasher
    std::vector<uint8_t> buffer;
    transcript_hasher<ppT> hasher(buffer);

    // initialize prover
    plonk_prover<ppT> prover;
    // compute proof
    plonk_proof<ppT> proof =
        prover.compute_proof(srs, witness, blind_scalars, hasher);

    // initialize verifier
    plonk_verifier<ppT> verifier;
    // verify proof
    bool b_valid_proof = verifier.verify_proof(proof, srs, hasher);
    ASSERT_TRUE(b_valid_proof);
    // assert that proof verification fails when the proof is
    // manipulated
    test_verify_invalid_proof(proof, srs, hasher);
}

TEST(TestPlonk, BLS12_381)
{
    test_plonk_srs<libff::bls12_381_pp>();
    test_plonk_prover_rounds<libff::bls12_381_pp>();
    test_plonk_prover<libff::bls12_381_pp>();
    test_plonk_verifier_steps<libff::bls12_381_pp>();
    test_plonk_verifier<libff::bls12_381_pp>();
}

} // namespace libsnark
