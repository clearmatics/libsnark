/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libsnark/zk_proof_systems/plonk/prover.hpp"
#include "libsnark/zk_proof_systems/plonk/tests/bls12_381_test_vector_transcript_hasher.hpp"
#include "libsnark/zk_proof_systems/plonk/verifier.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <gtest/gtest.h>
#include <libff/algebra/curves/public_params.hpp>

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
template<typename ppT, class transcript_hasher>
void test_verify_invalid_proof(
    const plonk_proof<ppT> &valid_proof,
    const srs<ppT> &srs,
    const std::vector<libff::Fr<ppT>> &PI_value_list,
    transcript_hasher &hasher)
{
    // initialize verifier
    plonk_verifier<ppT, transcript_hasher> verifier;
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
        hasher.buffer_clear();
        proof = valid_proof;
        G1_noise = libff::G1<ppT>::random_element();
        proof.W_polys_blinded_at_secret_g1[i] =
            proof.W_polys_blinded_at_secret_g1[i] + G1_noise;
        b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
        ASSERT_FALSE(b_accept);
    }
    // manipulate [z]_1
    hasher.buffer_clear();
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.z_poly_at_secret_g1 = proof.z_poly_at_secret_g1 + G1_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate [t_lo]_1, [t_mi]_1, [t_hi]_1
    for (size_t i = 0; i < valid_proof.t_poly_at_secret_g1.size(); ++i) {
        // re-initialize the manipulated proof
        hasher.buffer_clear();
        proof = valid_proof;
        G1_noise = libff::G1<ppT>::random_element();
        proof.t_poly_at_secret_g1[i] = proof.t_poly_at_secret_g1[i] + G1_noise;
        b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
        ASSERT_FALSE(b_accept);
    }
    // manipulate \bar{a}
    hasher.buffer_clear();
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.a_zeta = proof.a_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{b}
    hasher.buffer_clear();
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.b_zeta = proof.b_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{c}
    hasher.buffer_clear();
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.c_zeta = proof.c_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{S_sigma1}
    hasher.buffer_clear();
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.S_0_zeta = proof.S_0_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{S_sigma2}
    hasher.buffer_clear();
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.S_1_zeta = proof.S_1_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate \bar{z_w}
    hasher.buffer_clear();
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.z_poly_xomega_zeta = proof.z_poly_xomega_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate [W_zeta]_1
    hasher.buffer_clear();
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.W_zeta_at_secret = proof.W_zeta_at_secret + G1_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate [W_{zeta omega_roots}]_1
    hasher.buffer_clear();
    proof = valid_proof;
    G1_noise = libff::G1<ppT>::random_element();
    proof.W_zeta_omega_at_secret = proof.W_zeta_omega_at_secret + G1_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
    // manipulate r_zeta
    hasher.buffer_clear();
    proof = valid_proof;
    Fr_noise = libff::Fr<ppT>::random_element();
    proof.r_zeta = proof.r_zeta + Fr_noise;
    b_accept = verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_FALSE(b_accept);
}

template<typename ppT>
void test_plonk_compute_accumulator(
    const plonk_example &example,
    const libff::Fr<ppT> &beta,
    const libff::Fr<ppT> &gamma,
    const std::vector<libff::Fr<ppT>> &witness,
    const srs<ppT> &srs,
    std::shared_ptr<libfqfft::evaluation_domain<libff::Fr<ppT>>> domain)
{
    using Field = libff::Fr<ppT>;
    // A[0] = 1; ... A[i] = computed from (i-1)
    std::vector<Field> A_vector = plonk_compute_accumulator(
        srs.num_gates, beta, gamma, witness, srs.H_prime, srs.H_prime_permute);
    polynomial<Field> A_poly(srs.num_gates);
    plonk_interpolate_polynomial_from_points<Field>(A_vector, A_poly, domain);

    // initialize hard-coded values from example circuit
    printf("[%s:%d] A_poly\n", __FILE__, __LINE__);
    libff::print_vector(A_poly);
    ASSERT_EQ(A_poly, example.A_poly);
}

template<typename ppT, class transcript_hasher>
void test_plonk_prover_round_one(
    const plonk_example &example,
    const round_zero_out_t<ppT> &round_zero_out,
    const std::vector<libff::Fr<ppT>> &witness,
    const srs<ppT> &srs,
    std::shared_ptr<libfqfft::evaluation_domain<libff::Fr<ppT>>> domain,
    transcript_hasher &hasher)
{
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    round_one_out_t<ppT> round_one_out =
        plonk_prover<ppT, transcript_hasher>::round_one(
            round_zero_out, blind_scalars, witness, srs, domain, hasher);
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

template<typename ppT, class transcript_hasher>
void test_plonk_prover_round_two(
    const plonk_example &example,
    const libff::Fr<ppT> &beta,
    const libff::Fr<ppT> &gamma,
    const round_zero_out_t<ppT> &round_zero_out,
    const std::vector<libff::Fr<ppT>> &blind_scalars,
    const std::vector<libff::Fr<ppT>> &witness,
    const srs<ppT> &srs,
    std::shared_ptr<libfqfft::evaluation_domain<libff::Fr<ppT>>> domain,
    transcript_hasher &hasher)
{
    round_two_out_t<ppT> round_two_out =
        plonk_prover<ppT, transcript_hasher>::round_two(
            beta,
            gamma,
            round_zero_out,
            blind_scalars,
            witness,
            srs,
            domain,
            hasher);
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

template<typename ppT, class transcript_hasher>
void test_plonk_prover_round_three(
    const plonk_example &example,
    const libff::Fr<ppT> &alpha,
    const libff::Fr<ppT> &beta,
    const libff::Fr<ppT> &gamma,
    const round_zero_out_t<ppT> &round_zero_out,
    const round_one_out_t<ppT> &round_one_out,
    const round_two_out_t<ppT> &round_two_out,
    const std::vector<libff::Fr<ppT>> &witness,
    const srs<ppT> &srs,
    transcript_hasher &hasher)
{
    round_three_out_t<ppT> round_three_out =
        plonk_prover<ppT, transcript_hasher>::round_three(
            alpha,
            beta,
            gamma,
            round_zero_out,
            round_one_out,
            round_two_out,
            witness,
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

template<typename ppT, class transcript_hasher>
void test_plonk_prover_round_four(
    const plonk_example &example,
    const libff::Fr<ppT> &zeta,
    const round_one_out_t<ppT> &round_one_out,
    const round_three_out_t<ppT> &round_three_out,
    const srs<ppT> &srs,
    transcript_hasher &hasher)
{
    round_four_out_t<ppT> round_four_out =
        plonk_prover<ppT, transcript_hasher>::round_four(
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

template<typename ppT, class transcript_hasher>
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
    transcript_hasher &hasher)
{
    round_five_out_t<ppT> round_five_out =
        plonk_prover<ppT, transcript_hasher>::round_five(
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
template<typename ppT, class transcript_hasher> void test_plonk_prover_rounds()
{
    using Field = libff::Fr<ppT>;

    ppT::init_public_params();

    // load test vector values from example circuit
    plonk_example example;
    // random hidden element secret (toxic waste)
    Field secret = example.secret;
    // example witness
    std::vector<Field> witness = example.witness;

    // hard-coded values for the "random" blinding constants from
    // example circuit
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(example.num_gates);

    // prepare srs
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(
        usrs,
        example.gates_matrix,
        example.wire_permutation,
        example.PI_wire_indices);

    // initialize hasher
    transcript_hasher hasher;

    // Prover Round 0 (initialization)
    round_zero_out_t<ppT> round_zero_out =
        plonk_prover<ppT, transcript_hasher>::round_zero(srs);

    // --- Unit test Prover Round 1 ---
    // reset buffer at the start of the round (needed for testing only)
    printf("[%s:%d] Unit test Prover Round 1...\n", __FILE__, __LINE__);
    test_plonk_prover_round_one<ppT, transcript_hasher>(
        example, round_zero_out, witness, srs, domain, hasher);

    // --- Unit test Prover Round 2 ---
    // reset buffer at the start of the round (needed for testing only)
    printf("[%s:%d] Unit test Prover Round 2...\n", __FILE__, __LINE__);
    round_one_out_t<ppT> round_one_out =
        plonk_prover<ppT, transcript_hasher>::round_one(
            round_zero_out, blind_scalars, witness, srs, domain, hasher);
    // clear hash buffer
    hasher.buffer_clear();
    // add outputs from Round 1 to the hash buffer
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[a]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[b]);
    hasher.add_element(round_one_out.W_polys_blinded_at_secret_g1[c]);
    const libff::Fr<ppT> beta = hasher.get_hash();
    hasher.add_element(libff::Fr<ppT>::one());
    const libff::Fr<ppT> gamma = hasher.get_hash();
    test_plonk_prover_round_two<ppT, transcript_hasher>(
        example,
        beta,
        gamma,
        round_zero_out,
        blind_scalars,
        witness,
        srs,
        domain,
        hasher);

    // Unit test plonk_compute_accumulator
    test_plonk_compute_accumulator<ppT>(
        example, beta, gamma, witness, srs, domain);

    // --- Unit test Prover Round 3 ---
    // reset buffer at the start of the round (needed for testing only)
    printf("[%s:%d] Prover Round 3...\n", __FILE__, __LINE__);
    round_two_out_t<ppT> round_two_out =
        plonk_prover<ppT, transcript_hasher>::round_two(
            beta,
            gamma,
            round_zero_out,
            blind_scalars,
            witness,
            srs,
            domain,
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
    const libff::Fr<ppT> alpha = hasher.get_hash();
    test_plonk_prover_round_three<ppT, transcript_hasher>(
        example,
        alpha,
        beta,
        gamma,
        round_zero_out,
        round_one_out,
        round_two_out,
        witness,
        srs,
        hasher);

    // --- Unit test Prover Round 4 ---
    printf("[%s:%d] Prover Round 4...\n", __FILE__, __LINE__);
    round_three_out_t<ppT> round_three_out =
        plonk_prover<ppT, transcript_hasher>::round_three(
            alpha,
            beta,
            gamma,
            round_zero_out,
            round_one_out,
            round_two_out,
            witness,
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
    test_plonk_prover_round_four<ppT, transcript_hasher>(
        example, zeta, round_one_out, round_three_out, srs, hasher);

    // --- Unit test Prover Round 5 ---
    printf("[%s:%d] Unit test Prover Round 5...\n", __FILE__, __LINE__);
    round_four_out_t<ppT> round_four_out =
        plonk_prover<ppT, transcript_hasher>::round_four(
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
    test_plonk_prover_round_five<ppT, transcript_hasher>(
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
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(example.num_gates);

    // --- USRS ---
    // compute SRS = powers of secret times G1: 1*G1, secret^1*G1,
    // secret^2*G1, ... and secret times G2: 1*G2, secret^1*G2
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    // --- SRS ---
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(
        usrs,
        example.gates_matrix,
        example.wire_permutation,
        example.PI_wire_indices);
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
template<typename ppT, class transcript_hasher> void test_plonk_prover()
{
    using Field = libff::Fr<ppT>;

    ppT::init_public_params();
    // load test vector values from example circuit
    plonk_example example;
    // random hidden element secret (toxic waste)
    Field secret = example.secret;
    // example witness
    std::vector<Field> witness = example.witness;
    // hard-coded values for the "random" blinding constants from
    // example circuit
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(example.num_gates);

    // prepare srs
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(
        usrs,
        example.gates_matrix,
        example.wire_permutation,
        example.PI_wire_indices);

    // initialize hasher
    transcript_hasher hasher;

    // initialize prover
    plonk_prover<ppT, transcript_hasher> prover;
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
template<typename ppT, class transcript_hasher>
void test_plonk_verifier_preprocessed_input(
    const plonk_example &example, const srs<ppT> &srs)
{
    // compute verifier preprocessed input
    const verifier_preprocessed_input_t<ppT> preprocessed_input =
        plonk_verifier<ppT, transcript_hasher>::preprocessed_input(srs);

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

template<typename ppT, class transcript_hasher>
void test_plonk_verifier_step_five(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    std::shared_ptr<libfqfft::evaluation_domain<libff::Fr<ppT>>> domain)
{
    const step_five_out_t<ppT> step_five_out =
        plonk_verifier<ppT, transcript_hasher>::step_five(
            step_four_out, domain);
    printf("[%s:%d] zh_zeta ", __FILE__, __LINE__);
    step_five_out.zh_zeta.print();
    ASSERT_EQ(step_five_out.zh_zeta, example.zh_zeta);
}

template<typename ppT, class transcript_hasher>
void test_plonk_verifier_step_six(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const srs<ppT> &srs)
{
    const step_six_out_t<ppT> step_six_out =
        plonk_verifier<ppT, transcript_hasher>::step_six(step_four_out, srs);
    printf("L_0_zeta ");
    step_six_out.L_0_zeta.print();
    ASSERT_EQ(step_six_out.L_0_zeta, example.L_0_zeta);
}

template<typename ppT, class transcript_hasher>
void test_plonk_verifier_step_seven(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const std::vector<libff::Fr<ppT>> &PI_value_list,
    const srs<ppT> &srs)
{
    const step_seven_out_t<ppT> step_seven_out =
        plonk_verifier<ppT, transcript_hasher>::step_seven(
            step_four_out, PI_value_list, srs);
    printf("PI_zeta ");
    step_seven_out.PI_zeta.print();
    ASSERT_EQ(step_seven_out.PI_zeta, example.PI_zeta);
}

template<typename ppT, class transcript_hasher>
void test_plonk_verifier_step_eight(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const step_five_out_t<ppT> &step_five_out,
    const step_six_out_t<ppT> &step_six_out,
    const step_seven_out_t<ppT> &step_seven_out,
    const plonk_proof<ppT> &proof)
{
    const step_eight_out_t<ppT> step_eight_out =
        plonk_verifier<ppT, transcript_hasher>::step_eight(
            step_four_out, step_five_out, step_six_out, step_seven_out, proof);
    ASSERT_EQ(step_eight_out.r_prime_zeta, example.r_prime_zeta);
}

template<typename ppT, class transcript_hasher>
void test_plonk_verifier_step_nine(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const step_six_out_t<ppT> &step_six_out,
    const plonk_proof<ppT> &proof,
    const verifier_preprocessed_input_t<ppT> &preprocessed_input,
    const srs<ppT> &srs)
{
    step_nine_out_t<ppT> step_nine_out =
        plonk_verifier<ppT, transcript_hasher>::step_nine(
            step_four_out, step_six_out, proof, preprocessed_input, srs);
    step_nine_out.D1.print();
    libff::G1<ppT> D1_aff(step_nine_out.D1);
    D1_aff.to_affine_coordinates();
    ASSERT_EQ(D1_aff.X, example.D1[0]);
    ASSERT_EQ(D1_aff.Y, example.D1[1]);
}

template<typename ppT, class transcript_hasher>
void test_plonk_verifier_step_ten(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const step_nine_out_t<ppT> &step_nine_out,
    const plonk_proof<ppT> &proof,
    const verifier_preprocessed_input_t<ppT> &preprocessed_input,
    const srs<ppT> &srs)
{
    step_ten_out_t<ppT> step_ten_out =
        plonk_verifier<ppT, transcript_hasher>::step_ten(
            step_four_out, step_nine_out, proof, preprocessed_input, srs);
    printf("[%s:%d] F1\n", __FILE__, __LINE__);
    step_ten_out.F1.print();
    libff::G1<ppT> F1_aff(step_ten_out.F1);
    F1_aff.to_affine_coordinates();
    ASSERT_EQ(F1_aff.X, example.F1[0]);
    ASSERT_EQ(F1_aff.Y, example.F1[1]);
}

template<typename ppT, class transcript_hasher>
void test_plonk_verifier_step_eleven(
    const plonk_example &example,
    const step_four_out_t<ppT> &step_four_out,
    const step_eight_out_t<ppT> &step_eight_out,
    const plonk_proof<ppT> &proof)
{
    const step_eleven_out_t<ppT> step_eleven_out =
        plonk_verifier<ppT, transcript_hasher>::step_eleven(
            step_four_out, step_eight_out, proof);
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
template<typename ppT, class transcript_hasher> void test_plonk_verifier_steps()
{
    using Field = libff::Fr<ppT>;

    ppT::init_public_params();
    // load test vector values from example circuit
    plonk_example example;
    // random hidden element secret (toxic waste)
    Field secret = example.secret;
    // example witness
    std::vector<Field> witness = example.witness;
    // hard-coded values for the "random" blinding constants from
    // example circuit
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(example.num_gates);

    // prepare srs
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(
        usrs,
        example.gates_matrix,
        example.wire_permutation,
        example.PI_wire_indices);

    // initialize hasher
    transcript_hasher hasher;

    // initialize prover
    plonk_prover<ppT, transcript_hasher> prover;
    // compute proof
    plonk_proof<ppT> proof =
        prover.compute_proof(srs, witness, blind_scalars, hasher);

    // clear the hasher buffer in order to re-use the same
    // transcript_hasher object for the verifier
    hasher.buffer_clear();

    // Unit test verifier preprocessed input
    test_plonk_verifier_preprocessed_input<ppT, transcript_hasher>(
        example, srs);

    // prepare the list of PI values for the example circuit
    std::vector<Field> PI_value_list;
    for (size_t i = 0; i < example.PI_wire_indices.size(); i++) {
        Field PI_value = example.witness[example.PI_wire_indices[i]];
        PI_value_list.push_back(PI_value);
    }

    // compute step 4
    const step_four_out_t<ppT> step_four_out =
        plonk_verifier<ppT, transcript_hasher>::step_four(proof, hasher);

    // unit test verifier step 5
    test_plonk_verifier_step_five<ppT, transcript_hasher>(
        example, step_four_out, domain);

    // unit test verifier step 6
    test_plonk_verifier_step_six<ppT, transcript_hasher>(
        example, step_four_out, srs);

    // unit test verifier step 7
    test_plonk_verifier_step_seven<ppT, transcript_hasher>(
        example, step_four_out, PI_value_list, srs);

    // unit test verifier step 8
    const step_five_out_t<ppT> step_five_out =
        plonk_verifier<ppT, transcript_hasher>::step_five(
            step_four_out, domain);
    const step_six_out_t<ppT> step_six_out =
        plonk_verifier<ppT, transcript_hasher>::step_six(step_four_out, srs);
    const step_seven_out_t<ppT> step_seven_out =
        plonk_verifier<ppT, transcript_hasher>::step_seven(
            step_four_out, PI_value_list, srs);
    test_plonk_verifier_step_eight<ppT, transcript_hasher>(
        example,
        step_four_out,
        step_five_out,
        step_six_out,
        step_seven_out,
        proof);

    // unit test verifier step 9
    const verifier_preprocessed_input_t<ppT> preprocessed_input =
        plonk_verifier<ppT, transcript_hasher>::preprocessed_input(srs);
    test_plonk_verifier_step_nine<ppT, transcript_hasher>(
        example, step_four_out, step_six_out, proof, preprocessed_input, srs);

    // unit test verifier step 10
    step_nine_out_t<ppT> step_nine_out =
        plonk_verifier<ppT, transcript_hasher>::step_nine(
            step_four_out, step_six_out, proof, preprocessed_input, srs);
    test_plonk_verifier_step_ten<ppT, transcript_hasher>(
        example, step_four_out, step_nine_out, proof, preprocessed_input, srs);

    // unit test verifier step 11
    const step_eight_out_t<ppT> step_eight_out =
        plonk_verifier<ppT, transcript_hasher>::step_eight(
            step_four_out, step_five_out, step_six_out, step_seven_out, proof);
    test_plonk_verifier_step_eleven<ppT, transcript_hasher>(
        example, step_four_out, step_eight_out, proof);

    // unit test verifier pairing (step 12)
    step_ten_out_t<ppT> step_ten_out =
        plonk_verifier<ppT, transcript_hasher>::step_ten(
            step_four_out, step_nine_out, proof, preprocessed_input, srs);
    const step_eleven_out_t<ppT> step_eleven_out =
        plonk_verifier<ppT, transcript_hasher>::step_eleven(
            step_four_out, step_eight_out, proof);
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
template<typename ppT, class transcript_hasher> void test_plonk_verifier()
{
    using Field = libff::Fr<ppT>;

    ppT::init_public_params();
    // load test vector values from example circuit
    plonk_example example;
    // random hidden element secret (toxic waste)
    Field secret = example.secret;
    // example witness
    std::vector<Field> witness = example.witness;
    // hard-coded values for the "random" blinding constants from
    // example circuit
    std::vector<libff::Fr<ppT>> blind_scalars = example.prover_blind_scalars;
    // maximum degree of the encoded monomials in the usrs
    size_t max_degree = PLONK_MAX_DEGREE;

    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(example.num_gates);

    // prepare srs
    usrs<ppT> usrs = plonk_usrs_derive_from_secret<ppT>(secret, max_degree);
    srs<ppT> srs = plonk_srs_derive_from_usrs<ppT>(
        usrs,
        example.gates_matrix,
        example.wire_permutation,
        example.PI_wire_indices);

    // initialize hasher
    transcript_hasher hasher;

    // initialize prover
    plonk_prover<ppT, transcript_hasher> prover;
    // compute proof
    plonk_proof<ppT> proof =
        prover.compute_proof(srs, witness, blind_scalars, hasher);

    // clear the hasher buffer in order to re-use the same
    // transcript_hasher object for the verifier
    hasher.buffer_clear();

    // initialize verifier
    plonk_verifier<ppT, transcript_hasher> verifier;
    // prepare the list of PI values for the example circuit
    std::vector<Field> PI_value_list;
    for (size_t i = 0; i < example.PI_wire_indices.size(); i++) {
        Field PI_value = example.witness[example.PI_wire_indices[i]];
        PI_value_list.push_back(PI_value);
    }
    // verify proof
    bool b_valid_proof =
        verifier.verify_proof(proof, srs, PI_value_list, hasher);
    ASSERT_TRUE(b_valid_proof);

    // clear the hasher buffer in order to re-use the same
    // transcript_hasher object
    hasher.buffer_clear();
    // assert that proof verification fails when the proof is
    // manipulated
    test_verify_invalid_proof(proof, srs, PI_value_list, hasher);
}

/// \attention the example class is defined specifically for the BLS12-381
/// curve, so make sure we are using this curve
template<typename ppT> void test_plonk_gates_matrix_transpose()
{
    using Field = libff::Fr<ppT>;
    // load gates matrix from example circuit
    plonk_example example;
    std::vector<std::vector<Field>> gates_matrix_transpose =
        plonk_gates_matrix_transpose(example.gates_matrix);
    ASSERT_EQ(gates_matrix_transpose, example.gates_matrix_transpose);
}

template<typename ppT> void test_plonk_constants_k1_k2()
{
    using Field = libff::Fr<ppT>;
    Field k1, k2;
    // n = 2^s
    const size_t n = std::pow(2, k1.s);
    bool b_valid = false;
    // load k1,k2 from example circuit
    plonk_example example;
    k1 = example.k1;
    k2 = example.k2;
    b_valid = plonk_are_valid_constants_k1_k2(n, k1, k2);
    ASSERT_TRUE(b_valid);
    // check invalid k1,k2
    for (size_t i = 1; i <= example.num_gates; ++i) {
        size_t ipower = i;
        // invalid k2=k1*(omega^i)
        k1 = example.k1;
        k2 = k1 * (example.omega_base ^ ipower);
        b_valid = plonk_are_valid_constants_k1_k2(n, k1, k2);
        ASSERT_FALSE(b_valid);
        // invalid k1=k2*(omega^i)
        k2 = example.k2;
        k1 = k2 * (example.omega_base ^ ipower);
        b_valid = plonk_are_valid_constants_k1_k2(n, k1, k2);
        ASSERT_FALSE(b_valid);
        // invalid k1=omega^i
        k1 = (example.omega_base ^ ipower);
        k2 = example.k2;
        b_valid = plonk_are_valid_constants_k1_k2(n, k1, k2);
        ASSERT_FALSE(b_valid);
        // invalid k2=omega^i
        k1 = example.k1;
        k2 = (example.omega_base ^ ipower);
        b_valid = plonk_are_valid_constants_k1_k2(n, k1, k2);
        ASSERT_FALSE(b_valid);
    }
    // generate new k1,k2 and assert they are valid for a random
    // number of tests
    size_t ntests = 1UL << 10;
    for (size_t i = 0; i < ntests; ++i) {
        k1 = 0;
        k2 = 0;
        plonk_generate_constants_k1_k2(n, k1, k2);
        b_valid = plonk_are_valid_constants_k1_k2(n, k1, k2);
        // printf("k1 "); k1.print(); printf("k2 "); k2.print();
        ASSERT_TRUE(b_valid);
    }
}

TEST(TestPlonk, BLS12_381)
{
    test_plonk_srs<libff::bls12_381_pp>();
    test_plonk_prover_rounds<
        libff::bls12_381_pp,
        bls12_381_test_vector_transcript_hasher>();
    test_plonk_prover<
        libff::bls12_381_pp,
        bls12_381_test_vector_transcript_hasher>();
    test_plonk_verifier_steps<
        libff::bls12_381_pp,
        bls12_381_test_vector_transcript_hasher>();
    test_plonk_verifier<
        libff::bls12_381_pp,
        bls12_381_test_vector_transcript_hasher>();
    test_plonk_gates_matrix_transpose<libff::bls12_381_pp>();
    test_plonk_constants_k1_k2<libff::bls12_381_pp>();
}

} // namespace libsnark
