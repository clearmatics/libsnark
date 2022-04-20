/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libsnark/gadgetlib1/gadgets/pairing/bw6_761_bls12_377/bw6_761_pairing_params.hpp"
#include "libsnark/gadgetlib1/gadgets/verifiers/kzg10_batched_verifier_gadget.hpp"
#include "libsnark/gadgetlib1/gadgets/verifiers/kzg10_verifier_gadget.hpp"
#include "libsnark/polynomial_commitments/kzg10_batched.hpp"
#include "libsnark/polynomial_commitments/tests/polynomial_commitment_test_utils.hpp"
#include "libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp"

#include <gtest/gtest.h>
#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/algebra/curves/bw6_761/bw6_761_pp.hpp>
#include <libff/algebra/fields/field_utils.hpp>

using namespace libsnark;

static const size_t POLYNOMIAL_MAX_DEGREE = 9;
static const size_t POLYNOMIAL_SIZE = 5;

namespace
{

template<typename wppT, typename scheme>
void test_polynomial_commitment_verifier_gadget(
    const typename scheme::srs &srs,
    const typename scheme::commitment &C,
    const libff::Fr<other_curve<wppT>> &i,
    const libff::Fr<other_curve<wppT>> &evaluation,
    const typename kzg10<other_curve<wppT>>::evaluation_witness &eval_witness,
    const bool expected_result)
{
    using Field = libff::Fr<wppT>;
    using npp = other_curve<wppT>;

    // Ensure that the native implementation gives the expected result.
    ASSERT_EQ(
        expected_result,
        scheme::verify_evaluation(i, evaluation, srs, eval_witness, C));

    // Perform the equivalent check in an r1cs circuit
    protoboard<Field> pb;

    kzg10_srs_variable<wppT> srs_var(pb, POLYNOMIAL_MAX_DEGREE, "srs_var");
    kzg10_commitment_variable<wppT> C_var(pb, "C_var");
    pb_variable<Field> i_var;
    i_var.allocate(pb, "i_var");
    pb_variable<Field> poly_eval_var;
    poly_eval_var.allocate(pb, "poly_eval_var");
    kzg10_witness_variable<wppT> witness_var(pb, "witness_var");
    pb_variable<Field> result_var;
    result_var.allocate(pb, "result_var");

    kzg10_verifier_gadget<wppT> verifier_gadget(
        pb,
        srs_var,
        C_var,
        i_var,
        poly_eval_var,
        witness_var,
        result_var,
        "verifier_gadget");

    verifier_gadget.generate_r1cs_constraints();

    // i (the value at which the polynomial is evaluated) and poly_i (the
    // evaluation of the polynomial) must be converted from libff::Fr<npp> to
    // libff::Fr<wppT>.
    libff::Fr<wppT> wrapping_i;
    libff::Fr<wppT> wrapping_evaluation;
    libff::fp_from_fp(wrapping_i, i);
    libff::fp_from_fp(wrapping_evaluation, evaluation);

    srs_var.generate_r1cs_witness(srs);
    C_var.generate_r1cs_witness(C);
    pb.val(i_var) = wrapping_i;
    pb.val(poly_eval_var) = wrapping_evaluation;
    witness_var.generate_r1cs_witness(eval_witness);
    verifier_gadget.generate_r1cs_witness();

    // Check some members of verifier_gadget
    {
        const libff::G2<npp> i_in_G2_val = i * libff::G2<npp>::one();
        const libff::G2<npp> B_val = srs.alpha_g2 - i_in_G2_val;

        const libff::G1<npp> poly_eval_in_G1_val =
            evaluation * libff::G1<npp>::one();
        const libff::G1<npp> C_val = C - poly_eval_in_G1_val;

        ASSERT_EQ(i_in_G2_val, verifier_gadget.i_in_G2.get_element());
        ASSERT_EQ(B_val, verifier_gadget.B.get_element());
        ASSERT_EQ(
            poly_eval_in_G1_val, verifier_gadget.poly_eval_in_G1.get_element());
        ASSERT_EQ(C_val, verifier_gadget.C.get_element());
    }

    ASSERT_TRUE(pb.is_satisfied());
    ASSERT_EQ(
        expected_result ? Field::one() : Field::zero(), pb.val(result_var));

    // Test in zkproof
    const r1cs_gg_ppzksnark_keypair<wppT> keypair =
        r1cs_gg_ppzksnark_generator<wppT>(pb.get_constraint_system(), true);
    const r1cs_gg_ppzksnark_proof<wppT> proof = r1cs_gg_ppzksnark_prover<wppT>(
        keypair.pk, pb.primary_input(), pb.auxiliary_input(), true);
    ASSERT_TRUE(r1cs_gg_ppzksnark_verifier_strong_IC<wppT>(
        keypair.vk, pb.primary_input(), proof));
}

template<typename wppT> void test_kzg10_verifier_gadget()
{
    using npp = other_curve<wppT>;
    using scheme = kzg10<npp>;

    // SRS
    const typename scheme::srs srs = scheme::setup(POLYNOMIAL_MAX_DEGREE);

    // Generate polynomial and commitment
    const polynomial<libff::Fr<npp>> poly =
        gen_polynomial<npp>(POLYNOMIAL_SIZE);
    const typename scheme::commitment C = scheme::commit(srs, poly);

    // Evaluation and witness
    const libff::Fr<npp> i = libff::Fr<npp>::random_element();
    const libff::Fr<npp> evaluation = scheme::evaluate_polynomial(poly, i);
    const typename scheme::evaluation_witness eval_witness =
        scheme::create_evaluation_witness(poly, i, evaluation, srs);

    // Check evaluation and proof natively
    ASSERT_TRUE(scheme::verify_evaluation(i, evaluation, srs, eval_witness, C));

    test_polynomial_commitment_verifier_gadget<wppT, scheme>(
        srs, C, i, evaluation, eval_witness, true);

    // Test some failure cases:

    // Invalid cases
    {
        // Invalid evaluation point
        test_polynomial_commitment_verifier_gadget<wppT, scheme>(
            srs, C, i + 1, evaluation, eval_witness, false);
        // Invalid evaluation
        test_polynomial_commitment_verifier_gadget<wppT, scheme>(
            srs, C, i, evaluation + 1, eval_witness, false);
        // Invalid evaluation witness
        test_polynomial_commitment_verifier_gadget<wppT, scheme>(
            srs, C, i, evaluation, eval_witness + eval_witness, false);
        // Invalid commitment
        test_polynomial_commitment_verifier_gadget<wppT, scheme>(
            srs, C + C, i, evaluation, eval_witness, false);
    }
}

template<typename wppT, size_t num_entries>
void do_test_kzg10_batched_compute_gamma_eval_sum_gadget(
    const libff::Fr<other_curve<wppT>> &gamma,
    const std::vector<libff::Fr<other_curve<wppT>>> &evals,
    const libff::Fr<other_curve<wppT>> &expected_result)
{
    using Field = libff::Fr<wppT>;

    ASSERT_EQ(num_entries, evals.size());

    // Protoboard and constraints

    protoboard<Field> pb;

    pb_variable<Field> gamma_var = pb_variable_allocate<Field>(pb, "gamma_var");
    pb_variable_array<Field> gamma_powers_var;
    gamma_powers_var.allocate(pb, num_entries - 2, "gamma_powers_var");
    pb_variable_array<Field> evals_var;
    evals_var.allocate(pb, num_entries, "evals_var");
    pb_variable<Field> result_var;
    result_var.allocate(pb, "result_var");

    kzg10_batched_compute_gamma_eval_sum<wppT, num_entries>
        compute_gamma_evals_sum(
            pb,
            gamma_var,
            gamma_powers_var,
            evals_var,
            result_var,
            "compute_gamma_evals_sum");

    compute_gamma_evals_sum.generate_r1cs_constraints();

    // Witness

    Field wrapping_gamma;
    fp_from_fp(wrapping_gamma, gamma);
    pb.val(gamma_var) = wrapping_gamma;

    for (size_t i = 0; i < num_entries - 2; ++i) {
        Field wrapping_eval;
        fp_from_fp(wrapping_eval, gamma ^ (i + 2));
        pb.val(gamma_powers_var[i]) = wrapping_eval;
    }

    for (size_t i = 0; i < num_entries; ++i) {
        Field wrapping_eval;
        fp_from_fp(wrapping_eval, evals[i]);
        pb.val(evals_var[i]) = wrapping_eval;
    }
    compute_gamma_evals_sum.generate_r1cs_witness();

    // Check result value

    Field wrapping_expected_result;
    fp_from_fp(wrapping_expected_result, expected_result);

    ASSERT_TRUE(pb.is_satisfied());
    ASSERT_EQ(wrapping_expected_result, pb.val(result_var));

    // Test in proof

    const r1cs_gg_ppzksnark_keypair<wppT> keypair =
        r1cs_gg_ppzksnark_generator<wppT>(pb.get_constraint_system(), true);
    const r1cs_gg_ppzksnark_proof<wppT> proof = r1cs_gg_ppzksnark_prover<wppT>(
        keypair.pk, pb.primary_input(), pb.auxiliary_input(), true);
    ASSERT_TRUE(r1cs_gg_ppzksnark_verifier_strong_IC<wppT>(
        keypair.vk, pb.primary_input(), proof));
}

template<typename wppT> void test_kzg10_compute_gamma_evals_sum_gadget()
{
    using nField = libff::Fr<other_curve<wppT>>;

    const nField gamma("51");
    const std::vector<nField> evals{{"3"}, {"5"}, {"7"}, {"11"}};

    const nField r_3(3 + 51 * 5 + (51 * 51) * 7);
    do_test_kzg10_batched_compute_gamma_eval_sum_gadget<wppT, 3>(
        gamma, {evals[0], evals[1], evals[2]}, r_3);

    // 4-entry case
    const nField r_4(3 + 51 * 5 + (51 * 51) * 7 + (51 * 51 * 51) * 11);
    do_test_kzg10_batched_compute_gamma_eval_sum_gadget<wppT, 4>(
        gamma, evals, r_4);
}

template<typename wppT, size_t num_entries>
void do_test_kzg10_batched_commit_minus_eval_sum(
    const libff::Fr<other_curve<wppT>> &gamma,
    const std::vector<libff::Fr<other_curve<wppT>>> &evals,
    const std::vector<typename kzg10<other_curve<wppT>>::commitment> &cms,
    const libff::G1<other_curve<wppT>> &expect_result)
{
    using Field = libff::Fr<wppT>;

    ASSERT_EQ(num_entries, evals.size());
    ASSERT_EQ(num_entries, cms.size());

    // Protoboard and constraints

    protoboard<Field> pb;

    pb_variable<Field> gamma_var = pb_variable_allocate<Field>(pb, "gamma_var");
    std::vector<G1_variable<wppT>> cms_var =
        internal::allocate_variable_array<wppT, G1_variable<wppT>>(
            pb, num_entries, "cms_var");
    pb_variable_array<Field> evals_var;
    evals_var.allocate(pb, num_entries, "evals_var");
    G1_variable<wppT> result_var(pb, "result_var");

    kzg10_batched_compute_commit_minus_eval_sum<wppT, num_entries> compute_sum(
        pb, gamma_var, cms_var, evals_var, result_var, "compute_sum");

    compute_sum.generate_r1cs_constraints();

    // Witness

    Field wrapping_gamma;
    fp_from_fp(wrapping_gamma, gamma);
    pb.val(gamma_var) = wrapping_gamma;

    for (size_t i = 0; i < num_entries; ++i) {
        cms_var[i].generate_r1cs_witness(cms[i]);

        Field wrapping_eval;
        fp_from_fp(wrapping_eval, evals[i]);
        pb.val(evals_var[i]) = wrapping_eval;
    }
    compute_sum.generate_r1cs_witness();

    // Check result value

    ASSERT_TRUE(pb.is_satisfied());
    ASSERT_EQ(expect_result, result_var.get_element());

    // Test in proof

    const r1cs_gg_ppzksnark_keypair<wppT> keypair =
        r1cs_gg_ppzksnark_generator<wppT>(pb.get_constraint_system(), true);
    const r1cs_gg_ppzksnark_proof<wppT> proof = r1cs_gg_ppzksnark_prover<wppT>(
        keypair.pk, pb.primary_input(), pb.auxiliary_input(), true);
    ASSERT_TRUE(r1cs_gg_ppzksnark_verifier_strong_IC<wppT>(
        keypair.vk, pb.primary_input(), proof));
}

template<typename wppT> void test_kzg10_batched_commit_minus_eval_sum_gadget()
{
    using nField = libff::Fr<other_curve<wppT>>;
    using nG1 = libff::G1<other_curve<wppT>>;

    const nField gamma("51");
    const std::vector<nField> evals{{"3"}, {"5"}, {"7"}, {"11"}};
    const std::vector<nG1> cms{{
        nField("13") * nG1::one(),
        nField("17") * nG1::one(),
        nField("19") * nG1::one(),
        nField("23") * nG1::one(),
    }};

    // 2-entry case
    const nG1 result_2 = nField((13 - 3) + 51 * (17 - 5)) * nG1::one();
    do_test_kzg10_batched_commit_minus_eval_sum<wppT, 2>(
        gamma, {evals[0], evals[1]}, {cms[0], cms[1]}, result_2);

    // 3-entry case
    const nG1 result_3 =
        nField((13 - 3) + 51 * (17 - 5) + (51 * 51) * (19 - 7)) * nG1::one();
    do_test_kzg10_batched_commit_minus_eval_sum<wppT, 3>(
        gamma,
        {evals[0], evals[1], evals[2]},
        {cms[0], cms[1], cms[2]},
        result_3);

    // 4-entry case
    const nG1 result_4 = nField(
                             (13 - 3) + 51 * (17 - 5) + (51 * 51) * (19 - 7) +
                             (51 * 51 * 51) * (23 - 11)) *
                         nG1::one();
    do_test_kzg10_batched_commit_minus_eval_sum<wppT, 4>(
        gamma, evals, cms, result_4);
}

template<typename wppT, size_t num_polynomials_1, size_t num_polynomials_2>
void do_test_kzg10_batched_verifier_gadget(
    const libff::Fr<other_curve<wppT>> &z_1,
    const libff::Fr<other_curve<wppT>> &z_2,
    const typename kzg10_batched_2_point<other_curve<wppT>>::evaluations
        &evaluations,
    const typename kzg10<other_curve<wppT>>::srs &srs,
    const libff::Fr<other_curve<wppT>> &gamma_1,
    const libff::Fr<other_curve<wppT>> &gamma_2,
    const typename kzg10_batched_2_point<other_curve<wppT>>::evaluation_witness
        &eval_witness,
    const std::vector<typename kzg10<other_curve<wppT>>::commitment> &cm_1s,
    const std::vector<typename kzg10<other_curve<wppT>>::commitment> &cm_2s,
    const libff::Fr<other_curve<wppT>> &r,
    const bool expected_result)
{
    using Field = libff::Fr<wppT>;
    using npp = other_curve<wppT>;
    using nG1 = libff::G1<npp>;

    protoboard<Field> pb;

    // z_1 and z_2
    pb_variable<Field> z_1_var;
    z_1_var.allocate(pb, "z_1_var");
    pb_variable<Field> z_2_var;
    z_2_var.allocate(pb, "z_2_var");

    // Evaluations
    pb_variable_array<Field> poly_1_evals_var;
    poly_1_evals_var.allocate(pb, evaluations.s_1s.size(), "poly_1_evals_var");

    pb_variable_array<Field> poly_2_evals_var;
    poly_2_evals_var.allocate(pb, evaluations.s_2s.size(), "poly_2_evals_var");

    // srs
    kzg10_srs_variable<wppT> srs_var(pb, POLYNOMIAL_MAX_DEGREE, "srs_var");

    // Gammas
    pb_variable<Field> gamma_1_var;
    gamma_1_var.allocate(pb, "gamma_1_var");
    pb_variable<Field> gamma_2_var;
    gamma_2_var.allocate(pb, "gamma_2_var");

    // Witness
    kzg10_batched_witness_variable<wppT> eval_witness_var(
        pb, "eval_witness_var");

    // Commitments
    std::vector<kzg10_commitment_variable<wppT>> cm_1s_var;
    std::vector<kzg10_commitment_variable<wppT>> cm_2s_var;

    for (size_t i = 0; i < cm_1s.size(); ++i) {
        cm_1s_var.push_back(
            kzg10_commitment_variable<wppT>(pb, FMT("", "cm_1s_var[%zu]", i)));
    }

    for (size_t i = 0; i < cm_2s.size(); ++i) {
        cm_2s_var.push_back(
            kzg10_commitment_variable<wppT>(pb, FMT("", "cm_2s_var[%zu]", i)));
    }

    // r
    pb_variable<Field> r_var;
    r_var.allocate(pb, "r_var");

    // result
    pb_variable<Field> result_var;
    result_var.allocate(pb, "result_var");

    // Verifier gadget
    kzg10_batched_verifier_gadget<wppT, num_polynomials_1, num_polynomials_2>
        verifier_gadget(
            pb,
            z_1_var,
            z_2_var,
            poly_1_evals_var,
            poly_2_evals_var,
            srs_var,
            gamma_1_var,
            gamma_2_var,
            eval_witness_var,
            cm_1s_var,
            cm_2s_var,
            r_var,
            result_var,
            "verifier_gadget");

    verifier_gadget.generate_r1cs_constraints();

    // Field containers of nField elements
    Field wrapping_z_1;
    libff::fp_from_fp(wrapping_z_1, z_1);
    Field wrapping_z_2;
    libff::fp_from_fp(wrapping_z_2, z_2);
    std::vector<Field> wrapping_poly_1_evals(evaluations.s_1s.size());
    for (size_t i = 0; i < evaluations.s_1s.size(); ++i) {
        libff::fp_from_fp(wrapping_poly_1_evals[i], evaluations.s_1s[i]);
    }
    std::vector<Field> wrapping_poly_2_evals(evaluations.s_2s.size());
    for (size_t i = 0; i < evaluations.s_2s.size(); ++i) {
        libff::fp_from_fp(wrapping_poly_2_evals[i], evaluations.s_2s[i]);
    }
    Field wrapping_gamma_1;
    libff::fp_from_fp(wrapping_gamma_1, gamma_1);
    Field wrapping_gamma_2;
    libff::fp_from_fp(wrapping_gamma_2, gamma_2);
    Field wrapping_r;
    libff::fp_from_fp(wrapping_r, r);

    // Assign witnesses to all parameters
    pb.val(z_1_var) = wrapping_z_1;
    pb.val(z_2_var) = wrapping_z_2;
    poly_1_evals_var.fill_with_field_elements(pb, wrapping_poly_1_evals);
    poly_2_evals_var.fill_with_field_elements(pb, wrapping_poly_2_evals);
    srs_var.generate_r1cs_witness(srs);
    pb.val(gamma_1_var) = wrapping_gamma_1;
    pb.val(gamma_2_var) = wrapping_gamma_2;
    eval_witness_var.generate_r1cs_witness(eval_witness);
    for (size_t i = 0; i < cm_1s.size(); ++i) {
        cm_1s_var[i].generate_r1cs_witness(cm_1s[i]);
    }
    for (size_t i = 0; i < cm_2s.size(); ++i) {
        cm_2s_var[i].generate_r1cs_witness(cm_2s[i]);
    }
    pb.val(r_var) = wrapping_r;

    verifier_gadget.generate_r1cs_witness();

    // Check intermediate values

    if (expected_result) {
        const nG1 G_expect = internal::gamma_times_commit_minus_eval_sum<npp>(
            gamma_1, evaluations.s_1s, cm_1s);
        const nG1 H_expect = internal::gamma_times_commit_minus_eval_sum<npp>(
            gamma_2, evaluations.s_2s, cm_2s);
        const nG1 F_expect = G_expect + r * H_expect;

        const nG1 r_times_W_2_expect = r * eval_witness.W_2;
        const nG1 A_expect = eval_witness.W_1 + r_times_W_2_expect;

        const nG1 r_times_z_2_times_W_2_expect = z_2 * r_times_W_2_expect;
        const nG1 z_1_times_W_1_expect = z_1 * eval_witness.W_1;
        const nG1 F_plus_z_1_times_W_1_expect = F_expect + z_1_times_W_1_expect;
        const nG1 B_expect =
            F_plus_z_1_times_W_1_expect + r_times_z_2_times_W_2_expect;

        ASSERT_EQ(G_expect, verifier_gadget.G.get_element());
        ASSERT_EQ(H_expect, verifier_gadget.H.get_element());
        ASSERT_EQ(F_expect, verifier_gadget.F.get_element());
        ASSERT_EQ(
            r_times_W_2_expect, verifier_gadget.r_times_W_2.get_element());
        ASSERT_EQ(A_expect, verifier_gadget.A.get_element());
        ASSERT_EQ(
            r_times_z_2_times_W_2_expect,
            verifier_gadget.r_times_z_2_times_W_2.get_element());
        ASSERT_EQ(
            z_1_times_W_1_expect, verifier_gadget.z_1_times_W_1.get_element());
        ASSERT_EQ(
            F_plus_z_1_times_W_1_expect,
            verifier_gadget.F_plus_z_1_times_W_1.get_element());
        ASSERT_EQ(B_expect, verifier_gadget.B.get_element());
    }

    // Check the result.

    ASSERT_TRUE(pb.is_satisfied());
    ASSERT_EQ(
        expected_result ? Field::one() : Field::zero(), pb.val(result_var));

    // Test proof gen

    const r1cs_gg_ppzksnark_keypair<wppT> keypair =
        r1cs_gg_ppzksnark_generator<wppT>(pb.get_constraint_system(), true);
    const r1cs_gg_ppzksnark_proof<wppT> proof = r1cs_gg_ppzksnark_prover<wppT>(
        keypair.pk, pb.primary_input(), pb.auxiliary_input(), true);
    ASSERT_TRUE(r1cs_gg_ppzksnark_verifier_strong_IC<wppT>(
        keypair.vk, pb.primary_input(), proof));
}

template<typename wppT> void test_kzg10_batched_verifier_gadget()
{
    using npp = other_curve<wppT>;
    using scheme = kzg10<npp>;

    using nField = libff::Fr<npp>;

    // 2 polynomials to be evaluated at the first point, 3 at the second.
    const size_t num_polynomials_1 = 2;
    const size_t num_polynomials_2 = 3;

    // SRS
    const typename scheme::srs srs = scheme::setup(POLYNOMIAL_MAX_DEGREE);

    // Generate polynomials and commitment
    std::vector<polynomial<nField>> polynomials_1;
    std::vector<polynomial<nField>> polynomials_2;
    std::vector<typename kzg10<npp>::commitment> cm_1s;
    std::vector<typename kzg10<npp>::commitment> cm_2s;

    for (size_t i = 0; i < num_polynomials_1; ++i) {
        polynomials_1.push_back(gen_polynomial<npp>(POLYNOMIAL_SIZE));
        cm_1s.push_back(kzg10<npp>::commit(srs, polynomials_1.back()));
    }

    for (size_t i = 0; i < num_polynomials_2; ++i) {
        polynomials_2.push_back(gen_polynomial<npp>(POLYNOMIAL_SIZE));
        cm_2s.push_back(kzg10<npp>::commit(srs, polynomials_2.back()));
    }

    std::vector<typename kzg10<npp>::commitment> cm_1s_invalid{
        cm_1s[0] + cm_1s[0], cm_1s[1]};
    std::vector<typename kzg10<npp>::commitment> cm_2s_invalid{
        cm_2s[0] + cm_2s[0], cm_2s[1], cm_2s[2]};

    // Evaluations
    const nField z_1 = nField("13");
    const nField z_2 = nField("17");

    const typename kzg10_batched_2_point<npp>::evaluations evals =
        kzg10_batched_2_point<npp>::evaluate_polynomials(
            polynomials_1, polynomials_2, z_1, z_2);

    // Witness
    const nField gamma_1 = nField("3");
    const nField gamma_2 = nField("5");

    const typename kzg10_batched_2_point<npp>::evaluation_witness eval_witness =
        kzg10_batched_2_point<npp>::create_evaluation_witness(
            polynomials_1,
            polynomials_2,
            z_1,
            z_2,
            evals,
            srs,
            gamma_1,
            gamma_2);

    // Check evaluations and witness natively
    const nField r = nField("7");
    ASSERT_TRUE(kzg10_batched_2_point<npp>::verify_evaluations(
        z_1, z_2, evals, srs, gamma_1, gamma_2, eval_witness, cm_1s, cm_2s, r));

    // Test gadget in the positive case
    do_test_kzg10_batched_verifier_gadget<
        wppT,
        num_polynomials_1,
        num_polynomials_2>(
        z_1,
        z_2,
        evals,
        srs,
        gamma_1,
        gamma_2,
        eval_witness,
        cm_1s,
        cm_2s,
        r,
        true);

    // Test some failure cases:

    // Invalid cases
    {
        // Invalid evaluation point
        do_test_kzg10_batched_verifier_gadget<
            wppT,
            num_polynomials_1,
            num_polynomials_2>(
            z_1 + z_1,
            z_2,
            evals,
            srs,
            gamma_1,
            gamma_2,
            eval_witness,
            cm_1s,
            cm_2s,
            r,
            false);

        do_test_kzg10_batched_verifier_gadget<
            wppT,
            num_polynomials_1,
            num_polynomials_2>(
            z_1,
            z_2 + z_2,
            evals,
            srs,
            gamma_1,
            gamma_2,
            eval_witness,
            cm_1s,
            cm_2s,
            r,
            false);

        do_test_kzg10_batched_verifier_gadget<
            wppT,
            num_polynomials_1,
            num_polynomials_2>(
            z_1,
            z_2,
            evals,
            srs,
            gamma_1 + gamma_1,
            gamma_2,
            eval_witness,
            cm_1s,
            cm_2s,
            r,
            false);

        do_test_kzg10_batched_verifier_gadget<
            wppT,
            num_polynomials_1,
            num_polynomials_2>(
            z_1,
            z_2,
            evals,
            srs,
            gamma_1,
            gamma_2 + gamma_2,
            eval_witness,
            cm_1s,
            cm_2s,
            r,
            false);

        do_test_kzg10_batched_verifier_gadget<
            wppT,
            num_polynomials_1,
            num_polynomials_2>(
            z_1,
            z_2,
            evals,
            srs,
            gamma_1,
            gamma_2,
            eval_witness,
            cm_1s_invalid,
            cm_2s,
            r,
            false);

        do_test_kzg10_batched_verifier_gadget<
            wppT,
            num_polynomials_1,
            num_polynomials_2>(
            z_1,
            z_2,
            evals,
            srs,
            gamma_1,
            gamma_2,
            eval_witness,
            cm_1s,
            cm_2s_invalid,
            r,
            false);
    }
}

TEST(TestKZG10VerifierGadget, ValidEvaluation)
{
    test_kzg10_verifier_gadget<libff::bw6_761_pp>();
}

TEST(TestKZG10VerifierGadget, BatchedComputeGamma_EvalsSum)
{
    test_kzg10_compute_gamma_evals_sum_gadget<libff::bw6_761_pp>();
}

TEST(TestKZG10VerifierGadget, BatchedCommitMinusEvalSum)
{
    test_kzg10_batched_commit_minus_eval_sum_gadget<libff::bw6_761_pp>();
}

TEST(TestKZG10VerifierGadget, BatchedValidEvaluation)
{
    test_kzg10_batched_verifier_gadget<libff::bw6_761_pp>();
}

} // namespace

int main(int argc, char **argv)
{
    libff::bw6_761_pp::init_public_params();
    libff::bls12_377_pp::init_public_params();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
