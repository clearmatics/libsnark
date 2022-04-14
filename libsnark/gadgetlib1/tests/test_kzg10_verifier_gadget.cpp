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
void do_test_kzg10_batched_commit_minus_eval_sum(
    const libff::Fr<other_curve<wppT>> &gamma,
    const std::vector<libff::Fr<other_curve<wppT>>> &evals,
    const std::vector<typename kzg10<other_curve<wppT>>::commitment> &cms,
    const libff::G1<other_curve<wppT>> &r,
    const bool expected_result)
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

    if (expected_result) {
        ASSERT_TRUE(pb.is_satisfied());
        ASSERT_EQ(r, result_var.get_element());
    } else {
        ASSERT_NE(r, result_var.get_element());
    }

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
    const nG1 r_2 = nField((13 - 3) + 51 * (17 - 5)) * nG1::one();
    do_test_kzg10_batched_commit_minus_eval_sum<wppT, 2>(
        gamma, {evals[0], evals[1]}, {cms[0], cms[1]}, r_2, true);

    // 3-entry case
    const nG1 r_3 =
        nField((13 - 3) + 51 * (17 - 5) + (51 * 51) * (19 - 7)) * nG1::one();
    do_test_kzg10_batched_commit_minus_eval_sum<wppT, 3>(
        gamma,
        {evals[0], evals[1], evals[2]},
        {cms[0], cms[1], cms[2]},
        r_3,
        true);

    // 4-entry case
    const nG1 r_4 = nField(
                        (13 - 3) + 51 * (17 - 5) + (51 * 51) * (19 - 7) +
                        (51 * 51 * 51) * (23 - 11)) *
                    nG1::one();
    do_test_kzg10_batched_commit_minus_eval_sum<wppT, 4>(
        gamma, evals, cms, r_4, true);
}

TEST(TestKZG10VerifierGadget, ValidEvaluation)
{
    test_kzg10_verifier_gadget<libff::bw6_761_pp>();
}

TEST(TestKZG10VerifierGadget, BatchedCommitMinusEvalSum)
{
    test_kzg10_batched_commit_minus_eval_sum_gadget<libff::bw6_761_pp>();
}

} // namespace

int main(int argc, char **argv)
{
    libff::bw6_761_pp::init_public_params();
    libff::bls12_377_pp::init_public_params();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
