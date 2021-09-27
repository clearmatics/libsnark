/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libsnark/gadgetlib1/gadgets/pairing/bw6_761_bls12_377/bw6_761_pairing_params.hpp"
#include "libsnark/gadgetlib1/gadgets/verifiers/kzg10_verifier_gadget.hpp"
#include "libsnark/polynomial_commitments/kzg10.hpp"
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
    const evaluation_and_witness<other_curve<wppT>, scheme> &eval_witness,
    const bool expected_result)
{
    using Field = libff::Fr<wppT>;
    using npp = other_curve<wppT>;

    // Ensure that the native implementation gives the expected result.
    ASSERT_EQ(expected_result, scheme::verify_eval(srs, C, eval_witness));

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
    libff::Fr<wppT> wrapping_poly_i;
    libff::fp_from_fp(wrapping_i, eval_witness.i);
    libff::fp_from_fp(wrapping_poly_i, eval_witness.phi_i);

    srs_var.generate_r1cs_witness(srs);
    C_var.generate_r1cs_witness(C);
    pb.val(i_var) = wrapping_i;
    pb.val(poly_eval_var) = wrapping_poly_i;
    witness_var.generate_r1cs_witness(eval_witness.w);
    verifier_gadget.generate_r1cs_witness();

    // Check some members of verifier_gadget
    {
        const libff::G2<npp> i_in_G2_val =
            eval_witness.i * libff::G2<npp>::one();
        const libff::G2<npp> B_val = srs.alpha_g2 - i_in_G2_val;

        const libff::G1<npp> poly_eval_in_G1_val =
            eval_witness.phi_i * libff::G1<npp>::one();
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
    const evaluation_and_witness<npp, scheme> eval_witness =
        scheme::create_witness(srs, poly, i);

    // Check evaluation and proof natively
    ASSERT_TRUE(scheme::verify_eval(srs, C, eval_witness));

    test_polynomial_commitment_verifier_gadget<wppT, scheme>(
        srs, C, eval_witness, true);

    // Test some failure cases:

    // Invalid commitment
    test_polynomial_commitment_verifier_gadget<wppT, scheme>(
        srs, C + C, eval_witness, false);

    // Invalid evaluation point / polyomial evaluation.
    {
        evaluation_and_witness<npp, scheme> eval_witness_invalid = eval_witness;
        eval_witness_invalid.i += libff::Fr<npp>::one();
        test_polynomial_commitment_verifier_gadget<wppT, scheme>(
            srs, C, eval_witness_invalid, false);
    }

    // Invalid evaluation witness
    {
        evaluation_and_witness<npp, scheme> eval_witness_invalid = eval_witness;
        eval_witness_invalid.w =
            eval_witness_invalid.w + scheme::witness::one();
        test_polynomial_commitment_verifier_gadget<wppT, scheme>(
            srs, C, eval_witness_invalid, false);
    }
}

TEST(TestKZG10VerifierGadget, ValidEvaluation)
{
    test_kzg10_verifier_gadget<libff::bw6_761_pp>();
}

} // namespace

int main(int argc, char **argv)
{
    libff::bw6_761_pp::init_public_params();
    libff::bls12_377_pp::init_public_params();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
