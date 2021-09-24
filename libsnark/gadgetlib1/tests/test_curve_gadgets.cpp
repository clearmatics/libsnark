/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libsnark/gadgetlib1/gadgets/pairing/bw6_761_bls12_377/bw6_761_pairing_params.hpp"
#include "libsnark/gadgetlib1/gadgets/pairing/mnt/mnt_pairing_params.hpp"
#include "libsnark/gadgetlib1/gadgets/pairing/pairing_params.hpp"

#include <gtest/gtest.h>
#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/algebra/curves/bw6_761/bw6_761_pp.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>

using namespace libsnark;

using wpp = libff::bw6_761_pp;
using npp = other_curve<wpp>;

namespace
{

template<typename ppT>
void test_G2_checker_gadget(const std::string &annotation)
{
    protoboard<libff::Fr<ppT>> pb;
    G2_variable<ppT> g(pb, "g");
    G2_checker_gadget<ppT> g_check(pb, g, "g_check");
    g_check.generate_r1cs_constraints();

    printf("positive test\n");
    g.generate_r1cs_witness(libff::G2<other_curve<ppT>>::one());
    g_check.generate_r1cs_witness();
    assert(pb.is_satisfied());

    printf("negative test\n");
    g.generate_r1cs_witness(libff::G2<other_curve<ppT>>::zero());
    g_check.generate_r1cs_witness();
    assert(!pb.is_satisfied());

    printf(
        "number of constraints for G2 checker (Fr is %s)  = %zu\n",
        annotation.c_str(),
        pb.num_constraints());
}

TEST(TestCurveGadgets, G2Checker)
{
    test_G2_checker_gadget<libff::mnt4_pp>("mnt4");
    test_G2_checker_gadget<libff::mnt6_pp>("mnt6");
}

TEST(TestCurveGadgets, G1MulByConstScalar)
{
    // Compute inputs and results
    const libff::G1<npp> P_val = libff::Fr<npp>(13) * libff::G1<npp>::one();
    const libff::Fr<npp> scalar_val_a = libff::Fr<npp>(127);
    const libff::G1<npp> expect_result_val_a = scalar_val_a * P_val;
    const libff::Fr<npp> scalar_val_b = libff::Fr<npp>(122);
    const libff::G1<npp> expect_result_val_b = scalar_val_b * P_val;
    // Circuit
    protoboard<libff::Fr<wpp>> pb;
    G1_variable<wpp> P(pb, "P");
    G1_variable<wpp> result_a(pb, "result");
    G1_mul_by_const_scalar_gadget<wpp, libff::Fr<npp>::num_limbs> mul_gadget_a(
        pb, scalar_val_a.as_bigint(), P, result_a, "mul_gadget_a");
    G1_variable<wpp> result_b(pb, "result");
    G1_mul_by_const_scalar_gadget<wpp, libff::Fr<npp>::num_limbs> mul_gadget_b(
        pb, scalar_val_b.as_bigint(), P, result_b, "mul_gadget_b");

    mul_gadget_a.generate_r1cs_constraints();
    mul_gadget_b.generate_r1cs_constraints();

    P.generate_r1cs_witness(P_val);
    mul_gadget_a.generate_r1cs_witness();
    mul_gadget_b.generate_r1cs_witness();

    ASSERT_TRUE(pb.is_satisfied());

    const libff::G1<npp> result_a_val = g1_variable_get_element(result_a);
    ASSERT_EQ(expect_result_val_a, result_a_val);
    const libff::G1<npp> result_b_val = g1_variable_get_element(result_b);
    ASSERT_EQ(expect_result_val_b, result_b_val);
}

TEST(TestCurveGadgets, G1MulByConstScalarWithKnownResult)
{
    // Compute inputs and results
    const libff::G1<npp> P_val = libff::Fr<npp>(13) * libff::G1<npp>::one();
    const libff::G1<npp> Q_val = libff::Fr<npp>(12) * libff::G1<npp>::one();
    const libff::Fr<npp> scalar_val = libff::Fr<npp>(127);
    const libff::G1<npp> result_val = scalar_val * P_val;

    // Valid case
    {
        // Circuit
        protoboard<libff::Fr<wpp>> pb;
        G1_variable<wpp> P(pb, "P");
        G1_variable<wpp> result(pb, "result");
        G1_mul_by_const_scalar_gadget<wpp, libff::Fr<npp>::num_limbs>
            mul_gadget(pb, scalar_val.as_bigint(), P, result, "mul_gadget");

        mul_gadget.generate_r1cs_constraints();

        // Witness the input, gadget AND output
        P.generate_r1cs_witness(P_val);
        mul_gadget.generate_r1cs_witness();
        result.generate_r1cs_witness(result_val);
        ASSERT_TRUE(pb.is_satisfied());
    }

    // Invalid case. Use the gadget to ensure a specific value in the result,
    // by assigning the expected value after the gadget.
    {
        // Circuit
        protoboard<libff::Fr<wpp>> pb;
        G1_variable<wpp> P(pb, "P");
        G1_variable<wpp> result(pb, "result");
        G1_mul_by_const_scalar_gadget<wpp, libff::Fr<npp>::num_limbs>
            mul_gadget(pb, scalar_val.as_bigint(), P, result, "mul_gadget");

        mul_gadget.generate_r1cs_constraints();

        // Witness the input, gadget AND (invalid) output
        P.generate_r1cs_witness(Q_val);
        mul_gadget.generate_r1cs_witness();
        result.generate_r1cs_witness(result_val);
        ASSERT_FALSE(pb.is_satisfied());
    }
}

TEST(TestCurveGadgets, G2AddGadget)
{
    // Compute inputs and results
    const libff::G2<npp> A_val = libff::Fr<npp>(13) * libff::G2<npp>::one();
    const libff::G2<npp> B_val = libff::Fr<npp>(12) * libff::G2<npp>::one();
    const libff::G2<npp> expect_C_val =
        libff::Fr<npp>(12 + 13) * libff::G2<npp>::one();
    ASSERT_EQ(expect_C_val, A_val + B_val);

    protoboard<libff::Fr<wpp>> pb;
    G2_variable<wpp> A(pb, "A");
    G2_variable<wpp> B(pb, "B");
    G2_variable<wpp> C(pb, "C");
    G2_add_gadget<wpp> add_gadget(pb, A, B, C, "add_gadget");

    add_gadget.generate_r1cs_constraints();

    A.generate_r1cs_witness(A_val);
    B.generate_r1cs_witness(B_val);
    add_gadget.generate_r1cs_witness();

    const libff::G2<npp> C_val = g2_variable_get_element(C);
    ASSERT_TRUE(pb.is_satisfied());
    ASSERT_EQ(expect_C_val, C_val);
}

TEST(TestCurveGadgets, G2DblGadget)
{
    // Compute inputs and results
    const libff::G2<npp> A_val = libff::Fr<npp>(13) * libff::G2<npp>::one();
    const libff::G2<npp> expect_B_val =
        libff::Fr<npp>(13 + 13) * libff::G2<npp>::one();
    ASSERT_EQ(A_val.dbl(), expect_B_val);

    protoboard<libff::Fr<wpp>> pb;
    G2_variable<wpp> A(pb, "A");
    G2_variable<wpp> B(pb, "B");
    G2_dbl_gadget<wpp> dbl_gadget(pb, A, B, "dbl_gadget");

    dbl_gadget.generate_r1cs_constraints();

    A.generate_r1cs_witness(A_val);
    dbl_gadget.generate_r1cs_witness();

    const libff::G2<npp> B_val = g2_variable_get_element(B);
    ASSERT_TRUE(pb.is_satisfied());
    ASSERT_EQ(expect_B_val, B_val);
}

TEST(TestCurveGadgets, G2MulByConstScalar)
{
    // Compute inputs and results
    const libff::G2<npp> P_val = libff::Fr<npp>(13) * libff::G2<npp>::one();
    const libff::Fr<npp> scalar_val = libff::Fr<npp>(127);
    const libff::G2<npp> expect_result_val = scalar_val * P_val;

    // Circuit
    protoboard<libff::Fr<wpp>> pb;
    G2_variable<wpp> P(pb, "P");
    G2_variable<wpp> result(pb, "result");
    G2_mul_by_const_scalar_gadget<wpp, libff::Fr<npp>::num_limbs> mul_gadget(
        pb, scalar_val.as_bigint(), P, result, "mul_gadget");

    mul_gadget.generate_r1cs_constraints();

    P.generate_r1cs_witness(P_val);
    mul_gadget.generate_r1cs_witness();

    ASSERT_TRUE(pb.is_satisfied());

    const libff::G2<npp> result_val = g2_variable_get_element(result);
    ASSERT_EQ(expect_result_val, result_val);
}

TEST(TestCurveGadgets, G2MulByConstScalarWithKnownResult)
{
    // Compute inputs and results
    const libff::G2<npp> P_val = libff::Fr<npp>(13) * libff::G2<npp>::one();
    const libff::G2<npp> Q_val = libff::Fr<npp>(12) * libff::G2<npp>::one();
    const libff::Fr<npp> scalar_val = libff::Fr<npp>(127);
    const libff::G2<npp> result_val = scalar_val * P_val;

    // Valid case
    {
        // Circuit
        protoboard<libff::Fr<wpp>> pb;
        G2_variable<wpp> P(pb, "P");
        G2_variable<wpp> result(pb, "result");
        G2_mul_by_const_scalar_gadget<wpp, libff::Fr<npp>::num_limbs>
            mul_gadget(pb, scalar_val.as_bigint(), P, result, "mul_gadget");

        mul_gadget.generate_r1cs_constraints();

        // Witness the input and output
        result.generate_r1cs_witness(result_val);
        P.generate_r1cs_witness(P_val);
        mul_gadget.generate_r1cs_witness();
        result.generate_r1cs_witness(result_val);
        ASSERT_TRUE(pb.is_satisfied());
    }

    // Invalid case
    {
        // Circuit
        protoboard<libff::Fr<wpp>> pb;
        G2_variable<wpp> P(pb, "P");
        G2_variable<wpp> result(pb, "result");
        G2_mul_by_const_scalar_gadget<wpp, libff::Fr<npp>::num_limbs>
            mul_gadget(pb, scalar_val.as_bigint(), P, result, "mul_gadget");

        mul_gadget.generate_r1cs_constraints();

        // Witness the input and output
        result.generate_r1cs_witness(result_val);
        P.generate_r1cs_witness(Q_val);
        mul_gadget.generate_r1cs_witness();
        result.generate_r1cs_witness(result_val);
        ASSERT_FALSE(pb.is_satisfied());
    }
}

TEST(TestCurveGadgets, G2EqualityGadget)
{
    // Compute inputs and results
    const libff::G2<npp> P_val = libff::Fr<npp>(13) * libff::G2<npp>::one();
    const libff::G2<npp> Q_val = libff::Fr<npp>(12) * libff::G2<npp>::one();

    // Circuit
    protoboard<libff::Fr<wpp>> pb;
    G2_variable<wpp> P(pb, "P");
    G2_variable<wpp> Q(pb, "Q");
    G2_equality_gadget<wpp> equality_gadget(pb, P, Q, "equality_gadget");

    equality_gadget.generate_r1cs_constraints();

    // P == Q case
    P.generate_r1cs_witness(P_val);
    Q.generate_r1cs_witness(P_val);
    equality_gadget.generate_r1cs_witness();
    ASSERT_TRUE(pb.is_satisfied());

    // P != Q case
    P.generate_r1cs_witness(P_val);
    Q.generate_r1cs_witness(Q_val);
    equality_gadget.generate_r1cs_witness();
    ASSERT_FALSE(pb.is_satisfied());
}

} // namespace

int main(int argc, char **argv)
{
    libff::bls12_377_pp::init_public_params();
    libff::bw6_761_pp::init_public_params();
    libff::mnt4_pp::init_public_params();
    libff::mnt6_pp::init_public_params();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
