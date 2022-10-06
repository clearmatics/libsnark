/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <gtest/gtest.h>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/common/default_types/ec_pp.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_constants.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>
//#include
//<libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

using namespace libsnark;

template<typename ppT>
void test_pb_verify_circuit(protoboard<libff::Fr<ppT>> &pb)
{
    ASSERT_TRUE(pb.is_satisfied());
    const r1cs_gg_ppzksnark_keypair<ppT> keypair =
        r1cs_gg_ppzksnark_generator<ppT>(pb.get_constraint_system(), true);
    r1cs_primary_input<libff::Fr<ppT>> primary_input = pb.primary_input();
    r1cs_auxiliary_input<libff::Fr<ppT>> auxiliary_input = pb.auxiliary_input();
    r1cs_gg_ppzksnark_proof<ppT> proof = r1cs_gg_ppzksnark_prover(
        keypair.pk, primary_input, auxiliary_input, true);
    ASSERT_TRUE(r1cs_gg_ppzksnark_verifier_strong_IC<ppT>(
        keypair.vk, primary_input, proof));
}

template<typename FieldT>
void test_flystel_Q_gamma_prime_field_gadget(const size_t n)
{
    printf("testing flystel_power_two_gadget on all %zu bit strings\n", n);
    protoboard<FieldT> pb;
    pb_variable<FieldT> x;
    pb_variable<FieldT> y;

    // input
    x.allocate(pb, "x");
    // output
    y.allocate(pb, "y");

    // create gadget
    flystel_Q_gamma_prime_field_gadget<
        FieldT,
        FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR>
        d(pb, x, y, "d");
    // generate contraints
    d.generate_r1cs_constraints();
    // set input value
    pb.val(x) = 2;
    // generate witness for the given input
    d.generate_r1cs_witness();

    // the expected output is 13 for input 2
    ASSERT_EQ(pb.val(y), 13);
    ASSERT_TRUE(pb.is_satisfied());
    libff::print_time("flystel_power_two_gadget tests successful");
}

template<typename FieldT>
void test_flystel_Q_gamma_binary_field_gadge(const size_t n)
{
    printf("testing flystel_power_three_gadget on all %zu bit strings\n", n);

    protoboard<FieldT> pb;
    pb_variable<FieldT> x;
    pb_variable<FieldT> y;

    // input
    x.allocate(pb, "x");
    // output
    y.allocate(pb, "y");

    // create gadget
    flystel_Q_gamma_binary_field_gadget<
        FieldT,
        FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR>
        d(pb, x, y, "d");
    // generate contraints
    d.generate_r1cs_constraints();
    // set input value
    pb.val(x) = 2;
    // generate witness for the given input
    d.generate_r1cs_witness();

    // the expected output is 21 for input 2
    ASSERT_EQ(pb.val(y), 21);
    ASSERT_TRUE(pb.is_satisfied());

    libff::print_time("flystel_power_three_gadget tests successful");
}

template<typename FieldT> void test_flystel_E_power_five_gadget(const size_t n)
{
    printf("testing flystel_E_power_five_gadget on all %zu bit strings\n", n);

    protoboard<FieldT> pb;
    pb_variable<FieldT> x;
    pb_variable<FieldT> y;

    // input
    x.allocate(pb, "x");
    // output
    y.allocate(pb, "y");

    // create gadget
    flystel_E_power_five_gadget<FieldT> d(pb, x, y, "d");
    // generate contraints
    d.generate_r1cs_constraints();
    // set input value
    pb.val(x) = 2;
    // generate witness for the given input
    d.generate_r1cs_witness();

    // the expected output is 32 for input 2
    ASSERT_EQ(pb.val(y), 32);
    ASSERT_TRUE(pb.is_satisfied());

    libff::print_time("flystel_E_power_five_gadget tests successful");
}

template<typename FieldT> void test_flystel_E_root_five_gadget(const size_t n)
{
    printf("testing flystel_E_root_five_gadget on all %zu bit strings\n", n);

    protoboard<FieldT> pb;
    pb_variable<FieldT> x;
    pb_variable<FieldT> y;

    // input
    x.allocate(pb, "x");
    // output
    y.allocate(pb, "y");

    // create gadget
    flystel_E_root_five_gadget<FieldT> d(pb, x, y, "d");
    // generate contraints
    d.generate_r1cs_constraints();
    // set input value
    pb.val(x) = 22;
    // generate witness for the given input
    d.generate_r1cs_witness();

    // computed using Sage
    FieldT y_expected = FieldT("10357913779704000956629425810748166374506105653"
                               "828973721142406533896278368512");

    // the expected output is 32 for input 2
    ASSERT_EQ(pb.val(y), y_expected);
    ASSERT_TRUE(pb.is_satisfied());

    libff::print_time("flystel_E_root_five_gadget tests successful");
}

template<typename FieldT> void test_flystel_prime_field_gadget(const size_t n)
{
    printf("testing flystel_prime_field_gadget on all %zu bit strings\n", n);

    protoboard<FieldT> pb;

    // input
#if 0    
    pb_variable<FieldT> x0;
    pb_variable<FieldT> x1;
    x0.allocate(pb, "x0");
    x1.allocate(pb, "x1");
    // output
    pb_variable<FieldT> y0;
    pb_variable<FieldT> y1;
    y0.allocate(pb, "y0");
    y1.allocate(pb, "y1");
#endif
    linear_combination<FieldT> x0_lc;
    linear_combination<FieldT> x1_lc;

    pb_linear_combination<FieldT> x0(pb, x0_lc);
    pb_linear_combination<FieldT> x1(pb, x1_lc);

    linear_combination<FieldT> y0;
    linear_combination<FieldT> y1;

    flystel_prime_field_gadget<
        FieldT,
        FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR>
        d(pb, x0, x1, y0, y1, "flystel");

    // generate constraints
    d.generate_r1cs_constraints();

    //    const linear_combination<FieldT> x0 = 55;
    //    const linear_combination<FieldT> x1 = 3;
    //    pb.lc_val(x1_lc) = 3;
    pb.lc_val(x0) = 55;
    pb.lc_val(x1) = 3;

    // generate witness for the given input
    d.generate_r1cs_witness();

    FieldT x0_val = pb.lc_val(x0); // x0_lc.terms[0].coeff;
    FieldT x1_val = pb.lc_val(x1); // x1_lc.terms[0].coeff;

    // a0 = 23
    FieldT a0_expected = FieldT(23);
    // a1 = 22^{1/5}
    FieldT a1_expected =
        FieldT("10357913779704000956629425810748166374506105653"
               "828973721142406533896278368512");
    // a2 = 2 (3-a1)^2
    FieldT a2_expected =
        FieldT(2) * (FieldT(3) - a1_expected) * (FieldT(3) - a1_expected);
    // y0 = x0 - a0 + a2 = 22 + a2
    FieldT y0_expected = x0_val - a0_expected + a2_expected;
    // y1 = x1 - a1 = 3 - a1
    FieldT y1_expected = x1_val - a1_expected;

    //    std::vector<linear_term<FieldT>> terms;
    std::vector<FieldT> y0_assignment({x0_val, -a0_expected, a2_expected});
    std::vector<FieldT> y1_assignment({x1_val, -a1_expected});

    FieldT t1 = y0.evaluate(y0_assignment);
    printf("[%s:%d] y0.evaluate ", __FILE__, __LINE__);
    t1.print();
    printf("[%s:%d] y0_expected ", __FILE__, __LINE__);
    y0_expected.print();

    ASSERT_EQ(y0.evaluate(y0_assignment), y0_expected);
    ASSERT_EQ(y1.evaluate(y1_assignment), y1_expected);
    ASSERT_TRUE(pb.is_satisfied());
    libff::print_time("flystel_prime_field_gadget tests successful");
}

template<typename FieldT> void test_root_five()
{
    // alpha_inv =
    // 20974350070050476191779096203274386335076221000211055129041463479975432473805
    //    FieldT x = FieldT::random_element();
    //    FieldT y = power(x, 5);
    //    x.print();
    //    y.print();
    FieldT x = 5;
    FieldT x_mod_inv =
        FieldT("2097435007005047619177909620327438633507622100021"
               "1055129041463479975432473805");
    printf("Fr modulus   \n");
    x.mod.print();
    printf("x + x_mod_inv\n");
    FieldT z = x + x_mod_inv;
    z.print();
    printf("\n");
    x.print();
    x.inverse().print();
}

template<typename ppT> void test_bug()
{
    using FieldT = libff::Fr<ppT>;

    protoboard<FieldT> pb;
    pb_variable<FieldT> v1 = pb_variable_allocate(pb, "v1");
    pb_variable<FieldT> v2 = pb_variable_allocate(pb, "v2");
    pb_variable<FieldT> a0 = pb_variable_allocate(pb, "a0");
    pb_linear_combination<FieldT> x1;

    x1.assign(pb, v1 + v2);

    flystel_Q_gamma_prime_field_gadget<
        FieldT,
        FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR>
        d(pb, x1, a0, "flystel_Q_gamma");
    d.generate_r1cs_constraints();

    pb.val(v1) = FieldT(3);
    pb.val(v2) = FieldT(0);

    const FieldT expect_a0("23");

    d.generate_r1cs_witness();
    ASSERT_EQ(expect_a0, pb.val(a0));
    ASSERT_TRUE(pb.is_satisfied());

    //    test_pb_verify_circuit<ppT>(pb);
}

template<typename ppT> void test_bug_dt()
{
    using FieldT = libff::Fr<ppT>;

    // Circuit showing x_3 = beta * (x_1+x_2)^2 + gamma
    FieldT x1 = FieldT(7);
    FieldT x2 = FieldT(11);
    linear_combination<FieldT> lc(x1 + x2);

    protoboard<FieldT> pb;
    pb_variable<FieldT> x3 = pb_variable_allocate(pb, "x3");
    pb_linear_combination<FieldT> pb_lc; //(pb, lc);
    pb_lc.assign(pb, lc);

    flystel_Q_gamma_prime_field_gadget<FieldT, 2> d(
        pb, pb_lc, x3, "flystel_Q_gamma");
    d.generate_r1cs_constraints();

    // Expect x3 = 2 * (7+11)^2 + 5 = 653
    const FieldT expect_x3("653");

    d.generate_r1cs_witness();
    ASSERT_EQ(expect_x3, pb.val(x3));
    ASSERT_TRUE(pb.is_satisfied());

    //    test_pb_verify_circuit<ppT>(pb);
}

TEST(TestAnemoiGadget, TestBug) { test_bug<libff::bls12_381_pp>(); }
TEST(TestAnemoiGadget, TestBugDt) { test_bug_dt<libff::bls12_381_pp>(); }

int main(int argc, char **argv)
{
    libff::start_profiling();

    //    libff::default_ec_pp::init_public_params();
    //    using FieldT = libff::Fr<libff::default_ec_pp>;

    libff::bls12_381_pp::init_public_params();
    using ppT = libff::bls12_381_pp;
    using FieldT = libff::Fr<ppT>;

    // for BLS12-381
    // beta = g = first multiplicative generator = 7.
    // delta = g^(-1)
    // 14981678621464625851270783002338847382197300714436467949315331057125308909861
    // Fr modulus
    // 52435875175126190479447740508185965837690552500527637822603658699938581184513
#if 0
    FieldT a = FieldT(7);
    FieldT a_inv = a.inverse();
    assert((a * a_inv) == FieldT::one());
    printf("a_inv      ");
    a_inv.print();
    printf("\n");
    printf("Fr modulus ");
    a.mod.print();
    printf("\n");
#endif
#if 0    
    test_flystel_Q_gamma_prime_field_gadget<FieldT>(10);
    test_flystel_Q_gamma_binary_field_gadge<FieldT>(10);
    test_flystel_E_power_five_gadget<FieldT>(10);
    test_flystel_E_root_five_gadget<FieldT>(10);
#endif
    //    test_flystel_prime_field_gadget<FieldT>(10);
    //    test_bug<ppT>();
    test_bug_dt<ppT>();
    //    test_bug_two<FieldT>();
    //    test_bug_one<FieldT>();
    //    test_root_five<FieldT>();
    //    ::testing::InitGoogleTest(&argc, argv);
    //    return RUN_ALL_TESTS();
    return 0;
}
