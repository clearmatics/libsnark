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
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

using namespace libsnark;

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
    const linear_combination<FieldT> x0 = 55;
    const linear_combination<FieldT> x1 = 3;
    linear_combination<FieldT> y0;
    linear_combination<FieldT> y1;

    FieldT x0_val = x0.terms[0].coeff;
    FieldT x1_val = x1.terms[0].coeff;

    flystel_prime_field_gadget<
        FieldT,
        FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR>
        d(pb, x0, x1, y0, y1, "flystel");

    // generate constraints
    d.generate_r1cs_constraints();

    // generate witness for the given input
    d.generate_r1cs_witness();

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

    std::vector<FieldT> y0_assignment({x0_val, -a0_expected, a2_expected});
    std::vector<FieldT> y1_assignment({x1_val, -a1_expected});
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

template<typename FieldT> void test_bug()
{
    protoboard<FieldT> pb;
    //    pb_variable<FieldT> x1;
    linear_combination<FieldT> x1_lc;
    pb_linear_combination<FieldT> x1(pb, x1_lc);
    pb_variable<FieldT> a0;

    // input
    //    x1.allocate(pb, "x1");
    // output
    a0.allocate(pb, "a0");

    // create gadget
    flystel_Q_gamma_prime_field_gadget<
        FieldT,
        FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR>
        d(pb, x1, a0, "d");
    // generate contraints
    d.generate_r1cs_constraints();
    // set input value
    pb.lc_val(x1) = 3;
    // generate witness for the given input
    d.generate_r1cs_witness();

    // the expected output is 13 for input 2
    ASSERT_EQ(pb.val(a0), 23);
    ASSERT_TRUE(pb.is_satisfied());
}

template<typename FieldT> void test_bug_two()
{
    protoboard<FieldT> pb;
    pb_variable<FieldT> a0; // <-- allocate
    linear_combination<FieldT> x1 = 3;
    pb_linear_combination<FieldT> x1_lc(pb, x1); // <--- use assign

    // assert(!a)

    assert(x1_lc.is_variable == false);

    FieldT x1_lc_val = pb.lc_val(x1_lc);

    printf("[%s:%d] x1_lc print\n", __FILE__, __LINE__);
    //    pb.lc_val(x1_lc).print();
    x1_lc_val.print();

#if 1
    // create gadget
    flystel_Q_gamma_prime_field_gadget<
        FieldT,
        FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR>
        d(pb, x1_lc, a0, "d");

    // generate contraints
    d.generate_r1cs_constraints();

    // generate witness for the given input
    d.generate_r1cs_witness();
#endif
}

int main(void)
{
    libff::start_profiling();

    //    libff::default_ec_pp::init_public_params();
    //    using FieldT = libff::Fr<libff::default_ec_pp>;

    libff::bls12_381_pp::init_public_params();
    using FieldT = libff::Fr<libff::bls12_381_pp>;

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
    test_flystel_prime_field_gadget<FieldT>(10);
#endif
    test_bug<FieldT>();
    //    test_root_five<FieldT>();
}
