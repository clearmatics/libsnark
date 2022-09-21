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
#if 1
    test_flystel_Q_gamma_prime_field_gadget<FieldT>(10);
    test_flystel_Q_gamma_binary_field_gadge<FieldT>(10);
    test_flystel_E_power_five_gadget<FieldT>(10);
#endif
    //    //    test_flystel_power_two_gadget<libff::bls12_381_Fr>(10);
}
