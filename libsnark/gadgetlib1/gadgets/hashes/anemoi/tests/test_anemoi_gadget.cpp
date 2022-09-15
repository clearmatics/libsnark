/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <gtest/gtest.h>
//#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/common/default_types/ec_pp.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp> // VV
#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp> // VV

using namespace libsnark;

template<typename FieldT> void test_flystel_power_two_gadget(const size_t n)
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
    flystel_power_two_gadget<FieldT> d(pb, x, y, "d");
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

template<typename FieldT> void test_flystel_power_three_gadget(const size_t n)
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
    flystel_power_three_gadget<FieldT> d(pb, x, y, "d");
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

template<typename FieldT> void test_flystel_power_five_gadget(const size_t n)
{
    printf("testing flystel_power_five_gadget on all %zu bit strings\n", n);

    protoboard<FieldT> pb;
    pb_variable<FieldT> x;
    pb_variable<FieldT> y;

    // input
    x.allocate(pb, "x");
    // output
    y.allocate(pb, "y");

    // create gadget
    flystel_power_five_gadget<FieldT> d(pb, x, y, "d");
    // generate contraints
    d.generate_r1cs_constraints();
    // set input value
    pb.val(x) = 2;
    // generate witness for the given input
    d.generate_r1cs_witness();

    // the expected output is 32 for input 2
    ASSERT_EQ(pb.val(y), 32);
    ASSERT_TRUE(pb.is_satisfied());

    libff::print_time("flystel_power_five_gadget tests successful");
}

int main(void)
{
    libff::start_profiling();
    libff::default_ec_pp::init_public_params();
    test_flystel_power_two_gadget<libff::Fr<libff::default_ec_pp>>(10);
    test_flystel_power_three_gadget<libff::Fr<libff::default_ec_pp>>(10);
    test_flystel_power_five_gadget<libff::Fr<libff::default_ec_pp>>(10);
    //    test_flystel_power_two_gadget<libff::bls12_381_Fr>(10);
}
