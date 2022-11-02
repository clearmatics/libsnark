/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
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

template<typename ppT> void test_flystel_Q_gamma_prime_field_gadget()
{
    using FieldT = libff::Fr<ppT>;
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
    test_pb_verify_circuit<ppT>(pb);

    libff::print_time("flystel_power_two_gadget tests successful");
}

template<typename ppT> void test_flystel_Q_gamma_binary_field_gadge()
{
    using FieldT = libff::Fr<ppT>;

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
    test_pb_verify_circuit<ppT>(pb);

    libff::print_time("flystel_power_three_gadget tests successful");
}

template<typename ppT> void test_flystel_E_power_five_gadget()
{
    using FieldT = libff::Fr<ppT>;

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
    test_pb_verify_circuit<ppT>(pb);

    libff::print_time("flystel_E_power_five_gadget tests successful");
}

template<typename ppT> void test_flystel_E_root_five_gadget()
{
    using FieldT = libff::Fr<ppT>;

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
    test_pb_verify_circuit<ppT>(pb);

    libff::print_time("flystel_E_root_five_gadget tests successful");
}

template<typename ppT> void test_flystel_prime_field_gadget()
{
    using FieldT = libff::Fr<ppT>;

    protoboard<FieldT> pb;

    pb_variable<FieldT> x0 = pb_variable_allocate(pb, "x0");
    pb_variable<FieldT> x1 = pb_variable_allocate(pb, "x1");
    pb_variable<FieldT> y0 = pb_variable_allocate(pb, "y0");
    pb_variable<FieldT> y1 = pb_variable_allocate(pb, "y1");

    flystel_prime_field_gadget<
        FieldT,
        FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR>
        d(pb, x0, x1, y0, y1, "flystel");

    // generate constraints
    d.generate_r1cs_constraints();

    pb.val(x0) = FieldT(55);
    pb.val(x1) = FieldT(3);

    // generate witness for the given input
    d.generate_r1cs_witness();

    FieldT y0_expect = FieldT(34);
    FieldT y1_expect = FieldT(1);

    ASSERT_EQ(y0_expect, pb.val(y0));
    ASSERT_EQ(y1_expect, pb.val(y1));
    ASSERT_TRUE(pb.is_satisfied());

    test_pb_verify_circuit<ppT>(pb);

    libff::print_time("flystel_prime_field_gadget tests successful");
}

int main()
{
    libff::start_profiling();

    //    libff::default_ec_pp::init_public_params();
    //    using FieldT = libff::Fr<libff::default_ec_pp>;

    libff::bls12_381_pp::init_public_params();
    using ppT = libff::bls12_381_pp;

    test_flystel_Q_gamma_prime_field_gadget<ppT>();
    test_flystel_Q_gamma_binary_field_gadge<ppT>();
    test_flystel_E_power_five_gadget<ppT>();
    test_flystel_E_root_five_gadget<ppT>();
    test_flystel_prime_field_gadget<ppT>();
    return 0;
}
