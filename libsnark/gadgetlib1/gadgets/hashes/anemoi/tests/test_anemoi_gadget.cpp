/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/tests/anemoi_outputs.hpp"

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

class parameters_debug_bls12_381
{
public:
    using ppT = libff::bls12_381_pp;
    using FieldT = libff::Fr<ppT>;
    using BignumT = libff::bigint<FieldT::num_limbs>;
    static const bool b_prime_field = false;
    static constexpr size_t multiplicative_generator_g = 7;
    static constexpr size_t alpha = 5;
    static constexpr size_t beta = 2;
    static constexpr size_t gamma = 5;
    static constexpr size_t quad_exponent = 2;
    static const BignumT alpha_inv;
    static const BignumT delta;
};

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    parameters_debug_bls12_381::alpha_inv =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>(
            "209743500700504761917790962032743863350762210002110551290414634799"
            "75432473805");

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    parameters_debug_bls12_381::delta =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>("0");

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

template<typename ppT, class parameters>
void test_flystel_Q_gamma_prime_field_gadget()
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
    flystel_Q_prime_field_gadget<ppT> d(
        pb, parameters::beta, parameters::gamma, x, y, "d");

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

template<typename ppT, class parameters = anemoi_parameters<libff::Fr<ppT>>>
void test_flystel_Q_gamma_binary_field_gadget()
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
    flystel_Q_binary_field_gadget<ppT> d(
        pb, parameters::beta, parameters::gamma, x, y, "d");
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
    flystel_E_power_five_gadget<ppT> d(pb, x, y, "d");
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

template<typename ppT, class parameters = anemoi_parameters<libff::Fr<ppT>>>
void test_flystel_E_root_five_gadget()
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
    flystel_E_root_five_gadget<ppT, parameters> d(pb, x, y, "d");
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

template<typename ppT, class parameters = anemoi_parameters<libff::Fr<ppT>>>
void test_flystel_prime_field_gadget()
{
    using FieldT = libff::Fr<ppT>;

    protoboard<FieldT> pb;

    pb_variable<FieldT> x0 = pb_variable_allocate(pb, "x0");
    pb_variable<FieldT> x1 = pb_variable_allocate(pb, "x1");
    pb_variable<FieldT> y0 = pb_variable_allocate(pb, "y0");
    pb_variable<FieldT> y1 = pb_variable_allocate(pb, "y1");

    flystel_prime_field_gadget<ppT, parameters> d(
        pb, x0, x1, y0, y1, "flystel");

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

template<
    typename ppT,
    size_t NumStateColumns_L,
    class parameters = anemoi_parameters<libff::Fr<ppT>>>
void test_anemoi_permutation_round_prime_field_gadget()
{
    using FieldT = libff::Fr<ppT>;

    //        const size_t NumStateColumns_L = L;

    protoboard<FieldT> pb;
    std::vector<FieldT> C;
    std::vector<FieldT> D;

    pb_variable_array<FieldT> X_left;
    pb_variable_array<FieldT> X_right;
    pb_variable_array<FieldT> Y_left;
    pb_variable_array<FieldT> Y_right;

    X_left.allocate(pb, NumStateColumns_L, "left inputs");
    X_right.allocate(pb, NumStateColumns_L, "right inputs");

    Y_left.allocate(pb, NumStateColumns_L, "left outputs");
    Y_right.allocate(pb, NumStateColumns_L, "right outputs");

    // WARNINIG! test with zero constants. TODO add the original
    // constants of Anemoi
    for (size_t i = 0; i < NumStateColumns_L; i++) {
        C.push_back(FieldT(0));
        D.push_back(FieldT(0));
    }

    anemoi_permutation_round_prime_field_gadget<
        ppT,
        NumStateColumns_L,
        parameters>
        d(pb, C, D, X_left, X_right, Y_left, Y_right, "anemoi permutation");

    // generate constraints
    d.generate_r1cs_constraints();

    // Input values: X_left = 0,1,2...L-1 ; X_right = L, L+1, 2L-1
    for (size_t i = 0; i < NumStateColumns_L; i++) {
        pb.val(X_left[i]) = FieldT(i);
        pb.val(X_right[i]) = FieldT(NumStateColumns_L + i);
    }

    // generate witness for the given input
    d.generate_r1cs_witness();

    std::vector<FieldT> Y_expect =
        anemoi_expected_output_one_round(NumStateColumns_L);

    for (size_t i = 0; i < NumStateColumns_L; i++) {
        ASSERT_EQ(Y_expect[i], pb.val(Y_left[i]));
        ASSERT_EQ(Y_expect[NumStateColumns_L + i], pb.val(Y_right[i]));
    }

    ASSERT_TRUE(pb.is_satisfied());
    test_pb_verify_circuit<ppT>(pb);

    libff::print_time(
        "anemoi_permutation_round_prime_field_gadget tests successful");
}

template<typename ppT> void test_for_curve()
{
    // Execute all tests for the given curve.

    ppT::init_public_params();
    // Use debug parameters with small values for the small gadgets
    using parameters_debug = parameters_debug_bls12_381;
    test_flystel_Q_gamma_prime_field_gadget<ppT, parameters_debug>();
    test_flystel_Q_gamma_binary_field_gadget<ppT, parameters_debug>();
    test_flystel_E_power_five_gadget<ppT>();
    test_flystel_E_root_five_gadget<ppT, parameters_debug>();
    test_flystel_prime_field_gadget<ppT, parameters_debug>();
    // Use the original parameters for the full permutation
    using parameters = anemoi_parameters<ppT>;
    test_anemoi_permutation_round_prime_field_gadget<ppT, 1, parameters>();
    test_anemoi_permutation_round_prime_field_gadget<ppT, 2, parameters>();
    test_anemoi_permutation_round_prime_field_gadget<ppT, 3, parameters>();
    //    test_anemoi_permutation_round_prime_field_gadget<ppT, 4,
    //    parameters>();
}

TEST(TestAnemoiGadget, BLS12_381) { test_for_curve<libff::bls12_381_pp>(); }

int main(int argc, char **argv)
{
    libff::bls12_381_pp::init_public_params();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
