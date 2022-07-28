/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <gtest/gtest.h>
#include <libff/common/default_types/ec_pp.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp> // VV
#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp> // VV

using namespace libsnark;

template<typename FieldT> void test_anemoi_power_three_gadget(const size_t n)
{
    printf("testing anemoi_power_three_gadget on all %zu bit strings\n", n);

    protoboard<FieldT> pb;
    pb_variable<FieldT> x;
    pb_variable<FieldT> y;
    pb_variable<FieldT> alpha;
    pb_variable<FieldT> beta;

    // input
    x.allocate(pb, "x");
    // output
    y.allocate(pb, "y");
    // constants
    alpha.allocate(pb, "alpha");
    beta.allocate(pb, "beta");

    // create gadget
    anemoi_power_three_gadget<FieldT> d(pb, x, y, alpha, beta, "d");
    // generate contraints
    d.generate_r1cs_constraints();
    // witness values
    pb.val(x) = 2;
    pb.val(alpha) = 2;
    pb.val(beta) = 5;

    // generate witness
    d.generate_r1cs_witness();

    ASSERT_EQ(pb.val(y), 21);
    ASSERT_TRUE(pb.is_satisfied());

    libff::print_time("anemoi_power_three_gadget tests successful");
}

int main(void)
{
    libff::start_profiling();
    libff::default_ec_pp::init_public_params();
    test_anemoi_power_three_gadget<libff::Fr<libff::default_ec_pp>>(10);
}
