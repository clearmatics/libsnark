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
#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.hpp>

using namespace libsnark;

template<typename FieldT> void test_anemoi_power_three_gadget(const size_t n)
{
    printf("testing anemoi_power_three_gadget on all %zu bit strings\n", n);

    protoboard<FieldT> pb;
    pb_variable_array<FieldT> inputs;
    inputs.allocate(pb, n, "inputs");

    pb_variable<FieldT> output;
    output.allocate(pb, "output");

    anemoi_power_three_gadget<FieldT> d(pb, inputs, output, "d");
    d.generate_r1cs_constraints();

    for (size_t w = 0; w < 1ul << n; ++w) {
        for (size_t j = 0; j < n; ++j) {
            pb.val(inputs[j]) = FieldT((w & (1ul << j)) ? 1 : 0);
        }

        d.generate_r1cs_witness();

#ifdef DEBUG
        printf("positive test for %zu\n", w);
#endif
        ASSERT_EQ(pb.val(output), (w ? FieldT::one() : FieldT::zero()));
        ASSERT_TRUE(pb.is_satisfied());

#ifdef DEBUG
        printf("negative test for %zu\n", w);
#endif
        pb.val(output) = (w ? FieldT::zero() : FieldT::one());
        ASSERT_FALSE(pb.is_satisfied());
    }

    libff::print_time("anemoi_power_three_gadget tests successful");
}

int main(void)
{
    libff::start_profiling();
    libff::default_ec_pp::init_public_params();
    test_anemoi_power_three_gadget<libff::Fr<libff::default_ec_pp>>(10);
}
