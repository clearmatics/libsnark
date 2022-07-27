/** @file
 *****************************************************************************

 Implementation of interfaces for top-level Anemoi hash function gadgets.

 See anemoi_gadget.hpp .

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_

namespace libsnark
{

template<typename FieldT>
void anemoi_power_three_gadget<FieldT>::generate_r1cs_constraints()
{
    linear_combination<FieldT> sum = pb_sum(inputs);

    // inv * sum = output
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(inv, sum, output),
        FMT(this->annotation_prefix, " inv*sum=output"));

    // (1-output) * sum = 0
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(FieldT::one() - output, sum, FieldT::zero()),
        FMT(this->annotation_prefix, " (1-output)*sum=0"));
}

template<typename FieldT>
void anemoi_power_three_gadget<FieldT>::generate_r1cs_witness()
{
    FieldT sum = FieldT::zero();

    for (size_t i = 0; i < inputs.size(); ++i) {
        sum += this->pb.lc_val(inputs[i]);
    }

    if (sum.is_zero()) {
        this->pb.val(inv) = FieldT::zero();
        this->pb.val(output) = FieldT::zero();
    } else {
        this->pb.val(inv) = sum.inverse();
        this->pb.val(output) = FieldT::one();
    }
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_
