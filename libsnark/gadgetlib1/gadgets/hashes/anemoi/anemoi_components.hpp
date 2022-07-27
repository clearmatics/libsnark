/** @file
 *****************************************************************************

 Declaration of interfaces for top-level SHA256 gadgets.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_

#include <libsnark/common/data_structures/merkle_tree.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/hash_io.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_components.hpp>

namespace libsnark
{

/*
  the gadgets below are Fp specific:
  I * X = R
  (1-R) * X = 0

  if X = 0 then R = 0
  if X != 0 then R = 1 and I = X^{-1}
*/

/// Compute x^3
/// (old) Output is 0 iff the sum of inputs is 0.  Output is 1 otherwise.
template<typename FieldT>
class anemoi_power_three_gadget : public gadget<FieldT>
{
private:
    pb_variable<FieldT> inv;

public:
    const pb_linear_combination_array<FieldT> inputs;
    const pb_variable<FieldT> output;

    anemoi_power_three_gadget(
        protoboard<FieldT> &pb,
        const pb_linear_combination_array<FieldT> &inputs,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "")
        : gadget<FieldT>(pb, annotation_prefix), inputs(inputs), output(output)
    {
        assert(inputs.size() >= 1);
        inv.allocate(pb, FMT(this->annotation_prefix, " inv"));
    }

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

} // namespace libsnark

#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.tcc>

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
