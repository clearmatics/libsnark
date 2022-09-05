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

/// Compute y = alpha x^2 + beta
/// x: input
/// y: output
/// alpha, beta: constants
template<typename FieldT> class anemoi_power_two_gadget : public gadget<FieldT>
{
private:
    // constants
    FieldT alpha;
    FieldT beta;

public:
    // input/output
    const pb_variable<FieldT> input;
    const pb_variable<FieldT> output;

    anemoi_power_two_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Compute y = alpha x^3 + beta
/// x: input
/// y: output
/// alpha, beta: constants
template<typename FieldT>
class anemoi_power_three_gadget : public gadget<FieldT>
{
private:
    /// internal (i.e. intermediate) variable
    pb_variable<FieldT> internal;
    /// constants
    FieldT alpha;
    FieldT beta;

public:
    /// input/output
    const pb_variable<FieldT> input;
    const pb_variable<FieldT> output;

    anemoi_power_three_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

} // namespace libsnark

#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.tcc>

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
