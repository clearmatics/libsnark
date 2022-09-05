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

#define ANEMOI_BLS12_381_CONST_ALPHA 5
// equals to G1
#define ANEMOI_BLS12_381_CONST_BETA 2
// TODO: value by spec is 0
#define ANEMOI_BLS12_381_CONST_GAMMA 5
// TODO: value by spec is G1.inv()
#define ANEMOI_BLS12_381_CONST_DELTA 0

/// Compute y = const_a x^2 + const_b
/// x: input
/// y: output
/// const_a, const_b: constants
template<typename FieldT> class anemoi_power_two_gadget : public gadget<FieldT>
{
private:
    // constants
    const FieldT const_a;
    const FieldT const_b;

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

/// Compute y = const_a x^3 + const_b
/// x: input
/// y: output
/// const_a, const_b: constants
template<typename FieldT>
class anemoi_power_three_gadget : public gadget<FieldT>
{
private:
    /// internal (i.e. intermediate) variable
    pb_variable<FieldT> internal;
    /// constants
    const FieldT const_a;
    const FieldT const_b;

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
