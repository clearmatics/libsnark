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
    // x*x = a
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(x, x, a),
        FMT(this->annotation_prefix, " x*x=a"));
    // a*x = b
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(a, x, b),
        FMT(this->annotation_prefix, " a*x=b"));
    // b*alpha = c
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(b, alpha, c),
        FMT(this->annotation_prefix, " b*alpha=c"));
    // 1*(c+beta) = y
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(1, c + beta, y),
        FMT(this->annotation_prefix, " 1*(c+beta)=y"));
}

template<typename FieldT>
void anemoi_power_three_gadget<FieldT>::generate_r1cs_witness()
{
    // x*x = a
    this->pb.val(a) = this->pb.val(x) * this->pb.val(x);
    // a*x = b
    this->pb.val(b) = this->pb.val(a) * this->pb.val(x);
    // b*alpha = c
    this->pb.val(c) = this->pb.val(b) * this->pb.val(alpha);
    // 1*(c+beta) = y
    this->pb.val(y) = this->pb.val(c) + this->pb.val(beta);
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_
