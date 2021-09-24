/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_CURVE_SCALAR_MULTIPLICATION_TCC_
#define LIBSNARK_GADGETLIB1_GADGETS_CURVE_SCALAR_MULTIPLICATION_TCC_

#include "libsnark/gadgetlib1/gadgets/curves/scalar_multiplication.hpp"

namespace libsnark
{

template<
    typename groupT,
    typename groupVariableT,
    typename add_gadget,
    typename dbl_gadget,
    typename scalarT>
point_mul_by_const_scalar_gadget<
    groupT,
    groupVariableT,
    add_gadget,
    dbl_gadget,
    scalarT>::
    point_mul_by_const_scalar_gadget(
        protoboard<FieldT> &pb,
        const scalarT &scalar,
        const groupVariableT &P,
        const groupVariableT &result,
        const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix), _scalar(scalar), _result(result)
{
    const size_t last_bit = _scalar.num_bits() - 1;
    const groupVariableT *last_value = &P;

    // Temporary vector of intermediate variables. Reserve the maximum number
    // of possible entries to ensure no reallocation (i.e. last_value is always
    // valid).
    std::vector<groupVariableT> values;
    values.reserve(2 * last_bit);

    for (size_t i = last_bit - 1; i > 0; --i) {
        // Double
        values.emplace_back(pb, FMT(annotation_prefix, " value[%zu]", i));
        _dbl_gadgets.emplace_back(new dbl_gadget(
            pb,
            *last_value,
            values.back(),
            FMT(annotation_prefix, " double[%zu]", i)));
        last_value = &values.back();

        // Add
        if (_scalar.test_bit(i)) {
            values.emplace_back(pb, FMT(annotation_prefix, " value[%zu]", i));
            _add_gadgets.emplace_back(new add_gadget(
                pb,
                *last_value,
                P,
                values.back(),
                FMT(annotation_prefix, " add[%zu]", i)));
            last_value = &values.back();
        }
    }

    // Depending on the value of the final (lowest-order) bit, perform final
    // double or double-and-add into result.

    if (_scalar.test_bit(0)) {
        // Double
        values.emplace_back(pb, FMT(annotation_prefix, " value[0]"));
        _dbl_gadgets.emplace_back(new dbl_gadget(
            pb,
            *last_value,
            values.back(),
            FMT(annotation_prefix, " double[0]")));
        last_value = &values.back();

        // Add into result
        _add_gadgets.emplace_back(new add_gadget(
            pb, *last_value, P, result, FMT(annotation_prefix, " add[0]")));
    } else {
        // Double
        _dbl_gadgets.emplace_back(new dbl_gadget(
            pb, *last_value, result, FMT(annotation_prefix, " double[0]")));
    }
}

template<
    typename groupT,
    typename groupVariableT,
    typename add_gadget,
    typename dbl_gadget,
    typename scalarT>
void point_mul_by_const_scalar_gadget<
    groupT,
    groupVariableT,
    add_gadget,
    dbl_gadget,
    scalarT>::generate_r1cs_constraints()
{
    const size_t last_bit = _scalar.num_bits() - 1;
    size_t dbl_idx = 0;
    size_t add_idx = 0;
    for (ssize_t i = last_bit - 1; i >= 0; --i) {
        // Double gadget constraints
        _dbl_gadgets[dbl_idx++]->generate_r1cs_constraints();

        // Add gadget constraints
        if (_scalar.test_bit(i)) {
            _add_gadgets[add_idx++]->generate_r1cs_constraints();
        }
    }
}

template<
    typename groupT,
    typename groupVariableT,
    typename add_gadget,
    typename dbl_gadget,
    typename scalarT>
void point_mul_by_const_scalar_gadget<
    groupT,
    groupVariableT,
    add_gadget,
    dbl_gadget,
    scalarT>::generate_r1cs_witness()
{
    const size_t last_bit = _scalar.num_bits() - 1;
    size_t dbl_idx = 0;
    size_t add_idx = 0;
    for (ssize_t i = last_bit - 1; i >= 0; --i) {
        // Double gadget constraints
        _dbl_gadgets[dbl_idx++]->generate_r1cs_witness();

        // Add gadget constraints
        if (_scalar.test_bit(i)) {
            _add_gadgets[add_idx++]->generate_r1cs_witness();
        }
    }
}

template<
    typename groupT,
    typename groupVariableT,
    typename add_gadget,
    typename dbl_gadget,
    typename scalarT>
const groupVariableT &point_mul_by_const_scalar_gadget<
    groupT,
    groupVariableT,
    add_gadget,
    dbl_gadget,
    scalarT>::result() const
{
    return _result;
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_CURVE_SCALAR_MULTIPLICATION_TCC_
