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

template<typename ppT, typename groupT, typename groupVariableT>
variable_or_identity<ppT, groupT, groupVariableT>::variable_or_identity(
    protoboard<FieldT> &pb, const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix), value(pb, annotation_prefix)
{
    is_identity_var.allocate(pb, " is_identity");
    is_identity = pb_linear_combination<FieldT>(is_identity_var);
    generate_boolean_r1cs_constraint(
        pb, is_identity, FMT(annotation_prefix, " is_identity_is_bool"));
}

template<typename ppT, typename groupT, typename groupVariableT>
variable_or_identity<ppT, groupT, groupVariableT>::variable_or_identity(
    protoboard<FieldT> &pb,
    const groupT &P,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix), value(pb, P, annotation_prefix)
{
    is_identity.assign(pb, P.is_zero() ? FieldT::one() : FieldT::zero());
    is_identity.evaluate(pb);
}

template<typename ppT, typename groupT, typename groupVariableT>
void variable_or_identity<ppT, groupT, groupVariableT>::generate_r1cs_witness(
    const groupT &elt)
{
    const bool is_zero = elt.is_zero();
    value.generate_r1cs_witness(is_zero ? groupT::one() : elt);
    generate_r1cs_witness(is_zero);
}

template<typename ppT, typename groupT, typename groupVariableT>
void variable_or_identity<ppT, groupT, groupVariableT>::generate_r1cs_witness(
    bool is_zero)
{
    this->pb.val(is_identity_var) = is_zero ? FieldT::one() : FieldT::zero();
}

template<typename ppT, typename groupT, typename groupVariableT>
groupT variable_or_identity<ppT, groupT, groupVariableT>::get_element() const
{
    if (this->pb.lc_val(is_identity) == FieldT::one()) {
        return groupT::zero();
    }

    return value.get_element();
}

template<
    typename ppT,
    typename groupT,
    typename groupVariableT,
    typename variableSelectorT,
    typename addGadgetT>
add_variable_or_identity<
    ppT,
    groupT,
    groupVariableT,
    variableSelectorT,
    addGadgetT>::
    add_variable_or_identity(
        protoboard<FieldT> &pb,
        const variableOrIdentity &A,
        const variableOrIdentity &B,
        const variableOrIdentity &result,
        const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , add_result(pb, FMT(annotation_prefix, " add_result"))
    , add(pb, A.value, B.value, add_result, FMT(annotation_prefix, " add"))
    , A_not_identity_result(
          pb, FMT(annotation_prefix, " A_not_identity_result"))
    , selector_A_not_identity(
          pb,
          B.is_identity,
          add_result,
          A.value,
          A_not_identity_result,
          FMT(annotation_prefix, " selector_A_not_identity"))
    , selector_A(
          pb,
          A.is_identity,
          A_not_identity_result,
          B.value,
          result.value,
          FMT(annotation_prefix, " selector_A"))
    , result(result)
{
}

template<
    typename ppT,
    typename groupT,
    typename groupVariableT,
    typename variableSelectorT,
    typename addGadgetT>
void add_variable_or_identity<
    ppT,
    groupT,
    groupVariableT,
    variableSelectorT,
    addGadgetT>::generate_r1cs_constraints()
{
    add.generate_r1cs_constraints();
    selector_A_not_identity.generate_r1cs_constraints();
    selector_A.generate_r1cs_constraints();

    // result.is_identity = A.is_identity * B.is_identity
    // (no need to call generate_r1cs_constraints() on result itself)
    this->pb.add_r1cs_constraint(
        r1cs_constraint<libff::Fr<ppT>>(
            selector_A.selector,              // A.is_identity
            selector_A_not_identity.selector, // B.is_identity
            result.is_identity),
        FMT(this->annotation_prefix, " result.is_identity"));
}

template<
    typename ppT,
    typename groupT,
    typename groupVariableT,
    typename variableSelectorT,
    typename addGadgetT>
void add_variable_or_identity<
    ppT,
    groupT,
    groupVariableT,
    variableSelectorT,
    addGadgetT>::generate_r1cs_witness()
{
    add.generate_r1cs_witness();
    selector_A_not_identity.generate_r1cs_witness();

    // Generate result.value via the result of selector_A, and set
    // result.is_identity manually via result.generate_r1cs_witness().
    selector_A.generate_r1cs_witness();

    const libff::Fr<ppT> result_is_identity =
        this->pb.lc_val(selector_A.selector) *
        this->pb.lc_val(selector_A_not_identity.selector);
    result.generate_r1cs_witness(result_is_identity == libff::Fr<ppT>::one());
}

template<
    typename ppT,
    typename groupT,
    typename groupVariableT,
    typename dblGadgetT>
dbl_variable_or_identity<ppT, groupT, groupVariableT, dblGadgetT>::
    dbl_variable_or_identity(
        protoboard<FieldT> &pb,
        const variableOrIdentity &A,
        const variableOrIdentity &result,
        const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , A_is_identity(A.is_identity)
    , result(result)
    , double_gadget(
          pb, A.value, result.value, FMT(annotation_prefix, " double_gadget"))
{
}

template<
    typename ppT,
    typename groupT,
    typename groupVariableT,
    typename dblGadgetT>
void dbl_variable_or_identity<ppT, groupT, groupVariableT, dblGadgetT>::
    generate_r1cs_constraints()
{
    // It should be possible to do this by simply assigning A.is_identity to
    // result.is_identity, but result.is_identity has already been allocated at
    // this stage.

    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(
            A_is_identity, libff::Fr<ppT>::one(), result.is_identity),
        FMT(this->annotation_prefix, " result_is_identity_constraint"));
    double_gadget.generate_r1cs_constraints();
}

template<
    typename ppT,
    typename groupT,
    typename groupVariableT,
    typename dblGadgetT>
void dbl_variable_or_identity<ppT, groupT, groupVariableT, dblGadgetT>::
    generate_r1cs_witness()
{
    A_is_identity.evaluate(this->pb);
    this->pb.lc_val(result.is_identity) = this->pb.lc_val(A_is_identity);
    double_gadget.generate_r1cs_witness();
}

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
