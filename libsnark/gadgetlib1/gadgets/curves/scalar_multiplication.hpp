/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_CURVE_SCALAR_MULTIPLICATION_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_CURVE_SCALAR_MULTIPLICATION_HPP_

#include "libsnark/gadgetlib1/gadget.hpp"

#include <memory>

namespace libsnark
{

/// Wrapper for group variable, adding a flag that indicates whether the
/// variable is the identity.
///
/// For use in 2 possible configurations:
///
///   1) Construct with variable_or_identity(pb, annotation_prefix), which
///      internally allocates variables for the group element and flag.
///      Values can later be set via generate_r1cs_witness() calls.
///
///   2) Construct using a group element value, which defines values as linear
///      combinations (no variable allocation). generate_r1cs_witness() calls
///      should NOT be made to such objects.
template<typename ppT, typename groupT, typename groupVariableT>
class variable_or_identity : public gadget<libff::Fr<ppT>>
{
public:
    using FieldT = libff::Fr<ppT>;

    groupVariableT value;
    pb_linear_combination<FieldT> is_identity;

    // Internally allocates the variable and flag, to be set later via the
    // generate_r1cs_witness calls.
    variable_or_identity(
        protoboard<FieldT> &pb, const std::string &annotation_prefix);

    // Initialize the variable and flag with constants. Do not try to call
    // generate_r1cs_witness().
    variable_or_identity(
        protoboard<FieldT> &pb,
        const groupT &P,
        const std::string &annotation_prefix);

    /// Constrains the `is_identity` flag to be boolean. Do not call this if
    /// the setter can guarantee boolean-ness in some other way.
    void generate_r1cs_constraints();

    /// If elt == 0, this sets `value` to 1 and `is_identity` to true,
    /// otherwise it sets `value` to `elt` and `is_identity` to false.
    void generate_r1cs_witness(const groupT &elt);

    /// For the case where value has been set by other gadgets, and we simply
    /// wish to set whether this variable should be considered the identity.
    void generate_r1cs_witness(bool is_identity);

    groupT get_element() const;

protected:
    // Only allocated if necessary, in which case it should be set using
    // generate_r1cs_witness().
    pb_variable<FieldT> is_identity_var;
};

/// Wrap an add gadget, extending it to support adding variable_or_identity
/// elements. At most one element may be the identity, otherwise witness
/// generation fails due to the implementation of the underlying add gadgets.
template<
    typename ppT,
    typename groupT,
    typename variableT,
    typename variableSelectorT,
    typename addGadgetT>
class add_variable_or_identity : public gadget<libff::Fr<ppT>>
{
public:
    using FieldT = libff::Fr<ppT>;
    using variableOrIdentity = variable_or_identity<ppT, groupT, variableT>;

    // Recall the form of the variable selector gadget parameters:
    //   variableSelectorT(<flag>, <zero_case>, <one_case>)
    //
    // The `is_identity` flags of A, B and result are known to be boolean by
    // the constraints in variable_or_identity.
    //
    // The output from this gadget should therefore be be:
    //
    //   result.is_identity = A.is_identity * B.is_identity
    //
    //   result.value = variableSelectorT(  -- selector_A
    //       A.is_identity,
    //       variableSelectorT(             -- selector_A_not_identity
    //           B.is_identity,
    //           add_result,                   -- A != 0 && B != 0
    //           A.value),                     -- A != 0 && B == 0
    //       variableSelectorT(             -- selector_A_is_identity
    //           B.is_identity,
    //           B.value,                      -- A == 0 && B != 0
    //           groupT::one())                -- B == 0 && B == 0
    //
    // Note, this can be reduced to:
    //
    //   result.is_identity = A.is_identity * B.is_identity
    //
    //   result.value = variableSelectorT(  -- selector_A
    //       A.is_identity,
    //       variableSelectorT(             -- selector_A_not_identity
    //           B.is_identity,
    //           add_result,                   -- A != 0 && B != 0
    //           A.value),                     -- A != 0 && B == 0
    //       B.value)

    variableT add_result;
    addGadgetT add;

    // A_not_identity_result = variableSelectorT(
    //     B.is_identity, add_result, A.value)
    variableT A_not_identity_result;
    variableSelectorT selector_A_not_identity;

    // result = variableSelectorT(
    //     A.is_identity, A_not_identity_result, B.value)
    variableSelectorT selector_A;

    variableOrIdentity result;

    add_variable_or_identity(
        protoboard<FieldT> &pb,
        const variableOrIdentity &A,
        const variableOrIdentity &B,
        const variableOrIdentity &result,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();

    void generate_r1cs_witness();
};

/// Generic gadget to perform scalar multiplication of group variables. Used by
/// the individual group element implementations.
template<
    typename groupT,
    typename groupVariableT,
    typename add_gadget,
    typename dbl_gadget,
    typename scalarT>
class point_mul_by_const_scalar_gadget
    : public gadget<typename groupT::base_field>
{
public:
    using FieldT = typename groupT::base_field;

    const scalarT _scalar;
    const groupVariableT _result;
    std::vector<std::shared_ptr<add_gadget>> _add_gadgets;
    std::vector<std::shared_ptr<dbl_gadget>> _dbl_gadgets;

    point_mul_by_const_scalar_gadget(
        protoboard<FieldT> &pb,
        const scalarT &scalar,
        const groupVariableT &P,
        const groupVariableT &result,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
    const groupVariableT &result() const;
};

} // namespace libsnark

#include "libsnark/gadgetlib1/gadgets/curves/scalar_multiplication.tcc"

#endif // LIBSNARK_GADGETLIB1_GADGETS_CURVE_SCALAR_MULTIPLICATION_HPP_
