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

/*
class flystel_power_two_round_1_gadget : public flystel_power_two_gadget
{
flystel_power_two_round_1_gadget(
rotoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "")
  : flystel_power_two_gadget(..., ...ALPHA, ... BETA, ...)
  {
  }
}

class flystel_power_two_round_2_gadget : public flystel_power_two_gadget
{
flystel_power_two_round_2_gadget(
rotoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "")
  : flystel_power_two_gadget(..., flystel_constants_selector<ppT>::BETA, ...
flystel_constants_selector<ppT>::GAMMA, ...)
  {
  }
}
*/

template<typename FieldT>
flystel_power_two_gadget<FieldT>::flystel_power_two_gadget(
    protoboard<FieldT> &pb,
    const pb_variable<FieldT> &input,
    const pb_variable<FieldT> &output,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , const_a(FLYSTEL_BLS12_381_BETA)
    , const_b(FLYSTEL_BLS12_381_GAMMA)
    , input(input)
    , output(output)
{
}

// A R1CS constraint is a formal expression of the form
//
//                < A , X > * < B , X > = < C , X > ,
//
// where X = (x_0,x_1,...,x_m) is a vector of formal variables and
// A,B,C each consist of 1+m elements in <FieldT> and <A, X> = \sum_i
// (a_i x_i) is the dot product between vectors A and X. Equivalently,
// the vectors A,B,C are linear combinations of X.
//
// A R1CS constraint is used to construct a R1CS constraint system.
//
// See also class \ef r1cs_constraint

// R1CS constraints for the operation y = const_a x^2 + const_b with x =
// input, y = output. The latter is represented with one
// multiplication as
//
// (const_a x) * x = y-const_b
//
// for the variables vector X = (x0=1, x1=x, x2=y). This computation
// is represented with 1 R1CS constraint:
//
// < A , X > * < B , X > = < C , X >
//
// where A =(0, const_a, 0), B=(0, 1, 0) and C =(-const_b, 0, 1)
template<typename FieldT>
void flystel_power_two_gadget<FieldT>::generate_r1cs_constraints()
{
    // Constraint has the form:
    //   const_a * input^2 + const_b = output
    // which can be written as
    //   (const_a * input) * input = output - const_b
    this->pb.add_r1cs_constraint(
        {input * const_a, input, output - const_b},
        FMT(this->annotation_prefix, " const_a * x = y - const_b"));
}

// compute a witness y for a given input x for the computation y =
// const_a x^2 + const_b, where x=input, y=output
template<typename FieldT>
void flystel_power_two_gadget<FieldT>::generate_r1cs_witness()
{
    // y = const_a x^2 + const_b
    this->pb.val(output) =
        this->const_a * this->pb.val(input) * this->pb.val(input) +
        this->const_b;
}

template<typename FieldT>
flystel_Qi_power_two_gadget<FieldT>::flystel_Qi_power_two_gadget(
    protoboard<FieldT> &pb,
    const pb_variable<FieldT> &input,
    const pb_variable<FieldT> &output,
    const std::string &annotation_prefix)
    : flystel_power_two_gadget<FieldT>(
          pb, FLYSTEL_BLS12_381_BETA, FLYSTEL_BLS12_381_GAMMA, input, output)
{
}

template<typename FieldT>
flystel_Qf_power_two_gadget<FieldT>::flystel_Qf_power_two_gadget(
    protoboard<FieldT> &pb,
    const pb_variable<FieldT> &input,
    const pb_variable<FieldT> &output,
    const std::string &annotation_prefix)
    : flystel_power_two_gadget<FieldT>(
          pb, FLYSTEL_BLS12_381_BETA, FLYSTEL_BLS12_381_DELTA, input, output)
{
}

template<typename FieldT>
flystel_power_three_gadget<FieldT>::flystel_power_three_gadget(
    protoboard<FieldT> &pb,
    const pb_variable<FieldT> &input,
    const pb_variable<FieldT> &output,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , internal(pb_variable_allocate<FieldT>(
          pb, FMT(this->annotation_prefix, " internal")))
    , const_a(FLYSTEL_BLS12_381_BETA)
    , const_b(FLYSTEL_BLS12_381_GAMMA)
    , input(input)
    , output(output)
{
}

// R1CS constraints for the operation y = const_a x^3 + const_b with
// x=input, y=output. This operation is represented with two
// multiplications as y-const_b = ((const_a x * x) * x). Equivalently:
//
// const_a x1 * x1 = x2
// x2 * x1 = x3-const_b
//
// for the variables vector X = (x0=1, x1=input, x2=intermediate,
// x3=output). The above system is represented with 2 R1CS
// constraints resp.:
//
// < A0 , X > * < B0 , X > = < C0 , X > ,
// < A1 , X > * < B1 , X > = < C1 , X >
//
// where A0=(0, const_a, 0, 0), B0=(0, 1, 0, 0), C0=(0, 0, 1, 0) and
// A1=(0, 0, 1, 0), B1=(0, 1, 0, 0), C1=(-const_b, 0, 0, 1)
template<typename FieldT>
void flystel_power_three_gadget<FieldT>::generate_r1cs_constraints()
{
    // (const_a * input) * input = internal
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(const_a * input, input, internal),
        FMT(this->annotation_prefix, " const_a * x * x = x_square"));
    // internal * input = output - const_b
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(internal, input, output - const_b),
        FMT(this->annotation_prefix, " x_square * x = y - const_b"));
}

template<typename FieldT>
void flystel_power_three_gadget<FieldT>::generate_r1cs_witness()
{
    // x_internal = const_a x * x
    this->pb.val(internal) =
        (this->const_a * this->pb.val(input)) * this->pb.val(input);
    // y = const_a x^3 + const_b = x_internal * x + const_b
    this->pb.val(output) =
        this->pb.val(internal) * this->pb.val(input) + this->const_b;
}

template<typename FieldT>
flystel_power_five_gadget<FieldT>::flystel_power_five_gadget(
    protoboard<FieldT> &pb,
    const pb_variable<FieldT> &input,
    const pb_variable<FieldT> &output,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix), input(input), output(output)
{
    internal[0].allocate(this->pb, " internal 1");
    internal[1].allocate(this->pb, " internal 2");
}

template<typename FieldT>
flystel_Qi_power_three_gadget<FieldT>::flystel_Qi_power_three_gadget(
    protoboard<FieldT> &pb,
    const pb_variable<FieldT> &input,
    const pb_variable<FieldT> &output,
    const std::string &annotation_prefix)
    : flystel_power_three_gadget<FieldT>(
          pb, FLYSTEL_BLS12_381_BETA, FLYSTEL_BLS12_381_GAMMA, input, output)
{
}

template<typename FieldT>
flystel_Qf_power_three_gadget<FieldT>::flystel_Qf_power_three_gadget(
    protoboard<FieldT> &pb,
    const pb_variable<FieldT> &input,
    const pb_variable<FieldT> &output,
    const std::string &annotation_prefix)
    : flystel_power_three_gadget<FieldT>(
          pb, FLYSTEL_BLS12_381_BETA, FLYSTEL_BLS12_381_GAMMA, input, output)
{
}

// R1CS constraints for the operation y = x^5 with x=input,
// y=output. This operation is represented with three multiplications
// as y = (temp * temp * x), temp = x * x. Equivalently:
//
// x1 * x1 = x2
// x2 * x2 = x3
// x1 * x3 = x4
//
// for the variables vector X = (x0=1, x1=input, x2=internal,
// x3=internal, x4=output). The above system is represented with 3
// R1CS constraints resp.:
//
// < A0 , X > * < B0 , X > = < C0 , X > ,
// < A1 , X > * < B1 , X > = < C1 , X > ,
// < A2 , X > * < B2 , X > = < C2 , X >
//
// where A0=(01000), B0=(01000), C0=(00100); A1=(00100), B0=(00100),
// C0=(00010) and A2=(01000), B2=(00010), C2=(00001)
template<typename FieldT>
void flystel_power_five_gadget<FieldT>::generate_r1cs_constraints()
{
    // x1*x1 = x2
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(input, input, internal[0]),
        FMT(this->annotation_prefix, " x * x = x^2"));
    // x2*x2 = x3
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(internal[0], internal[0], internal[1]),
        FMT(this->annotation_prefix, " x^2 * x^2 = x^4"));
    // x1*x3 = x4
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(input, internal[1], output),
        FMT(this->annotation_prefix, " x^1 * x^4 = x^5"));
}

template<typename FieldT>
void flystel_power_five_gadget<FieldT>::generate_r1cs_witness()
{
    // x2 = x1 * x1
    this->pb.val(internal[0]) = (this->pb.val(input)) * this->pb.val(input);
    // x3 = x2 * x2
    this->pb.val(internal[1]) =
        (this->pb.val(internal[0])) * this->pb.val(internal[0]);
    // y = x1 * x3
    this->pb.val(output) = this->pb.val(input) * this->pb.val(internal[1]);
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_
