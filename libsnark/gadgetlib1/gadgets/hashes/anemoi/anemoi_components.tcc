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

// R1CS constraints for the operation y = const_a x^2 + const_b with x = input,
// y = output. This operation is realized by the components \ref
// flystel_Q_gamma_prime_field_gadget and \ref
// flystel_Q_delta_prime_field_gadget
//
// The operation is represented with one multiplication as
//
// (const_a x) * x = y-const_b
//
// for the variables vector X = (x0=1, x1=x, x2=y). This computation
// is represented with 1 R1CS constraint:
//
// < A , X > * < B , X > = < C , X >
//
// where A =(0, const_a, 0), B=(0, 1, 0) and C =(-const_b, 0, 1)

template<typename FieldT, size_t generator>
flystel_Q_gamma_prime_field_gadget<FieldT, generator>::
    flystel_Q_gamma_prime_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
#ifdef FLYSTEL_DEBUG
    , beta(DEBUG_FLYSTEL_BETA)
    , gamma(DEBUG_FLYSTEL_GAMMA)
#else
    , beta(FieldT(generator))
    , gamma(FieldT(0))
#endif // #ifdef FLYSTEL_DEBUG
    , input(input)
    , output(output)
{
}

template<typename FieldT, size_t generator>
void flystel_Q_gamma_prime_field_gadget<FieldT, generator>::
    generate_r1cs_constraints()
{
    // Constraint has the form:
    //   beta * input^2 + gamma = output
    // which can be written as
    //   (beta * input) * input = output - gamma
    this->pb.add_r1cs_constraint(
        {input * beta, input, output - gamma},
        FMT(this->annotation_prefix, " beta * x = y - gamma"));
}

// compute a witness y for a given input x for the computation y =
// beta x^2 + gamma, where x=input, y=output
template<typename FieldT, size_t generator>
void flystel_Q_gamma_prime_field_gadget<FieldT, generator>::
    generate_r1cs_witness()
{
    // y = beta x^2 + gamma
    this->pb.val(output) =
        this->beta * this->pb.val(input) * this->pb.val(input) + this->gamma;
}

template<typename FieldT, size_t generator>
flystel_Q_delta_prime_field_gadget<FieldT, generator>::
    flystel_Q_delta_prime_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
#ifdef FLYSTEL_DEBUG
    , beta(DEBUG_FLYSTEL_BETA)
    , delta(DEBUG_FLYSTEL_DELTA)
#elif
    , beta(FieldT(generator))
    , delta(FieldT(generator).inverse())
#endif // #ifdef FLYSTEL_DEBUG
    , input(input)
    , output(output)
{
}

template<typename FieldT, size_t generator>
void flystel_Q_delta_prime_field_gadget<FieldT, generator>::
    generate_r1cs_constraints()
{
    // Constraint has the form:
    //   beta * input^2 + delta = output
    // which can be written as
    //   (beta * input) * input = output - delta
    this->pb.add_r1cs_constraint(
        {input * beta, input, output - delta},
        FMT(this->annotation_prefix, " beta * x = y - delta"));
}

// compute a witness y for a given input x for the computation y =
// beta x^2 + delta, where x=input, y=output
template<typename FieldT, size_t generator>
void flystel_Q_delta_prime_field_gadget<FieldT, generator>::
    generate_r1cs_witness()
{
    // y = beta x^2 + delta
    this->pb.val(output) =
        this->beta * this->pb.val(input) * this->pb.val(input) + this->delta;
}

// R1CS constraints for the operation y = beta x^3 + gamma with
// x=input, y=output. This operation is represented with two
// multiplications as y-gamma = ((beta x * x) * x). Equivalently:
//
// beta x1 * x1 = x2
// x2 * x1 = x3-gamma
//
// for the variables vector X = (x0=1, x1=input, x2=intermediate,
// x3=output). The above system is represented with 2 R1CS
// constraints resp.:
//
// < A0 , X > * < B0 , X > = < C0 , X > ,
// < A1 , X > * < B1 , X > = < C1 , X >
//
// where A0=(0, beta, 0, 0), B0=(0, 1, 0, 0), C0=(0, 0, 1, 0) and
// A1=(0, 0, 1, 0), B1=(0, 1, 0, 0), C1=(-gamma, 0, 0, 1)
template<typename FieldT, size_t generator>
flystel_Q_gamma_binary_field_gadget<FieldT, generator>::
    flystel_Q_gamma_binary_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , internal(pb_variable_allocate<FieldT>(
          pb, FMT(this->annotation_prefix, " internal")))
#ifdef FLYSTEL_DEBUG
    , beta(DEBUG_FLYSTEL_BETA)
    , gamma(DEBUG_FLYSTEL_GAMMA)
#elif
    , beta(FieldT(generator))
    , gamma(FieldT(0))
#endif // #ifdef FLYSTEL_DEBUG
    , input(input)
    , output(output)
{
}

template<typename FieldT, size_t generator>
void flystel_Q_gamma_binary_field_gadget<FieldT, generator>::
    generate_r1cs_constraints()
{
    // (beta * input) * input = internal
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(beta * input, input, internal),
        FMT(this->annotation_prefix, " beta * x * x = x_square"));
    // internal * input = output - gamma
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(internal, input, output - gamma),
        FMT(this->annotation_prefix, " x_square * x = y - gamma"));
}

template<typename FieldT, size_t generator>
void flystel_Q_gamma_binary_field_gadget<FieldT, generator>::
    generate_r1cs_witness()
{
    // x_internal = beta x * x
    this->pb.val(internal) =
        (this->beta * this->pb.val(input)) * this->pb.val(input);
    // y = beta x^3 + gamma = x_internal * x + gamma
    this->pb.val(output) =
        this->pb.val(internal) * this->pb.val(input) + this->gamma;
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
flystel_E_power_five_gadget<FieldT>::flystel_E_power_five_gadget(
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
void flystel_E_power_five_gadget<FieldT>::generate_r1cs_constraints()
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
void flystel_E_power_five_gadget<FieldT>::generate_r1cs_witness()
{
    // x2 = x1 * x1
    this->pb.val(internal[0]) = (this->pb.val(input)) * this->pb.val(input);
    // x3 = x2 * x2
    this->pb.val(internal[1]) =
        (this->pb.val(internal[0])) * this->pb.val(internal[0]);
    // y = x1 * x3
    this->pb.val(output) = this->pb.val(input) * this->pb.val(internal[1]);
}

template<typename FieldT, size_t generator>
flystel_closed_prime_field_gadget<FieldT, generator>::
    flystel_closed_prime_field_gadget(
        protoboard<FieldT> &pb,
        const std::array<pb_variable<FieldT>, 2> &input,
        const std::array<pb_variable<FieldT>, 2> &output,
        const std::string &annotation_prefix)
    : flystel_Q_gamma_prime_field_gadget<FieldT, generator>(
          pb, input[0], internal[0])
    , flystel_Q_delta_prime_field_gadget<FieldT, generator>(
          pb, input[1], internal[2])
    , flystel_E_power_five_gadget<FieldT>(pb, input[0] - input[1], internal[1])
{
    internal[0].allocate(this->pb, " v3");
    internal[1].allocate(this->pb, " v4");
    internal[2].allocate(this->pb, " v5");
}

// R1CS constraints for the operation
//
// y0 = Q_gamma(x0) + power_five(x0-x1)
// y1 = Q_delta(x1) + power_five(x0-x1)
//
// x0=input[0], x1=input[1], y0=output[0], y1=output[1].
//
// The function generates the constraints for the three gadgets:
// Q_gamma, Q_delta, power_five by calling their corresponding
// generate_r1cs_constraints() methods
template<typename FieldT, size_t generator>
void flystel_closed_prime_field_gadget<FieldT, generator>::
    generate_r1cs_constraints()
{
    Q_gamma.generate_r1cs_constraints();
    Q_delta.generate_r1cs_constraints();
    power_five.generate_r1cs_constraints();
}

template<typename FieldT, size_t generator>
void flystel_closed_prime_field_gadget<FieldT, generator>::
    generate_r1cs_witness()
{
    Q_gamma.generate_r1cs_witness();
    Q_delta.generate_r1cs_witness();
    power_five.generate_r1cs_witness();
}

template<typename FieldT, size_t generator, size_t NumStateColumns_L>
anemoi_permutation_round_prime_field_gadget<
    FieldT,
    generator,
    NumStateColumns_L>::
    anemoi_permutation_round_prime_field_gadget(
        protoboard<FieldT> &pb,
        std::array<pb_variable<FieldT>, 2 * NumStateColumns_L> &input,
        std::array<pb_variable<FieldT>, 2 * NumStateColumns_L> &output,
        std::string &annotation_prefix)
    : input(input), output(output)
{
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_
