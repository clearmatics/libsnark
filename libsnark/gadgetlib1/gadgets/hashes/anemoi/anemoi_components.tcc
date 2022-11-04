/** @file
 *****************************************************************************

 Implementation of interfaces for top-level Anemoi hash function gadgets.

 See anemoi_gadget.hpp .

 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_

namespace libsnark
{

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

template<typename ppT>
flystel_Q_prime_field_gadget<ppT>::flystel_Q_prime_field_gadget(
    protoboard<libff::Fr<ppT>> &pb,
    const libff::Fr<ppT> A,
    const libff::Fr<ppT> B,
    const linear_combination<libff::Fr<ppT>> &input,
    const pb_variable<libff::Fr<ppT>> &output,
    const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , A(A)
    , B(B)
    , input(input)
    , output(output)
{
}

template<typename ppT>
void flystel_Q_prime_field_gadget<ppT>::generate_r1cs_constraints()
{
    // Constraint has the form:
    //   A * input^2 + B = output
    // which can be written as
    //   (A * input) * input = output - B
    this->pb.add_r1cs_constraint(
        {input * A, input, output - B},
        FMT(this->annotation_prefix, " A * x = y - B"));
}

// compute a witness y for a given input x for the computation y =
// A x^2 + B
template<typename ppT>
void flystel_Q_prime_field_gadget<ppT>::generate_r1cs_witness()
{
    using FieldT = libff::Fr<ppT>;
    const FieldT input_value =
        input.evaluate(this->pb.full_variable_assignment());
    // y = A x^2 + B
    this->pb.val(output) = A * input_value * input_value + B;
}

// R1CS constraints for the operation y = beta x^3 + gamma with
// x=input, y=output, represented as y-gamma = ((beta x * x) *
// x) or equivalently:
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

template<typename ppT>
flystel_Q_binary_field_gadget<ppT>::flystel_Q_binary_field_gadget(
    protoboard<libff::Fr<ppT>> &pb,
    const libff::Fr<ppT> A,
    const libff::Fr<ppT> B,
    const linear_combination<libff::Fr<ppT>> &input,
    const pb_variable<libff::Fr<ppT>> &output,
    const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , internal(pb_variable_allocate<libff::Fr<ppT>>(
          pb, FMT(this->annotation_prefix, " internal")))
    , A(A)
    , B(B)
    , input(input)
    , output(output)
{
}

template<typename ppT>
void flystel_Q_binary_field_gadget<ppT>::generate_r1cs_constraints()
{
    // (A * input) * input = internal
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(A * input, input, internal),
        FMT(this->annotation_prefix, " A * x * x = x_square"));
    // internal * input = output - B
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(internal, input, output - B),
        FMT(this->annotation_prefix, " x_square * x = y - B"));
}

template<typename ppT>
void flystel_Q_binary_field_gadget<ppT>::generate_r1cs_witness()
{
    const FieldT input_value =
        input.evaluate(this->pb.full_variable_assignment());
    // x_internal = A x * x
    this->pb.val(internal) = (A * input_value) * input_value;
    // y = A x^3 + B = x_internal * x + B
    this->pb.val(output) = this->pb.val(internal) * input_value + B;
}

// R1CS constraints for the operation y = x^5 with x=input,
// y=output. This operation is represented with three multiplications
// as y = (temp * temp * x), temp = x * x. Equivalently:
//
// x  * x  = a0
// a0 * a0 = a1
// x  * a1 =  y
//
// for the variables vector X = (x0=1, x1=x, x2=a0, x3=a1, x4=y). The
// above system is represented with 3 R1CS constraints resp.:
//
// < A0 , X > * < B0 , X > = < C0 , X > ,
// < A1 , X > * < B1 , X > = < C1 , X > ,
// < A2 , X > * < B2 , X > = < C2 , X >
//
// where A0=(01000), B0=(01000), C0=(00100); A1=(00100), B0=(00100),
// C0=(00010) and A2=(01000), B2=(00010), C2=(00001)
template<typename ppT>
flystel_E_power_five_gadget<ppT>::flystel_E_power_five_gadget(
    protoboard<libff::Fr<ppT>> &pb,
    const linear_combination<libff::Fr<ppT>> &input,
    const pb_variable<libff::Fr<ppT>> &output,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , a0(pb_variable_allocate(pb, FMT(annotation_prefix, " a0")))
    , a1(pb_variable_allocate(pb, FMT(annotation_prefix, " a1")))
    , input(input)
    , output(output)
{
}

template<typename ppT>
void flystel_E_power_five_gadget<ppT>::generate_r1cs_constraints()
{
    // x1*x1 = x2
    this->pb.add_r1cs_constraint(
        r1cs_constraint<libff::Fr<ppT>>(input, input, a0),
        FMT(this->annotation_prefix, " x * x = x^2"));
    // x2*x2 = x3
    this->pb.add_r1cs_constraint(
        r1cs_constraint<libff::Fr<ppT>>(a0, a0, a1),
        FMT(this->annotation_prefix, " x^2 * x^2 = x^4"));
    // x1*x3 = x4
    this->pb.add_r1cs_constraint(
        r1cs_constraint<libff::Fr<ppT>>(input, a1, output),
        FMT(this->annotation_prefix, " x^1 * x^4 = x^5"));
}

template<typename ppT>
void flystel_E_power_five_gadget<ppT>::generate_r1cs_witness()
{
    const libff::Fr<ppT> input_value =
        input.evaluate(this->pb.full_variable_assignment());
    // x2 = x1 * x1
    this->pb.val(a0) = (input_value)*input_value;
    // x3 = x2 * x2
    this->pb.val(a1) = (this->pb.val(a0)) * this->pb.val(a0);
    // y = x1 * x3
    this->pb.val(output) = input_value * this->pb.val(a1);
}

// R1CS constraints for the operation y = x^1/5 with x=input,
// y=output. The constraints are computed using the equivalent
// operation y^5=x (\see flystel_E_power_five_gadget). This operation
// is represented with three multiplications as x = (temp * temp * y),
// temp = y * y. Equivalently:
//
// y  *  y = a0
// a0 * a0 = a1
// y  * a1 = x
//
// for the variables vector X = (x0=1, x1=x, x2=a0, x3=a1, x4=y). The
// above system is represented with 3 R1CS constraints resp.:
//
// < A0 , X > * < B0 , X > = < C0 , X > ,
// < A1 , X > * < B1 , X > = < C1 , X > ,
// < A2 , X > * < B2 , X > = < C2 , X >
//
// where A0=(01000), B0=(01000), C0=(00100); A1=(00100), B0=(00100),
// C0=(00010) and A2=(01000), B2=(00010), C2=(00001)
template<typename ppT, class parameters>
flystel_E_root_five_gadget<ppT, parameters>::flystel_E_root_five_gadget(
    protoboard<libff::Fr<ppT>> &pb,
    const linear_combination<libff::Fr<ppT>> &input,
    const pb_variable<libff::Fr<ppT>> &output,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , a0(pb_variable_allocate(pb, FMT(annotation_prefix, " a0")))
    , a1(pb_variable_allocate(pb, FMT(annotation_prefix, " a1")))
    , input(input)
    , output(output)
{
}

template<typename ppT, class parameters>
void flystel_E_root_five_gadget<ppT, parameters>::generate_r1cs_constraints()
{
    using FieldT = libff::Fr<ppT>;
    // y1*y1 = y2
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(output, output, a0),
        FMT(this->annotation_prefix, " y * y = y^2"));
    // y2*y2 = y3
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(a0, a0, a1),
        FMT(this->annotation_prefix, " y^2 * y^2 = y^4"));
    // y1*y3 = y4
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(output, a1, input),
        FMT(this->annotation_prefix, " y * y^4 = y^5"));
}

template<typename ppT, class parameters>
void flystel_E_root_five_gadget<ppT, parameters>::generate_r1cs_witness()
{
    const libff::Fr<ppT> input_value =
        input.evaluate(this->pb.full_variable_assignment());
    libff::Fr<ppT> x = input_value; // this->pb.lc_val(input);
    libff::Fr<ppT> y = power(x, parameters::alpha_inv);

    // x2 = x1 * x1
    this->pb.val(a0) = y * y;
    // x3 = x2 * x2
    this->pb.val(a1) = (this->pb.val(a0)) * this->pb.val(a0);
    // y = x1 * x3
    this->pb.val(output) = y;
}

template<typename ppT, class parameters>
flystel_prime_field_gadget<ppT, parameters>::flystel_prime_field_gadget(
    protoboard<libff::Fr<ppT>> &pb,
    const linear_combination<libff::Fr<ppT>> &x0,
    const linear_combination<libff::Fr<ppT>> &x1,
    const pb_variable<libff::Fr<ppT>> &y0,
    const pb_variable<libff::Fr<ppT>> &y1,
    const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , a0(pb_variable_allocate(pb, FMT(annotation_prefix, " a0")))
    , a1(pb_variable_allocate(pb, FMT(annotation_prefix, " a1")))
    , a2(pb_variable_allocate(pb, FMT(annotation_prefix, " a2")))
    , input_x0(x0)
    , input_x1(x1)
    , output_y0(y0)
    , output_y1(y1)
    , Q_gamma(
          pb,
          parameters::beta,
          parameters::gamma,
          x1,
          a0,
          FMT(annotation_prefix, " Q_gamma"))
    , Q_delta(
          pb,
          parameters::beta,
          parameters::delta,
          x1 - a1,
          a2,
          FMT(annotation_prefix, " Q_delta"))
    , E_root_five(pb, x0 - a0, a1, FMT(annotation_prefix, " E_root_five"))
{
    static_assert((parameters::alpha == 5), "Parameter alpha must be 5");
}

template<typename ppT, class parameters>
void flystel_prime_field_gadget<ppT, parameters>::generate_r1cs_constraints()
{
    Q_gamma.generate_r1cs_constraints();
    Q_delta.generate_r1cs_constraints();
    E_root_five.generate_r1cs_constraints();
}

template<typename ppT, class parameters>
void flystel_prime_field_gadget<ppT, parameters>::generate_r1cs_witness()
{
    Q_gamma.generate_r1cs_witness();
    E_root_five.generate_r1cs_witness();
    Q_delta.generate_r1cs_witness();

    const FieldT input_x0_value =
        input_x0.evaluate(this->pb.full_variable_assignment());
    const FieldT input_x1_value =
        input_x1.evaluate(this->pb.full_variable_assignment());

    this->pb.lc_val(output_y0) =
        input_x0_value - this->pb.val(a0) + this->pb.val(a2);
    this->pb.lc_val(output_y1) = input_x1_value - this->pb.val(a1);
}

template<typename FieldT, size_t NumStateColumns_L>
std::array<std::array<FieldT, NumStateColumns_L>, NumStateColumns_L>
anemoi_permutation_mds(const FieldT g)
{
    static_assert(
        (NumStateColumns_L == 2) || (NumStateColumns_L == 3) ||
            (NumStateColumns_L == 4),
        "NumStateColumns_L must be 2,3 or 4");

    std::array<std::array<FieldT, NumStateColumns_L>, NumStateColumns_L> M;
    const FieldT g2 = g * g;
    if (NumStateColumns_L == 2) {
        M = {{1, g}, {g, g2 + 1}};
        return M;
    }
    if (NumStateColumns_L == 3) {
        M = {{g + 1, 1, g + 1}, {1, 1, g}, {g, 1, 1}};
        return M;
    }
    if (NumStateColumns_L == 4) {
        M = {
            {1, 1 + g, g, g},
            {g2, g + g2, 1 + g, 1 + 2 * g},
            {g2, g2, 1, 1 + g},
            {1 + g, 1 + 2 * g, g, 1 + g}};
        return M;
    }
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_
