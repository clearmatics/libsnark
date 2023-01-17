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

// Fast matrix-vector multiplication algorithm for Anemoi MDS layer with \ell =
// 1,2 for inputs of type "linear combination of FieldT elements"
template<typename ppT>
std::vector<linear_combination<libff::Fr<ppT>>> anemoi_fast_multiply_mds_2x2(
    const std::vector<linear_combination<libff::Fr<ppT>>> X_input,
    const libff::Fr<ppT> g)
{
    if (!(X_input.size() == 2)) {
        throw std::invalid_argument("input vector must be of length 2");
    }
    std::vector<linear_combination<libff::Fr<ppT>>> X = X_input;
    X[0] = X[0] + (g * X[1]);
    X[1] = X[1] + (g * X[0]);
    return X;
}

// Fast matrix-vector multiplication algorithm for Anemoi MDS layer with \ell
// = 3 for inputs of type "linear combination of FieldT elements". From Figure 6
// of [DL18](https://tosc.iacr.org/index.php/ToSC/article/view/888).
template<typename ppT>
std::vector<linear_combination<libff::Fr<ppT>>> anemoi_fast_multiply_mds_3x3(
    const std::vector<linear_combination<libff::Fr<ppT>>> X_input,
    const libff::Fr<ppT> g)
{
    if (!(X_input.size() == 3)) {
        throw std::invalid_argument("input vector must be of length 3");
    }
    std::vector<linear_combination<libff::Fr<ppT>>> X = X_input;
    linear_combination<libff::Fr<ppT>> t = X[0] + (g * X[2]);
    X[2] = X[2] + X[1];
    X[2] = X[2] + (g * X[0]);
    X[0] = t + X[2];
    X[1] = X[1] + t;
    return X;
}

// Fast matrix-vector multiplication algorithm for Anemoi MDS layer with \ell
// = 4 for inputs of type "linear combination of FieldT elements". Figure 8 of
// [DL18](https://tosc.iacr.org/index.php/ToSC/article/view/888).
template<typename ppT>
std::vector<linear_combination<libff::Fr<ppT>>> anemoi_fast_multiply_mds_4x4(
    const std::vector<linear_combination<libff::Fr<ppT>>> X_input,
    const libff::Fr<ppT> g)
{
    if (!(X_input.size() == 4)) {
        throw std::invalid_argument("input vector must be of length 4");
    }
    std::vector<linear_combination<libff::Fr<ppT>>> X = X_input;
    X[0] = X[0] + X[1];
    X[2] = X[2] + X[3];
    X[3] = X[3] + (g * X[0]);
    X[1] = g * (X[1] + X[2]);
    X[0] = X[0] + X[1];
    X[2] = X[2] + (g * X[3]);
    X[1] = X[1] + X[2];
    X[3] = X[3] + X[0];
    return X;
}

// multiply matrix by a vector of elements of type "linear combination of FieldT
// elements"
template<typename ppT, size_t NumStateColumns_L>
std::vector<linear_combination<libff::Fr<ppT>>> anemoi_fast_multiply_mds(
    const std::vector<linear_combination<libff::Fr<ppT>>> X,
    const libff::Fr<ppT> g)
{
    static_assert(
        (NumStateColumns_L == 1) || (NumStateColumns_L == 2) ||
            (NumStateColumns_L == 3) || (NumStateColumns_L == 4),
        "NumStateColumns_L must be 2,3 or 4");
    if (!(X.size() == NumStateColumns_L)) {
        throw std::invalid_argument("invalid length of input vector");
    }

    std::vector<linear_combination<libff::Fr<ppT>>> Y;
    if (NumStateColumns_L == 2) {
        Y = anemoi_fast_multiply_mds_2x2<ppT>(X, g);
    }
    if (NumStateColumns_L == 3) {
        Y = anemoi_fast_multiply_mds_3x3<ppT>(X, g);
    }
    if (NumStateColumns_L == 4) {
        Y = anemoi_fast_multiply_mds_4x4<ppT>(X, g);
    }
    return Y;
}

// multiply matrix by a vector of elements of type "linear combination of FieldT
// elements"
template<typename ppT, size_t NumStateColumns_L>
std::vector<linear_combination<libff::Fr<ppT>>> anemoi_fast_multiply_mds(
    const pb_linear_combination_array<libff::Fr<ppT>> X, const libff::Fr<ppT> g)
{
    static_assert(
        (NumStateColumns_L == 2) || (NumStateColumns_L == 3) ||
            (NumStateColumns_L == 4),
        "NumStateColumns_L must be 2,3 or 4");
    if (!(X.size() == NumStateColumns_L)) {
        throw std::invalid_argument("invalid length of input vector");
    }

    std::vector<linear_combination<libff::Fr<ppT>>> Y;
    if (NumStateColumns_L == 2) {
        Y = anemoi_fast_multiply_mds_2x2<ppT>(X, g);
    }
    if (NumStateColumns_L == 3) {
        Y = anemoi_fast_multiply_mds_3x3<ppT>(X, g);
    }
    if (NumStateColumns_L == 4) {
        Y = anemoi_fast_multiply_mds_4x4<ppT>(X, g);
    }
    return Y;
}

// rotate left by 1 a vector of elements of type "linear combination of FieldT
// elements": (x1_0 x1_1 ... x1_{L-1}) -> (x1_1 ... x1_{L-1} x_0)
template<typename ppT>
std::vector<linear_combination<libff::Fr<ppT>>> anemoi_vector_left_rotate_by_one(
    const std::vector<linear_combination<libff::Fr<ppT>>> X_input)
{
    if (!((X_input.size() == 2) || (X_input.size() == 3) ||
          (X_input.size() == 4))) {
        throw std::invalid_argument("invalid length of input vector");
    }
    std::vector<linear_combination<libff::Fr<ppT>>> X = X_input;
    rotate(X.begin(), X.begin() + 1, X.end());
    return X;
}

// rotate left by 1 a vector of elements of type "linear combination of FieldT
// elements": (x1_0 x1_1 ... x1_{L-1}) -> (x1_1 ... x1_{L-1} x_0)
template<typename ppT>
pb_linear_combination_array<libff::Fr<ppT>> anemoi_vector_left_rotate_by_one(
    const pb_linear_combination_array<libff::Fr<ppT>> X_input)
{
    if (!((X_input.size() == 2) || (X_input.size() == 3) ||
          (X_input.size() == 4))) {
        throw std::invalid_argument("invalid length of input vector");
    }
    pb_linear_combination_array<libff::Fr<ppT>> X = X_input;
    rotate(X.begin(), X.begin() + 1, X.end());
    return X;
}

template<typename ppT, size_t NumStateColumns_L, class parameters>
anemoi_permutation_round_prime_field_gadget<
    ppT,
    NumStateColumns_L,
    parameters>::
    anemoi_permutation_round_prime_field_gadget(
        protoboard<libff::Fr<ppT>> &pb,
        const std::vector<FieldT> &C,
        const std::vector<FieldT> &D,
        const pb_linear_combination_array<FieldT> &X_left,
        const pb_linear_combination_array<FieldT> &X_right,
        const pb_variable_array<FieldT> &Y_left,
        const pb_variable_array<FieldT> &Y_right,
        const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , C_const(C)
    , D_const(D)
    , X_left_input(X_left)
    , X_right_input(X_right)
    , Y_left_output(Y_left)
    , Y_right_output(Y_right)
{
    const libff::Fr<ppT> g = parameters::multiplicative_generator_g;
    const size_t ncols = NumStateColumns_L;

    // temporary variables (Z_left, Z_right) modified in-place during
    // the computation from (X_left, X_right) to (Y_left, Y_right)
    std::vector<linear_combination<libff::Fr<ppT>>> Z_left;
    std::vector<linear_combination<libff::Fr<ppT>>> Z_right;

    // add constants Z_left[i]+=C[i], Z_right[i]+=D[i]
    for (size_t i = 0; i < ncols; i++) {
        Z_left.push_back(X_left[i] + C[i]);
        Z_right.push_back(X_right[i] + D[i]);
    }

    if (ncols > 1) {
        M_matrix = anemoi_permutation_mds<ppT, ncols>::permutation_mds(g);
    } else { // ncols == 1
        // the MDS matrix for a state with 1 column (L=1) is the same as
        // for a state with 2 columns (L=2)
        M_matrix = anemoi_permutation_mds<ppT, 2>::permutation_mds(g);
    }

    // multiply by matrix M
    if (ncols > 1) {
        // l > 1:
        // Z_left  = (zL_0 zL_1 ... zL_{l-1})
        // Z_right = (zR_0 zR_1 ... zR_{l-1})
        // Z_left = M Z_left
        // Z_right = M (Z_right <<< 1)
        // where (Z_right <<< 1) = (zR_1 ... zR_{l-1} zR_0)
        Z_left = anemoi_fast_multiply_mds<ppT, NumStateColumns_L>(Z_left, g);
        std::vector<linear_combination<libff::Fr<ppT>>> Z_right_lrot =
            anemoi_vector_left_rotate_by_one<ppT>(Z_right);
        Z_right =
            anemoi_fast_multiply_mds<ppT, NumStateColumns_L>(Z_right_lrot, g);
    } else { // ncols == 1
        // l = 1:
        // Z_left = zL_0
        // Z_right = zR_0
        // Z = Z_left || Z_right
        // Z = M Z
        // Z_left  = Z[0]
        // Z_right = Z[1]
        assert(Z_left.size() == 1);
        assert(Z_right.size() == 1);
        std::vector<linear_combination<libff::Fr<ppT>>> Z;
        Z.push_back(Z_left[0]);
        Z.push_back(Z_right[0]);
        // for L=1 still calling multiply routine for L=2
        Z = anemoi_fast_multiply_mds<ppT, 2>(Z, g);
        assert(Z.size() == 2);
        Z_left.clear();
        Z_right.clear();
        Z_left.push_back(Z[0]);
        Z_right.push_back(Z[1]);
    }

    // apply layer of L Flystel S-boxes
    for (size_t i = 0; i < ncols; i++) {
        Flystel.emplace_back(flystel_prime_field_gadget<ppT, parameters>(
            pb,
            Z_left[i],
            Z_right[i],
            Y_left[i],
            Y_right[i],
            FMT(this->annotation_prefix, " Flystel[%zu]", i)));
    }
}

template<typename ppT, size_t NumStateColumns_L, class parameters>
void anemoi_permutation_round_prime_field_gadget<
    ppT,
    NumStateColumns_L,
    parameters>::generate_r1cs_constraints()
{
    for (size_t i = 0; i < NumStateColumns_L; i++) {
        Flystel[i].generate_r1cs_constraints();
    }
}

template<typename ppT, size_t NumStateColumns_L, class parameters>
void anemoi_permutation_round_prime_field_gadget<
    ppT,
    NumStateColumns_L,
    parameters>::generate_r1cs_witness()
{
    for (size_t i = 0; i < NumStateColumns_L; i++) {
        Flystel[i].generate_r1cs_witness();
    }

    for (size_t i = 0; i < NumStateColumns_L; i++) {
        this->pb.val(Y_left_output[i]) = this->pb.val(Flystel[i].output_y0);
        this->pb.val(Y_right_output[i]) = this->pb.val(Flystel[i].output_y1);
    }
}

// TODO: consdier applying the following changes to all
// anemoi_permutation_mds::permutation_mds functions in order to
// remove the input g parameter:
//
// - extract the ppT part from the anemoi_parameters class
//
// - use the ppT part from the anemoi_parameters class to implicitly
//   get the value of g as
//   anemoi_parameters<ppT>::multiplicative_generator_g;
//
// - remove the input parameter const libff::Fr<ppT> g from all
//   permutation_mds functions and extract g as above
//
// see: https://github.com/clearmatics/libsnark/pull/102#discussion_r1071444422
template<typename ppT>
std::array<std::array<libff::Fr<ppT>, 2>, 2> anemoi_permutation_mds<ppT, 2>::
    permutation_mds(const libff::Fr<ppT> g)
{
    using FieldT = libff::Fr<ppT>;
    const FieldT g2 = g * g;
    anemoi_mds_matrix_t M = {{{1, g}, {g, g2 + 1}}};
    return M;
}

template<typename ppT>
std::array<std::array<libff::Fr<ppT>, 3>, 3> anemoi_permutation_mds<ppT, 3>::
    permutation_mds(const libff::Fr<ppT> g)
{
    anemoi_mds_matrix_t M = {{{g + 1, 1, g + 1}, {1, 1, g}, {g, 1, 1}}};
    return M;
}

template<typename ppT>
std::array<std::array<libff::Fr<ppT>, 4>, 4> anemoi_permutation_mds<ppT, 4>::
    permutation_mds(const libff::Fr<ppT> g)
{
    using FieldT = libff::Fr<ppT>;
    const FieldT g2 = g * g;
    anemoi_mds_matrix_t M = {
        {{1, g + 1, g, g},
         {g2, g + g2, g + 1, g + g + 1},
         {g2, g2, 1, g + 1},
         {g + 1, g + g + 1, g, g + 1}}};
    return M;
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_
