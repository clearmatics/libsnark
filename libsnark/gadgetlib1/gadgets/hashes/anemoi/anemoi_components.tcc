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
anemoi_power_two_gadget<FieldT>::anemoi_power_two_gadget(
    protoboard<FieldT> &pb,
    const pb_variable<FieldT> &input,
    const pb_variable<FieldT> &output,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , const_a(ANEMOI_BLS12_381_CONST_BETA)
    , const_b(ANEMOI_BLS12_381_CONST_GAMMA)
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
void anemoi_power_two_gadget<FieldT>::generate_r1cs_constraints()
{
    // X variables (input, output, intermediate) x_0=1, x_1, x_2
    std::vector<pb_variable<FieldT>> X{ONE, input, output};
    // A constants a_0, a_1, a_2
    std::vector<FieldT> A{0, const_a, 0};
    // B constants b_0, b_1, b_2
    std::vector<FieldT> B{0, FieldT(1), 0};
    // C constants c_0, c_1, c_2
    std::vector<FieldT> C{-const_b, 0, FieldT(1)};
    // < A , X >
    std::vector<linear_term<FieldT>> A_lc_terms{
        {X[0], A[0]}, {X[1], A[1]}, {X[2], A[2]}};
    linear_combination<FieldT> A_lc(A_lc_terms);
    // < B , X >
    std::vector<linear_term<FieldT>> B_lc_terms{
        {X[0], B[0]}, {X[1], B[1]}, {X[2], B[2]}};
    linear_combination<FieldT> B_lc(B_lc_terms);
    // < C , X >
    std::vector<linear_term<FieldT>> C_lc_terms{
        {X[0], C[0]}, {X[1], C[1]}, {X[2], C[2]}};
    linear_combination<FieldT> C_lc(C_lc_terms);
    // < A , X > * < B , X > = < C , X >
    this->pb.add_r1cs_constraint(
        r1cs_constraint<FieldT>(A_lc, B_lc, C_lc),
        FMT(this->annotation_prefix, " A*B=C"));
}

// compute a witness y for a given input x for the computation y =
// const_a x^2 + const_b, where x=input, y=output
template<typename FieldT>
void anemoi_power_two_gadget<FieldT>::generate_r1cs_witness()
{
    // y = const_a x^2 + const_b
    this->pb.val(output) =
        this->const_a * this->pb.val(input) * this->pb.val(input) +
        this->const_b;
}

template<typename FieldT>
anemoi_power_three_gadget<FieldT>::anemoi_power_three_gadget(
    protoboard<FieldT> &pb,
    const pb_variable<FieldT> &input,
    const pb_variable<FieldT> &output,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , internal(pb_variable_allocate<FieldT>(
          pb, FMT(this->annotation_prefix, " internal")))
    , const_a(ANEMOI_BLS12_381_CONST_BETA)
    , const_b(ANEMOI_BLS12_381_CONST_GAMMA)
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
void anemoi_power_three_gadget<FieldT>::generate_r1cs_constraints()
{
    // X variables (input, output, internal) x_0=1, x_1, x_2, x_3
    std::vector<pb_variable<FieldT>> X{ONE, input, internal, output};
    // A constants
    std::vector<std::vector<FieldT>> A{{0, const_a, 0, 0}, {0, 0, 1, 0}};
    // B constants
    std::vector<std::vector<FieldT>> B{{0, 1, 0, 0}, {0, 1, 0, 0}};
    // C constants
    std::vector<std::vector<FieldT>> C{{0, 0, 1, 0}, {-const_b, 0, 0, 1}};
    // add j-th R1CS constraint
    for (size_t j = 0; j < 2; ++j) {
        std::vector<linear_term<FieldT>> A_lc_terms, B_lc_terms, C_lc_terms;
        for (size_t i = 0; i < X.size(); ++i) {
            A_lc_terms.push_back({X[i], A[j][i]});
            B_lc_terms.push_back({X[i], B[j][i]});
            C_lc_terms.push_back({X[i], C[j][i]});
        }
        // < Aj , X >
        linear_combination<FieldT> A_lc(A_lc_terms);
        // < Bj , X >
        linear_combination<FieldT> B_lc(B_lc_terms);
        // < Cj , X >
        linear_combination<FieldT> C_lc(C_lc_terms);
        // < Aj , X > * < Bj , X > = < Cj , X >
        this->pb.add_r1cs_constraint(
            r1cs_constraint<FieldT>(A_lc, B_lc, C_lc),
            FMT(this->annotation_prefix, " Aj*Bj=Cj"));
    }
}

template<typename FieldT>
void anemoi_power_three_gadget<FieldT>::generate_r1cs_witness()
{
    // x_internal = const_a x * x
    this->pb.val(internal) =
        (this->const_a * this->pb.val(input)) * this->pb.val(input);
    // y = const_a x^3 + const_b = x_internal * x + const_b
    this->pb.val(output) =
        this->pb.val(internal) * this->pb.val(input) + this->const_b;
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_TCC_
