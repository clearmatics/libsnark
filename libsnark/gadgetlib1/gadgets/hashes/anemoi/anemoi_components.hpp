/** @file
 *****************************************************************************

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_

#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>

namespace libsnark
{

#define FLYSTEL_DEBUG
#define FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR 7
// alpha constant used in the Flystel E transformation. note that
// alpha is independent of the choice of the curve, but 1/alpha (see
// below) depends on the curve (mores specifically -- on the modulus r of
// its scalar field Fr)
#define FLYSTEL_ALPHA_FIVE 5
// the mapping f(x)=x^a=y: x,y \in Fr applied in the Flystel E
// transformation (where a is alpha) is invertible if 1/a exists. then
// f^-1(y)=y^1/a=x. 1/a exists if gcd(a,r-1)=1 where r is the modulus
// of Fr. 1/a can be found with the extended Euclidean algorithm which
// finds u,v s.t. ua+v(r-1)=1 mod (r-1)=ua and so u=1/a. parameter
// FLYSTEL_ALPHA_FIVE_INVERSE gives the value of u=1/a for a=5 for the
// curve BLS12-381 precomputed using the Sage command
// inverse_mod(alpha, r-1). TODO: write a function anemoi_parameters()
// specialized by ppT that loads the precomputed constants (including
// alpha and the multiplicative subgroup generator
// FLYSTEL_MULTIPLICATIVE_SUBGROUP_GENERATOR) for any curve
#define FLYSTEL_ALPHA_FIVE_INVERSE                                             \
    "209743500700504761917790962032743863350762210002110551290414634799754324" \
    "73805"

// original constants by specification
// for BLS12-381
// beta = g = first multiplicative generator = 7.
// delta = g^(-1) =
// 14981678621464625851270783002338847382197300714436467949315331057125308909861
// gamma = 0

// constants used for debug
#define DEBUG_FLYSTEL_ALPHA 5
#define DEBUG_FLYSTEL_BETA 2
#define DEBUG_FLYSTEL_GAMMA 5
#define DEBUG_FLYSTEL_DELTA 0

/// Flystel Q_gamma function for prime fields:
/// Qf(x) = beta x^2 + gamma
/// x: input
/// y: output
template<typename FieldT, size_t generator>
class flystel_Q_gamma_prime_field_gadget : public gadget<FieldT>
{
private:
    // constants
    const FieldT beta;
    const FieldT gamma;

public:
    // input/output
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_gamma_prime_field_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Q_delta function for prime fields:
/// Qf(x) = beta x^2 + delta
/// x: input
/// y: output
template<typename FieldT, size_t generator>
class flystel_Q_delta_prime_field_gadget : public gadget<FieldT>
{
private:
    // constants
    const FieldT beta;
    const FieldT delta;

public:
    // input/output
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_delta_prime_field_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Q_gamma function for binary fields:
/// Qi(x) = beta x^3 + gamma
///
/// Compute y = beta x^3 + gamma
/// x: input
/// y: output
/// beta, gamma: constants
template<typename FieldT, size_t generator>
class flystel_Q_gamma_binary_field_gadget : public gadget<FieldT>
{
private:
    /// internal (i.e. intermediate) variable
    pb_variable<FieldT> internal;
    /// constants
    const FieldT beta;
    const FieldT gamma;

public:
    /// input/output
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_gamma_binary_field_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Q_delta function for binary fields:
/// Qi(x) = beta x^3 + delta
///
/// Compute y = beta x^3 + delta
/// x: input
/// y: output
/// beta, delta: constants
template<typename FieldT, size_t generator>
class flystel_Q_delta_binary_field_gadget : public gadget<FieldT>
{
private:
    /// internal (i.e. intermediate) variable
    pb_variable<FieldT> internal;
    /// constants
    const FieldT beta;
    const FieldT delta;

public:
    /// input/output
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_delta_binary_field_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Compute y = x^5
/// x: input
/// y: output
template<typename FieldT>
class flystel_E_power_five_gadget : public gadget<FieldT>
{
private:
    // internal (i.e. intermediate) variables
    pb_variable<FieldT> a0;
    pb_variable<FieldT> a1;

public:
    /// input/output
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_E_power_five_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Compute y = x^1/5
/// x: input
/// y: output
template<typename FieldT>
class flystel_E_root_five_gadget : public gadget<FieldT>
{
private:
    // internal (i.e. intermediate) variables
    pb_variable<FieldT> a0;
    pb_variable<FieldT> a1;

public:
    /// input/output
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_E_root_five_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Anemoi closed Flystel component for fields of prime characteristic
///
/// x0,x1: input (x,y in the paper)
/// y0,y1: output (u,v in the paper)
///
/// The component performs the following computation:
///
/// a0 = (beta x1^2 + gamma) = Q_gamma(x1)
/// a1 = (x0 - a0)^{1/alpha} = E_root_five(x0-a0)
/// a2 = beta (x1-a1)^2 + delta = Q_delta(x1-a1)
/// y0 = x0 - a0 + a2
/// y1 = x1 - a1
///
/// \note: in the paper (x0,x1)->(y0,y1) is denoted with (x,y)->(u,v)
template<typename FieldT, size_t generator>
class flystel_prime_field_gadget : public gadget<FieldT>
{
private:
    // internal (i.e. intermediate) variables
    pb_variable<FieldT> a0;
    pb_variable<FieldT> a1;
    pb_variable<FieldT> a2;

public:
    const linear_combination<FieldT> input_x0;
    const linear_combination<FieldT> input_x1;
    const pb_variable<FieldT> output_y0;
    const pb_variable<FieldT> output_y1;

    flystel_Q_gamma_prime_field_gadget<FieldT, generator> Q_gamma;
    flystel_Q_delta_prime_field_gadget<FieldT, generator> Q_delta;
    flystel_E_root_five_gadget<FieldT> E_root_five;

    flystel_prime_field_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &x0,
        const linear_combination<FieldT> &x1,
        const pb_variable<FieldT> &y0,
        const pb_variable<FieldT> &y1,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

// get the MDS matrix from the number of columns 2,3 or 4
template<typename FieldT, size_t NumStateColumns_L>
std::array<std::array<FieldT, NumStateColumns_L>, NumStateColumns_L>
anemoi_permutation_mds(const FieldT g);

/// One round of the Anemoi permutation mapping (Fr)^{2l} -> (Fr)^{2l}
///
/// NumStateColumns_L : l parameter - number of columns in the
///                     state. can be 1,2,3,4. each column is composed of 2
///                     elements in F_r. One Flystel Sbox accepts 1 column as
///                     input. There are l Flystel-s in 1 round of the
///                     Anemoi permutation applied in parallel.
///
/// x0,x1: input
/// y0,y1: output
///
// template<typename FieldT, FieldT beta, FieldT gamma, FieldT delta>
template<typename FieldT, size_t generator, size_t NumStateColumns_L>
class anemoi_permutation_round_prime_field_gadget : public gadget<FieldT>
{
private:
    // array of C round constants
    std::array<FieldT, NumStateColumns_L> c_const;
    // array of D round constants
    std::array<FieldT, NumStateColumns_L> d_const;
    // matrix M
    std::array<std::array<FieldT, NumStateColumns_L>, NumStateColumns_L> M;
    // array of Flystel S-boxes
    std::array<flystel_prime_field_gadget<FieldT, generator>, NumStateColumns_L>
        flystel;

public:
    std::array<pb_variable<FieldT>, 2 * NumStateColumns_L> input;
    std::array<pb_variable<FieldT>, 2 * NumStateColumns_L> output;

    anemoi_permutation_round_prime_field_gadget(
        protoboard<FieldT> &pb,
        std::array<pb_variable<FieldT>, 2 * NumStateColumns_L> &input,
        std::array<pb_variable<FieldT>, 2 * NumStateColumns_L> &output,
        std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

} // namespace libsnark

#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.tcc>

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
