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

// original constants by specification
// for BLS12-381
// beta = g = first multiplicative generator = 7.
// delta = g^(-1) =
// 14981678621464625851270783002338847382197300714436467949315331057125308909861
// gamma = 0

// constants used for debug
#define DEBUG_FLYSTEL_ALPHA FLYSTEL_ALPHA
#define DEBUG_FLYSTEL_BETA 2
#define DEBUG_FLYSTEL_GAMMA 5
#define DEBUG_FLYSTEL_DELTA 0

/// Flystel Qf function for prime fields:
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
    const pb_variable<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_gamma_prime_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Qf function for prime fields:
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
    const pb_variable<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_delta_prime_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Qi function for binary fields:
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
    const pb_variable<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_gamma_binary_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
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
    /// internal (i.e. intermediate) variable: x2,x3
    std::array<pb_variable<FieldT>, 2> internal;

public:
    /// input/output: x1,x4
    const pb_variable<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_E_power_five_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Anemoi closed Flystel component for fields of prime characteristic
///
/// x0,x1: input (y,v in the paper)
/// y0,y1: output (x,u in the paper)
///
/// The component performs the following computation:
///
/// y0 = (beta x0^2 + gamma) + (x0-x1)^5
/// y1 = (beta x1^2 + delta) + (x0-x1)^5
///
/// Using Q_gamma, Q_delta and power_five gadgets the above is
/// equivalent to
///
/// y0 = Q_gamma(x0) + power_five(x0-x1)
/// y1 = Q_delta(x1) + power_five(x0-x1)
///
/// \note: in the paper (x0,x1)->(y0,y1) is denoted with (y,v)->(x,u)
// template<typename FieldT, FieldT beta, FieldT gamma, DeildT delta>
template<typename FieldT, size_t generator>
class flystel_closed_prime_field_gadget : public gadget<FieldT>
{
private:
    // internal (i.e. intermediate) variables: v3,v4,v5
    std::array<pb_variable<FieldT>, 4> internal;

public:
    // (v1,v2)=(x0,x1)
    std::array<pb_variable<FieldT>, 2> input;
    // (v7,v8)=(y0,y1)
    std::array<pb_variable<FieldT>, 2> output;

    flystel_Q_gamma_prime_field_gadget<FieldT, generator> Q_gamma;
    flystel_Q_delta_prime_field_gadget<FieldT, generator> Q_delta;
    flystel_E_power_five_gadget<FieldT> power_five;

    flystel_closed_prime_field_gadget(
        protoboard<FieldT> &pb,
        const std::array<pb_variable<FieldT>, 2> &input,
        const std::array<pb_variable<FieldT>, 2> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

// get the MDS matrix from the number of columns 2,3 or 4
template<typename FieldT, size_t NumStateColumns_L>
std::array<std::array<FieldT, NumStateColumns_L>, NumStateColumns_L>
anemoi_permutation_get_mds(const FieldT g);

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
    std::array<
        flystel_closed_prime_field_gadget<FieldT, generator>,
        NumStateColumns_L>
        flystel;

public:
    std::array<pb_variable<FieldT>, 2 * NumStateColumns_L> input;
    std::array<pb_variable<FieldT>, 2 * NumStateColumns_L> output;

    anemoi_permutation_round_prime_field_gadget(
        protoboard<FieldT> &pb,
        std::array<pb_variable<FieldT>, (2 * NumStateColumns_L)> &input,
        std::array<pb_variable<FieldT>, (2 * NumStateColumns_L)> &output,
        std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

} // namespace libsnark

#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.tcc>

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
