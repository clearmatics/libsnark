/** @file
 *****************************************************************************

 Declaration of interfaces for top-level SHA256 gadgets.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_

#include <libsnark/common/data_structures/merkle_tree.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/hash_io.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_components.hpp>

namespace libsnark
{

// TODO: template-ize the following constants

#define FLYSTEL_BLS12_381_ALPHA 5
// equals to G1
#define FLYSTEL_BLS12_381_BETA 2
// TODO: value by spec is 0
#define FLYSTEL_BLS12_381_GAMMA 5
// TODO: value by spec is G1.inv()
#define FLYSTEL_BLS12_381_DELTA 0

// --- Prime fields ---

/// Compute y = const_a x^2 + const_b
/// x: input
/// y: output
/// const_a, const_b: constants
template<typename FieldT> class flystel_power_two_gadget : public gadget<FieldT>
{
private:
    // constants
    const FieldT const_a;
    const FieldT const_b;

public:
    // input/output
    const pb_variable<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_power_two_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Qi function for prime fields:
/// Qi(x) = beta x^2 + gamma
template<typename FieldT>
class flystel_Q_gamma_prime_field_gadget
    : public flystel_power_two_gadget<FieldT>
{
    flystel_Q_gamma_prime_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");
};

/// Flystel Qf function for prime fields:
/// Qf(x) = beta x^2 + delta
template<typename FieldT>
class flystel_Q_delta_prime_field_gadget
    : public flystel_power_two_gadget<FieldT>
{
    flystel_Q_delta_prime_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");
};

/// Compute y = x^5
/// x: input
/// y: output
template<typename FieldT>
class flystel_power_five_gadget : public gadget<FieldT>
{
private:
    /// internal (i.e. intermediate) variable: x2,x3
    std::array<pb_variable<FieldT>, 2> internal;

public:
    /// input/output: x1,x4
    const pb_variable<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_power_five_gadget(
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
template<typename FieldT>
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

    flystel_Q_gamma_prime_field_gadget<FieldT> Q_gamma;
    flystel_Q_delta_prime_field_gadget<FieldT> Q_delta;
    flystel_power_five_gadget<FieldT> power_five;

    flystel_closed_prime_field_gadget(
        protoboard<FieldT> &pb,
        const std::array<pb_variable<FieldT>, 2> &input,
        const std::array<pb_variable<FieldT>, 2> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// One round of the Anemoi permutation mapping (Fr)^{2l} -> (Fr)^{2l}
///
/// NumStateColumns_L : l parameter - number of columns in the
///                     state. each column is composed of 2 elements
///                     in F_r. One Flystel Sbox accepts 1 column as
///                     input. There are l Flystel-s in 1 round of the
///                     Anemoi permutation applied in parallel.
///
/// x0,x1: input
/// y0,y1: output
///
// template<typename FieldT, FieldT beta, FieldT gamma, FieldT delta>
template<typename FieldT, size_t NumStateColumns_L>
class anemoi_permutation_round_prime_field_gadget : public gadget<FieldT>
{
private:
    // internal (i.e. intermediate) variables: v3,v4,v5
    std::array<pb_variable<FieldT>, 4> internal;

public:
    // (v1,v2)=(x0,x1)
    std::array<pb_variable<FieldT>, 2> input;
    // (v7,v8)=(y0,y1)
    std::array<pb_variable<FieldT>, 2> output;

    flystel_Q_gamma_prime_field_gadget<FieldT> Q_gamma;
    flystel_Q_delta_prime_field_gadget<FieldT> Q_delta;
    flystel_power_five_gadget<FieldT> power_five;

    anemoi_permutation_round_prime_field_gadget(
        protoboard<FieldT> &pb,
        const std::array<pb_variable<FieldT>, 2> &input,
        const std::array<pb_variable<FieldT>, 2> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

// --- Binary fields ---

/// Compute y = const_a x^3 + const_b
/// x: input
/// y: output
/// const_a, const_b: constants
template<typename FieldT>
class flystel_power_three_gadget : public gadget<FieldT>
{
private:
    /// internal (i.e. intermediate) variable
    pb_variable<FieldT> internal;
    /// constants
    const FieldT const_a;
    const FieldT const_b;

public:
    /// input/output
    const pb_variable<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_power_three_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Qi function for binary fields:
/// Qi(x) = beta x^3 + gamma
template<typename FieldT>
class flystel_Q_gamma_binary_field_gadget
    : public flystel_power_three_gadget<FieldT>
{
    flystel_Q_gamma_binary_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");
};

/// Flystel Qf function for binary fields:
/// Qf(x) = beta x^3 + delta
template<typename FieldT>
class flystel_Q_delta_binary_field_gadget
    : public flystel_power_three_gadget<FieldT>
{
    flystel_Q_delta_binary_field_gadget(
        protoboard<FieldT> &pb,
        const pb_variable<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix = "");
};

} // namespace libsnark

#include <libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.tcc>

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
