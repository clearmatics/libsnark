/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_

#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_parameters.hpp"

#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>

/// Implementation of the Anenoi arithmetization-oriented hash function
///
/// Reference:
/// - \[BBCPSVW22]:
///   Title: "New Design Techniques for Efficient
///   Arithmetization-Oriented Hash Functions: Anemoi Permutations and
///   Jive Compression Mode", Clemence Bouvier, Pierre Briaud, Pyrros
///   Chaidos, Leo Perrin, Robin Salen, Vesselin Velichkov, Danny
///   Willems, Cryptology ePrint Archive, Report 2022/840, 2019,
///   <https://eprint.iacr.org/2022/840>

namespace libsnark
{

/// Combined gadget for the Flystel Q-functions for prime fields Q(x) = A x^2 +
/// B:
/// Q_gamma(x) = beta x^2 + gamma: A = beta, B = gamma
/// Q_delta(x) = beta x^2 + delta: A = beta, B = delta
template<typename ppT>
class flystel_Q_prime_field_gadget : public gadget<libff::Fr<ppT>>
{
    using FieldT = libff::Fr<ppT>;

private:
    const FieldT A;
    const FieldT B;

public:
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_prime_field_gadget(
        protoboard<FieldT> &pb,
        const FieldT A,
        const FieldT B,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Combined gadget for the Flystel Q-functions for binary fields Q(x) = A x^3 +
/// B:
/// Q_gamma(x) = beta x^3 + gamma: A = beta, B = gamma
/// Q_delta(x) = beta x^3 + delta: A = beta, B = delta
template<typename ppT>
class flystel_Q_binary_field_gadget : public gadget<libff::Fr<ppT>>
{
    using FieldT = libff::Fr<ppT>;

private:
    const pb_variable<FieldT> internal;
    const FieldT A;
    const FieldT B;

public:
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_binary_field_gadget(
        protoboard<FieldT> &pb,
        const FieldT A,
        const FieldT B,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Compute y = x^5
template<typename ppT>
class flystel_E_power_five_gadget : public gadget<libff::Fr<ppT>>
{
    using FieldT = libff::Fr<ppT>;

private:
    // internal (i.e. intermediate) variables
    const pb_variable<FieldT> a0;
    const pb_variable<FieldT> a1;

public:
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_E_power_five_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Compute y = x^1/5, x=input, y=output/result
template<typename ppT, class parameters = anemoi_parameters<libff::Fr<ppT>>>
class flystel_E_root_five_gadget : public gadget<libff::Fr<ppT>>
{
    using FieldT = libff::Fr<ppT>;

private:
    // internal (i.e. intermediate) variables
    const pb_variable<FieldT> a0;
    const pb_variable<FieldT> a1;

public:
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_E_root_five_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix);

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
/// \note: in [BBCPSVW22] (x0,x1)->(y0,y1) is denoted with (x,y)->(u,v)
template<typename ppT, class parameters = anemoi_parameters<libff::Fr<ppT>>>
class flystel_prime_field_gadget : public gadget<libff::Fr<ppT>>
{
    using FieldT = libff::Fr<ppT>;

private:
    // internal (i.e. intermediate) variables
    const pb_variable<FieldT> a0;
    const pb_variable<FieldT> a1;
    const pb_variable<FieldT> a2;

public:
    const linear_combination<FieldT> input_x0;
    const linear_combination<FieldT> input_x1;
    const pb_variable<FieldT> output_y0;
    const pb_variable<FieldT> output_y1;

    flystel_Q_prime_field_gadget<ppT> Q_gamma;
    flystel_Q_prime_field_gadget<ppT> Q_delta;
    flystel_E_root_five_gadget<ppT, parameters> E_root_five;

    flystel_prime_field_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &x0,
        const linear_combination<FieldT> &x1,
        const pb_variable<FieldT> &y0,
        const pb_variable<FieldT> &y1,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

// get the MDS matrix from the number of columns 2,3 or 4
template<typename FieldT, size_t NumStateColumns_L>
std::array<std::array<FieldT, NumStateColumns_L>, NumStateColumns_L>
anemoi_permutation_mds(const FieldT g);

} // namespace libsnark

#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_components.tcc"

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_COMPONENTS_HPP_
