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

/// Flystel Q_gamma function for prime fields:
/// Qf(x) = beta x^2 + gamma
template<typename ppT, class parameters = anemoi_parameters<libff::Fr<ppT>>>
class flystel_Q_gamma_prime_field_gadget : public gadget<libff::Fr<ppT>>
{
    using FieldT = libff::Fr<ppT>;

private:
    const FieldT beta;
    const FieldT gamma;

public:
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_gamma_prime_field_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Q_delta function for prime fields:
/// Qf(x) = beta x^2 + delta
template<typename ppT, class parameters = anemoi_parameters<libff::Fr<ppT>>>
class flystel_Q_delta_prime_field_gadget : public gadget<libff::Fr<ppT>>
{
    using FieldT = libff::Fr<ppT>;

private:
    const FieldT beta;
    const FieldT delta;

public:
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_delta_prime_field_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Q_gamma function for binary fields:
/// Qi(x) = beta x^3 + gamma
template<typename ppT, class parameters = anemoi_parameters<libff::Fr<ppT>>>
class flystel_Q_gamma_binary_field_gadget : public gadget<libff::Fr<ppT>>
{
    using FieldT = libff::Fr<ppT>;

private:
    const pb_variable<FieldT> internal;
    const FieldT beta;
    const FieldT gamma;

public:
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_gamma_binary_field_gadget(
        protoboard<FieldT> &pb,
        const linear_combination<FieldT> &input,
        const pb_variable<FieldT> &output,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Flystel Q_delta function for binary fields:
/// Qi(x) = beta x^3 + delta
template<typename ppT, class parameters = anemoi_parameters<libff::Fr<ppT>>>
class flystel_Q_delta_binary_field_gadget : public gadget<libff::Fr<ppT>>
{
    using FieldT = libff::Fr<ppT>;

private:
    const pb_variable<FieldT> internal;
    const FieldT beta;
    const FieldT delta;

public:
    const linear_combination<FieldT> input;
    const pb_variable<FieldT> output;

    flystel_Q_delta_binary_field_gadget(
        protoboard<FieldT> &pb,
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

    flystel_Q_gamma_prime_field_gadget<ppT, parameters> Q_gamma;
    flystel_Q_delta_prime_field_gadget<ppT, parameters> Q_delta;
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
