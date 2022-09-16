/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_CIRCUIT_HPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_CIRCUIT_HPP_

#include "libsnark/zk_proof_systems/plonk/utils.hpp"

#include <libff/algebra/curves/public_params.hpp>

/// Declaration of Common Preprocessed Input data structures for a
/// specific arithmetic circuit.
///
/// This includes:
/// - class for common preprocessed input
///
/// Reference:
/// - \[GWC19]:
///   Title: "Plonk: Permutations over lagrange-bases for oecumenical
///   noninteractive arguments of knowledge", Ariel Gabizon, Zachary
///   J. Williamson, and Oana Ciobotaru, Cryptology ePrint Archive,
///   Report 2019/953, 2019, <https://eprint.iacr.org/2019/953>

namespace libsnark
{

/// Plonk circuit-specific data
template<typename ppT> struct circuit_t {
    using Field = libff::Fr<ppT>;

    /// number of gates / constraints
    size_t num_gates;

    /// number of selector polynomials (q-polynomials) (= 5 in the
    /// vanilla Plonk proposal [GWC19])
    size_t num_qpolys;

    /// Public input polynomial
    polynomial<Field> PI_poly;

    /// Circuit selector polynomials (Q-polynomials)
    std::vector<polynomial<Field>> Q_polys;

    /// Permutation polynomials S_sigma_1, S_sigma_2, S_sigma_2 (see
    /// [GWC19], Sect. 8.1)
    std::vector<polynomial<Field>> S_polys;

    /// omega[0] are the n roots of unity, omega[1] are omega[0]*k1,
    /// omega[2] are omega[0]*k2
    std::vector<std::vector<Field>> omega_roots;

    /// H_gen contains the generators of H, k1 H and k2 H in one place
    /// ie. omega, omega_k1 and omega_k2
    std::vector<Field> H_gen;

    /// H_gen permuted according to the wire permutation
    std::vector<Field> H_gen_permute;

    /// constants for H, k1 H, k2 H
    libff::Fr<ppT> k1;
    libff::Fr<ppT> k2;

    circuit_t(
        size_t num_gates,
        size_t num_qpolys,
        polynomial<Field> &&PI_poly,
        std::vector<polynomial<Field>> &&Q_polys,
        std::vector<polynomial<Field>> &&S_polys,
        std::vector<std::vector<Field>> &&omega_roots,
        std::vector<Field> &&H_gen,
        std::vector<Field> &&H_gen_permute,
        libff::Fr<ppT> &&k1,
        libff::Fr<ppT> &&k2);
};

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/circuit.tcc"

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_CIRCUIT_HPP_
