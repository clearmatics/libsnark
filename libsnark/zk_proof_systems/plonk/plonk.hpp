/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef PLONK_PPZKSNARK_HPP_
#define PLONK_PPZKSNARK_HPP_

#include <libff/algebra/curves/public_params.hpp>
#include <libsnark/zk_proof_systems/plonk/tests/example.hpp>
#include <memory>

// enable debug checks. in particular enable comparisons to test
// vector values.
#define DEBUG
//#undef DEBUG

/// Declaration of common untility functions for ppzkSNARK proof
/// system Plonk. The implementation instantiates the protocol of
/// PlonK \[GWC19],
///
/// References:
/// - \[GWC19]:
///   Title: "Plonk: Permutations over lagrange-bases for oecumenical
///   noninteractive arguments of knowledge", Ariel Gabizon, Zachary
///   J. Williamson, and Oana Ciobotaru, Cryptology ePrint Archive,
///   Report 2019/953, 2019, <https://eprint.iacr.org/2019/953>

// maximum degree of the encoded monomials in the usrs
const size_t MAX_DEGREE = 254;

// number of generators for H, Hk1, Hk2
const size_t NUM_HGEN = 3;

namespace libsnark
{

enum W_polys_id { a, b, c };
enum Q_polys_id { L, R, M, O, C };
enum t_polys_id { lo, mid, hi };
enum omega_id { base, base_k1, base_k2 };

template<typename FieldT> using polynomial = std::vector<FieldT>;

/***************************** Main algorithms *******************************/

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/plonk.tcc"

#endif // PLONK_PPZKSNARK_HPP_
