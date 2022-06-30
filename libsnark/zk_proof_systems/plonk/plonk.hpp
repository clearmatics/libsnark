/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_PLONK_HPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_PLONK_HPP_

#include "libsnark/zk_proof_systems/plonk/tests/example.hpp"

#include <libff/algebra/curves/public_params.hpp>
#include <memory>

// enable debug checks. in particular enable comparisons to test
// vector values.
#define DEBUG_PLONK
//#undef DEBUG_PLONK

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

/// maximum degree of the encoded monomials in the usrs
const size_t MAX_DEGREE = 254;

/// number of generators for H, Hk1, Hk2
const size_t NUM_HGEN = 3;

namespace libsnark
{

enum W_polys_id { a, b, c };
enum Q_polys_id { L, R, M, O, C };
enum t_polys_id { lo, mid, hi };
enum omega_id { base, base_k1, base_k2 };

template<typename FieldT> using polynomial = std::vector<FieldT>;

/***************************** Main algorithms *******************************/

/// print the elements of a vector
template<typename FieldT> void print_vector(const std::vector<FieldT> &v);

/// Compute the Lagrange basis polynomials for interpolating sets of
/// n points
///
/// INPUT:
///
/// \param[in] npoints - number of points
///
/// OUTPUT:
///
/// \param[out] L[0..n-1][0..n-1]: Lagrange basis over the n roots of
///             unity omega_0, ..., omega_{n-1} i.e. L[omega_i] = [a0,
///             a1, ..., a_{n-1}] is a vector representing the
///             coefficients of the i-th Lagrange polynomial L_i(x) =
///             a0+a1x+a2x^2+..+a_{n-1}x^{n-1} s.t. L_i(x=omega_i)=1
///             and L_i(x\neq{omega_i)})=0
///
/// Note: uses libfqfft iFFT for the interpolation
template<typename FieldT>
void plonk_compute_lagrange_basis(
    const size_t npoints, std::vector<polynomial<FieldT>> &L);

/// Interpolate a polynomial from a set of points through inverse FFT
///
/// INPUT:
///
/// \param[in] f_points[0..n-1]: a set of points (0,y0), (1,y1),
///            ... (n-1,y_{n-1}) s.t. y0=f_points[0], y1=f_points[1],
///            ... which we want to interpolate as a polynomial
///
/// OUTPUT:
///
/// \param[out] f_poly[0..n-1]: the coefficients [a0, a1, ..., a_{n-1}]
///             of the polynomial f(x) interpolating the set of points
///             f_points. For example if f_poly[0..n-1] = [a0, a1, ...,
///             a_{n-1}] then this represents the polynomial f(x) =
///             a0+a1x+a1x^2+...+a_{n-1}x^{n-1} such that
///             f(omega_i)=f_points[i], where omega_0, ..., omega_{n-1}
///             are the n roots of unity.
///
/// \note uses libfqfft iFFT for the interpolation
template<typename FieldT>
void plonk_interpolate_polynomial_from_points(
    const std::vector<FieldT> &f_points, polynomial<FieldT> &f_poly);

/// Interpolate a polynomial from a set of points over Lagrange basis
///
/// INPUT:
///
/// \param[in] f_points[0..n-1]: a set of points (0,y0), (1,y1),
///            ... (n-1,y_{n-1}) s.t. y0=f_points[0], y1=f_points[1],
///            ... which we want to interpolate as a polynomial
///
/// \param[in] L[0..n-1][0..n-1]: Lagrange basis over the n roots of
///            unity omega_0, ..., omega_{n-1} i.e. L[omega_i] = [a0,
///            a1, ..., a_{n-1}] is a vector representing the
///            coefficients of the i-th Lagrange polynomial L_i(x) =
///            a0+a1x+a2x^2+..+a_{n-1}x^{n-1} s.t. L_i(x=omega_i)=1
///            and L_i(x\neq{omega_i)})=0
///
/// OUTPUT:
///
/// \param[out] f_poly[0..n-1]: the coefficients [a0, a1, ...,
///             a_{n-1}] of the polynomial f(x) interpolating the set
///             of points f_points over the Lagrange basis. For
///             example if f_poly[0..n-1] = [a0, a1, ..., a_{n-1}]
///             then this represents the polynomial f(x) =
///             \sum^{n-1}_{i=0} f_vec[i] * L[i] =
///             a0+a1x+a1x^2+...+a_{n-1}x^{n-1} such that
///             f(omega_i)=f_vec[i].
template<typename FieldT>
void plonk_interpolate_over_lagrange_basis(
    const std::vector<FieldT> &f_points,
    const std::vector<polynomial<FieldT>> &L,
    polynomial<FieldT> &f_poly);

/// Compute the selector polynomials of the given circuit (also
/// called here "q-polynomials"). See Sect. 8.1.  The matrix
/// gates_matrix_transpose has 5 rows, each corresponding to the
/// values L, R, M, O and C for each gate; the number of columns is
/// equal to the number of gates. L_basis is the Lagrange basis.
template<typename FieldT>
void plonk_compute_selector_polynomials(
    const std::vector<std::vector<FieldT>> &gates_matrix_transpose,
    std::vector<polynomial<FieldT>> &Q_polys);

/// output: omega[0] are the n roots of unity, omega[1] are
/// omega[0]*k1, omega[2] are omega[0]*k2; k1 H is a coset of H with
/// generator omega_k1 distinct from H; k2 H is a coset of H with
/// generator omega_k2, distinct from H and k1 H
template<typename FieldT>
void plonk_compute_roots_of_unity_omega(
    const FieldT k1, const FieldT k2, std::vector<std::vector<FieldT>> &omega);

/// copy the roots of unity omega[base], omega[k1] and omega[k2] into
/// a single vector H_gen representing the multiplicative subgroups
/// H, k1 H and k2 H (see [GWC19] Sect. 8); \see
/// plonk_compute_roots_of_unity_omega
template<typename FieldT>
void plonk_roots_of_unity_omega_to_subgroup_H(
    const std::vector<std::vector<FieldT>> &omega, std::vector<FieldT> &H_gen);

/// permute the multiplicative subgroup H according to the wire
/// permutation: (see [GWC19] Sect. 8), \see
/// plonk_compute_roots_of_unity_omega, \see
/// plonk_roots_of_unity_omega_to_subgroup_H
template<typename FieldT>
void plonk_permute_subgroup_H(
    const std::vector<FieldT> &H_gen,
    const std::vector<size_t> &wire_permutation,
    std::vector<FieldT> &H_gen_permute);

/// compute the permutation polynomials S_sigma_1, S_sigma_2,
/// S_sigma_2 (see [GWC19], Sect. 8.1)
template<typename FieldT>
void plonk_compute_permutation_polynomials(
    const std::vector<FieldT> &H_gen_permute,
    const size_t num_gates,
    std::vector<polynomial<FieldT>> &S_polys);

// A wrapper for multi_exp_method_BDLO12_signed() dot-product a
// vector of group elements in G1 (curve points) with a vector of
// scalar elements in Fr
template<typename ppT>
libff::G1<ppT> plonk_multi_exp_G1(
    const std::vector<libff::G1<ppT>> &curve_points,
    const std::vector<libff::Fr<ppT>> &scalar_elements);

/// Evaluate a polynomial F at the encrypted secret input
/// \secret^i*G_1 ie. compute f(\secret)*G1 = [f(\secret)]_i
///
/// INPUT
///
/// \param[in] secret_powers_g1: \secret^i*G1:
///            0\le{i}<max_degree(Q[j]): 0\le{j}<n f_poly: a
///            polynomial
///
/// OUTPUT
///
/// \param[out] [f_poly(\secret)]_1, 0\le 1<n : the "encrypted"
///             evaluation of the polynomial f_poly in the secret
///             parameter \secret (the toxic waste) multiplied by the
///             group generator G_1 i.e. compute f_poly(\secret)*G_1
template<typename ppT>
libff::G1<ppT> plonk_evaluate_poly_at_secret_G1(
    const std::vector<libff::G1<ppT>> &secret_powers_g1,
    const polynomial<libff::Fr<ppT>> &f_poly);

/// Evaluate a list of polynomials in the encrypted secret input: see
/// plonk_evaluate_poly_at_secret_G1
template<typename ppT>
void plonk_evaluate_polys_at_secret_G1(
    const std::vector<libff::G1<ppT>> &secret_powers_g1,
    const std::vector<polynomial<libff::Fr<ppT>>> &Q_polys,
    std::vector<libff::G1<ppT>> &Q_polys_at_secret_g1);

/// Compute the factors in the product of the permutation polynomial
/// z(X) in Prover Round 2. Note that accumulator A[0]=1 and A[i],
/// i>0 is computed from values at i-1 for witness[i-1], H_gen[i-1],
/// H_gen_permute[i-1]m etc.
template<typename FieldT>
FieldT plonk_compute_accumulator_factor(
    const size_t i,
    const size_t n, // nconstraimts
    const FieldT beta,
    const FieldT gamma,
    const std::vector<FieldT> &witness,
    const std::vector<FieldT> &H_gen, // H, Hk1, Hk2
    const std::vector<FieldT> &H_gen_permute,
    const std::vector<FieldT> &A);

/// A: accumulatro vector
template<typename FieldT>
void plonk_compute_accumulator(
    const size_t n, // num_gates
    const FieldT beta,
    const FieldT gamma,
    const std::vector<FieldT> &witness,
    const std::vector<FieldT> &H_gen, // H, Hk1, Hk2
    const std::vector<FieldT> &H_gen_permute,
    std::vector<FieldT> &A);

/// check that the input is an element of the field
template<typename FieldT> bool check_field_element(const FieldT x);

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/plonk.tcc"

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_PLONK_HPP_
