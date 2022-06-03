/** @file
*****************************************************************************
Implementation of interfaces for a ppzkSNARK for Plonk.

See plonk.hpp .

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_TCC_
#define PLONK_PPZKSNARK_TCC_

namespace libsnark
{

template<typename FieldT> void print_vector(std::vector<FieldT> v)
{
    for (size_t i = 0; i < v.size(); ++i) {
        printf("[%2d]: ", (int)i);
        v[i].print();
    }
}

//
// Compute the Lagrange basis polynomails for interpolating sets of
// n points
//
// INPUT:
//
// - npoints - number of points
//
// OUTPUT:
//
// - L[0..n-1][0..n-1]: Lagrange basis over the n roots of unity
//                      omega_0, ..., omega_{n-1} i.e. L[omega_i] =
//                      [a0, a1, ..., a_{n-1}] is a vector
//                      representing the coefficients of the i-th
//                      Lagrange polynomial L_i(x) =
//                      a0+a1x+a2x^2+..+a_{n-1}x^{n-1}
//                      s.t. L_i(x=omega_i)=1 and
//                      L_i(x\neq{omega_i)})=0
//
// Note: uses libfqfft iFFT for the interpolation
//
template<typename FieldT>
void plonk_compute_lagrange_basis(
    size_t npoints, std::vector<polynomial<FieldT>> &L)
{
    assert(L.size() != 0);
    assert(L.size() == L[0].size());
    assert(L.size() == npoints);

    std::shared_ptr<libfqfft::evaluation_domain<FieldT>> domain =
        libfqfft::get_evaluation_domain<FieldT>(npoints);
    for (size_t i = 0; i < npoints; ++i) {
        polynomial<FieldT> u(npoints, FieldT(0));
        u[i] = FieldT(1);
        // compute i-th Lagrange basis vector via inverse FFT
        domain->iFFT(u);
        L[i] = u;
    }
}

//
// Interpolate a set of points as a polynomial over Lagrange basis
//
// INPUT:
//
// - f_points[0..n-1]: a set of points (0,y0), (1,y1),
//                     ... (n-1,y_{n-1}) s.t. y0=f_points[0],
//                     y1=f_points[1], ... which we want to interpolate
//                     as a polynomial
//
// - L[0..n-1][0..n-1]: Lagrange basis over the n roots of unity
//                      omega_0, ..., omega_{n-1} i.e. L[omega_i] =
//                      [a0, a1, ..., a_{n-1}] is a vector
//                      representing the coefficients of the i-th
//                      Lagrange polynomial L_i(x) =
//                      a0+a1x+a2x^2+..+a_{n-1}x^{n-1}
//                      s.t. L_i(x=omega_i)=1 and
//                      L_i(x\neq{omega_i)})=0
//
// OUTPUT:
//
// - f_poly[0..n-1]: the coefficients [a0, a1, ..., a_{n-1}] of the
//                   polynomial f(x) interpolating the set of points
//                   f_points over the Lagrange basis. For example if
//                   f_poly[0..n-1] = [a0, a1, ..., a_{n-1}] then
//                   this represents the polynomial f(x) =
//                   \sum^{n-1}_{i=0} f_vec[i] * L[i] =
//                   a0+a1x+a1x^2+...+a_{n-1}x^{n-1} such that
//                   f(omega_i)=f_vec[i].
//
template<typename FieldT>
void plonk_interpolate_over_lagrange_basis(
    std::vector<FieldT> f_points,
    polynomial<FieldT> &f_poly,
    std::vector<polynomial<FieldT>> L)
{
    assert(L.size() != 0);
    assert(L.size() == L[0].size());
    assert(L.size() == f_points.size());

    size_t num_gates = L.size();

    // the num_gates components of f_poly. each components is
    // a polynomial that is the multiplication of the i-th element of
    // f_points to the i-th polynomial in the lagrange basis L[i]:
    // f_poly_component[i] = f_points[i] * L[i]
    std::vector<polynomial<FieldT>> f_poly_component(num_gates);
    for (size_t i = 0; i < num_gates; ++i) {
        // represent the scalar f_points[i] as an all-zero vector with
        // only the first element set to f_points[i]. this is done in
        // order to use libfqfft multiplication function as a function
        // for multiply vector by scalar.
        std::vector<FieldT> f_points_coeff_i(num_gates, FieldT(0));
        f_points_coeff_i[0] = f_points[i];
        // f_poly_component[i] = f_points[i] * L[i]
        libfqfft::_polynomial_multiplication<FieldT>(
            f_poly_component[i], f_points_coeff_i, L[i]);
    }
    std::fill(f_poly.begin(), f_poly.end(), FieldT(0));
    // f_poly[i] = \sum_i (f_points[i] * L[i])
    for (size_t i = 0; i < num_gates; ++i) {
        // f_poly[i] = f_poly[i] + f_poly_component[i];
        libfqfft::_polynomial_addition<FieldT>(
            f_poly, f_poly, f_poly_component[i]);
    }
}

// Compute the selector polynomials of the given circuit (also
// called here "q-polynomials"). See Sect. 8.1.  The matrix
// gates_matrix_transpose has 5 rows, each corresponding to the
// values L, R, M, O and C for each gate; the number of columns is
// equal to the number of gates. L_basis is the Lagrange basis.
template<typename FieldT>
void plonk_compute_selector_polynomials(
    std::vector<std::vector<FieldT>> gates_matrix_transpose,
    std::vector<polynomial<FieldT>> &Q_polys,
    std::vector<polynomial<FieldT>> L_basis)
{
    assert(gates_matrix_transpose.size() == Q_polys.size());
    assert(gates_matrix_transpose[0].size() == Q_polys[0].size());
    assert(gates_matrix_transpose[0].size() == L_basis.size());
    size_t num_qpolys = gates_matrix_transpose.size();
    for (size_t i = 0; i < num_qpolys; ++i) {
        std::vector<FieldT> q_vec = gates_matrix_transpose[i];
        plonk_interpolate_over_lagrange_basis<FieldT>(
            q_vec, Q_polys[i], L_basis);
    }
};

template<typename FieldT>
void plonk_compute_public_input_polynomial(
    std::vector<FieldT> PI_points,
    polynomial<FieldT> &PI_poly,
    std::vector<polynomial<FieldT>> L_basis)
{
    assert(PI_points.size() == L_basis.size());
    plonk_interpolate_over_lagrange_basis<FieldT>(PI_points, PI_poly, L_basis);
};

// output: omega[0] are the n roots of unity, omega[1] are
// omega[0]*k1, omega[2] are omega[0]*k2; k1 H is a coset of H with
// generator omega_k1 distinct from H; k2 H is a coset of H with
// generator omega_k2, distinct from H and k1 H
template<typename FieldT>
void plonk_compute_roots_of_unity_omega(
    const FieldT k1, const FieldT k2, std::vector<std::vector<FieldT>> &omega)
{
    assert(omega.size() == NUM_HGEN);
    assert(omega[0].size() > 0);
    size_t num_gates = omega[0].size();
    // Get the n-th root of unity omega in Fq (n=8 is the number of
    // constraints in the example). omega is a generator of the
    // multiplicative subgroup H.  Example (2**32)-th primitive root
    // of unity in the base field Fq of bls12-381 i.e. such that
    // omega_base**(2**32) = 1. The bls12-381 prime q is such that any
    // power of 2 divides (q-1). In particular 2**32|(q-1) and so the
    // 2**32-th root of unity exists.
    FieldT omega_base = libff::get_root_of_unity<FieldT>(num_gates);
#ifdef DEBUG
    // assert that omega_base is a 2^32-th root of unity in Fq
    omega_base.print();
    FieldT temp = libff::power(omega_base, libff::bigint<1>(std::pow(2, 32)));
    assert(temp == 1);
#endif // #ifdef DEBUG
    for (int i = 0; i < (int)num_gates; ++i) {
        FieldT omega_i = libff::power(omega_base, libff::bigint<1>(i));
        omega[base][i] = omega_i;
        FieldT omega_k1_i = omega[base][i] * k1;
        FieldT omega_k2_i = omega[base][i] * k2;
        omega[base_k1][i] = omega_k1_i;
        omega[base_k2][i] = omega_k2_i;
    }
}

// copy the roots of unity omega[base], omega[k1] and omega[k2] into
// a single vector H_gen representing the multiplicative subgroups
// H, k1 H and k2 H (see [GWC19] Sect. 8); \see
// plonk_compute_roots_of_unity_omega
template<typename FieldT>
void plonk_roots_of_unity_omega_to_subgroup_H(
    std::vector<std::vector<FieldT>> omega, std::vector<FieldT> &H_gen)
{
    assert(omega.size() == NUM_HGEN);
    assert(omega[0].size() > 0);
    std::copy(omega[base].begin(), omega[base].end(), back_inserter(H_gen));
    std::copy(
        omega[base_k1].begin(), omega[base_k1].end(), back_inserter(H_gen));
    std::copy(
        omega[base_k2].begin(), omega[base_k2].end(), back_inserter(H_gen));
}

// permute the multiplicative subgroup H accordingto the wire
// permutation: (see [GWC19] Sect. 8), \see
// plonk_compute_roots_of_unity_omega, \see
// plonk_roots_of_unity_omega_to_subgroup_H
template<typename FieldT>
void plonk_permute_subgroup_H(
    std::vector<FieldT> &H_gen_permute,
    std::vector<FieldT> H_gen,
    std::vector<size_t> wire_permutation)
{
    assert(H_gen.size() > 0);
    for (size_t i = 0; i < H_gen.size(); ++i) {
        H_gen_permute[i] = H_gen[wire_permutation[i] - 1];
    }
}
// compute the permutation polynomials S_sigma_1, S_sigma_2,
// S_sigma_2 (see [GWC19], Sect. 8.1)
template<typename FieldT>
void plonk_compute_permutation_polynomials(
    std::vector<polynomial<FieldT>> &S_polys,
    std::vector<FieldT> H_gen_permute,
    std::vector<polynomial<FieldT>> L_basis,
    size_t num_gates)
{
    assert(S_polys.size() == NUM_HGEN);
    assert(S_polys[0].size() == num_gates);
    assert(H_gen_permute.size() == (NUM_HGEN * num_gates));
    for (size_t i = 0; i < NUM_HGEN; ++i) {
        typename std::vector<FieldT>::iterator begin =
            H_gen_permute.begin() + (i * num_gates);
        typename std::vector<FieldT>::iterator end =
            H_gen_permute.begin() + (i * num_gates) + (num_gates);
        std::vector<FieldT> S_points(begin, end);
        plonk_interpolate_over_lagrange_basis<FieldT>(
            S_points, S_polys[i], L_basis);
    }
}

// A wrapper for multi_exp_method_BDLO12_signed() for multiplying a
// single group element in G1 (a point on the curve) to a scalar
// element in Fr
template<typename ppT>
libff::G1<ppT> plonk_exp_G1(
    libff::G1<ppT> curve_point, libff::Fr<ppT> scalar_element)
{
    std::vector<libff::Fr<ppT>> scalar{scalar_element};
    std::vector<libff::G1<ppT>> point{curve_point};
    const size_t chunks = 1;
    libff::G1<ppT> product = libff::multi_exp<
        libff::G1<ppT>,
        libff::Fr<ppT>,
        libff::multi_exp_method_BDLO12_signed>(
        point.begin(), point.end(), scalar.begin(), scalar.end(), chunks);
    return product;
}

// A wrapper for multi_exp_method_BDLO12_signed() dot-product a
// vector of group elements in G1 (curve points) with a vector of
// scalar elements in Fr
template<typename ppT>
libff::G1<ppT> plonk_multi_exp_G1(
    std::vector<libff::G1<ppT>> curve_points,
    std::vector<libff::Fr<ppT>> scalar_elements)
{
    assert(curve_points.size() == scalar_elements.size());
    const size_t chunks = 1;
    libff::G1<ppT> product = libff::multi_exp<
        libff::G1<ppT>,
        libff::Fr<ppT>,
        libff::multi_exp_method_BDLO12_signed>(
        curve_points.begin(),
        curve_points.end(),
        scalar_elements.begin(),
        scalar_elements.end(),
        chunks);
    return product;
}

//
// Evaluate a polynomial F at the encrypted secret input
// \secret^i*G_1 ie. compute f(\secret)*G1 = [f(\secret)]_i
//
// INPUT
//
// secret_powers_g1: \secret^i*G1: 0\le{i}<max_degree(Q[j]): 0\le{j}<n
// f_poly: a polynomial
//
// OUTPUT
//
// [f_poly(\secret)]_1, 0\le 1<n : the "encrypted" evaluation of
// the polynomial f_poly in the secret parameter \secret (the
// toxic waste) multiplied by the group generator G_1 i.e. compute
// f_poly(\secret)*G_1
//
template<typename ppT>
libff::G1<ppT> plonk_evaluate_poly_at_secret_G1(
    std::vector<libff::G1<ppT>> secret_powers_g1,
    polynomial<libff::Fr<ppT>> f_poly)
{
    assert(f_poly.size() <= secret_powers_g1.size());
    const size_t chunks = 1;
    const size_t num_coefficients = f_poly.size();
    libff::G1<ppT> f_poly_at_secret_g1 = libff::multi_exp<
        libff::G1<ppT>,
        libff::Fr<ppT>,
        libff::multi_exp_method_BDLO12_signed>(
        secret_powers_g1.begin(),
        secret_powers_g1.begin() + num_coefficients,
        f_poly.begin(),
        f_poly.end(),
        chunks);
    return f_poly_at_secret_g1;
}

// Evaluate a list of polynomials in the encrypted secret input: see
// plonk_evaluate_poly_at_secret_G1
template<typename ppT>
void plonk_evaluate_polys_at_secret_G1(
    std::vector<libff::G1<ppT>> secret_powers_g1,
    std::vector<polynomial<libff::Fr<ppT>>> Q_polys,
    std::vector<libff::G1<ppT>> &Q_polys_at_secret_g1)
{
    assert(secret_powers_g1.size());
    assert(Q_polys.size());
    assert(Q_polys.size() <= secret_powers_g1.size());
    const size_t npolys = Q_polys.size();
    for (size_t i = 0; i < npolys; ++i) {
        Q_polys_at_secret_g1[i] =
            plonk_evaluate_poly_at_secret_G1<ppT>(secret_powers_g1, Q_polys[i]);
    }
}

// Compute the factors in the product of the permutation polynomial
// z(X) in Prover Round 2. Note that accumulator A[0]=1 and A[i],
// i>0 is computed from values at i-1 for witness[i-1], H_gen[i-1],
// H_gen_permute[i-1]m etc.
template<typename FieldT>
FieldT plonk_compute_accumulator_factor(
    size_t i,
    size_t n, // nconstraimts
    FieldT beta,
    FieldT gamma,
    std::vector<FieldT> witness,
    std::vector<FieldT> H_gen, // H, Hk1, Hk2
    std::vector<FieldT> H_gen_permute,
    std::vector<FieldT> &A // accumulatro vector
)
{
    assert(n);
    assert((i >= 0) && (i < n));
    assert(witness.size() == (3 * n));
    assert(H_gen.size() == (3 * n));
    assert(H_gen_permute.size() == (3 * n));
    assert(A.size() == n);
    FieldT res = FieldT(1);
    if (i > 0) {
        FieldT nom_1 = witness[i - 1] + (beta * H_gen[i - 1]) + gamma;
        FieldT den_1 = witness[i - 1] + (beta * H_gen_permute[i - 1]) + gamma;

        FieldT nom_2 = witness[n + i - 1] + (beta * H_gen[n + i - 1]) + gamma;
        FieldT den_2 =
            witness[n + i - 1] + (beta * H_gen_permute[n + i - 1]) + gamma;

        FieldT nom_3 =
            witness[2 * n + i - 1] + (beta * H_gen[2 * n + i - 1]) + gamma;
        FieldT den_3 = witness[2 * n + i - 1] +
                       (beta * H_gen_permute[2 * n + i - 1]) + gamma;

        FieldT nom = nom_1 * nom_2 * nom_3;
        FieldT den = den_1 * den_2 * den_3;

        res = nom * den.inverse() * A[i - 1];
    }
    return res;
}

template<typename FieldT>
void plonk_compute_accumulator(
    size_t n, // num_gates
    FieldT beta,
    FieldT gamma,
    std::vector<FieldT> witness,
    std::vector<FieldT> H_gen, // H, Hk1, Hk2
    std::vector<FieldT> H_gen_permute,
    std::vector<FieldT> &A // accumulatro vector
)
{
    assert(n);
    assert(witness.size() == (3 * n));
    assert(H_gen.size() == (3 * n));
    assert(H_gen_permute.size() == (3 * n));
    assert(A.size() == n);
    for (size_t i = 0; i < n; ++i) {
        A[i] = plonk_compute_accumulator_factor(
            i, n, beta, gamma, witness, H_gen, H_gen_permute, A);
    }
}

//
// Note: the following function was copied from
// libff/algebra/curves/tests/test_groups.cpp TODO: maybe add the
// declaration in a header to use it directly from here. Uded in the
// Verifier code.
//
// Check that a group element satisfies the curve equation Y^2 = X^3
// + a X + b
//
template<typename GroupT> bool check_curve_equation(GroupT P)
{
    P.to_affine_coordinates();
    using Fq = typename std::decay<decltype(P.X)>::type;
    Fq LHS = (P.Y * P.Y);
    Fq RHS = ((P.X * P.X * P.X) + (GroupT::coeff_a * P.X) + GroupT::coeff_b);
    bool b_equal = (LHS == RHS);
    return b_equal;
}

//
// check that the input is an elemnet of the field
// Warning: under development
//
template<typename FieldT> bool check_field_element(FieldT x)
{
    bool b_valid = (typeid(x) == typeid(FieldT));
    return b_valid;
}

} // namespace libsnark

#endif // PLONK_PPZKSNARK_TCC_
