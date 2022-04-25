/** @file
 *****************************************************************************
 Test program that exercises the Plonk protocol (first setup, then
 prover, then verifier) on a synthetic R1CS instance.

 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <gtest/gtest.h>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>
#include <libff/algebra/curves/curve_serialization.hpp>

#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>
#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>

#define DEBUG 1
const size_t MAX_DEGREE = 254;

/*

example circuit

P(x) = x**3 + x + 5 = 3

circuit has 6 gates + 2 dummy gates (to make power of 2)

gates / constraints

        L   R    O
1 mul   x * x  = v1
2 mul  v1 * x  = v2
3 add  v2 + x  = v3
4 con5  1 * 5  = 5
5 pin   1 * 35 = 35
6 add  v3 + 5  = 35
7 dum   /    /    /
8 dum   /    /    /

wire polynomials

w_L = [ x, v1, v2,  1,  1, v3,  /,  /] = [a1, a2, a3, a4, a5, a6, a7, a8] = a
w_R = [ x,  x,  x,  5, 35,  5,  /,  /] = [b1, b2, b3, b4, b5, b6, b7, b8] = b
w_O = [v1, v2, v3,  5, 35, 35,  /,  /] = [c1, c2, c3, c4, c5, c6, c7, c8] = c

wires = [a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7, c8]
index = [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
perm  = [ 9, 17, 18,  5,  4, 19,  7,  8, 10, 11,  1, 14, 21, 20, 15, 16,  2,  3,  6, 12, 22, 13, 23, 24]

witness

x = 3 => v1 = 9, v2 = 27, v3 = 30

w_L = a = [ 3,  9, 27,  1,  1, 30,  0,  0]
w_R = b = [ 3,  3,  3,  5, 35,  5,  0,  0]
w_O = c = [ 9, 27, 30,  5, 35, 35,  0,  0]

W = w_L + w_R + w_O

encoding of plonk gates

add  = [1, 1,-1, 0, 0]
mul  = [0, 0,-1, 1, 0]
con5 = [0, 1, 0, 0, 5]
pi   = [0, 1, 0, 0, 0]
dum  = [0, 0, 0, 0, 0]

gates matrix

(q_L * a) + (q_R * b) + (q_O * c) + (q_M * a * b) + (q_C) = 0

gates                  q_L q_R q_O q_M q_C
1 mul   x * x  = v1 : [  0,  0, -1,  1,  0]
2 mul  v1 * x  = v2 : [  0,  0, -1,  1,  0]
3 add  v2 + x  = v3 : [  1,  1, -1,  0,  0]
4 con5  1 * 5  = 5  : [  0,  1,  0,  0, -5]
5 pin   1 * 35 = 35 : [  0,  1,  0,  0,  0]
6 add  v2 + 5  = 35 : [  1,  1, -1,  0,  0]
7 dum   /    /    / : [  0,  0,  0,  0,  0]
8 dum   /    /    / : [  0,  0,  0,  0,  0]

q_L = [ 0,  0,  1,  0,  0,  1,  0,  0]
q_R = [ 0,  0,  1,  1,  1,  1,  0,  0]
q_O = [-1, -1, -1,  0,  0,  0,  0,  0]
q_M = [ 1,  1,  0,  0,  0, -1,  0,  0]
q_C = [ 0,  0,  0, -1,  0,  0,  0,  0]

setup

input: 
- gate matrix
- permutation vector
- public input

output: 
- SRS
- gate polynomials Q: q_L, q_R, q_O, q_M, q_C
- public input (PI) polynomial
- permutation polynomials S_1, S_2, S_3
- permutation precomputation: id_doman, perm_domain

*/
  
namespace libsnark
{

  template<typename FieldT> using polynomial = std::vector<FieldT>; // kzg10.hpp

  //  void compute_qpolynomials();
  //  void compute_wire_permutation();

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
				    size_t npoints,
				    std::vector<polynomial<FieldT>>& L
				    )
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
					     std::vector<polynomial<FieldT>> L,
					     std::vector<FieldT> f_points,
					     polynomial<FieldT> &f_poly								       
					     )
  {
    assert(L.size() != 0);
    assert(L.size() == L[0].size());
    assert(L.size() == f_points.size());

    size_t nconstraints = L.size();
    
    // the nconstraints components of f_poly. each components is
    // a polynomial that is the multiplication of the i-th element of
    // f_points to the i-th polynomial in the lagrange basis L[i]:
    // f_poly_component[i] = f_points[i] * L[i]
    std::vector<polynomial<FieldT>> f_poly_component(nconstraints);
    for (size_t i = 0; i < nconstraints; ++i) {
      // represent the scalar f_points[i] as an all-zero vector with
      // only the first element set to f_points[i]. this is done in
      // order to use libfqfft multiplication function as a function
      // for multiply vector by scalar.
      std::vector<FieldT> f_points_coeff_i(nconstraints, FieldT(0));
      f_points_coeff_i[0] = f_points[i];
      // f_poly_component[i] = f_points[i] * L[i]
      libfqfft::_polynomial_multiplication<FieldT>(f_poly_component[i], f_points_coeff_i, L[i]);
    }
    std::fill(f_poly.begin(), f_poly.end(), FieldT(0));
    // f_poly[i] = \sum_i (f_points[i] * L[i])
    for (size_t i = 0; i < nconstraints; ++i) {
      // f_poly[i] = f_poly[i] + f_poly_component[i];
      libfqfft::_polynomial_addition<FieldT>(f_poly, f_poly, f_poly_component[i]);
    }      

  }

  //
  // Evaluate a polynomial F at the encrypted secret input
  // \secret^i*G_1 ie. compute f(\secret)*G1 = [f(\secret)]_i
  //
  // INPUT
  //
  // secret_powers_g1: \secret^i*G1: 0\le{i}<max_degree(Q[j]): 0\le{j}<n
  // Q_polys[0..n-1]: a set of n polynomials
  //
  // OUTPUT
  //
  // [Q_polys[i](\secret)]_1, 0\le 1<n : the "encrypted" evaluation of
  // the polynomials Q_polys[i] in the secret parameter \secret (the
  // toxic waste) multiplied by the group generator G_1 i.e. compute
  // Q_polys[i](\secret)*G_1
  //
  template<typename ppT>
  libff::G1<ppT> plonk_evaluate_poly_at_secret_G1(
						  std::vector<libff::G1<ppT>> secret_powers_g1,
						  polynomial<libff::Fr<ppT>> f_poly
						  )
  {
    assert(f_poly.size() <= secret_powers_g1.size());
    const size_t chunks = 1;
    const size_t num_coefficients = f_poly.size();
    libff::G1<ppT> f_poly_at_secret_g1 = 
      libff::multi_exp<
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
					 std::vector<libff::G1<ppT>>& Q_polys_at_secret_g1
					 )
  {
    assert(secret_powers_g1.size());
    assert(Q_polys.size());
    assert(Q_polys.size() <= secret_powers_g1.size());
    const size_t npolys = Q_polys.size();
    for (size_t i = 0; i < npolys; ++i) {
      Q_polys_at_secret_g1[i] = plonk_evaluate_poly_at_secret_G1<ppT>(secret_powers_g1, Q_polys[i]);
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
					  std::vector<FieldT>& A // accumulatro vector
					  )
  {
    assert(n);
    assert((i >= 0) && (i < n));
    assert(witness.size() == (3*n));
    assert(H_gen.size() == (3*n));
    assert(H_gen_permute.size() == (3*n));
    assert(A.size() == n);    
    FieldT res = FieldT(1);
    if(i > 0) {
      FieldT nom_1 = witness[i-1]         + (beta * H_gen[i-1])             + gamma;
      FieldT den_1 = witness[i-1]         + (beta * H_gen_permute[i-1])     + gamma;

      FieldT nom_2 = witness[n + i-1]     + (beta * H_gen[n+i-1])           + gamma;
      FieldT den_2 = witness[n + i-1]     + (beta * H_gen_permute[n+i-1])   + gamma;

      FieldT nom_3 = witness[2 * n + i-1] + (beta * H_gen[2*n+i-1])         + gamma;
      FieldT den_3 = witness[2 * n + i-1] + (beta * H_gen_permute[2*n+i-1]) + gamma;

      FieldT nom = nom_1 * nom_2 * nom_3;
      FieldT den = den_1 * den_2 * den_3;

      res = nom * den.inverse() * A[i-1];
    }
    return res;
  }

  template<typename FieldT>
  void plonk_compute_accumulator(
				 size_t n, // nconstraints
				 FieldT beta,
				 FieldT gamma,
				 std::vector<FieldT> witness,
				 std::vector<FieldT> H_gen, // H, Hk1, Hk2
				 std::vector<FieldT> H_gen_permute,
				 std::vector<FieldT>& A // accumulatro vector
				 )
  {
    assert(n);
    assert(witness.size() == (3*n));
    assert(H_gen.size() == (3*n));
    assert(H_gen_permute.size() == (3*n));
    assert(A.size() == n);    
    for (size_t i = 0; i < n; ++i) {
      A[i] = plonk_compute_accumulator_factor(i, n, beta, gamma, witness, H_gen, H_gen_permute, A);
    }
  }

  //
  // P(x) plonk_evaluate_poly_at_point(P, x)
  //
  // evaluate a polynomial P given as list of coefficients (p0, p1,
  // ..., p_{n-1}) at a given point t using the lagrange basis L_i:
  // 0 <= i <= deg(P) ie. compute P(t) as:
  //
  // P(t) = p0 L0(t) + p1 L1(t) + ... + p_{n-1} L_{n-1}(t)
  //
  template<typename FieldT>
  FieldT plonk_evaluate_poly_at_point(
				      polynomial<FieldT> P,
				      FieldT t,
				      size_t n // nconstraints = size of domain
				      )
  {
    // make sure that the size of the domain is equal (can it be
    // less?) to the degree of the polyomial (to be able to evaluate
    // the polynomial on this domain)
    assert(P.size() == n);
    std::shared_ptr<libfqfft::evaluation_domain<FieldT>> domain = libfqfft::get_evaluation_domain<FieldT>(n);
    // evaluate all Lagrange polynomials at point x
    std::vector<FieldT> Lt = domain->evaluate_all_lagrange_polynomials(t);
    FieldT res = FieldT(0);    
    // P(t) = p0 L0(t) + p1 L1(t) + ... + p_{n-1} L_{n-1}(t)
    for (size_t i = 0; i < P.size(); ++i) {
      res = res + (P[i] * Lt[i]);
    }
    return res;
  }

  template<typename ppT> void test_plonk()
  {
    // Execute all tests for the given curve.
    ppT::init_public_params();
    
    using Field = libff::Fr<ppT>;
    //    using Group = libff::G1<ppT>;
    //    libff::G1<ppT> G1 = libff::G1<ppT>::G1_one;
    // libff::G2<ppT> G2 = libff::G2<ppT>::G2_one;
    // Setting G2 equal to G1 (Type1 Bilinear group) for simplicity
    // and in order to match the test vectors which are produced under
    // this setting
    //    libff::G1<ppT> G2 = libff::G1<ppT>::G1_one;

    enum W_polys_id{a, b, c};
    enum Q_polys_id{L, R, M, O, C};
    enum t_polys_id{lo, mid, hi};
    enum omega_id{base, base_k1, base_k2}; // roots of unity

    // --- SETUP ---

    printf("[%s:%d] Start setup...\n", __FILE__, __LINE__);
    
    // number of gates / constraints. we have 6 gates for the example
    // circuit + 2 dummy gates to make it a power of 2 (for the fft)
    const int nconstraints = 8;

    // number of q-polynomials
    const int nqpoly = 5;
#ifdef DEBUG
    // ensure nconstraints is power of 2
    bool b_is_power2 = ((nconstraints & (nconstraints - 1)) == 0);
    assert(b_is_power2);
    // ensure that nconstraints is not 0
    assert(nconstraints);
#endif // #ifdef DEBUG 

    // hard-coded gates matrix for the example circuit
    // P(x) = x**3 + x + 5 = 3
    // Each column is a q-vector
    std::vector<std::vector<Field>>gates_matrix
      {
       // q_L     q_R        q_O         q_M       q_C
       {Field(0), Field(0), -Field("1"), Field(1),  Field(0)},   // mul
       {Field(0), Field(0), -Field("1"), Field(1),  Field(0)},   // mul
       {Field(1), Field(1), -Field("1"), Field(0),  Field(0)},   // add
       {Field(0), Field(1),  Field(0),   Field(0), -Field("5")}, // con5
       {Field(0), Field(1),  Field(0),   Field(0),  Field(0)},   // PI
       {Field(1), Field(1), -Field("1"), Field(0),  Field(0)},   // add
       {Field(0), Field(0),  Field(0),   Field(0),  Field(0)},   // dummy
       {Field(0), Field(0),  Field(0),   Field(0),  Field(0)},   // dummy
      };

    // Transposed gates matrix: each row is a q-vector WARN: rows 2
    // q_O and 3 q_M are swapped ti match the Plonk_Py test vectors
    // implementation (reason unclear)
    std::vector<std::vector<Field>>gates_matrix_transpose
      {
       //  mul          mul          add          con5       PI           add        dum        dum
       {   Field(0),    Field(0),    Field(1),    Field(0),  Field(0),    Field(1),  Field(0),  Field(0)}, // q_L
       {   Field(0),    Field(0),    Field(1),    Field(1),  Field(1),    Field(1),  Field(0),  Field(0)}, // q_R
       {   Field(1),    Field(1),    Field(0),    Field(0),  Field(0),    Field(0),  Field(0),  Field(0)}, // q_M
       {-Field("1"), -Field("1"), -Field("1"),    Field(0),  Field(0), -Field("1"),  Field(0),  Field(0)}, // q_O
       {   Field(0),    Field(0),    Field(0), -Field("5"),  Field(0),    Field(0),  Field(0),  Field(0)}, // q_C
      };

    // witness values
    // w_L = a = [ 3,  9, 27,  1,  1, 30,  0,  0]
    // w_R = b = [ 3,  3,  3,  5, 35,  5,  0,  0]
    // w_O = c = [ 9, 27, 30,  5, 35, 35,  0,  0]
    // W = w_L + w_R + w_O
    std::vector<Field> witness
      {
       3,  9, 27,  1,  1, 30,  0,  0, // w_L 
       3,  3,  3,  5, 35,  5,  0,  0, // w_R
       9, 27, 30,  5, 35, 35,  0,  0  // w_O
      };
    
    // output from plonk_compute_permutation()
    std::vector<int> wire_permutation{9, 17, 18, 5, 4, 19, 7, 8, 10, 11, 1, 14, 21, 20, 15, 16, 2, 3, 6, 12, 22, 13, 23, 24};

    // Get the n-th root of unity omega in Fq (n=8 is the number of
    // constraints in the example). omega is a generator of the
    // multiplicative subgroup H.  Example (2**32)-th primitive root
    // of unity in the base field Fq of bls12-381 i.e. such that
    // omega_base**(2**32) = 1. The bls12-381 prime q is such that any
    // power of 2 divides (q-1). In particular 2**32|(q-1) and so the
    // 2**32-th root of unity exists.
    Field omega_base = libff::get_root_of_unity<Field>(nconstraints);
    omega_base.print();
    assert(omega_base == Field("23674694431658770659612952115660802947967373701506253797663184111817857449850"));
    
#ifdef DEBUG
    // assert that omega_base is a 2^32-th root of unity in Fq
    Field temp = libff::power(omega_base, libff::bigint<1>(std::pow(2,32)));
    assert(temp == 1);
#endif // #ifdef DEBUG

    // number of generators for H, Hk1, Hk2
    int ngen = 3;
    
    // output from plonk_compute_roots_of_unity omega[0] are the n
    // roots of unity, omega[1] are omega[0]*k1, omega[2] are
    // omega[0]*k2
    std::vector<std::vector<Field>> omega(ngen);    
    for (int i = 0; i < nconstraints; ++i) {
      Field omega_i = libff::power(omega_base, libff::bigint<1>(i));
      omega[base].push_back(omega_i);
    }

#ifdef DEBUG
    for (int i = 0; i < nconstraints; ++i) {
      printf("w^%d: ", i);
      omega[base][i].print();
    }
    // check that omega^8 = 1 i.e. omega is a generator of the
    // multiplicative subgroup H of Fq of order 'nconstraints'
    Field omega_temp = libff::power(omega_base, libff::bigint<1>(nconstraints));
    printf("w^%d: ", nconstraints);
    omega_temp.print();
    assert(omega_temp == 1);
#endif // #ifdef DEBUG
    
    // We represent the constraints q_L, q_R, q_O, q_M, q_C and the
    // witness w_L, w_R, w_O as polynomials in the roots of unity
    // e.g. f_{q_L}(omega_i) = q_L[i], 0\le{i}<8
    
    // output from plonk_compute_lagrange_basis
    std::vector<polynomial<Field>> L_basis(nconstraints, polynomial<Field>(nconstraints));
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain = libfqfft::get_evaluation_domain<Field>(nconstraints);
    plonk_compute_lagrange_basis<Field>(nconstraints, L_basis);

    // public input (PI)
    Field PI_value = Field(35);
    // index of the row of the PI in the non-transposed gates_matrix 
    int PI_index = 4;
    assert(nconstraints == omega[base].size());
    // compute public input (PI) polynomial
    std::vector<Field> PI_points(nconstraints, Field(0));
    PI_points[PI_index] = Field(-PI_value);
    polynomial<Field> PI_poly;
    plonk_interpolate_over_lagrange_basis<Field>(L_basis, PI_points, PI_poly);
#ifdef DEBUG
    printf("[%s:%d] PI_poly\n", __FILE__, __LINE__);
    print_vector(PI_poly);
#endif // #ifdef DEBUG

#if 1 // DEBUG
    // test Lagrange polynomials
    //template<typename FieldT>
    //FieldT evaluate_lagrange_polynomial(const size_t &m, const std::vector<FieldT> &domain, const FieldT &t, const size_t &idx);
#endif // #if 1 // DEBUG
    
    // output from plonk_compute_constraints_polynomials() compute the q-polynomials from the
    // (transposed) gates matrix over the Lagrange basis q_poly =
    // \sum_i q[i] * L[i] where q[i] is a coefficient (a Field
    // element) and L[i] is a polynomial with Field coeffs
    std::vector<polynomial<Field>> Q_polys(nqpoly, polynomial<Field>(nconstraints));
    for (size_t i = 0; i < nqpoly; ++i) {
      std::vector<Field> q_vec = gates_matrix_transpose[i];
      plonk_interpolate_over_lagrange_basis<Field>(L_basis, q_vec, Q_polys[i]);
    }

#if 1 // DEBUG
    for (int i = 0; i < nqpoly; ++i) {
      printf("\n[%s:%d] Q_polys[%2d]\n", __FILE__, __LINE__, i);
      print_vector(Q_polys[i]);
    }
#endif // #if 1 // DEBUG
    
    // Generate domains on which to evaluate the witness
    // polynomials. k can be random, but we fix it for debug to match
    // against the test vector values
    Field k1 = Field("7069874114745813936829552608791213902061117400356596714713673571023200548519"); // k1
    // similarly, k2 can be random, but we fix it to match the test
    // vectors
    Field k2 = libff::power(k1, libff::bigint<1>(2)); // k2 = k1^2
#if 1 // DEBUG
    printf("[%s:%d] k1 ", __FILE__, __LINE__);
    k1.print();    
    printf("[%s:%d] k2 ", __FILE__, __LINE__);
    k2.print();    
#endif // #if 1 // DEBUG

    // k1 H is a coset of H with generator omega_k1 distinct from H k2
    // H is a coset of H with generator omega_k2, distinct from H and
    // k1 H
    for (int i = 0; i < nconstraints; ++i) {
      Field omega_k1_i = omega[base][i] * k1;
      Field omega_k2_i = omega[base][i] * k2;
      omega[base_k1].push_back(omega_k1_i);
      omega[base_k2].push_back(omega_k2_i);
    }
    
    // H_gen contains the generators of H, k1 H and K2 H in one place
    // ie. omega, omega_k1 and omega_k2
    std::vector<Field> H_gen;
    std::copy(omega[base].begin(), omega[base].end(), back_inserter(H_gen));
    std::copy(omega[base_k1].begin(), omega[base_k1].end(), back_inserter(H_gen));
    std::copy(omega[base_k2].begin(), omega[base_k2].end(), back_inserter(H_gen));
    printf("[%s:%d] H_gen\n", __FILE__, __LINE__);
    print_vector(H_gen);

    // permute H_gen according to the wire permutation
    std::vector<Field> H_gen_permute(ngen*nconstraints, Field(0));
    for (int i = 0; i < ngen*nconstraints; ++i) {
      printf("[%s:%d] i %2d -> %2d, \n", __FILE__, __LINE__, i, wire_permutation[i]-1);
      H_gen_permute[i] = H_gen[wire_permutation[i]-1];
    }
    printf("[%s:%d] H_gen_permute\n", __FILE__, __LINE__);
    print_vector(H_gen_permute);

    std::vector<polynomial<Field>> S_polys(ngen, polynomial<Field>(nconstraints));
    for (int i = 0; i < ngen; ++i) {
      typename std::vector<Field>::iterator begin = H_gen_permute.begin()+(i*nconstraints);
      typename std::vector<Field>::iterator end = H_gen_permute.begin()+(i*nconstraints)+(nconstraints);
      std::vector<Field> S_points(begin, end);
      plonk_interpolate_over_lagrange_basis<Field>(L_basis, S_points, S_polys[i]);
    }

#if 1 // DEBUG
    for (int i = 0; i < ngen; ++i) {
      printf("[%s:%d] S_polys[%d]\n", __FILE__, __LINE__, i);
      print_vector(S_polys[i]);
    }
#endif // #if 1 // DEBUG

    // random hidden element secret (toxic waste). we fix it to a
    // constant in order to match against the test vectors
    Field secret = Field("13778279493383315901513166932749987230291710199728570152123261818328463629146");
#if 1 // DEBUG
    printf("[%s:%d] secret ", __FILE__, __LINE__);
    secret.print();
#endif // #if 1 // DEBUG
    
    // compute powers of secret * G1: 1*G1, secret*G1, secret^2*G1, ...
    const libff::bigint<Field::num_limbs> secret_bigint = secret.as_bigint();
    const size_t window_size = std::max(
					libff::wnaf_opt_window_size<libff::G1<ppT>>(secret_bigint.num_bits()),
					1ul);
    const std::vector<long> naf =
      libff::find_wnaf<Field::num_limbs>(window_size, secret_bigint);    
    std::vector<libff::G1<ppT>> secret_powers_g1;
    secret_powers_g1.reserve(nconstraints + 3);
    libff::G1<ppT> secret_i_g1 = libff::G1<ppT>::one();
    secret_powers_g1.push_back(secret_i_g1);
    for (size_t i = 1; i < (nconstraints + 3); ++i) {
      // secret^i * G1
      secret_i_g1 = libff::fixed_window_wnaf_exp<libff::G1<ppT>>(
								window_size, secret_i_g1, naf);
      secret_powers_g1.push_back(secret_i_g1);
    }    

#if 1 // DEBUG
    for (int i = 0; i < nconstraints + 3; ++i) {
      printf("secret_power[%2d] ", i);
      secret_powers_g1[i].print();
    }
#endif // #if 1 // DEBUG

    //    std::vector<polynomial<Field>> Q_polys(nqpoly, polynomial<Field>(nconstraints));
    // Q_polys[i] is polynomial<Field>

    //    const size_t num_coefficients = phi.size();
    //    const libff::Fr<ppT> i
    //    const Field phi_i = libfqfft::evaluate_polynomial(num_coefficients, phi, i);

    // Verifier precomputation:
    std::vector<libff::G1<ppT>> Q_polys_at_secret_g1(Q_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(secret_powers_g1, Q_polys, Q_polys_at_secret_g1);
    std::vector<libff::G1<ppT>> S_polys_at_secret_g1(S_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(secret_powers_g1, S_polys, S_polys_at_secret_g1);
    
#if 1 // DEBUG
    for (int i = 0; i < (int)Q_polys.size(); ++i) {
      printf("Q_polys_at_secret_G1[%d] \n", i);
      Q_polys_at_secret_g1[i].print();
    }
    for (int i = 0; i < (int)S_polys.size(); ++i) {
      printf("S_polys_at_secret_G1[%d] \n", i);
      S_polys_at_secret_g1[i].print();
    }
#endif // #if 1 // DEBUG

    // --- PROVER ---

    printf("[%s:%d] Prover preparation...\n", __FILE__, __LINE__);
    
    int nwitness = 3;
    std::vector<polynomial<Field>> W_polys(nwitness, polynomial<Field>(nconstraints));
    for (int i = 0; i < nwitness; ++i) {
      typename std::vector<Field>::iterator begin = witness.begin()+(i*nconstraints);
      typename std::vector<Field>::iterator end = witness.begin()+(i*nconstraints)+(nconstraints);
      std::vector<Field> W_points(begin, end);
      plonk_interpolate_over_lagrange_basis<Field>(L_basis, W_points, W_polys[i]);
    }

#if 1 // DEBUG
    for (int i = 0; i < nwitness; ++i) {
      printf("[%s:%d] W_polys[%d]\n", __FILE__, __LINE__, i);
      print_vector(W_polys[i]);
    }
#endif // #if 1 // DEBUG


    // Evaluate polynomials over 2*n points
    Field omega2_base = libff::get_root_of_unity<Field>(2*nconstraints);
    omega2_base.print();
    
    std::vector<Field> omega2;    
    for (int i = 0; i < (2*nconstraints); ++i) {
      Field omega2_i = libff::power(omega2_base, libff::bigint<1>(i));
      omega2.push_back(omega2_i);
    }
    
#ifdef NDEBUG
    for (int i = 0; i < (2*nconstraints); ++i) {
      printf("w^%d: ", i);
      omega2[i].print();
    }
    // check that omega2^8 = 1 i.e. omega2 is a generator of the
    // multiplicative subgroup H of Fq of order 'nconstraints'
    Field omega2_temp = libff::power(omega2_base, libff::bigint<1>(2*nconstraints));
    printf("w^%2d: ", 2*nconstraints);
    omega2_temp.print();
    assert(omega2_temp == 1);
#endif // #ifdef DEBUG

    // Evaluate polynomials over 8*n points
    Field omega3_base = libff::get_root_of_unity<Field>(8*nconstraints);
    omega3_base.print();
    
    std::vector<Field> omega3;    
    for (int i = 0; i < (8*nconstraints); ++i) {
      Field omega3_i = libff::power(omega3_base, libff::bigint<1>(i));
      omega3.push_back(omega3_i);
    }
    
#ifdef NDEBUG
    for (int i = 0; i < (8*nconstraints); ++i) {
      printf("w^%2d: ", i);
      omega3[i].print();
    }
    // check that omega3^8 = 1 i.e. omega3 is a generator of the
    // multiplicative subgroup H of Fq of order 'nconstraints'
    Field omega3_temp = libff::power(omega3_base, libff::bigint<1>(8*nconstraints));
    printf("w^%d: ", 8*nconstraints);
    omega3_temp.print();
    assert(omega3_temp == 1);
#endif // #ifdef DEBUG

    //    Field t = Field(omega_base);
    //    Field a;
    //    a = domain->compute_vanishing_polynomial(t);
    //    printf("[%s:%d] a ", __FILE__, __LINE__);
    //    a.print();

    printf("[%s:%d] Prover Round 1...\n", __FILE__, __LINE__);

    // vanishing polynomial zh_poly(X) = x^n-1. vanishes on all n roots of
    // unity omega
    std::vector<Field> zh_poly(nconstraints+1, Field(0));
    zh_poly[0] = Field(-1);
    zh_poly[nconstraints] = Field(1);
    printf("[%s:%d] Vanishing polynomial\n", __FILE__, __LINE__);
    print_vector(zh_poly);

    // blinding scalars b1, b2, ..., b9. random but fixed to match the
    // python test vectors
    //    size_t nscalars = 9;
    std::vector<Field> rand_scalars
      {
       Field("8063396892870388055806370369789704857755116044327394765020751373651916505604"),
       Field("6827026430120056597679453111370682306316948363643792417785314991392447377909"),
       Field("20903799562102387073359556962112335230020378309596765284572286455779329747315"),
       Field("27445824854335787523979734401573136947589999159092723101543900479804718923773"),
       Field("5216447975508541021290757380485235885597475407449078898274648785082871817687"),
       Field("48720156268681305476740160454555587712391923971264646500554531948708930944069"),
       Field("2131891516651518828698089707163191982683101187444208330829153527689737950718"),
       Field("47548532878000795436471885496554996210469829388180983864669623532585348412472"),
       Field("2534599719872500160900817853393315885420633320379105447254598708515031311667")
      };
#ifdef DEBUG
    printf("[%s:%d] rand_scalars\n", __FILE__, __LINE__);
    print_vector(rand_scalars);
#endif // #ifdef DEBUG

    std::vector<std::vector<Field>> blind_polys
      {
       {rand_scalars[1], rand_scalars[0]}, // b1 + b0 X
       {rand_scalars[3], rand_scalars[2]}, // b3 + b2 X
       {rand_scalars[5], rand_scalars[4]}  // b5 + b4 X
      };
    
    // a_poly = blind_polys[0] * zh_poly + W_polys[0]
    std::vector<std::vector<Field>> W_polys_blinded(nwitness);
    for (int i = 0; i < nwitness; ++i) {
      libfqfft::_polynomial_multiplication<Field>(W_polys_blinded[i], blind_polys[i], zh_poly);
      libfqfft::_polynomial_addition<Field>(W_polys_blinded[i], W_polys_blinded[i], W_polys[i]);
    }
#ifdef DEBUG
    for (int i = 0; i < nwitness; ++i) {
      printf("[%s:%d] W_polys_blinded[%d]\n", __FILE__, __LINE__, i);
      print_vector(W_polys_blinded[i]);
    }
#endif // #ifdef DEBUG

    std::vector<libff::G1<ppT>> W_polys_blinded_at_secret_g1(W_polys_blinded.size());
    printf("[%s:%d] poly %d %d %d secret %d\n", __FILE__, __LINE__,
	   (int)W_polys_blinded[0].size(),
	   (int)W_polys_blinded[1].size(),
	   (int)W_polys_blinded[2].size(),
	   (int)secret_powers_g1.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(secret_powers_g1, W_polys_blinded, W_polys_blinded_at_secret_g1);
#ifdef DEBUG
    for (int i = 0; i < nwitness; ++i) {
      printf("W_polys_at_secret_g1[%d]\n", i);
      W_polys_blinded_at_secret_g1[i].print();
    }
#endif // #ifdef DEBUG

    printf("[%s:%d] Prover Round 2...\n", __FILE__, __LINE__);

    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field beta = Field("3710899868510394644410941212967766116886736137326022751891187938298987182388");
    Field gamma = Field("11037930384083194587907709665332116843267274045828802249545114995763715746939");

    // compute permutation polynomial

    // blinding polynomial
    std::vector<Field> z1_blind_poly{rand_scalars[8], rand_scalars[7], rand_scalars[6]}; // b8 + b7 X + b6 X^2
    // multiply by the vanishing polynomial: z1 = z1 * zh_poly
    libfqfft::_polynomial_multiplication<Field>(z1_blind_poly, z1_blind_poly, zh_poly);
#ifdef DEBUG
    printf("[%s:%d] z1_blind_poly * zh_poly\n", __FILE__, __LINE__);
    print_vector(z1_blind_poly);
#endif // #ifdef DEBUG

    // A[0] = 1; ... A[i] = computed from (i-1)
    std::vector<Field> A_vector(nconstraints, Field(0));
    plonk_compute_accumulator(nconstraints, beta, gamma, witness, H_gen, H_gen_permute, A_vector);
#ifdef DEBUG
    for (int i = 0; i < nconstraints; ++i) {
      printf("A[%d] ", i);
      A_vector[i].print();
    }
#endif // #ifdef DEBUG

    polynomial<Field> A_poly(nconstraints);
    plonk_interpolate_over_lagrange_basis<Field>(L_basis, A_vector, A_poly);
#ifdef DEBUG
    printf("[%s:%d] A_poly\n", __FILE__, __LINE__);
    print_vector(A_poly);
#endif // #ifdef DEBUG

    // add blinding polynomial z_1 to the accumulatir polynomial A_poly
    polynomial<Field> z_poly;
    libfqfft::_polynomial_addition<Field>(z_poly, z1_blind_poly, A_poly);
#ifdef DEBUG
    printf("[%s:%d] z_poly\n", __FILE__, __LINE__);
    print_vector(z_poly);
#endif // #ifdef DEBUG
    
    libff::G1<ppT> z_poly_at_secret_g1 = plonk_evaluate_poly_at_secret_G1<ppT>(secret_powers_g1, z_poly);
#ifdef DEBUG
    printf("[%s:%d] z_poly_at_secret_g1\n", __FILE__, __LINE__);
    z_poly_at_secret_g1.print();
#endif // #ifdef DEBUG


    // -- Verified up to here with the Plonk-Py implementation ---
    
    printf("[%s:%d] Prover Round 3...\n", __FILE__, __LINE__);
    // W_polys_blinded = [a_poly, b_poly, c_poly]

    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field alpha = Field("37979978999274723893071781986484838492111162341880360022719385400306128734648");
    
    // plonk_poly_shifted    
    // Computing the polynomial z(x*w) i.e. z(x) shifted by w where
    // w=omega is the base root of unity and z is z_poly. we do this
    // by multiplying the coefficients of z by w
    std::vector<Field> z_poly_xomega(z_poly.size(), Field(0));
    for (size_t i = 0; i < z_poly.size(); ++i) {
      z_poly_xomega[i] = z_poly[i] * omega[base][1];
    }

#ifdef DEBUG
    printf("[%s:%d] z_poly_xomega\n", __FILE__, __LINE__);
    print_vector(z_poly_xomega);
#endif // #ifdef DEBUG

    // start computation of polynomial t(X) in round 3. we break t
    // into 4 parts which we compute separately. each of the 4 parts
    // is multiplied by 1/zh_poly in the paper
    std::vector<polynomial<Field>> t_part(4);

    // --- Computation of t_part[0]
    
    // a(x)b(x)q_M(x)
    polynomial<Field> abqM;
    libfqfft::_polynomial_multiplication<Field>(abqM, W_polys[a], W_polys[b]);
    libfqfft::_polynomial_multiplication<Field>(abqM, abqM, Q_polys[M]);    
    // a(x)q_L(x)
    polynomial<Field> aqL;
    libfqfft::_polynomial_multiplication<Field>(aqL, W_polys[a], Q_polys[L]);
    // b(x)q_R(x)
    polynomial<Field> bqR;
    libfqfft::_polynomial_multiplication<Field>(bqR, W_polys[b], Q_polys[R]);
    // c(x)q_O(x)
    polynomial<Field> cqO;
    libfqfft::_polynomial_multiplication<Field>(cqO, W_polys[c], Q_polys[O]);
    // t_part[0](x) = a(x)b(x)q_M(x) + a(x)q_L(x) + b(x)q_R(x) + c(x)q_O(x) + PI(x) + q_C(x)
    libfqfft::_polynomial_addition<Field>(t_part[0], abqM, aqL);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], bqR);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], cqO);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], PI_poly);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], Q_polys[C]);
    
    // --- Computation of t_part[1]
    
    // X*beta as polynomial in X
    std::vector<polynomial<Field>> xbeta_poly
      {
       {Field(0), beta}, // X*beta
       {Field(0), beta*k1}, // X*beta*k1
       {Field(0), beta*k2} // X*beta*k2
      };    
    // represent gamma as polynomial in X, needed for prover Round 3
    polynomial<Field> gamma_poly{gamma}; // gamma
    // represent alpha as polynomial in X, needed for prover Round 3
    polynomial<Field> alpha_poly{alpha}; // alpha

    // a(x) + beta*x + gamma 
    polynomial<Field> a_xbeta_gamma;
    libfqfft::_polynomial_addition<Field>(a_xbeta_gamma, W_polys[a], xbeta_poly[base]);    
    libfqfft::_polynomial_addition<Field>(a_xbeta_gamma, a_xbeta_gamma, gamma_poly);    
    // b(x) + beta_k1*x + gamma 
    polynomial<Field> b_xbeta_gamma_k1;
    libfqfft::_polynomial_addition<Field>(b_xbeta_gamma_k1, W_polys[b], xbeta_poly[base_k1]);
    libfqfft::_polynomial_addition<Field>(b_xbeta_gamma_k1, b_xbeta_gamma_k1, gamma_poly);
    // c(x) + beta_k1*x + gamma 
    polynomial<Field> c_xbeta_gamma_k2;
    libfqfft::_polynomial_addition<Field>(c_xbeta_gamma_k2, W_polys[c], xbeta_poly[base_k2]);
    libfqfft::_polynomial_addition<Field>(c_xbeta_gamma_k2, c_xbeta_gamma_k2, gamma_poly);
    // t_part[1] = (a(x) + beta*x + gamma)*(b(x) + beta_k1*x +
    // gamma)*(c(x) + beta_k1*x + gamma)*z(x)*alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[1], a_xbeta_gamma, b_xbeta_gamma_k1);
    libfqfft::_polynomial_multiplication<Field>(t_part[1], t_part[1], c_xbeta_gamma_k2);
    libfqfft::_polynomial_multiplication<Field>(t_part[1], t_part[1], z_poly);
    libfqfft::_polynomial_multiplication<Field>(t_part[1], t_part[1], alpha_poly);

    // --- Computation of t_part[2]
    
    // represent beta as polynomial in X, needed for prover Round 3
    polynomial<Field> beta_poly{beta};
    // S*beta as polynomial
    // S_sigma1(x)*beta, S_sigma2(x)*beta, S_sigma3(x)*beta
    std::vector<polynomial<Field>> sbeta_poly(ngen);
    for (int i = 0; i < ngen; ++i) {
      libfqfft::_polynomial_multiplication<Field>(sbeta_poly[i], S_polys[i], beta_poly);
    }
    // a(x) + S_sigma1(x)*beta + gamma
    polynomial<Field> a_sbeta_gamma;
    libfqfft::_polynomial_addition<Field>(a_sbeta_gamma, W_polys[a], sbeta_poly[base]);    
    libfqfft::_polynomial_addition<Field>(a_sbeta_gamma, a_sbeta_gamma, gamma_poly);    
    // b(x) + S_sigma2(x)*beta + gamma
    polynomial<Field> b_sbeta_gamma_k1;
    libfqfft::_polynomial_addition<Field>(b_sbeta_gamma_k1, W_polys[b], sbeta_poly[base_k1]);    
    libfqfft::_polynomial_addition<Field>(b_sbeta_gamma_k1, b_sbeta_gamma_k1, gamma_poly);    
    // b(x) + S_sigma2(x)*beta + gamma
    polynomial<Field> c_sbeta_gamma_k2;
    libfqfft::_polynomial_addition<Field>(c_sbeta_gamma_k2, W_polys[c], sbeta_poly[base_k2]);    
    libfqfft::_polynomial_addition<Field>(c_sbeta_gamma_k2, c_sbeta_gamma_k2, gamma_poly);    
    // t_part[2] = (a(x) + S_sigma1(x)*beta + gamma)*(b(x) +
    // S_sigma2(x)*beta + gamma)*(b(x) + S_sigma2(x)*beta +
    // gamma)*z(x*omega)*alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[2], a_sbeta_gamma, b_sbeta_gamma_k1);
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], c_sbeta_gamma_k2);
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], z_poly_xomega);
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], alpha_poly);
    // -t_part[2]
    polynomial<Field> neg_one_poly = {-Field("1")};
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], neg_one_poly);
    
    // --- Computation of t_part[3]

    // z(x) - 1
    polynomial<Field> z_neg_one;
    libfqfft::_polynomial_addition<Field>(z_neg_one, z_poly, neg_one_poly);
    // (z(x)-1) * L_1(x)
    libfqfft::_polynomial_multiplication<Field>(t_part[3], z_neg_one, L_basis[0]);
    // (z(x)-1) * L_1(x) * alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[3], t_part[3], alpha_poly);
    // (z(x)-1) * L_1(x) * alpha * alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[3], t_part[3], alpha_poly);


    // --- computation of t(x)

    // t(x) = (t[0] + t[1] + (-t[2]) + t[3]) / zh(x)
    polynomial<Field> t_poly{Field(0)};
    libfqfft::_polynomial_addition<Field>(t_poly, t_poly, t_part[1]);
    //    libfqfft::_polynomial_addition<Field>(t_poly, t_poly, t_part[2]);
    //    libfqfft::_polynomial_addition<Field>(t_poly, t_poly, t_part[2]);
    //    libfqfft::_polynomial_addition<Field>(t_poly, t_poly, t_part[3]);
#ifdef DEBUG
    printf("[%s:%d] t_poly\n", __FILE__, __LINE__);
    print_vector(t_poly);
#endif // #ifdef DEBUG

    //    t(x) = t(x) / zh(x): A/B = (Q, R) st. A = (Q * B) + R.
    polynomial<Field> remainder;
    libfqfft::_polynomial_division(t_poly, remainder, t_poly, zh_poly);
#ifdef DEBUG
    printf("[%s:%d] remainder\n", __FILE__, __LINE__);
    print_vector(remainder);
#endif // #ifdef DEBUG
    //    assert(libfqfft::_is_zero(remainder));

    /**
     * Perform the standard Euclidean Division algorithm.
     * Input: Polynomial A, Polynomial B, where A / B
     * Output: Polynomial Q, Polynomial R, such that A = (Q * B) + R.
     */
    //    template<typename FieldT>
    //      void _polynomial_division(std::vector<FieldT> &q, std::vector<FieldT> &r, const std::vector<FieldT> &a, const std::vector<FieldT> &b);
   
    // end 

    printf("[%s:%d] Test OK\n", __FILE__, __LINE__);
  }

  // output from plonk_compute_permutation_polynomials() compute the
  // S polynomials
  
  TEST(TestPlonk, BLS12_381)
  {
    test_plonk<libff::bls12_381_pp>();
  }

} // namespace libsnark
