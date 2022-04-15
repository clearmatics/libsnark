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
//#include <libsnark/polynomial_Commitments/kzg10.hpp>
//#include <libsnark/polynomial_Commitments/tests/polynomial_Commitment_test_utils.hpp>
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
- permutation polynomials S_sigma_1, S_sigma_2, S_sigma_3
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
  // Compute the Lagrange basis polynomails for interpolating sets of n points
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
      libfqfft::_polynomial_multiplication(f_poly_component[i], f_points_coeff_i, L[i]);
    }
    std::fill(f_poly.begin(), f_poly.end(), FieldT(0));
    // f_poly[i] = \sum_i (f_points[i] * L[i])
    for (size_t i = 0; i < nconstraints; ++i) {
      // f_poly[i] = f_poly[i] + f_poly_component[i];
      libfqfft::_polynomial_addition<FieldT>(f_poly, f_poly, f_poly_component[i]);
    }      

  }


  template<typename ppT> void test_plonk()
  {
    // Execute all tests for the given curve.
    ppT::init_public_params();
    
    using Field = libff::Fr<ppT>;
    using Group = libff::G1<ppT>;
    libff::G1<ppT> G1 = libff::G1<ppT>::G1_one;
    // libff::G2<ppT> G2 = libff::G2<ppT>::G2_one;
    // Setting G2 equal to G1 (Type1 Bilinear group) for simplicity
    // and in order to match the test vectors which are produced under
    // this setting
    libff::G1<ppT> G2 = libff::G1<ppT>::G1_one;

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

    // output from plonk_compute_permutation()
    std::vector<int> wire_permutation{9, 17, 18, 5, 4, 19, 7, 8, 10, 11, 1, 14, 21, 20, 15, 16, 2, 3, 6, 12, 22, 13, 23, 24};

    // public input (PI)
    Field public_input = Field(35);

    // index of the row of the PI in the non-transposed gates_matrix 
    int public_input_index = 4;

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

    // output from plonk_compute_roots_of_unity
    std::vector<Field> omega;    
    for (int i = 0; i < nconstraints; ++i) {
      Field omega_i = libff::power(omega_base, libff::bigint<1>(i));
      omega.push_back(omega_i);
    }

#ifdef DEBUG
    for (int i = 0; i < nconstraints; ++i) {
      printf("w^%d: ", i);
      omega[i].print();
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
    std::vector<polynomial<Field>> L(nconstraints, polynomial<Field>(nconstraints));
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain = libfqfft::get_evaluation_domain<Field>(nconstraints);
    plonk_compute_lagrange_basis<Field>(nconstraints, L);

#if 1 // DEBUG
    // test Lagrange polynomials
    //template<typename FieldT>
    //FieldT evaluate_lagrange_polynomial(const size_t &m, const std::vector<FieldT> &domain, const FieldT &t, const size_t &idx);
#endif // #if 1 // DEBUG
    
    // output from plonk_compute_constraints_polynomials() compute the q-polynomials from the
    // (transposed) gates matrix over the Lagrange basis q_poly =
    // \sum_i q[i] * L[i] where q[i] is a coefficient (a Field
    // element) and L[i] is a polynomial with Field coeffs
    std::vector<polynomial<Field>> Q(nqpoly, polynomial<Field>(nconstraints));
    for (size_t i = 0; i < nqpoly; ++i) {
      std::vector<Field> q_vec = gates_matrix_transpose[i];
      plonk_interpolate_over_lagrange_basis<Field>(L, q_vec, Q[i]);
    }

#if 1 // DEBUG
    for (int i = 0; i < nqpoly; ++i) {
      printf("\n[%s:%d] Q[%2d]\n", __FILE__, __LINE__, i);
      print_vector(Q[i]);
    }
#endif // #if 1 // DEBUG
    
    // Generate domains on which to evaluate the witness polynomials k
    // can be random, but we fix it for debug to match against the
    // test vector values
    Field k = Field("7069874114745813936829552608791213902061117400356596714713673571023200548519");
#if 1 // DEBUG
    printf("[%s:%d] k ", __FILE__, __LINE__);
    k.print();
#endif // #if 1 // DEBUG

    // k1 H is a coset of H with generator omega_k1 distinct from H
    std::vector<Field> omega_k1;    
    for (int i = 0; i < nconstraints; ++i) {
      Field omega_k1_i = omega[i] * k;
      omega_k1.push_back(omega_k1_i);
    }
    
    // k2 H is a coset of H with generator omega_k2, distinct from H
    // and k1 H
    std::vector<Field> omega_k2;    
    for (int i = 0; i < nconstraints; ++i) {
      Field omega_k2_i = omega[i] * libff::power(k, libff::bigint<1>(2));
      omega_k2.push_back(omega_k2_i);
    }

    // sigma contains the generators of H, k1 H and K2 H in one place
    // ie. omega, omega_k1 and omega_k2
    std::vector<Field> sigma;
    std::copy(omega.begin(), omega.end(), back_inserter(sigma));
    std::copy(omega_k1.begin(), omega_k1.end(), back_inserter(sigma));
    std::copy(omega_k2.begin(), omega_k2.end(), back_inserter(sigma));
    printf("[%s:%d] sigma\n", __FILE__, __LINE__);
    print_vector(sigma);

    // number of sigma polynomials S_sigma
    int nsigma = 3;
    // permute sigma according to the wire permutation
    std::vector<Field> sigma_star(nsigma*nconstraints, Field(0));
    for (int i = 0; i < nsigma*nconstraints; ++i) {
      printf("[%s:%d] i %2d -> %2d, \n", __FILE__, __LINE__, i, wire_permutation[i]-1);
      sigma_star[i] = sigma[wire_permutation[i]-1];
    }
    printf("[%s:%d] sigma_star\n", __FILE__, __LINE__);
    print_vector(sigma_star);

    std::vector<polynomial<Field>> S_sigma_poly(nsigma, polynomial<Field>(nconstraints));
    for (int i = 0; i < nsigma; ++i) {
      typename std::vector<Field>::iterator begin = sigma_star.begin()+(i*nconstraints);
      typename std::vector<Field>::iterator end = sigma_star.begin()+(i*nconstraints)+(nconstraints);
      std::vector<Field> S_sigma_points(begin, end);
      plonk_interpolate_over_lagrange_basis<Field>(L, S_sigma_points, S_sigma_poly[i]);
    }

#if 1 // DEBUG
    for (int i = 0; i < nsigma; ++i) {
      printf("[%s:%d] S_sigma_poly[%d]\n", __FILE__, __LINE__, i);
      print_vector(S_sigma_poly[i]);
    }
#endif // #if 1 // DEBUG

    // random hidden element alpha (toxic waste). we fix it to a
    // constant in order to match against the test vectors
    Field alpha = Field("13778279493383315901513166932749987230291710199728570152123261818328463629146");
#if 1 // DEBUG
    printf("[%s:%d] alpha ", __FILE__, __LINE__);
    alpha.print();
#endif // #if 1 // DEBUG
    
    printf("[%s:%d] G1 \n", __FILE__, __LINE__);
    G1.print();
    G1.print_coordinates();
    
    //    printf("[%s:%d] G2 \n", __FILE__, __LINE__);
    //    libff::G2<ppT>::G2_one.print();
    //    libff::G2<ppT>::G2_one.print_coordinates();

#if 1
    //    libff::G1<ppT> x = Field("2") * G1;
    libff::G1<ppT> x = alpha * G1;
    printf("[%s:%d] x * G1\n", __FILE__, __LINE__);
    x.print();
    x.print_coordinates();
#endif    

#if 0    
    libff::G1<ppT> y = alpha * libff::G1<ppT>::one();
    printf("[%s:%d]     alpha * G1 affine\n", __FILE__, __LINE__);
    y.print(); // print affine coordinates
    printf("[%s:%d]     alpha * G1 projec\n", __FILE__, __LINE__);
    y.print_coordinates(); // print projective coordinates
    //    std::ostringstream ss;
    libff::group_write<libff::encoding_json, libff::form_plain, libff::compression_off>(y, std::cout);
    //    printf("[%s:%d] %s\n", __FILE__, __LINE__, ss.str());
    //    libff::G1<ppT> z = y.to_special();
    // perform scalar multiplication: alpha * G1 where alpha \in Fr, G1=(x,y): x,y \in Fq
#endif
    
    // compute powers of alpha: 1, alpha, alpha^2, ...
    const libff::bigint<Field::num_limbs> alpha_bigint = alpha.as_bigint();
    const size_t window_size = std::max(
					libff::wnaf_opt_window_size<libff::G1<ppT>>(alpha_bigint.num_bits()),
					1ul);
    const std::vector<long> naf =
      libff::find_wnaf<Field::num_limbs>(window_size, alpha_bigint);
      
    std::vector<libff::G1<ppT>> alpha_powers_g1;
    alpha_powers_g1.reserve(nconstraints + 3);
    libff::G1<ppT> alpha_i_g1 = libff::G1<ppT>::one();
    alpha_powers_g1.push_back(alpha_i_g1);
    for (size_t i = 1; i < (nconstraints + 3); ++i) {
      // alpha^i * G1
      alpha_i_g1 = libff::fixed_window_wnaf_exp<libff::G1<ppT>>(
								window_size, alpha_i_g1, naf);
      alpha_powers_g1.push_back(alpha_i_g1);
    }
    
    //      alpha_i_g1 = libff::fixed_window_wnaf_exp<libff::G1<ppT>>(window_size, alpha_i_g1, naf);
    /*
    std::vector<libff::G1<ppT>> alpha_powers_g1;
    alpha_powers_g1.reserve(nconstraints + 3);
    libff::G1<ppT> alpha_i_g1 = libff::G1<ppT>::one();
    alpha_powers_g1.push_back(alpha_i_g1);
    for (size_t i = 1; i < nconstraints + 3; ++i) {
      Field omega_k2_i = omega[i] * libff::power(k, libff::bigint<1>(2));
      
        alpha_i_g1 =
	  libff::fixed_window_wnaf_exp<libff::G1<ppT>>(
            window_size, alpha_i_g1, naf);
        alpha_powers_g1.push_back(alpha_i_g1);
    }
    */
#if 0    
    std::vector<libff::G1<ppT>> alpha_powers_g1;
    alpha_powers_g1.reserve(nconstraints + 3);
    for (int i = 0; i < nconstraints + 3; ++i) {
      Field alpha_i = libff::power(alpha, libff::bigint<1>(i));
      libff::G1<ppT> alpha_powers_i = alpha_i * G1;
      alpha_powers_g1.push_back(alpha_powers_i);
    }
#endif
#if 1 // DEBUG
    for (int i = 0; i < nconstraints + 3; ++i) {
      printf("alpha_power[%2d] ", i);
      alpha_powers_g1[i].print();
    }
#endif // #if 1 // DEBUG

    
    printf("[%s:%d] Test OK\n", __FILE__, __LINE__);
  }

  // output from plonk_compute_permutation_polynomials() compute the
  // S_sigma polynomials
  
  TEST(TestPlonk, BLS12_381)
  {
    test_plonk<libff::bls12_381_pp>();
  }

} // namespace libsnark
