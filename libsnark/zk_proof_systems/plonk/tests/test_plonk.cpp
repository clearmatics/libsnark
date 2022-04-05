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

#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>
//#include <libfqfft/polynomial_arithmetic/basic_operations.hpp>
//#include <libfqfft/evaluation_domain/domains/arithmetic_sequence_domain.hpp>
#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>
//#include <libfqfft/evaluation_domain/domains/extended_radix2_domain.hpp>
//#include <libfqfft/evaluation_domain/domains/geometric_sequence_domain.hpp>
//#include <libfqfft/evaluation_domain/domains/step_radix2_domain.hpp>
//#include <libfqfft/polynomial_arithmetic/naive_evaluate.hpp>
//#include <libfqfft/evaluation_domain/domains/basic_radix2_domain_aux.hpp>

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

  template<typename FieldT>  void print_vector(std::vector<FieldT> v)
  {
    for (size_t i = 0; i < v.size(); ++i) {
      printf("v[%d]: ", (int)i);
      v[i].print();
    }
  }

  template<typename Field> void multiply_polynomial_by_scalar(Field a, polynomial<Field> &v)
  {
    transform(v.begin(), v.end(), v.begin(), [a](Field &c){ return c*a; });
  }
  
  template<typename Field> polynomial<Field> add_polynomials(const polynomial<Field> a, const polynomial<Field> b)
  {
    assert(a.size() == b.size());

    std::vector<Field> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
		   std::back_inserter(result), std::plus<Field>());
    return result;
  }

/*
  Below we make use of pseudocode from [CLRS 2n Ed, pp. 864].
  Also, note that it's the caller's responsibility to multiply by 1/N.
*/
  template<typename FieldT, typename GroupT>
  void _basic_serial_radix2_FFT(std::vector<GroupT> &a, const FieldT &omega, bool b_inverse)
  {
    const size_t n = a.size(), logn = log2(n);
    if (n != (1u << logn)) {
      printf("[%s:%d] Must be power of 2. Terminating...\n", __FILE__, __LINE__);
      assert(0);
    }

    /* swapping in place (from Storer's book) */
    for (size_t k = 0; k < n; ++k)
      {
        const size_t rk = libff::bitreverse(k, logn);
        if (k < rk)
	  std::swap(a[k], a[rk]);
      }

    size_t m = 1; // invariant: m = 2^{s-1}
    for (size_t s = 1; s <= logn; ++s)
      {
        // w_m is 2^s-th root of unity now
        const FieldT w_m = omega^(n/(2*m));

        asm volatile  ("/* pre-inner */");
        for (size_t k = 0; k < n; k += 2*m)
	  {
            FieldT w = FieldT::one();
            for (size_t j = 0; j < m; ++j)
	      {
                const GroupT t = w * a[k+j+m];
                a[k+j+m] = a[k+j] - t;
                a[k+j] = a[k+j] + t;
                w *= w_m;
	      }
	  }
        asm volatile ("/* post-inner */");
        m *= 2;
      }

    if(b_inverse) {
      const FieldT sconst = FieldT(a.size()).inverse();
      for (size_t i = 0; i < a.size(); ++i)
	{
	  a[i] *= sconst;
	}
    }
  }
  
  
  template<typename ppT> void test_plonk()
  {
    // Execute all tests for the given curve.
    ppT::init_public_params();
    
    using Field = libff::Fr<ppT>;

    // number of gates / constraints. we have 6 gates for the example
    // circuit + 2 dummy gates to make it a power of 2 (for the fft)
    const size_t nconstraints = 8;

    // number of q-polynomials
    const size_t nqpoly = 5;
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

    // Transposed gates matrix: each row is a q-vector
    std::vector<std::vector<Field>>gates_matrix_transpose
      {
       {   Field(0),    Field(0),    Field(1),    Field(0),  Field(0),    Field(1),  Field(0),  Field(0)}, // q_L
       {   Field(0),    Field(0),    Field(1),    Field(1),  Field(1),    Field(1),  Field(0),  Field(0)}, // q_R
       {-Field("1"), -Field("1"), -Field("1"),    Field(0),  Field(0), -Field("1"),  Field(0),  Field(0)}, // q_O
       {   Field(1),    Field(1),    Field(0),    Field(0),  Field(0),    Field(0),  Field(0),  Field(0)}, // q_M
       {   Field(0),    Field(0),    Field(0), -Field("5"),  Field(0),    Field(0),  Field(0),  Field(0)}, // q_C
      };

    // output from compute_permutation()
    std::vector<int> wire_permutation{9, 17, 18, 5, 4, 19, 7, 8, 10, 11, 1, 14, 21, 20, 15, 16, 2, 3, 6, 12, 22, 13, 23, 24};

    Field public_input = Field(35);

    // index of the row in the gates_matrix
    int public_input_index = 4;

    // output from compute_omega_base()
    // Example (2**32)-th primitive root of unity in the base field Fq
    // of bls12-381 i.e. such that omega_base**(2**32) = 1. The
    // bls12-381 prime q is such that any power of 2 divides (q-1). In
    // particular 2**32|(q-1) and so the 2**32-th root of unity
    // exists.
    Field omega_base = Field("5398635374615924329006194896463915232623626471712994625313653909674492148212");
#ifdef DEBUG
    // assert that omega_base is a 2^32-th root of unity in Fq
    Field temp = libff::power(omega_base, libff::bigint<1>(std::pow(2,32)));
    assert(temp == 1);
#endif // #ifdef DEBUG

    // Generate the n-th root of unity omega in Fq (n=8 is the number
    // of constraints in the example) using omega_base. omega is a
    // generator of the multiplicative subgroup H.
    const libff::bigint<1> n = std::pow(2,32) / nconstraints;
    n.print();
    assert(n == 536870912); 
    Field omega = libff::power(omega_base, n);
    omega.print();
    assert(omega == Field("8685283084174350996472453922654922162880456818468779543064782192722679779374"));

    // output from compute_roots_of_unity
    std::vector<Field> roots;    
    for (size_t i = 0; i < nconstraints; ++i) {
      Field omega_i = libff::power(omega, libff::bigint<1>(i));
      roots.push_back(omega_i);
    }

#ifdef DEBUG
    for (size_t i = 0; i < nconstraints; ++i) {
      printf("w^%d: ", i);
      roots[i].print();
    }
    // check that omega^8 = 1 i.e. omega is a generator of the
    // multiplicative subgroup H of Fq of order 'nconstraints'
    Field omega_temp = libff::power(omega, libff::bigint<1>(nconstraints));
    printf("w^%d: ", nconstraints);
    omega_temp.print();
    assert(omega_temp == 1);
#endif // #ifdef DEBUG

    // We want to represent the constraints q_L, q_R, q_O, q_M, q_C and
    // the witness w_L, w_R, w_O as polynomials in the roots of unity
    // e.g. f_{q_L}(omega_i) = q_L[i], 0\le{i}<8
    //    std::vector<Field> f = { 3, 4, 5, 9 };
    //    std::vector<Field> f = {1, 0, 0, 0, 0, 0, 0, 0};
    polynomial<Field> f{   Field(1),    Field(0),    Field(0),    Field(0),  Field(0),    Field(0),  Field(0),  Field(0)};
    const size_t m = nconstraints;
    //    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain;
    //    domain.reset(new libfqfft::basic_radix2_domain<Field>(m));
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain = libfqfft::get_evaluation_domain<Field>(m);

    // ---

    // compute Lagrange basis
    std::vector<polynomial<Field>> L(nconstraints);
    for (size_t i = 0; i < nconstraints; ++i) {
      polynomial<Field> u(nconstraints, Field(0));
      u[i] = Field(1);
      printf("[%s:%d] u[%d]\n", __FILE__, __LINE__, i);
      print_vector(u);      
      domain->iFFT(u);
      //      _basic_serial_radix2_FFT(u, omega.inverse(), true);
      // i-th Lagrange basis vector
      L[i] = u;
      printf("[%s:%d] L[%d]\n", __FILE__, __LINE__, i);
      print_vector(L[i]);      
    }
    
    // output from compute_qpoly() compute the q-polynomials from the
    // (transposed) gates matrix over the Lagrange basis q_poly =
    // \sum_i q[i] * L[i] where q[i] is a coefficient (a Field
    // element) and L[i] is a polynomial with Field coeffs
#if 1    
    std::vector<polynomial<Field>> Q(nqpoly, polynomial<Field>(nconstraints));
    for (size_t i = 0; i < nqpoly; ++i) {
      std::vector<Field> q_vec = gates_matrix_transpose[i];
      // a set of ncontraints polynomials. each polynomial is the
      // multiplication of the j-th element of q_vec to the j-th
      // polynomial in the lagrange basis
      std::vector<polynomial<Field>> q_poly(nconstraints);
      for (size_t j = 0; j < nconstraints; ++j) {
	// initialize the j-th slice of the i-th q-polynimial to the
	// j-th vector in the Lagrange basis
	q_poly[j] = L[j];
	// q[j] * L[j]: q[j] - scalar, L[j] - polynomial
	multiply_polynomial_by_scalar(q_vec[j], q_poly[j]);      
      }
      std::fill(Q[i].begin(), Q[i].end(), Field(0));
      for (size_t j = 0; j < nconstraints; ++j) {
	// Q[i] += q_poly[j];
	Q[i] = add_polynomials(Q[i], q_poly[j]);
      }      
    }
#endif    

    for (size_t i = 0; i < nqpoly; ++i) {
      printf("\n[%s:%d] Q[%2d]\n", __FILE__, __LINE__, i);
      print_vector(Q[i]);
    }

    Field omega_tmp = libff::get_root_of_unity<Field>(nconstraints);
    printf("[%s:%d] Lifqfft omega_tmp ", __FILE__, __LINE__);
    omega_tmp.print();
    printf("[%s:%d] PlonkPy omega     ", __FILE__, __LINE__);
    roots[1].print();
    
    for (size_t i = 0; i < nconstraints+1; ++i) {
      Field omega_i_tmp = libff::power(omega_tmp, libff::bigint<1>(i));
      printf("[%s:%d] omega_tmp[%d] ", __FILE__, __LINE__, i);
      omega_i_tmp.print();
    }
    
    for (size_t i = 0; i < nconstraints; ++i) {
      printf("w^%d: ", i);
      roots[i].print();
    }
    
    //    assert(0);
    printf("[%s:%d] Test OK\n", __FILE__, __LINE__);
  }

  TEST(TestPlonk, BLS12_381)
  {
    test_plonk<libff::bls12_381_pp>();
  }

} // namespace libsnark
