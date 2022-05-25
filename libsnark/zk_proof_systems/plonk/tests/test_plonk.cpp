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
#include <gtest/gtest.h>

#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>
#include <libff/algebra/curves/curve_serialization.hpp>

#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>
#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>
#include <libfqfft/polynomial_arithmetic/naive_evaluate.hpp>

#include <libsnark/zk_proof_systems/plonk/plonk.hpp>
#include <libsnark/zk_proof_systems/plonk/tests/plonk_example.hpp>

#define DEBUG 1
// maximum polynomial degree (resp. maximum number of gates)
const size_t MAX_DEGREE = 254;

namespace libsnark
{

  // Compute the Lagrange interpolation of the witness
  // points w_i. The result is the W_polys polynomials.
  // INPUT
  // - witness: witness values
  // - common preprocessed input
  // OUTPUT
  // - W_polys: Lagrange interpolation of the witness values
  template<typename ppT>
  void plonk_prover_compute_witness_polys(
					  std::vector<polynomial<libff::Fr<ppT>>>& W_polys,
					  const std::vector<libff::Fr<ppT>> witness,
					  const common_preprocessed_input<ppT> common_input
					  )
  {
    using Field = libff::Fr<ppT>;
    
    int nwitness = 3;
    for (int i = 0; i < nwitness; ++i) {
      typename std::vector<Field>::const_iterator begin = witness.begin()+(i*common_input.num_gates);
      typename std::vector<Field>::const_iterator end = witness.begin()+(i*common_input.num_gates)+(common_input.num_gates);
      std::vector<Field> W_points(begin, end);
      plonk_interpolate_over_lagrange_basis<Field>(W_points, W_polys[i], common_input.L_basis);
    }

  };

  
  template<typename ppT>
  plonk_proof<ppT> plonk_prover(
				const srs<ppT> srs,
				const common_preprocessed_input<ppT> common_input
				)
  {
    // return class
    plonk_proof<ppT> proof;
    
    using Field = libff::Fr<ppT>;
    // initialize hard-coded values from example circuit
    plonk_example<ppT> example;
    
    std::vector<Field> witness = example.witness;
    Field k1 = example.k1;
    Field k2 = example.k2;
    int num_hgen = NUM_HGEN;
    int nwitness = 3;

    printf("[%s:%d] Prover Round 1...\n", __FILE__, __LINE__);
#if 1 // prover round 1
    
    std::vector<polynomial<Field>> W_polys(nwitness, polynomial<Field>(common_input.num_gates));
    plonk_prover_compute_witness_polys<ppT>(W_polys, example.witness, common_input);
    
#if 1 // DEBUG
    printf("[%s:%d] Output from Prover setup\n", __FILE__, __LINE__);
    for (int i = 0; i < nwitness; ++i) {
      printf("[%s:%d] W_polys[%d]\n", __FILE__, __LINE__, i);
      print_vector(W_polys[i]);
      assert(W_polys[i] == example.W_polys[i]);
    }
#endif // #if 1 // DEBUG
    
    // vanishing polynomial zh_poly(X) = x^n-1. vanishes on all n roots of
    // unity common_input.omega_roots
    std::vector<Field> zh_poly(common_input.num_gates+1, Field(0));
    zh_poly[0] = Field(-1);
    zh_poly[common_input.num_gates] = Field(1);
    printf("[%s:%d] Vanishing polynomial\n", __FILE__, __LINE__);
    print_vector(zh_poly);

    // blinding scalars b1, b2, ..., b9. random but fixed to match the
    // python test vectors
    std::vector<Field> blind_scalars = example.prover_blind_scalars;
#ifdef DEBUG
    printf("[%s:%d] blind_scalars\n", __FILE__, __LINE__);
    print_vector(blind_scalars);
#endif // #ifdef DEBUG

    std::vector<std::vector<Field>> blind_polys
      {
       {blind_scalars[1], blind_scalars[0]}, // b1 + b0 X
       {blind_scalars[3], blind_scalars[2]}, // b3 + b2 X
       {blind_scalars[5], blind_scalars[4]}  // b5 + b4 X
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
	   (int)srs.secret_powers_g1.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(srs.secret_powers_g1, W_polys_blinded, W_polys_blinded_at_secret_g1);
#ifdef DEBUG
    printf("[%s:%d] Output from Round 1\n", __FILE__, __LINE__);
    for (int i = 0; i < nwitness; ++i) {
      printf("W_polys_at_secret_g1[%d]\n", i);
      W_polys_blinded_at_secret_g1[i].print();
    }
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 1

    printf("[%s:%d] Prover Round 2...\n", __FILE__, __LINE__);
#if 1 // prover round 2
    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field beta = example.beta;
    Field gamma = example.gamma;

    // compute permutation polynomial

    // blinding polynomial
    std::vector<Field> z1_blind_poly{blind_scalars[8], blind_scalars[7], blind_scalars[6]}; // b8 + b7 X + b6 X^2
    // multiply by the vanishing polynomial: z1 = z1 * zh_poly
    libfqfft::_polynomial_multiplication<Field>(z1_blind_poly, z1_blind_poly, zh_poly);
#ifdef DEBUG
    printf("[%s:%d] z1_blind_poly * zh_poly\n", __FILE__, __LINE__);
    print_vector(z1_blind_poly);
#endif // #ifdef DEBUG

    // A[0] = 1; ... A[i] = computed from (i-1)
    std::vector<Field> A_vector(common_input.num_gates, Field(0));
    plonk_compute_accumulator(common_input.num_gates, beta, gamma, witness, common_input.H_gen, common_input.H_gen_permute, A_vector);
#ifdef DEBUG
    for (int i = 0; i < (int)common_input.num_gates; ++i) {
      printf("A[%d] ", i);
      A_vector[i].print();
    }
#endif // #ifdef DEBUG

    polynomial<Field> A_poly(common_input.num_gates);
    plonk_interpolate_over_lagrange_basis<Field>(A_vector, A_poly, common_input.L_basis);
#ifdef DEBUG
    printf("[%s:%d] A_poly\n", __FILE__, __LINE__);
    print_vector(A_poly);
#endif // #ifdef DEBUG

    // add blinding polynomial z_1 to the accumulator polynomial A_poly
    polynomial<Field> z_poly;
    libfqfft::_polynomial_addition<Field>(z_poly, z1_blind_poly, A_poly);
#ifdef DEBUG
    printf("[%s:%d] z_poly\n", __FILE__, __LINE__);
    print_vector(z_poly);
#endif // #ifdef DEBUG
    
    libff::G1<ppT> z_poly_at_secret_g1 = plonk_evaluate_poly_at_secret_G1<ppT>(srs.secret_powers_g1, z_poly);
#ifdef DEBUG
    printf("[%s:%d] Output from Round 2\n", __FILE__, __LINE__);
    printf("[%s:%d] z_poly_at_secret_g1\n", __FILE__, __LINE__);
    z_poly_at_secret_g1.print();
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 2

    printf("[%s:%d] Prover Round 3...\n", __FILE__, __LINE__);
#if 1 // prover round 3
    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field alpha = example.alpha;
    
    // plonk_poly_shifted    
    // Computing the polynomial z(x*w) i.e. z(x) shifted by w where
    // w=common_input.omega_roots is the base root of unity and z is z_poly. we do this
    // by multiplying the coefficients of z by w
    std::vector<Field> z_poly_xomega_roots(z_poly.size(), Field(0));
    for (size_t i = 0; i < z_poly.size(); ++i) {
      // common_input.omega_roots^i
      //      z_poly_xomega_roots[i] = z_poly[i] * common_input.omega_roots[base][1]; // !!!! <----- (common_input.omega_roots[base][i]**i)
      Field omega_roots_i = libff::power(common_input.omega_roots[base][1], libff::bigint<1>(i));
      z_poly_xomega_roots[i] = z_poly[i] * omega_roots_i;
    }

#ifdef DEBUG
    printf("[%s:%d] z_poly_xomega_roots\n", __FILE__, __LINE__);
    print_vector(z_poly_xomega_roots);
#endif // #ifdef DEBUG

    // start computation of polynomial t(X) in round 3. we break t
    // into 4 parts which we compute separately. each of the 4 parts
    // is multiplied by 1/zh_poly in the paper
    std::vector<polynomial<Field>> t_part(4);

    // --- Computation of t_part[0]

    // a(x)b(x)q_M(x)
    polynomial<Field> abqM;
    libfqfft::_polynomial_multiplication<Field>(abqM, W_polys_blinded[a], W_polys_blinded[b]);
    libfqfft::_polynomial_multiplication<Field>(abqM, abqM, common_input.Q_polys[M]);    
    // a(x)q_L(x)
    polynomial<Field> aqL;
    libfqfft::_polynomial_multiplication<Field>(aqL, W_polys_blinded[a], common_input.Q_polys[L]);
    // b(x)q_R(x)
    polynomial<Field> bqR;
    libfqfft::_polynomial_multiplication<Field>(bqR, W_polys_blinded[b], common_input.Q_polys[R]);
    // c(x)q_O(x)
    polynomial<Field> cqO;
    libfqfft::_polynomial_multiplication<Field>(cqO, W_polys_blinded[c], common_input.Q_polys[O]);
    // t_part[0](x) = a(x)b(x)q_M(x) + a(x)q_L(x) + b(x)q_R(x) + c(x)q_O(x) + PI(x) + q_C(x)
    polynomial<Field> poly_null{Field(0)};
    libfqfft::_polynomial_addition<Field>(t_part[0], poly_null, abqM);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], aqL);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], bqR);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], cqO);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], common_input.PI_poly);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], common_input.Q_polys[C]);
    
#ifdef NDEBUG
    printf("[%s:%d] W_polys_blinded[a]\n", __FILE__, __LINE__);
    print_vector(W_polys_blinded[a]);
    printf("[%s:%d] W_polys_blinded[b]\n", __FILE__, __LINE__);
    print_vector(W_polys_blinded[b]);
    printf("[%s:%d] common_input.Q_polys[M]\n", __FILE__, __LINE__);
    print_vector(common_input.Q_polys[M]);
    printf("[%s:%d] abqM\n", __FILE__, __LINE__);
    print_vector(abqM);    
    printf("[%s:%d] t_part[0]\n", __FILE__, __LINE__);
    print_vector(t_part[0]);    
    assert(0);    
#endif // #ifdef DEBUG
    
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
    libfqfft::_polynomial_addition<Field>(a_xbeta_gamma, W_polys_blinded[a], xbeta_poly[base]);    
    libfqfft::_polynomial_addition<Field>(a_xbeta_gamma, a_xbeta_gamma, gamma_poly);    
    // b(x) + beta_k1*x + gamma 
    polynomial<Field> b_xbeta_gamma_k1;
    libfqfft::_polynomial_addition<Field>(b_xbeta_gamma_k1, W_polys_blinded[b], xbeta_poly[base_k1]);
    libfqfft::_polynomial_addition<Field>(b_xbeta_gamma_k1, b_xbeta_gamma_k1, gamma_poly);
    // c(x) + beta_k1*x + gamma 
    polynomial<Field> c_xbeta_gamma_k2;
    libfqfft::_polynomial_addition<Field>(c_xbeta_gamma_k2, W_polys_blinded[c], xbeta_poly[base_k2]);
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
    std::vector<polynomial<Field>> sbeta_poly(num_hgen);
    for (int i = 0; i < num_hgen; ++i) {
      libfqfft::_polynomial_multiplication<Field>(sbeta_poly[i], common_input.S_polys[i], beta_poly);
    }
    // a(x) + S_sigma1(x)*beta + gamma
    polynomial<Field> a_sbeta_gamma;
    libfqfft::_polynomial_addition<Field>(a_sbeta_gamma, W_polys_blinded[a], sbeta_poly[base]);    
    libfqfft::_polynomial_addition<Field>(a_sbeta_gamma, a_sbeta_gamma, gamma_poly);    
    // b(x) + S_sigma2(x)*beta + gamma
    polynomial<Field> b_sbeta_gamma_k1;
    libfqfft::_polynomial_addition<Field>(b_sbeta_gamma_k1, W_polys_blinded[b], sbeta_poly[base_k1]);    
    libfqfft::_polynomial_addition<Field>(b_sbeta_gamma_k1, b_sbeta_gamma_k1, gamma_poly);    
    // b(x) + S_sigma2(x)*beta + gamma
    polynomial<Field> c_sbeta_gamma_k2;
    libfqfft::_polynomial_addition<Field>(c_sbeta_gamma_k2, W_polys_blinded[c], sbeta_poly[base_k2]);    
    libfqfft::_polynomial_addition<Field>(c_sbeta_gamma_k2, c_sbeta_gamma_k2, gamma_poly);    
    // t_part[2] = (a(x) + S_sigma1(x)*beta + gamma)*(b(x) +
    // S_sigma2(x)*beta + gamma)*(b(x) + S_sigma2(x)*beta +
    // gamma)*z(x*common_input.omega_roots)*alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[2], a_sbeta_gamma, b_sbeta_gamma_k1);
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], c_sbeta_gamma_k2);
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], z_poly_xomega_roots);
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], alpha_poly);
    // -t_part[2]
    polynomial<Field> neg_one_poly = {-Field("1")};
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], neg_one_poly);
    
    // --- Computation of t_part[3]

    // z(x) - 1
    polynomial<Field> z_neg_one;
    libfqfft::_polynomial_addition<Field>(z_neg_one, z_poly, neg_one_poly);
    // (z(x)-1) * L_1(x)
    libfqfft::_polynomial_multiplication<Field>(t_part[3], z_neg_one, common_input.L_basis[0]);
    // (z(x)-1) * L_1(x) * alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[3], t_part[3], alpha_poly);
    // (z(x)-1) * L_1(x) * alpha * alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[3], t_part[3], alpha_poly);


    // --- computation of t(x)

    // t(x) = (t[0] + t[1] + (-t[2]) + t[3]) / zh(x)
    polynomial<Field> t_poly_long{Field(0)};
    libfqfft::_polynomial_addition<Field>(t_poly_long, t_poly_long, t_part[0]);
#if 1 // DEBUG   
    libfqfft::_polynomial_addition<Field>(t_poly_long, t_poly_long, t_part[1]);
    libfqfft::_polynomial_addition<Field>(t_poly_long, t_poly_long, t_part[2]);
    libfqfft::_polynomial_addition<Field>(t_poly_long, t_poly_long, t_part[3]);
#endif // #if 0 // DEBUG   
    //    t(x) = t(x) / zh(x): A/B = (Q, R) st. A = (Q * B) + R.
    polynomial<Field> remainder;
    libfqfft::_polynomial_division(t_poly_long, remainder, t_poly_long, zh_poly);
#ifdef DEBUG
    printf("[%s:%d] t_poly_long\n", __FILE__, __LINE__);
    print_vector(t_poly_long);
#endif // #ifdef DEBUG
#ifdef DEBUG
    printf("[%s:%d] remainder\n", __FILE__, __LINE__);
    print_vector(remainder);
#endif // #ifdef DEBUG
    assert(libfqfft::_is_zero(remainder));

    // break t_poly_long into three parts: lo, mid, hi, each of degree 7
    // note: (common_input.num_gates+3) is the length of the CRS = (common_input.num_gates+2) powers of G1 + 1 power of G2
    std::vector<polynomial<Field>> t_poly(num_hgen);
    for (int i = 0; i < num_hgen; ++i) {
      typename std::vector<Field>::iterator begin = t_poly_long.begin()+(i*(common_input.num_gates+2));
      typename std::vector<Field>::iterator end = t_poly_long.begin()+(i*(common_input.num_gates+2))+(common_input.num_gates+2);
      std::vector<Field> tmp(begin, end);
      t_poly[i] = tmp;
    }
#ifdef DEBUG
    for (int i = 0; i < num_hgen; ++i) {
      printf("[%s:%d] t_poly[%d]\n", __FILE__, __LINE__, i);
      print_vector(t_poly[i]);
    }
#endif // #ifdef DEBUG
    // evaluate each part of t_poly in the secret input
    std::vector<libff::G1<ppT>> t_poly_at_secret_g1(num_hgen);
    for (int i = 0; i < num_hgen; ++i) {
      t_poly_at_secret_g1[i] = plonk_evaluate_poly_at_secret_G1<ppT>(srs.secret_powers_g1, t_poly[i]);
    }
#ifdef DEBUG
    printf("[%s:%d] Output from Round 3\n", __FILE__, __LINE__);
    // verify the output from Round 3 to the test vectors. test
    // vectors obtained from the Plonk Python reference implementation
    // (used for debug)
    for (int i = 0; i < num_hgen; ++i) {
      printf("[%s:%d] t_poly_at_secret_g1[%d]\n", __FILE__, __LINE__, i);
      t_poly_at_secret_g1[i].print();
      libff::G1<ppT> t_poly_at_secret_g1_i(t_poly_at_secret_g1[i]);
      t_poly_at_secret_g1_i.to_affine_coordinates();
      assert(t_poly_at_secret_g1_i.X == example.t_poly_at_secret_g1[i][0]);
      assert(t_poly_at_secret_g1_i.Y == example.t_poly_at_secret_g1[i][1]);
    }
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 3
    
    printf("[%s:%d] Prover Round 4...\n", __FILE__, __LINE__);
#if 1 // prover round 4
    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field zeta = example.zeta;
#ifdef DEBUG
    printf("[%s:%d] zeta\n", __FILE__, __LINE__);
    zeta.print();
#endif // #ifdef DEBUG
    Field a_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates + 2, W_polys_blinded[a], zeta);
    Field b_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates + 2, W_polys_blinded[b], zeta);
    Field c_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates + 2, W_polys_blinded[c], zeta);
    Field S_0_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates, common_input.S_polys[0], zeta);
    Field S_1_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates, common_input.S_polys[1], zeta);
    Field t_zeta = libfqfft::evaluate_polynomial<Field>(t_poly_long.size(), t_poly_long, zeta);
    Field z_poly_xomega_zeta = libfqfft::evaluate_polynomial<Field>(z_poly_xomega_roots.size(), z_poly_xomega_roots, zeta);
    
#ifdef DEBUG
    printf("[%s:%d] Output from Round 4\n", __FILE__, __LINE__);
    printf("a_zeta ");
    a_zeta.print();
    assert(a_zeta == example.a_zeta);
    printf("b_zeta ");
    b_zeta.print();
    assert(b_zeta == example.b_zeta);
    printf("c_zeta ");
    c_zeta.print();
    assert(c_zeta == example.c_zeta);
    printf("S_0_zeta ");
    S_0_zeta.print();
    assert(S_0_zeta == example.S_0_zeta);
    printf("S_1_zeta ");
    S_1_zeta.print();
    assert(S_1_zeta == example.S_1_zeta);
    printf("t_zeta ");
    t_zeta.print();
    assert(t_zeta == example.t_zeta);
    printf("z_poly_xomega_zeta ");
    z_poly_xomega_zeta.print();
    assert(z_poly_xomega_zeta == example.z_poly_xomega_zeta);
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 4
    
    printf("[%s:%d] Prover Round 5...\n", __FILE__, __LINE__);
#if 1 // prover round 5
    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field nu = example.nu;
    
    // compute linerisation polynomial r in five parts
    std::vector<polynomial<Field>> r_part(5);

    // --- Computation of r_part[0]
    
    // represent values as constant term polynomials in orderto use
    // the functions in the libfqfft library on polynomials
    polynomial<Field> a_zeta_poly{a_zeta}; 
    polynomial<Field> b_zeta_poly{b_zeta}; 
    polynomial<Field> c_zeta_poly{c_zeta}; 
    // a(z)b(z)q_M(x)
    polynomial<Field> abqM_zeta;
    libfqfft::_polynomial_multiplication<Field>(abqM_zeta, common_input.Q_polys[M], a_zeta_poly);
    libfqfft::_polynomial_multiplication<Field>(abqM_zeta, abqM_zeta, b_zeta_poly);    
    // a(z)q_L(x)
    polynomial<Field> aqL_zeta;
    libfqfft::_polynomial_multiplication<Field>(aqL_zeta, common_input.Q_polys[L], a_zeta_poly);
    // b(z)q_R(x)
    polynomial<Field> bqR_zeta;
    libfqfft::_polynomial_multiplication<Field>(bqR_zeta, common_input.Q_polys[R], b_zeta_poly);
    // c(z)q_O(x)
    polynomial<Field> cqO_zeta;
    libfqfft::_polynomial_multiplication<Field>(cqO_zeta, common_input.Q_polys[O], c_zeta_poly);
    // a(z)b(z)q_M(x) + a(z)q_L(x) + b(z)q_R(x) + c(z)q_O(x) + q_C(x)
    libfqfft::_polynomial_addition<Field>(r_part[0], poly_null, abqM_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], aqL_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], bqR_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], cqO_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], common_input.Q_polys[C]);

    // --- Computation of r_part[1]

    polynomial<Field> r1_const_poly
      {
       (a_zeta + (beta * zeta) + gamma) *
       (b_zeta + (beta * k1 * zeta) + gamma) *
       (c_zeta + (beta * k2 * zeta) + gamma) * alpha
      };
    libfqfft::_polynomial_multiplication<Field>(r_part[1], r1_const_poly, z_poly);

    // --- Computation of r_part[2]
    
    polynomial<Field> r2_const_poly
      {
       (a_zeta + (beta * S_0_zeta) + gamma) *
       (b_zeta + (beta * S_1_zeta) + gamma) *
       (alpha * beta * z_poly_xomega_zeta)
      };
    libfqfft::_polynomial_multiplication<Field>(r_part[2], r2_const_poly, common_input.S_polys[2]);
    // -r_part[2]
    libfqfft::_polynomial_multiplication<Field>(r_part[2], r_part[2], neg_one_poly);
    
    // --- Computation of r_part[3]
    
    //     r3 = accumulator_poly_ext3 * eval_poly(L_1, [zeta])[0] * alpha ** 2
    polynomial<Field> L_0_zeta_poly{libfqfft::evaluate_polynomial<Field>(common_input.L_basis[0].size(), common_input.L_basis[0], zeta)};
    polynomial<Field> alpha_power2_poly{libff::power(alpha, libff::bigint<1>(2))};
    libfqfft::_polynomial_multiplication<Field>(r_part[3], z_poly, L_0_zeta_poly);
    libfqfft::_polynomial_multiplication<Field>(r_part[3], r_part[3], alpha_power2_poly);

    // --- Computation of r_poly = (r0+r1-r2+r3)

    //
    // Note: here the reference Python implementation differs from the
    // paper where:
    //
    // r(x) = r(x) - zh(zeta) (t_lo(x) + zeta^n t_mid(x) + zeta^2n t_hi(x))
    //
    // In the reference implementation, the missing term is added in
    // the computation of the W_zeta(x) polynomial
    //
    polynomial<Field> r_poly; 
    libfqfft::_polynomial_addition<Field>(r_poly, poly_null, r_part[0]);
    libfqfft::_polynomial_addition<Field>(r_poly, r_poly, r_part[1]);
    libfqfft::_polynomial_addition<Field>(r_poly, r_poly, r_part[2]);
    libfqfft::_polynomial_addition<Field>(r_poly, r_poly, r_part[3]);
    
#if DEBUG    
    //    printf("abqM_zeta\n");
    //    print_vector(abqM_zeta);
    printf("[%s:%d] r_part[0]\n", __FILE__, __LINE__);
    print_vector(r_part[0]);
    printf("[%s:%d] r_part[1]\n", __FILE__, __LINE__);
    print_vector(r_part[1]);
    printf("[%s:%d] r_part[2]\n", __FILE__, __LINE__);
    print_vector(r_part[2]);
    printf("[%s:%d] r_part[3]\n", __FILE__, __LINE__);
    print_vector(r_part[3]);
    printf("[%s:%d] r_poly\n", __FILE__, __LINE__);
    print_vector(r_poly);
#endif // #if DEBUG

    // Evaluate the r-polynomial at zeta. Note: in the reference
    // implementation, r_zeta is added to the pi-SNARK proof. In the
    // paper this is omitted, which makes the proof shorter at the
    // epxense of a slightly heavier computation on the verifier's
    // side 
    Field r_zeta = libfqfft::evaluate_polynomial<Field>(r_poly.size(), r_poly, zeta);
#ifdef DEBUG
    printf("r_zeta ");
    r_zeta.print();
#endif // #ifdef DEBUG    
    assert(r_zeta == example.r_zeta);

    // W_zeta polynomial is of degree 6 in the random element nu and
    // hence has 7 terms
    std::vector<polynomial<Field>> W_zeta_part(7);

    // --- compute W_zeta_part[0]
    
    // t_lo(x)
    polynomial<Field> t_lo{t_poly[lo]};
    // t_mid(x) * zeta^(n+2)
    polynomial<Field> t_mid_zeta_n;
    polynomial<Field> zeta_powern_poly{libff::power(zeta, libff::bigint<1>(common_input.num_gates+2))};
    libfqfft::_polynomial_multiplication<Field>(t_mid_zeta_n, t_poly[mid], zeta_powern_poly);
    // t_hi(x) * zeta^(2(n+1))
    polynomial<Field> t_hi_zeta_2n;
    polynomial<Field> zeta_power2n_poly{libff::power(zeta, libff::bigint<1>(2*(common_input.num_gates+2)))};
    libfqfft::_polynomial_multiplication<Field>(t_hi_zeta_2n, t_poly[hi], zeta_power2n_poly);
    // -t_zeta as constant term polynomial
    polynomial<Field> t_zeta_poly{-t_zeta};
    // t_lo(x) + (t_mid(x) * zeta^n) + (t_hi(x) * zeta^2n) + t_zeta_poly
    libfqfft::_polynomial_addition<Field>(W_zeta_part[0], poly_null, t_lo);
    libfqfft::_polynomial_addition<Field>(W_zeta_part[0], W_zeta_part[0], t_mid_zeta_n);
    libfqfft::_polynomial_addition<Field>(W_zeta_part[0], W_zeta_part[0], t_hi_zeta_2n);
    libfqfft::_polynomial_addition<Field>(W_zeta_part[0], W_zeta_part[0], t_zeta_poly);    
    
    // --- compute W_zeta_part[1]

    // -r_zeta as constant term polynomial
    polynomial<Field> r_zeta_poly{-r_zeta};
    // r(x) - r_zeta
    polynomial<Field> r_sub_rzeta;
    libfqfft::_polynomial_addition<Field>(r_sub_rzeta, r_poly, r_zeta_poly);
    // (r(x) - r_zeta) * nu
    polynomial<Field> nu_poly{nu};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[1], r_sub_rzeta, nu_poly);
    
    // --- compute W_zeta_part[2]
    
    // -a_zeta as constant term polynomial
    polynomial<Field> a_zeta_poly_neg;
    libfqfft::_polynomial_multiplication<Field>(a_zeta_poly_neg, a_zeta_poly, neg_one_poly);
    // a(x) - a_zeta
    polynomial<Field> a_sub_azeta;
    libfqfft::_polynomial_addition<Field>(a_sub_azeta, W_polys_blinded[a], a_zeta_poly_neg);
    // (a(x) - a_zeta) * nu^2
    Field nu2 = libff::power(nu, libff::bigint<1>(2));
    polynomial<Field> nu2_poly{nu2};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[2], a_sub_azeta, nu2_poly);
    
    // -b_zeta as constant term polynomial
    polynomial<Field> b_zeta_poly_neg;
    libfqfft::_polynomial_multiplication<Field>(b_zeta_poly_neg, b_zeta_poly, neg_one_poly);
    // (b(x) - b_zeta)
    polynomial<Field> b_sub_bzeta;
    libfqfft::_polynomial_addition<Field>(b_sub_bzeta, W_polys_blinded[b], b_zeta_poly_neg);
    // (b(x) - b_zeta) * nu^3
    Field nu3 = libff::power(nu, libff::bigint<1>(3));
    polynomial<Field> nu3_poly{nu3};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[3], b_sub_bzeta, nu3_poly);

    // -c_zeta as constant term polynomial
    polynomial<Field> c_zeta_poly_neg;
    libfqfft::_polynomial_multiplication<Field>(c_zeta_poly_neg, c_zeta_poly, neg_one_poly);
    // (c(x) - c_zeta)
    polynomial<Field> c_sub_czeta;
    libfqfft::_polynomial_addition<Field>(c_sub_czeta, W_polys_blinded[c], c_zeta_poly_neg);
    // (c(x) - c_zeta) * nu^4
    Field nu4 = libff::power(nu, libff::bigint<1>(4));
    polynomial<Field> nu4_poly{nu4};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[4], c_sub_czeta, nu4_poly);

    // -S_0_zeta as constant term polynomial
    polynomial<Field> S_0_zeta_poly_neg{-S_0_zeta};
    //    libfqfft::_polynomial_multiplication<Field>(S_0_zeta_poly_neg, S_0_zeta_poly, neg_one_poly);
    // (S0(x) - S_0_zeta)
    polynomial<Field> S0_sub_szeta;
    libfqfft::_polynomial_addition<Field>(S0_sub_szeta, common_input.S_polys[0], S_0_zeta_poly_neg);
    // (S0(x) - S_0_zeta) * nu^5
    Field nu5 = libff::power(nu, libff::bigint<1>(5));
    polynomial<Field> nu5_poly{nu5};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[5], S0_sub_szeta, nu5_poly);

    // -S_1_zeta as constant term polynomial
    polynomial<Field> S_1_zeta_poly_neg{-S_1_zeta};
    //    libfqfft::_polynomial_multiplication<Field>(S_1_zeta_poly_neg, S_1_zeta_poly, neg_one_poly);
    // (S1(x) - S_1_zeta)
    polynomial<Field> S1_sub_szeta;
    libfqfft::_polynomial_addition<Field>(S1_sub_szeta, common_input.S_polys[1], S_1_zeta_poly_neg);
    // (S1(x) - S_1_zeta) * nu^6
    Field nu6 = libff::power(nu, libff::bigint<1>(6));
    polynomial<Field> nu6_poly{nu6};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[6], S1_sub_szeta, nu6_poly);

    // compute full zeta polynomial W_zeta = \sum W_zeta_part[i]
    int nzeta = 7;
    polynomial<Field> W_zeta(poly_null);
    for (int i = 0; i < nzeta; ++i) {
      libfqfft::_polynomial_addition<Field>(W_zeta, W_zeta, W_zeta_part[i]);
    }

    // compute 1/(X-zeta) * W_zeta
    polynomial<Field> x_sub_zeta_poly{-zeta, Field(1)};
    libfqfft::_polynomial_division(W_zeta, remainder, W_zeta, x_sub_zeta_poly);
#ifdef DEBUG
    printf("W_zeta\n");
    print_vector(W_zeta);
#endif // #ifdef DEBUG
    assert(libfqfft::_is_zero(remainder));

    //    polynomial<Field> z_poly;
    
    // Compute opening proof:
    // W_zeta_omega = z(X) - z(zeta*common_input.omega_roots) / X - (zeta*common_input.omega_roots)
    polynomial<Field> W_zeta_omega{poly_null};
    
    // -z(zeta*common_input.omega_roots)
    polynomial<Field> z_poly_xomega_zeta_neg{-z_poly_xomega_zeta};
    // z(X) - z(zeta*common_input.omega_roots) 
    libfqfft::_polynomial_addition<Field>(W_zeta_omega, z_poly, z_poly_xomega_zeta_neg);    
    // -zeta*common_input.omega_roots; common_input.omega_roots[base][1] = common_input.omega_roots_base
    polynomial<Field> x_sub_zeta_omega_roots{-(zeta*common_input.omega_roots[base][1]), Field(1)};
    
    // z(X) - z(zeta*common_input.omega_roots) / X - (zeta*common_input.omega_roots)
    libfqfft::_polynomial_division(W_zeta_omega, remainder, W_zeta_omega, x_sub_zeta_omega_roots);
    assert(libfqfft::_is_zero(remainder));

#ifdef DEBUG
    printf("W_zeta_part[0]\n");
    print_vector(W_zeta_part[0]);
    printf("W_zeta_part[1]\n");
    print_vector(W_zeta_part[1]);
    printf("W_zeta_part[2]\n");
    print_vector(W_zeta_part[2]);
    printf("W_zeta_part[3]\n");
    print_vector(W_zeta_part[3]);
    printf("W_zeta_part[4]\n");
    print_vector(W_zeta_part[4]);
    printf("W_zeta_part[5]\n");
    print_vector(W_zeta_part[5]);
    printf("W_zeta_part[6]\n");
    print_vector(W_zeta_part[6]);
    printf("W_zeta\n");
    print_vector(W_zeta);
    printf("W_zeta_omega\n");
    print_vector(W_zeta_omega);
#endif // #ifdef DEBUG
      
    assert(W_zeta == example.W_zeta);
    assert(W_zeta_omega == example.W_zeta_omega);
    
    // Evaluate polynomials W_zeta and W_zeta_omega at the seceret
    // input
    libff::G1<ppT> W_zeta_at_secret =
      plonk_evaluate_poly_at_secret_G1<ppT>(srs.secret_powers_g1, W_zeta);    
    libff::G1<ppT> W_zeta_omega_at_secret =
      plonk_evaluate_poly_at_secret_G1<ppT>(srs.secret_powers_g1, W_zeta_omega);
    
#ifdef DEBUG
    printf("[%s:%d] Outputs from Prover round 5\n", __FILE__, __LINE__);
    
    printf("r_zeta ");
    r_zeta.print();
    
    printf("[%s:%d] W_zeta_at_secret \n", __FILE__, __LINE__);
    W_zeta_at_secret.print();
    libff::G1<ppT> W_zeta_at_secret_aff(W_zeta_at_secret);
    W_zeta_at_secret_aff.to_affine_coordinates();
    assert(W_zeta_at_secret_aff.X == example.W_zeta_at_secret[0]);
    assert(W_zeta_at_secret_aff.Y == example.W_zeta_at_secret[1]);
    
    printf("[%s:%d] W_zeta_omega_at_secret \n", __FILE__, __LINE__);
    W_zeta_omega_at_secret.print();
    libff::G1<ppT> W_zeta_omega_at_secret_aff(W_zeta_omega_at_secret);
    W_zeta_omega_at_secret_aff.to_affine_coordinates();
    assert(W_zeta_omega_at_secret_aff.X == example.W_zeta_omega_at_secret[0]);
    assert(W_zeta_omega_at_secret_aff.Y == example.W_zeta_omega_at_secret[1]);
#endif // #ifdef DEBUG

    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field u = example.u;
#endif // #if 1 // prover round 5
    
    return proof;
  }

  
  template<typename ppT> void test_plonk()
  {
    // Execute all tests for the given curve.
    ppT::init_public_params();

    using Field = libff::Fr<ppT>;
    //    using BaseField = libff::Fq<ppT>;

    // initialize hard-coded values from example circuit
    plonk_example<ppT> example;
    
    // --- SETUP ---
    printf("[%s:%d] Start setup...\n", __FILE__, __LINE__);

    // common preprocessed input    
    common_preprocessed_input<ppT> common_input; 
    plonk_setup<ppT>(common_input);

    // random hidden element secret (toxic waste). we fix it to a
    // constant in order to match against the test vectors
    Field secret = example.secret;
#ifdef DEBUG
    printf("[%s:%d] secret ", __FILE__, __LINE__);
    secret.print();
#endif // #ifdef DEBUG

    // --- SRS ---
    printf("[%s:%d] Prepare SRS...\n", __FILE__, __LINE__);    
    // compute powers of secret times G1: 1*G1, secret^1*G1,
    // secret^2*G1, ... and secret times G2: 1*G2, secret^1*G2
    srs<ppT> srs = plonk_derive_srs_from_secret<ppT>(secret, common_input.num_gates);    
    // Compare against reference test values
#ifdef DEBUG
    for (int i = 0; i < (int)common_input.num_gates + 3; ++i) {
      printf("secret_power_G1[%2d] ", i);
      srs.secret_powers_g1[i].print();
      // test from generator
      libff::G1<ppT> srs_secret_powers_g1_i(srs.secret_powers_g1[i]);
      srs_secret_powers_g1_i.to_affine_coordinates();
      assert(srs_secret_powers_g1_i.X == example.secret_powers_g1[i][0]);
      assert(srs_secret_powers_g1_i.Y == example.secret_powers_g1[i][1]);
    }
    for (int i = 0; i < 2; ++i) {
      printf("secret_power_G2[%2d] ", i);
      srs.secret_powers_g2[i].print();
    }
#endif // #ifdef DEBUG
    
    plonk_proof<ppT> proof = plonk_prover(srs, common_input);

    // --- PROVER ---
    
    printf("[%s:%d] Prover preparation...\n", __FILE__, __LINE__);
#if 1 // prover preparation    
    // Generate domains on which to evaluate the witness
    // polynomials. k1,k2 can be random, but we fix them for debug to
    // match against the test vector values
    Field k1 = example.k1;
    Field k2 = example.k2;
    // number of generators for H, Hk1, Hk2
    int num_hgen = NUM_HGEN;
    
    // witness values
    std::vector<Field> witness = example.witness;
    
    int nwitness = 3;
    std::vector<polynomial<Field>> W_polys(nwitness, polynomial<Field>(common_input.num_gates));
    for (int i = 0; i < nwitness; ++i) {
      typename std::vector<Field>::iterator begin = witness.begin()+(i*common_input.num_gates);
      typename std::vector<Field>::iterator end = witness.begin()+(i*common_input.num_gates)+(common_input.num_gates);
      std::vector<Field> W_points(begin, end);
      plonk_interpolate_over_lagrange_basis<Field>(W_points, W_polys[i], common_input.L_basis);
    }

#if 1 // DEBUG
    for (int i = 0; i < nwitness; ++i) {
      printf("[%s:%d] W_polys[%d]\n", __FILE__, __LINE__, i);
      print_vector(W_polys[i]);
    }
#endif // #if 1 // DEBUG
#endif // #if 1 // prover preparation    
    
    printf("[%s:%d] Prover Round 1...\n", __FILE__, __LINE__);
#if 1 // prover round 1
    // vanishing polynomial zh_poly(X) = x^n-1. vanishes on all n roots of
    // unity common_input.omega_roots
    std::vector<Field> zh_poly(common_input.num_gates+1, Field(0));
    zh_poly[0] = Field(-1);
    zh_poly[common_input.num_gates] = Field(1);
    printf("[%s:%d] Vanishing polynomial\n", __FILE__, __LINE__);
    print_vector(zh_poly);

    // blinding scalars b1, b2, ..., b9. random but fixed to match the
    // python test vectors
    std::vector<Field> blind_scalars = example.prover_blind_scalars;
#ifdef DEBUG
    printf("[%s:%d] blind_scalars\n", __FILE__, __LINE__);
    print_vector(blind_scalars);
#endif // #ifdef DEBUG

    std::vector<std::vector<Field>> blind_polys
      {
       {blind_scalars[1], blind_scalars[0]}, // b1 + b0 X
       {blind_scalars[3], blind_scalars[2]}, // b3 + b2 X
       {blind_scalars[5], blind_scalars[4]}  // b5 + b4 X
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
	   (int)srs.secret_powers_g1.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(srs.secret_powers_g1, W_polys_blinded, W_polys_blinded_at_secret_g1);
#ifdef DEBUG
    printf("[%s:%d] Output from Round 1\n", __FILE__, __LINE__);
    for (int i = 0; i < nwitness; ++i) {
      printf("W_polys_at_secret_g1[%d]\n", i);
      W_polys_blinded_at_secret_g1[i].print();
    }
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 1
    
    printf("[%s:%d] Prover Round 2...\n", __FILE__, __LINE__);
#if 1 // prover round 2
    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field beta = example.beta;
    Field gamma = example.gamma;

    // compute permutation polynomial

    // blinding polynomial
    std::vector<Field> z1_blind_poly{blind_scalars[8], blind_scalars[7], blind_scalars[6]}; // b8 + b7 X + b6 X^2
    // multiply by the vanishing polynomial: z1 = z1 * zh_poly
    libfqfft::_polynomial_multiplication<Field>(z1_blind_poly, z1_blind_poly, zh_poly);
#ifdef DEBUG
    printf("[%s:%d] z1_blind_poly * zh_poly\n", __FILE__, __LINE__);
    print_vector(z1_blind_poly);
#endif // #ifdef DEBUG

    // A[0] = 1; ... A[i] = computed from (i-1)
    std::vector<Field> A_vector(common_input.num_gates, Field(0));
    plonk_compute_accumulator(common_input.num_gates, beta, gamma, witness, common_input.H_gen, common_input.H_gen_permute, A_vector);
#ifdef DEBUG
    for (int i = 0; i < (int)common_input.num_gates; ++i) {
      printf("A[%d] ", i);
      A_vector[i].print();
    }
#endif // #ifdef DEBUG

    polynomial<Field> A_poly(common_input.num_gates);
    plonk_interpolate_over_lagrange_basis<Field>(A_vector, A_poly, common_input.L_basis);
#ifdef DEBUG
    printf("[%s:%d] A_poly\n", __FILE__, __LINE__);
    print_vector(A_poly);
#endif // #ifdef DEBUG

    // add blinding polynomial z_1 to the accumulator polynomial A_poly
    polynomial<Field> z_poly;
    libfqfft::_polynomial_addition<Field>(z_poly, z1_blind_poly, A_poly);
#ifdef DEBUG
    printf("[%s:%d] z_poly\n", __FILE__, __LINE__);
    print_vector(z_poly);
#endif // #ifdef DEBUG
    
    libff::G1<ppT> z_poly_at_secret_g1 = plonk_evaluate_poly_at_secret_G1<ppT>(srs.secret_powers_g1, z_poly);
#ifdef DEBUG
    printf("[%s:%d] Output from Round 2\n", __FILE__, __LINE__);
    printf("[%s:%d] z_poly_at_secret_g1\n", __FILE__, __LINE__);
    z_poly_at_secret_g1.print();
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 2

    printf("[%s:%d] Prover Round 3...\n", __FILE__, __LINE__);
#if 1 // prover round 3
    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field alpha = example.alpha;
    
    // plonk_poly_shifted    
    // Computing the polynomial z(x*w) i.e. z(x) shifted by w where
    // w=common_input.omega_roots is the base root of unity and z is z_poly. we do this
    // by multiplying the coefficients of z by w
    std::vector<Field> z_poly_xomega_roots(z_poly.size(), Field(0));
    for (size_t i = 0; i < z_poly.size(); ++i) {
      // common_input.omega_roots^i
      //      z_poly_xomega_roots[i] = z_poly[i] * common_input.omega_roots[base][1]; // !!!! <----- (common_input.omega_roots[base][i]**i)
      Field omega_roots_i = libff::power(common_input.omega_roots[base][1], libff::bigint<1>(i));
      z_poly_xomega_roots[i] = z_poly[i] * omega_roots_i;
    }

#ifdef DEBUG
    printf("[%s:%d] z_poly_xomega_roots\n", __FILE__, __LINE__);
    print_vector(z_poly_xomega_roots);
#endif // #ifdef DEBUG

    // start computation of polynomial t(X) in round 3. we break t
    // into 4 parts which we compute separately. each of the 4 parts
    // is multiplied by 1/zh_poly in the paper
    std::vector<polynomial<Field>> t_part(4);

    // --- Computation of t_part[0]

    // a(x)b(x)q_M(x)
    polynomial<Field> abqM;
    libfqfft::_polynomial_multiplication<Field>(abqM, W_polys_blinded[a], W_polys_blinded[b]);
    libfqfft::_polynomial_multiplication<Field>(abqM, abqM, common_input.Q_polys[M]);    
    // a(x)q_L(x)
    polynomial<Field> aqL;
    libfqfft::_polynomial_multiplication<Field>(aqL, W_polys_blinded[a], common_input.Q_polys[L]);
    // b(x)q_R(x)
    polynomial<Field> bqR;
    libfqfft::_polynomial_multiplication<Field>(bqR, W_polys_blinded[b], common_input.Q_polys[R]);
    // c(x)q_O(x)
    polynomial<Field> cqO;
    libfqfft::_polynomial_multiplication<Field>(cqO, W_polys_blinded[c], common_input.Q_polys[O]);
    // t_part[0](x) = a(x)b(x)q_M(x) + a(x)q_L(x) + b(x)q_R(x) + c(x)q_O(x) + PI(x) + q_C(x)
    polynomial<Field> poly_null{Field(0)};
    libfqfft::_polynomial_addition<Field>(t_part[0], poly_null, abqM);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], aqL);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], bqR);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], cqO);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], common_input.PI_poly);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], common_input.Q_polys[C]);
    
#ifdef NDEBUG
    printf("[%s:%d] W_polys_blinded[a]\n", __FILE__, __LINE__);
    print_vector(W_polys_blinded[a]);
    printf("[%s:%d] W_polys_blinded[b]\n", __FILE__, __LINE__);
    print_vector(W_polys_blinded[b]);
    printf("[%s:%d] common_input.Q_polys[M]\n", __FILE__, __LINE__);
    print_vector(common_input.Q_polys[M]);
    printf("[%s:%d] abqM\n", __FILE__, __LINE__);
    print_vector(abqM);    
    printf("[%s:%d] t_part[0]\n", __FILE__, __LINE__);
    print_vector(t_part[0]);    
    assert(0);    
#endif // #ifdef DEBUG
    
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
    libfqfft::_polynomial_addition<Field>(a_xbeta_gamma, W_polys_blinded[a], xbeta_poly[base]);    
    libfqfft::_polynomial_addition<Field>(a_xbeta_gamma, a_xbeta_gamma, gamma_poly);    
    // b(x) + beta_k1*x + gamma 
    polynomial<Field> b_xbeta_gamma_k1;
    libfqfft::_polynomial_addition<Field>(b_xbeta_gamma_k1, W_polys_blinded[b], xbeta_poly[base_k1]);
    libfqfft::_polynomial_addition<Field>(b_xbeta_gamma_k1, b_xbeta_gamma_k1, gamma_poly);
    // c(x) + beta_k1*x + gamma 
    polynomial<Field> c_xbeta_gamma_k2;
    libfqfft::_polynomial_addition<Field>(c_xbeta_gamma_k2, W_polys_blinded[c], xbeta_poly[base_k2]);
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
    std::vector<polynomial<Field>> sbeta_poly(num_hgen);
    for (int i = 0; i < num_hgen; ++i) {
      libfqfft::_polynomial_multiplication<Field>(sbeta_poly[i], common_input.S_polys[i], beta_poly);
    }
    // a(x) + S_sigma1(x)*beta + gamma
    polynomial<Field> a_sbeta_gamma;
    libfqfft::_polynomial_addition<Field>(a_sbeta_gamma, W_polys_blinded[a], sbeta_poly[base]);    
    libfqfft::_polynomial_addition<Field>(a_sbeta_gamma, a_sbeta_gamma, gamma_poly);    
    // b(x) + S_sigma2(x)*beta + gamma
    polynomial<Field> b_sbeta_gamma_k1;
    libfqfft::_polynomial_addition<Field>(b_sbeta_gamma_k1, W_polys_blinded[b], sbeta_poly[base_k1]);    
    libfqfft::_polynomial_addition<Field>(b_sbeta_gamma_k1, b_sbeta_gamma_k1, gamma_poly);    
    // b(x) + S_sigma2(x)*beta + gamma
    polynomial<Field> c_sbeta_gamma_k2;
    libfqfft::_polynomial_addition<Field>(c_sbeta_gamma_k2, W_polys_blinded[c], sbeta_poly[base_k2]);    
    libfqfft::_polynomial_addition<Field>(c_sbeta_gamma_k2, c_sbeta_gamma_k2, gamma_poly);    
    // t_part[2] = (a(x) + S_sigma1(x)*beta + gamma)*(b(x) +
    // S_sigma2(x)*beta + gamma)*(b(x) + S_sigma2(x)*beta +
    // gamma)*z(x*common_input.omega_roots)*alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[2], a_sbeta_gamma, b_sbeta_gamma_k1);
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], c_sbeta_gamma_k2);
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], z_poly_xomega_roots);
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], alpha_poly);
    // -t_part[2]
    polynomial<Field> neg_one_poly = {-Field("1")};
    libfqfft::_polynomial_multiplication<Field>(t_part[2], t_part[2], neg_one_poly);
    
    // --- Computation of t_part[3]

    // z(x) - 1
    polynomial<Field> z_neg_one;
    libfqfft::_polynomial_addition<Field>(z_neg_one, z_poly, neg_one_poly);
    // (z(x)-1) * L_1(x)
    libfqfft::_polynomial_multiplication<Field>(t_part[3], z_neg_one, common_input.L_basis[0]);
    // (z(x)-1) * L_1(x) * alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[3], t_part[3], alpha_poly);
    // (z(x)-1) * L_1(x) * alpha * alpha
    libfqfft::_polynomial_multiplication<Field>(t_part[3], t_part[3], alpha_poly);


    // --- computation of t(x)

    // t(x) = (t[0] + t[1] + (-t[2]) + t[3]) / zh(x)
    polynomial<Field> t_poly_long{Field(0)};
    libfqfft::_polynomial_addition<Field>(t_poly_long, t_poly_long, t_part[0]);
#if 1 // DEBUG   
    libfqfft::_polynomial_addition<Field>(t_poly_long, t_poly_long, t_part[1]);
    libfqfft::_polynomial_addition<Field>(t_poly_long, t_poly_long, t_part[2]);
    libfqfft::_polynomial_addition<Field>(t_poly_long, t_poly_long, t_part[3]);
#endif // #if 0 // DEBUG   
    //    t(x) = t(x) / zh(x): A/B = (Q, R) st. A = (Q * B) + R.
    polynomial<Field> remainder;
    libfqfft::_polynomial_division(t_poly_long, remainder, t_poly_long, zh_poly);
#ifdef DEBUG
    printf("[%s:%d] t_poly_long\n", __FILE__, __LINE__);
    print_vector(t_poly_long);
#endif // #ifdef DEBUG
#ifdef DEBUG
    printf("[%s:%d] remainder\n", __FILE__, __LINE__);
    print_vector(remainder);
#endif // #ifdef DEBUG
    assert(libfqfft::_is_zero(remainder));

    // break t_poly_long into three parts: lo, mid, hi, each of degree 7
    // note: (common_input.num_gates+3) is the length of the CRS = (common_input.num_gates+2) powers of G1 + 1 power of G2
    std::vector<polynomial<Field>> t_poly(num_hgen);
    for (int i = 0; i < num_hgen; ++i) {
      typename std::vector<Field>::iterator begin = t_poly_long.begin()+(i*(common_input.num_gates+2));
      typename std::vector<Field>::iterator end = t_poly_long.begin()+(i*(common_input.num_gates+2))+(common_input.num_gates+2);
      std::vector<Field> tmp(begin, end);
      t_poly[i] = tmp;
    }
#ifdef DEBUG
    for (int i = 0; i < num_hgen; ++i) {
      printf("[%s:%d] t_poly[%d]\n", __FILE__, __LINE__, i);
      print_vector(t_poly[i]);
    }
#endif // #ifdef DEBUG
    // evaluate each part of t_poly in the secret input
    std::vector<libff::G1<ppT>> t_poly_at_secret_g1(num_hgen);
    for (int i = 0; i < num_hgen; ++i) {
      t_poly_at_secret_g1[i] = plonk_evaluate_poly_at_secret_G1<ppT>(srs.secret_powers_g1, t_poly[i]);
    }
#ifdef DEBUG
    printf("[%s:%d] Output from Round 3\n", __FILE__, __LINE__);
    // verify the output from Round 3 to the test vectors. test
    // vectors obtained from the Plonk Python reference implementation
    // (used for debug)
    for (int i = 0; i < num_hgen; ++i) {
      printf("[%s:%d] t_poly_at_secret_g1[%d]\n", __FILE__, __LINE__, i);
      t_poly_at_secret_g1[i].print();
      libff::G1<ppT> t_poly_at_secret_g1_i(t_poly_at_secret_g1[i]);
      t_poly_at_secret_g1_i.to_affine_coordinates();
      assert(t_poly_at_secret_g1_i.X == example.t_poly_at_secret_g1[i][0]);
      assert(t_poly_at_secret_g1_i.Y == example.t_poly_at_secret_g1[i][1]);
    }
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 3

    printf("[%s:%d] Prover Round 4...\n", __FILE__, __LINE__);
#if 1 // prover round 4
    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field zeta = example.zeta;
#ifdef DEBUG
    printf("[%s:%d] zeta\n", __FILE__, __LINE__);
    zeta.print();
#endif // #ifdef DEBUG
    Field a_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates + 2, W_polys_blinded[a], zeta);
    Field b_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates + 2, W_polys_blinded[b], zeta);
    Field c_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates + 2, W_polys_blinded[c], zeta);
    Field S_0_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates, common_input.S_polys[0], zeta);
    Field S_1_zeta = libfqfft::evaluate_polynomial<Field>(common_input.num_gates, common_input.S_polys[1], zeta);
    Field t_zeta = libfqfft::evaluate_polynomial<Field>(t_poly_long.size(), t_poly_long, zeta);
    Field z_poly_xomega_zeta = libfqfft::evaluate_polynomial<Field>(z_poly_xomega_roots.size(), z_poly_xomega_roots, zeta);
    
#ifdef DEBUG
    printf("[%s:%d] Output from Round 4\n", __FILE__, __LINE__);
    printf("a_zeta ");
    a_zeta.print();
    assert(a_zeta == example.a_zeta);
    printf("b_zeta ");
    b_zeta.print();
    assert(b_zeta == example.b_zeta);
    printf("c_zeta ");
    c_zeta.print();
    assert(c_zeta == example.c_zeta);
    printf("S_0_zeta ");
    S_0_zeta.print();
    assert(S_0_zeta == example.S_0_zeta);
    printf("S_1_zeta ");
    S_1_zeta.print();
    assert(S_1_zeta == example.S_1_zeta);
    printf("t_zeta ");
    t_zeta.print();
    assert(t_zeta == example.t_zeta);
    printf("z_poly_xomega_zeta ");
    z_poly_xomega_zeta.print();
    assert(z_poly_xomega_zeta == example.z_poly_xomega_zeta);
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 4
    
    printf("[%s:%d] Prover Round 5...\n", __FILE__, __LINE__);
#if 1 // prover round 5
    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field nu = example.nu;
    
    // compute linerisation polynomial r in five parts
    std::vector<polynomial<Field>> r_part(5);

    // --- Computation of r_part[0]
    
    // represent values as constant term polynomials in orderto use
    // the functions in the libfqfft library on polynomials
    polynomial<Field> a_zeta_poly{a_zeta}; 
    polynomial<Field> b_zeta_poly{b_zeta}; 
    polynomial<Field> c_zeta_poly{c_zeta}; 
    // a(z)b(z)q_M(x)
    polynomial<Field> abqM_zeta;
    libfqfft::_polynomial_multiplication<Field>(abqM_zeta, common_input.Q_polys[M], a_zeta_poly);
    libfqfft::_polynomial_multiplication<Field>(abqM_zeta, abqM_zeta, b_zeta_poly);    
    // a(z)q_L(x)
    polynomial<Field> aqL_zeta;
    libfqfft::_polynomial_multiplication<Field>(aqL_zeta, common_input.Q_polys[L], a_zeta_poly);
    // b(z)q_R(x)
    polynomial<Field> bqR_zeta;
    libfqfft::_polynomial_multiplication<Field>(bqR_zeta, common_input.Q_polys[R], b_zeta_poly);
    // c(z)q_O(x)
    polynomial<Field> cqO_zeta;
    libfqfft::_polynomial_multiplication<Field>(cqO_zeta, common_input.Q_polys[O], c_zeta_poly);
    // a(z)b(z)q_M(x) + a(z)q_L(x) + b(z)q_R(x) + c(z)q_O(x) + q_C(x)
    libfqfft::_polynomial_addition<Field>(r_part[0], poly_null, abqM_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], aqL_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], bqR_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], cqO_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], common_input.Q_polys[C]);

    // --- Computation of r_part[1]

    polynomial<Field> r1_const_poly
      {
       (a_zeta + (beta * zeta) + gamma) *
       (b_zeta + (beta * k1 * zeta) + gamma) *
       (c_zeta + (beta * k2 * zeta) + gamma) * alpha
      };
    libfqfft::_polynomial_multiplication<Field>(r_part[1], r1_const_poly, z_poly);

    // --- Computation of r_part[2]
    
    polynomial<Field> r2_const_poly
      {
       (a_zeta + (beta * S_0_zeta) + gamma) *
       (b_zeta + (beta * S_1_zeta) + gamma) *
       (alpha * beta * z_poly_xomega_zeta)
      };
    libfqfft::_polynomial_multiplication<Field>(r_part[2], r2_const_poly, common_input.S_polys[2]);
    // -r_part[2]
    libfqfft::_polynomial_multiplication<Field>(r_part[2], r_part[2], neg_one_poly);
    
    // --- Computation of r_part[3]
    
    //     r3 = accumulator_poly_ext3 * eval_poly(L_1, [zeta])[0] * alpha ** 2
    polynomial<Field> L_0_zeta_poly{libfqfft::evaluate_polynomial<Field>(common_input.L_basis[0].size(), common_input.L_basis[0], zeta)};
    polynomial<Field> alpha_power2_poly{libff::power(alpha, libff::bigint<1>(2))};
    libfqfft::_polynomial_multiplication<Field>(r_part[3], z_poly, L_0_zeta_poly);
    libfqfft::_polynomial_multiplication<Field>(r_part[3], r_part[3], alpha_power2_poly);

    // --- Computation of r_poly = (r0+r1-r2+r3)

    //
    // Note: here the reference Python implementation differs from the
    // paper where:
    //
    // r(x) = r(x) - zh(zeta) (t_lo(x) + zeta^n t_mid(x) + zeta^2n t_hi(x))
    //
    // In the reference implementation, the missing term is added in
    // the computation of the W_zeta(x) polynomial
    //
    polynomial<Field> r_poly; 
    libfqfft::_polynomial_addition<Field>(r_poly, poly_null, r_part[0]);
    libfqfft::_polynomial_addition<Field>(r_poly, r_poly, r_part[1]);
    libfqfft::_polynomial_addition<Field>(r_poly, r_poly, r_part[2]);
    libfqfft::_polynomial_addition<Field>(r_poly, r_poly, r_part[3]);
    
#if DEBUG    
    //    printf("abqM_zeta\n");
    //    print_vector(abqM_zeta);
    printf("[%s:%d] r_part[0]\n", __FILE__, __LINE__);
    print_vector(r_part[0]);
    printf("[%s:%d] r_part[1]\n", __FILE__, __LINE__);
    print_vector(r_part[1]);
    printf("[%s:%d] r_part[2]\n", __FILE__, __LINE__);
    print_vector(r_part[2]);
    printf("[%s:%d] r_part[3]\n", __FILE__, __LINE__);
    print_vector(r_part[3]);
    printf("[%s:%d] r_poly\n", __FILE__, __LINE__);
    print_vector(r_poly);
#endif // #if DEBUG

    // Evaluate the r-polynomial at zeta. Note: in the reference
    // implementation, r_zeta is added to the pi-SNARK proof. In the
    // paper this is omitted, which makes the proof shorter at the
    // epxense of a slightly heavier computation on the verifier's
    // side 
    Field r_zeta = libfqfft::evaluate_polynomial<Field>(r_poly.size(), r_poly, zeta);
#ifdef DEBUG
    printf("r_zeta ");
    r_zeta.print();
#endif // #ifdef DEBUG    
    assert(r_zeta == example.r_zeta);

    // W_zeta polynomial is of degree 6 in the random element nu and
    // hence has 7 terms
    std::vector<polynomial<Field>> W_zeta_part(7);

    // --- compute W_zeta_part[0]
    
    // t_lo(x)
    polynomial<Field> t_lo{t_poly[lo]};
    // t_mid(x) * zeta^(n+2)
    polynomial<Field> t_mid_zeta_n;
    polynomial<Field> zeta_powern_poly{libff::power(zeta, libff::bigint<1>(common_input.num_gates+2))};
    libfqfft::_polynomial_multiplication<Field>(t_mid_zeta_n, t_poly[mid], zeta_powern_poly);
    // t_hi(x) * zeta^(2(n+1))
    polynomial<Field> t_hi_zeta_2n;
    polynomial<Field> zeta_power2n_poly{libff::power(zeta, libff::bigint<1>(2*(common_input.num_gates+2)))};
    libfqfft::_polynomial_multiplication<Field>(t_hi_zeta_2n, t_poly[hi], zeta_power2n_poly);
    // -t_zeta as constant term polynomial
    polynomial<Field> t_zeta_poly{-t_zeta};
    // t_lo(x) + (t_mid(x) * zeta^n) + (t_hi(x) * zeta^2n) + t_zeta_poly
    libfqfft::_polynomial_addition<Field>(W_zeta_part[0], poly_null, t_lo);
    libfqfft::_polynomial_addition<Field>(W_zeta_part[0], W_zeta_part[0], t_mid_zeta_n);
    libfqfft::_polynomial_addition<Field>(W_zeta_part[0], W_zeta_part[0], t_hi_zeta_2n);
    libfqfft::_polynomial_addition<Field>(W_zeta_part[0], W_zeta_part[0], t_zeta_poly);    
    
    // --- compute W_zeta_part[1]

    // -r_zeta as constant term polynomial
    polynomial<Field> r_zeta_poly{-r_zeta};
    // r(x) - r_zeta
    polynomial<Field> r_sub_rzeta;
    libfqfft::_polynomial_addition<Field>(r_sub_rzeta, r_poly, r_zeta_poly);
    // (r(x) - r_zeta) * nu
    polynomial<Field> nu_poly{nu};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[1], r_sub_rzeta, nu_poly);
    
    // --- compute W_zeta_part[2]
    
    // -a_zeta as constant term polynomial
    polynomial<Field> a_zeta_poly_neg;
    libfqfft::_polynomial_multiplication<Field>(a_zeta_poly_neg, a_zeta_poly, neg_one_poly);
    // a(x) - a_zeta
    polynomial<Field> a_sub_azeta;
    libfqfft::_polynomial_addition<Field>(a_sub_azeta, W_polys_blinded[a], a_zeta_poly_neg);
    // (a(x) - a_zeta) * nu^2
    Field nu2 = libff::power(nu, libff::bigint<1>(2));
    polynomial<Field> nu2_poly{nu2};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[2], a_sub_azeta, nu2_poly);
    
    // -b_zeta as constant term polynomial
    polynomial<Field> b_zeta_poly_neg;
    libfqfft::_polynomial_multiplication<Field>(b_zeta_poly_neg, b_zeta_poly, neg_one_poly);
    // (b(x) - b_zeta)
    polynomial<Field> b_sub_bzeta;
    libfqfft::_polynomial_addition<Field>(b_sub_bzeta, W_polys_blinded[b], b_zeta_poly_neg);
    // (b(x) - b_zeta) * nu^3
    Field nu3 = libff::power(nu, libff::bigint<1>(3));
    polynomial<Field> nu3_poly{nu3};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[3], b_sub_bzeta, nu3_poly);

    // -c_zeta as constant term polynomial
    polynomial<Field> c_zeta_poly_neg;
    libfqfft::_polynomial_multiplication<Field>(c_zeta_poly_neg, c_zeta_poly, neg_one_poly);
    // (c(x) - c_zeta)
    polynomial<Field> c_sub_czeta;
    libfqfft::_polynomial_addition<Field>(c_sub_czeta, W_polys_blinded[c], c_zeta_poly_neg);
    // (c(x) - c_zeta) * nu^4
    Field nu4 = libff::power(nu, libff::bigint<1>(4));
    polynomial<Field> nu4_poly{nu4};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[4], c_sub_czeta, nu4_poly);

    // -S_0_zeta as constant term polynomial
    polynomial<Field> S_0_zeta_poly_neg{-S_0_zeta};
    //    libfqfft::_polynomial_multiplication<Field>(S_0_zeta_poly_neg, S_0_zeta_poly, neg_one_poly);
    // (S0(x) - S_0_zeta)
    polynomial<Field> S0_sub_szeta;
    libfqfft::_polynomial_addition<Field>(S0_sub_szeta, common_input.S_polys[0], S_0_zeta_poly_neg);
    // (S0(x) - S_0_zeta) * nu^5
    Field nu5 = libff::power(nu, libff::bigint<1>(5));
    polynomial<Field> nu5_poly{nu5};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[5], S0_sub_szeta, nu5_poly);

    // -S_1_zeta as constant term polynomial
    polynomial<Field> S_1_zeta_poly_neg{-S_1_zeta};
    //    libfqfft::_polynomial_multiplication<Field>(S_1_zeta_poly_neg, S_1_zeta_poly, neg_one_poly);
    // (S1(x) - S_1_zeta)
    polynomial<Field> S1_sub_szeta;
    libfqfft::_polynomial_addition<Field>(S1_sub_szeta, common_input.S_polys[1], S_1_zeta_poly_neg);
    // (S1(x) - S_1_zeta) * nu^6
    Field nu6 = libff::power(nu, libff::bigint<1>(6));
    polynomial<Field> nu6_poly{nu6};
    libfqfft::_polynomial_multiplication<Field>(W_zeta_part[6], S1_sub_szeta, nu6_poly);

    // compute full zeta polynomial W_zeta = \sum W_zeta_part[i]
    int nzeta = 7;
    polynomial<Field> W_zeta(poly_null);
    for (int i = 0; i < nzeta; ++i) {
      libfqfft::_polynomial_addition<Field>(W_zeta, W_zeta, W_zeta_part[i]);
    }

    // compute 1/(X-zeta) * W_zeta
    polynomial<Field> x_sub_zeta_poly{-zeta, Field(1)};
    libfqfft::_polynomial_division(W_zeta, remainder, W_zeta, x_sub_zeta_poly);
#ifdef DEBUG
    printf("W_zeta\n");
    print_vector(W_zeta);
#endif // #ifdef DEBUG
    assert(libfqfft::_is_zero(remainder));

    //    polynomial<Field> z_poly;
    
    // Compute opening proof:
    // W_zeta_omega = z(X) - z(zeta*common_input.omega_roots) / X - (zeta*common_input.omega_roots)
    polynomial<Field> W_zeta_omega{poly_null};
    
    // -z(zeta*common_input.omega_roots)
    polynomial<Field> z_poly_xomega_zeta_neg{-z_poly_xomega_zeta};
    // z(X) - z(zeta*common_input.omega_roots) 
    libfqfft::_polynomial_addition<Field>(W_zeta_omega, z_poly, z_poly_xomega_zeta_neg);    
    // -zeta*common_input.omega_roots; common_input.omega_roots[base][1] = common_input.omega_roots_base
    polynomial<Field> x_sub_zeta_omega_roots{-(zeta*common_input.omega_roots[base][1]), Field(1)};
    
    // z(X) - z(zeta*common_input.omega_roots) / X - (zeta*common_input.omega_roots)
    libfqfft::_polynomial_division(W_zeta_omega, remainder, W_zeta_omega, x_sub_zeta_omega_roots);
    assert(libfqfft::_is_zero(remainder));

#ifdef DEBUG
    printf("W_zeta_part[0]\n");
    print_vector(W_zeta_part[0]);
    printf("W_zeta_part[1]\n");
    print_vector(W_zeta_part[1]);
    printf("W_zeta_part[2]\n");
    print_vector(W_zeta_part[2]);
    printf("W_zeta_part[3]\n");
    print_vector(W_zeta_part[3]);
    printf("W_zeta_part[4]\n");
    print_vector(W_zeta_part[4]);
    printf("W_zeta_part[5]\n");
    print_vector(W_zeta_part[5]);
    printf("W_zeta_part[6]\n");
    print_vector(W_zeta_part[6]);
    printf("W_zeta\n");
    print_vector(W_zeta);
    printf("W_zeta_omega\n");
    print_vector(W_zeta_omega);
#endif // #ifdef DEBUG
      
    assert(W_zeta == example.W_zeta);
    assert(W_zeta_omega == example.W_zeta_omega);
    
    // Evaluate polynomials W_zeta and W_zeta_omega at the seceret
    // input
    libff::G1<ppT> W_zeta_at_secret =
      plonk_evaluate_poly_at_secret_G1<ppT>(srs.secret_powers_g1, W_zeta);    
    libff::G1<ppT> W_zeta_omega_at_secret =
      plonk_evaluate_poly_at_secret_G1<ppT>(srs.secret_powers_g1, W_zeta_omega);
    
#ifdef DEBUG
    printf("[%s:%d] Outputs from Prover round 5\n", __FILE__, __LINE__);
    
    printf("r_zeta ");
    r_zeta.print();
    
    printf("[%s:%d] W_zeta_at_secret \n", __FILE__, __LINE__);
    W_zeta_at_secret.print();
    libff::G1<ppT> W_zeta_at_secret_aff(W_zeta_at_secret);
    W_zeta_at_secret_aff.to_affine_coordinates();
    assert(W_zeta_at_secret_aff.X == example.W_zeta_at_secret[0]);
    assert(W_zeta_at_secret_aff.Y == example.W_zeta_at_secret[1]);
    
    printf("[%s:%d] W_zeta_omega_at_secret \n", __FILE__, __LINE__);
    W_zeta_omega_at_secret.print();
    libff::G1<ppT> W_zeta_omega_at_secret_aff(W_zeta_omega_at_secret);
    W_zeta_omega_at_secret_aff.to_affine_coordinates();
    assert(W_zeta_omega_at_secret_aff.X == example.W_zeta_omega_at_secret[0]);
    assert(W_zeta_omega_at_secret_aff.Y == example.W_zeta_omega_at_secret[1]);
#endif // #ifdef DEBUG

    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    Field u = example.u;
#endif // #if 1 // prover round 5

    // --- VERIFIER ---

    //
    // SNARK proof
    //
    // Pi ([a]_1, [b]_1, [c]_1, [z]_1,
    //     [t_lo]_1, [t_mi]_1, [t_hi]_1,
    //     \bar{a}, \bar{b}, \bar{c},
    //     \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
    //     [W_zeta]_1, [W_{zeta common_input.omega_roots}]_1
    //     r_zeta (*))
    //
    // (*) Note: in the reference Python implementation, r_zeta (the
    // evaluation of the r(X) polynomial at zeta from Prover round 5)
    // is added to the pi-SNARK proof. In the paper this is omitted,
    // which seems to make the proof shorter by 1 element at the
    // epxense of a slightly heavier computation on the verifier's
    // side. Here we follow the reference implementation to make sure
    // we match the test values. TODO: once all test vectors are
    // verified, we may remove r_zeta from the proof to be fully
    // compliant with the paper.
    //
    // Mapping code-to-paper quantities
    //
    // - W_polys_blinded_at_secret_g1[a, b, c]: [a]_1, [b]_1, [c]_1 (from Round 1)
    // - z_poly_at_secret_g1: [z]_1 (from Round 2)
    // - t_poly_at_secret_g1[lo, mi, hi]: [t_lo]_1, [t_mi]_1, [t_hi]_1 (from Round 3)
    // - a_zeta, b_zeta, c_zeta, S_0_zeta, S_1_zeta, z_poly_xomega_zeta: \bar{a}, \bar{b}, \bar{c}, \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w} (from Round 4)
    // - W_zeta_at_secret, W_zeta_omega_at_secret: [W_zeta]_1, [W_{zeta common_input.omega_roots}]_1 (from Round 5)
    //
    // Verifier preprocessed input
    //
    //  - Q_polys_at_secret_g1[0..3]: [q_M]_1, [q_L]_1, [q_R]_1, [q_O]_1
    //  - S_polys_at_secret_g1[0..2]: [S_sigma1]_1, [S_sigma2]_1, [S_sigma3]_1
    //  - srs.secret_powers_g2[1]: [secret]_2 = secret^1 * G2
    //
    // Public input polynomial
    //
    // common_input.PI_poly: w_i, 0\le{i}<l<n
    //    
    printf("[%s:%d] Verifier preparation...\n", __FILE__, __LINE__);
    
    // Verifier precomputation:
    std::vector<libff::G1<ppT>> Q_polys_at_secret_g1(common_input.Q_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(srs.secret_powers_g1, common_input.Q_polys, Q_polys_at_secret_g1);
    std::vector<libff::G1<ppT>> S_polys_at_secret_g1(common_input.S_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(srs.secret_powers_g1, common_input.S_polys, S_polys_at_secret_g1);

    // secret * G2
    
#ifdef DEBUG
    for (int i = 0; i < (int)common_input.Q_polys.size(); ++i) {
      printf("common_input.Q_polys_at_secret_G1[%d] \n", i);
      Q_polys_at_secret_g1[i].print();      
      libff::G1<ppT> Q_poly_at_secret_g1_i(Q_polys_at_secret_g1[i]);
      Q_poly_at_secret_g1_i.to_affine_coordinates();      
      assert(Q_poly_at_secret_g1_i.X == example.Q_polys_at_secret_g1[i][0]);
      assert(Q_poly_at_secret_g1_i.Y == example.Q_polys_at_secret_g1[i][1]);
    }
    for (int i = 0; i < (int)common_input.S_polys.size(); ++i) {
      printf("S_polys_at_secret_G1[%d] \n", i);
      S_polys_at_secret_g1[i].print();
      libff::G1<ppT> S_poly_at_secret_g1_i(S_polys_at_secret_g1[i]);
      S_poly_at_secret_g1_i.to_affine_coordinates();      
      assert(S_poly_at_secret_g1_i.X == example.S_polys_at_secret_g1[i][0]);
      assert(S_poly_at_secret_g1_i.Y == example.S_polys_at_secret_g1[i][1]);
    }
#endif // #ifdef DEBUG

    Field alpha_power2 = libff::power(alpha, libff::bigint<1>(2));
    
    bool b_valid = false;

    // --- Verifier Step 1: validate that elements belong to group G1
    b_valid = check_curve_equation<libff::G1<ppT>>(W_polys_blinded_at_secret_g1[a]);
    assert(b_valid);
    b_valid = check_curve_equation<libff::G1<ppT>>(W_polys_blinded_at_secret_g1[b]);
    assert(b_valid);
    b_valid = check_curve_equation<libff::G1<ppT>>(W_polys_blinded_at_secret_g1[c]);
    assert(b_valid);
    b_valid = check_curve_equation<libff::G1<ppT>>(z_poly_at_secret_g1);
    assert(b_valid);
    b_valid = check_curve_equation<libff::G1<ppT>>(t_poly_at_secret_g1[lo]);
    assert(b_valid);
    b_valid = check_curve_equation<libff::G1<ppT>>(t_poly_at_secret_g1[mid]);
    assert(b_valid);
    b_valid = check_curve_equation<libff::G1<ppT>>(t_poly_at_secret_g1[hi]);
    assert(b_valid);
    b_valid = check_curve_equation<libff::G1<ppT>>(W_zeta_at_secret);
    assert(b_valid);
    b_valid = check_curve_equation<libff::G1<ppT>>(W_zeta_omega_at_secret);
    assert(b_valid);

    // --- Verifier Step 2: validate that elements belong to scalar field Fr
    b_valid = check_field_element<Field>(a_zeta);    
    assert(b_valid);
    b_valid = check_field_element<Field>(b_zeta);    
    assert(b_valid);
    b_valid = check_field_element<Field>(c_zeta);    
    assert(b_valid);
    b_valid = check_field_element<Field>(S_0_zeta);    
    assert(b_valid);
    b_valid = check_field_element<Field>(S_1_zeta);    
    assert(b_valid);
    b_valid = check_field_element<Field>(z_poly_xomega_zeta);    
    assert(b_valid);
    b_valid = check_field_element<Field>(r_zeta);    
    assert(b_valid);
    
    // --- Verifier Step 3: validate that the public input belongs to scalar field Fr
    assert(common_input.PI_poly.size() <= common_input.num_gates);
    for (int i = 0; i < (int)common_input.PI_poly.size(); ++i) {
      b_valid = check_field_element<Field>(common_input.PI_poly[i]);    
      assert(b_valid);      
    }    

    // --- Verifier Step 4: compute challenges hashed transcript as in
    //     prover description, from the common inputs, public input,
    //     and elements of pi-SNARK .  TODO: fixed to the test vectors for now
    beta = example.beta;
    gamma = example.gamma;
    alpha = example.alpha;
    zeta = example.zeta;
    nu = example.nu;
    u = example.u;

    // --- Verifier Step 5: compute zero polynomial evaluation
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain = libfqfft::get_evaluation_domain<Field>(common_input.num_gates);
    Field zh_zeta = domain->compute_vanishing_polynomial(zeta);
    printf("[%s:%d] zh_zeta ", __FILE__, __LINE__);
    zh_zeta.print();
    assert(zh_zeta == example.zh_zeta);

    // --- Verifier Step 6: Compute Lagrange polynomial evaluation L1(zeta)
    // Note: the paper counts the L-polynomials from 1; we count from 0
    Field L_0_zeta = libfqfft::evaluate_polynomial<Field>(common_input.L_basis[0].size(), common_input.L_basis[0], zeta);   
#ifdef DEBUG
    printf("L_0_zeta ");
    L_0_zeta.print();
#endif // #ifdef DEBUG    
    assert(L_0_zeta == example.L_0_zeta);

    // --- Verifier Step 7: compute public input polynomial evaluation PI(zeta)
    Field PI_zeta = libfqfft::evaluate_polynomial<Field>(common_input.PI_poly.size(), common_input.PI_poly, zeta);   
#ifdef DEBUG
    printf("PI_zeta ");
    PI_zeta.print();
#endif // #ifdef DEBUG    
    assert(PI_zeta == example.PI_zeta);
    
    // --- Verifier Step 8: compute quotient polynomial evaluation
    // r'(zeta) = r(zeta) - r0, where r0 is a constant term Note:
    // follows the Python reference implementation, which slightly
    // deviates from the paper due to the presence of the r_zeta term
    // in the proof (not present in the paper.  In particular, the
    // reference code computes and uses r'(zeta) in step 8, while the
    // paper uses r0. In addition, the reference code divides r'(zeta)
    // by the vanishing polynomial at zeta zh_zeta, while the paper
    // does not do that (see also Step 9).

    // compute polynomial r'(zeta) = r(zeta) - r_0
    std::vector<Field> r_prime_parts(5);
    r_prime_parts[0] = r_zeta + PI_zeta;
    r_prime_parts[1] = (a_zeta + (beta * S_0_zeta) + gamma);
    r_prime_parts[2] = (b_zeta + (beta * S_1_zeta) + gamma);
    r_prime_parts[3] = (c_zeta + gamma) * z_poly_xomega_zeta * alpha;
    r_prime_parts[4] = (L_0_zeta * alpha_power2);
    Field r_prime_zeta = (r_prime_parts[0] - (r_prime_parts[1] * r_prime_parts[2] * r_prime_parts[3]) - r_prime_parts[4]) * zh_zeta.inverse();    
#ifdef DEBUG
    printf("r_prime_parts[%d] ", 0);
    r_prime_parts[0].print();
    printf("r_prime_parts[%d] ", 1);
    r_prime_parts[1].print();
    printf("r_prime_parts[%d] ", 2);
    r_prime_parts[2].print();
    printf("r_prime_parts[%d] ", 3);
    r_prime_parts[3].print();
    printf("r_prime_parts[%d] ", 4);
    r_prime_parts[4].print();
    printf("r_prime_zeta     ");
    r_prime_zeta.print();
#endif // #ifdef DEBUG    
    assert(r_prime_zeta == example.r_prime_zeta);

    // --- Verifier Step 9: compute first part of batched polynomial commitment [D]_1

    // Note: the reference implemention differs from the paper -- it
    // does not add the following term to D1, but to F1 (Step 10):
    // -Zh(zeta)([t_lo]_1 + zeta^n [t_mid]_1 + zeta^2n
    // [t_hi]_1). Instead ([t_lo]_1 + zeta^n [t_mid]_1 + zeta^2n
    // [t_hi]_1) is added to F1 in Step 10 and the multiplication by
    // Zh(zeta) is accounted for by dividing by Zh(zeta) of r'(zeta)
    // in Step 8.

    // D1 is computed in 3 parts
    std::vector<libff::G1<ppT>> D1_part(3);

    // compute D1_part[0]:    
    // (a_bar b_bar [q_M]_1 + a_bar [q_L]_1 + b_bar [q_R]_1 + c_bar [q_O]_1 + [q_C]_1) nu
    // Note: the paper omits the final multiplication by nu    
    std::vector<libff::G1<ppT>> curve_points_9
      {
       Q_polys_at_secret_g1[M],
       Q_polys_at_secret_g1[L],
       Q_polys_at_secret_g1[R],
       Q_polys_at_secret_g1[O],
       Q_polys_at_secret_g1[C]
      };
    std::vector<libff::Fr<ppT>> scalar_elements_9
     {
       a_zeta * b_zeta * nu,
       a_zeta * nu,
       b_zeta * nu,
       c_zeta * nu,
       nu
      };
    D1_part[0] = plonk_multi_exp_G1<ppT>(curve_points_9, scalar_elements_9);

    // compute D1_part[1]:
    // ((a_bar + beta zeta + gamma)(b_bar + beta k1 zeta + gamma)(c_bar + beta k2 zeta + gamma) alpha + L1(zeta) alpha^2 + u) [z]_1
    Field D1_part1_scalar = 
      (a_zeta + (beta * zeta) + gamma) *
      (b_zeta + (beta * k1 * zeta) + gamma) *
      (c_zeta + (beta * k2 * zeta) + gamma) * alpha * nu
      + L_0_zeta * alpha_power2 * nu + u;
    D1_part[1] = plonk_exp_G1<ppT>(z_poly_at_secret_g1, D1_part1_scalar);

    // compute D1_part[2]:
    // (a_bar + beta s_sigma1_bar + gamma)(b_bar + beta s_sigma2_bar + gamma)alpha beta z_common_input.omega_roots_bar [s_sigma3]_1
    Field D1_part2_scalar = 
      ((a_zeta + (beta * S_0_zeta) + gamma) *
       (b_zeta + (beta * S_1_zeta) + gamma) *
       alpha * beta * z_poly_xomega_zeta * nu) * Field(-1);
    D1_part[2] = plonk_exp_G1<ppT>(S_polys_at_secret_g1[2], D1_part2_scalar);

    // Compute D1 = D1_part[0] + D1_part[1] + D1_part[2]
    libff::G1<ppT> D1 = D1_part[0] + D1_part[1] + D1_part[2];
    
#ifdef DEBUG    
    printf("[%s:%d] D1_part[%d]\n", __FILE__, __LINE__, 0);
    D1_part[0].print();
    printf("[%s:%d] D1_part[%d]\n", __FILE__, __LINE__, 1);
    D1_part[1].print();
    printf("[%s:%d] D1_part[%d]\n", __FILE__, __LINE__, 2);
    D1_part[2].print();
    printf("[%s:%d] D1\n", __FILE__, __LINE__);
    D1.print();
    libff::G1<ppT> D1_aff(D1);
    D1_aff.to_affine_coordinates();      
    assert(D1_aff.X == example.D1[0]);
    assert(D1_aff.Y == example.D1[1]);
#endif // #ifdef DEBUG    
    
    // --- Verifier Step 10: compute full batched polynomial
    // commitment [F]_1 [F]_1 = [D]_1 + v [a]_1 + v^2 [b]_1 + v^3
    // [c]_1 + v^4 [s_sigma_1]_1 + v^5 [s_sigma_2]_1 Note: to [F]_1
    // the erefernce code also adds the term ([t_lo]_1 + zeta^n
    // [t_mid]_1 + zeta^2n [t_hi]_1) which is addedto [D]_1 in the
    // paper (see commenst to Steps 8,9)

    Field zeta_power_n = libff::power(zeta, libff::bigint<1>(common_input.num_gates+2));
    Field zeta_power_2n = libff::power(zeta, libff::bigint<1>(2*(common_input.num_gates+2)));
    std::vector<Field> nu_power(7);
    for (size_t i = 0; i < nu_power.size(); ++i) {
      nu_power[i] = libff::power(nu, libff::bigint<1>(i));
    }    
    std::vector<libff::G1<ppT>> curve_points_10
      {
       t_poly_at_secret_g1[lo], // nu^0
       t_poly_at_secret_g1[mid], // nu^0
       t_poly_at_secret_g1[hi], // nu^0
       D1, // nu^1
       W_polys_blinded_at_secret_g1[a], // nu^2
       W_polys_blinded_at_secret_g1[b], // nu^3
       W_polys_blinded_at_secret_g1[c], // nu^4
       S_polys_at_secret_g1[0], // nu^5
       S_polys_at_secret_g1[1] // nu^6
      };
    std::vector<libff::Fr<ppT>> scalar_elements_10
      {
       Field(1), 
       zeta_power_n,  // zeta^(n+2), 
       zeta_power_2n, // zeta^(2*(n+2)), 
       Field(1), 
       nu_power[2], // nu^2, 
       nu_power[3], // nu^3, 
       nu_power[4], // nu^4, 
       nu_power[5], // nu^5, 
       nu_power[6]  // nu^6
      };
    libff::G1<ppT> F1 = plonk_multi_exp_G1<ppT>(curve_points_10, scalar_elements_10);

#ifdef DEBUG    
    printf("[%s:%d] F1\n", __FILE__, __LINE__);
    F1.print();
    libff::G1<ppT> F1_aff(F1);
    F1_aff.to_affine_coordinates();      
    assert(F1_aff.X == example.F1[0]);
    assert(F1_aff.Y == example.F1[1]);
#endif // #ifdef DEBUG    

    // --- Verifier Step 11: compute group-encoded batch evaluation [E]_1

    std::vector<libff::G1<ppT>> curve_points_11
      {
       libff::G1<ppT>::one()
      };
    std::vector<libff::Fr<ppT>> scalar_elements_11
      {
       r_prime_zeta +
       nu_power[1] * r_zeta + // v^1
       nu_power[2] * a_zeta + // v^2
       nu_power[3] * b_zeta + // v^3
       nu_power[4] * c_zeta + // v^4
       nu_power[5] * S_0_zeta + // v^5
       nu_power[6] * S_1_zeta + // v^6
       u * z_poly_xomega_zeta
      };
    libff::G1<ppT> E1 = plonk_multi_exp_G1<ppT>(curve_points_11, scalar_elements_11);
    
#ifdef DEBUG    
    printf("[%s:%d] G\n", __FILE__, __LINE__);
    libff::G1<ppT>::one().print();
    printf("[%s:%d] E1\n", __FILE__, __LINE__);
    E1.print();
    libff::G1<ppT> E1_aff(E1);
    E1_aff.to_affine_coordinates();      
    assert(E1_aff.X == example.E1[0]);
    assert(E1_aff.Y == example.E1[1]);
#endif // #ifdef DEBUG    

    // --- Verifier Step 12: batch validate all evaluations

    // Check the following equality
    //
    // e( [W_zeta]_1 + u [W_{zeta common_input.omega_roots}]_1, [x]_2 ) *
    // e( -zeta [W_zeta ]_1 - u zeta common_input.omega_roots [W_{zeta common_input.omega_roots}]_1 - [F]_1 + [E]_1, [1]_2 )
    // = Field(1)
    //
    // Denoted as: 
    // e(first_lhs, second_lhs) * e(first_rhs, second_rhs) = 1
    //

    // add random element (noise) to the opening polynomials to check
    // that the pairing fails
#if 0
    libff::G1<ppT> noise = libff::G1<ppT>::random_element();
#else    
    libff::G1<ppT> noise = libff::G1<ppT>::zero();
#endif    
    
    std::vector<libff::G1<ppT>> curve_points_lhs
      {
       W_zeta_at_secret + noise,
       W_zeta_omega_at_secret
      };
    std::vector<libff::Fr<ppT>> scalar_elements_lhs
      {
       Field(1),
       u
      };
    libff::G1<ppT> pairing_first_lhs = plonk_multi_exp_G1<ppT>(curve_points_lhs, scalar_elements_lhs);
    libff::G2<ppT> pairing_second_lhs = srs.secret_powers_g2[1];
    
    std::vector<libff::G1<ppT>> curve_points_rhs
      {
       W_zeta_at_secret,
       W_zeta_omega_at_secret,
       F1,
       E1
      };
    std::vector<libff::Fr<ppT>> scalar_elements_rhs
      {
       // Warning! raise to the power of -1 to check e() * e()^-1 = 1
       Field(-1) * zeta,
       Field(-1) * u * zeta * common_input.omega_roots[base][1],
       Field(-1) * Field(1),
       Field(-1) * Field(-1)
      };
   
    libff::G1<ppT> pairing_first_rhs = plonk_multi_exp_G1<ppT>(curve_points_rhs, scalar_elements_rhs);
    libff::G2<ppT> pairing_second_rhs = srs.secret_powers_g2[0];

#ifdef DEBUG    
    printf("[%s:%d] pairing_first_lhs\n", __FILE__, __LINE__);
    pairing_first_lhs.print();
    libff::G1<ppT> pairing_first_lhs_aff(pairing_first_lhs);
    pairing_first_lhs_aff.to_affine_coordinates();      
    assert(pairing_first_lhs_aff.X == example.pairing_first_lhs[0]);
    assert(pairing_first_lhs_aff.Y == example.pairing_first_lhs[1]);

    printf("[%s:%d] pairing_first_rhs\n", __FILE__, __LINE__);
    pairing_first_rhs.print();
    libff::G1<ppT> pairing_first_rhs_aff(pairing_first_rhs);
    pairing_first_rhs_aff.to_affine_coordinates();
    assert(pairing_first_rhs_aff.X == example.pairing_first_rhs[0]);
    assert(pairing_first_rhs_aff.Y == example.pairing_first_rhs[1]);
#endif // #ifdef DEBUG    

    const libff::G1_precomp<ppT> _A = ppT::precompute_G1(pairing_first_lhs);
    const libff::G2_precomp<ppT> _B = ppT::precompute_G2(pairing_second_lhs);
    const libff::G1_precomp<ppT> _C = ppT::precompute_G1(pairing_first_rhs);
    const libff::G2_precomp<ppT> _D = ppT::precompute_G2(pairing_second_rhs);
    const libff::Fqk<ppT> miller_result =
        ppT::double_miller_loop(_A, _B, _C, _D);
    const libff::GT<ppT> result = ppT::final_exponentiation(miller_result);
    assert(result == libff::GT<ppT>::one());
    
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
