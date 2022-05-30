/** @file
*****************************************************************************

Implementation of Verifier interfaces for a ppzkSNARK for Plonk.

See verifier.hpp .

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_VERIFIER_TCC_
#define PLONK_PPZKSNARK_VERIFIER_TCC_

namespace libsnark
{
  //
  // SNARK proof
  //
  // Pi ([a]_1, [b]_1, [c]_1, [z]_1,
  //     [t_lo]_1, [t_mi]_1, [t_hi]_1,
  //     \bar{a}, \bar{b}, \bar{c},
  //     \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
  //     [W_zeta]_1, [W_{zeta omega}]_1
  //     r_zeta (*))
  //
  // (*) Note: in the reference Python implementation, r_zeta (the
  // evaluation of the linearlization polynomial r(X) at zeta from
  // Prover round 5) is added to the pi-SNARK proof. In the paper
  // this is omitted, which seems to make the proof shorter by 1
  // element at the epxense of a slightly heavier computation on the
  // verifier's side. Here we follow the reference implementation to
  // make sure we match the test values. TODO: once all test vectors
  // are verified, we may remove r_zeta from the proof to be fully
  // compliant with the paper.
  //
  // Mapping code-to-paper quantities
  //
  // - W_polys_blinded_at_secret_g1[a, b, c]: [a]_1, [b]_1, [c]_1 (from Round 1)
  // - z_poly_at_secret_g1: [z]_1 (from Round 2)
  // - t_poly_at_secret_g1[lo, mi, hi]: [t_lo]_1, [t_mi]_1, [t_hi]_1 (from Round 3)
  // - a_zeta, b_zeta, c_zeta, S_0_zeta, S_1_zeta, z_poly_xomega_zeta: \bar{a}, \bar{b}, \bar{c}, \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w} (from Round 4)
  // - W_zeta_at_secret, W_zeta_omega_at_secret: [W_zeta]_1, [W_{zeta omega}]_1 (from Round 5)
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
  
  // Verifier precomputation
  template<typename ppT>
  void plonk_verifier<ppT>::preprocessed_input(
					       const srs<ppT> srs,
					       const common_preprocessed_input<ppT> common_input
					       )
  {
    this->Q_polys_at_secret_g1.resize(common_input.Q_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(srs.secret_powers_g1, common_input.Q_polys, Q_polys_at_secret_g1);
    
    this->S_polys_at_secret_g1.resize(common_input.S_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(srs.secret_powers_g1, common_input.S_polys, S_polys_at_secret_g1);

  }
  
  // Verifier Step 1: validate that elements belong to group G1
  //
  // WARNING! This validation MUST be done by the caller. Empty
  // function here for consistency with the description in [GWC19]
  template<typename ppT>
  void plonk_verifier<ppT>::step_one(
				     plonk_proof<ppT> proof
				     )
  {
  }
  
  // Verifier Step 2: validate that elements belong to scalar field Fr
  //
  // WARNING! This validation MUST be done by the caller. Empty
  // function here for consistency with the description in [GWC19]
  template<typename ppT>
  void plonk_verifier<ppT>::step_two(
				     plonk_proof<ppT> proof
				     )
  {
  }
  
  // Verifier Step 3: validate that the public input belongs to scalar field Fr
  //
  // WARNING! This validation MUST be done by the caller. Empty
  // function here for consistency with the description in [GWC19]
  template<typename ppT>
  void plonk_verifier<ppT>::step_three(
				       const common_preprocessed_input<ppT> common_input
				       )
  {
  }
  
  // Verifier Step 4: compute challenges hashed transcript as in
  // prover description, from the common inputs, public input, and
  // elements of pi-SNARK .  TODO: fixed to the test vectors for now
  template<typename ppT>
  void plonk_verifier<ppT>::step_four()
  {
    // load challenges from example for debug
#ifdef DEBUG
    plonk_example<ppT> example;
#endif // #ifdef DEBUG
    this->beta = example.beta;
    this->gamma = example.gamma;
    this->alpha = example.alpha;
    this->zeta = example.zeta;
    this->nu = example.nu;
    this->u = example.u;
  }
  
  // Verifier Step 5: compute zero polynomial evaluation
  template<typename ppT>
  void plonk_verifier<ppT>::step_five(
				      const common_preprocessed_input<ppT> common_input
				      )
  {
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain = libfqfft::get_evaluation_domain<Field>(common_input.num_gates);
    this->zh_zeta = domain->compute_vanishing_polynomial(this->zeta);
  }
  
  // Verifier Step 6: Compute Lagrange polynomial evaluation L1(zeta)
  // Note: the paper counts the L-polynomials from 1; we count from 0
  template<typename ppT>
  void plonk_verifier<ppT>::step_six(
				     const common_preprocessed_input<ppT> common_input
				     )
  {
    this->L_0_zeta = libfqfft::evaluate_polynomial<Field>(common_input.L_basis[0].size(), common_input.L_basis[0], zeta);   
  }
  
  // Verifier Step 7: compute public input polynomial evaluation
  // PI(zeta)
  template<typename ppT>
  void plonk_verifier<ppT>::step_seven(
				       const common_preprocessed_input<ppT> common_input
				       )
  {
    this->PI_zeta = libfqfft::evaluate_polynomial<Field>(common_input.PI_poly.size(), common_input.PI_poly, zeta);   
  }
  
  // Verifier Step 8: compute quotient polynomial evaluation r'(zeta)
  // = r(zeta) - r0, where r0 is a constant term Note: follows the
  // Python reference implementation, which slightly deviates from the
  // paper due to the presence of the r_zeta term in the proof (not
  // present in the paper.  In particular, the reference code computes
  // and uses r'(zeta) in step 8, while the paper uses r0. In
  // addition, the reference code divides r'(zeta) by the vanishing
  // polynomial at zeta zh_zeta, while the paper does not do that (see
  // also Step 9).
  template<typename ppT>
  void plonk_verifier<ppT>::step_eight(
				       const plonk_proof<ppT> proof
				       )
  {
    Field alpha_power2 = libff::power(alpha, libff::bigint<1>(2));
    
    // compute polynomial r'(zeta) = r(zeta) - r_0
    std::vector<Field> r_prime_parts(5);
    r_prime_parts[0] = proof.r_zeta + this->PI_zeta;
    r_prime_parts[1] = (proof.a_zeta + (this->beta * proof.S_0_zeta) + this->gamma);
    r_prime_parts[2] = (proof.b_zeta + (this->beta * proof.S_1_zeta) + this->gamma);
    r_prime_parts[3] = (proof.c_zeta + this->gamma) * proof.z_poly_xomega_zeta * this->alpha;
    r_prime_parts[4] = (this->L_0_zeta * alpha_power2);
    this->r_prime_zeta = (r_prime_parts[0] - (r_prime_parts[1] * r_prime_parts[2] * r_prime_parts[3]) - r_prime_parts[4]) * this->zh_zeta.inverse();    
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
  }
  
  // Verifier Step 9: compute first part of batched polynomial
  // commitment [D]_1 Note: the reference implemention differs from
  // the paper -- it does not add the following term to D1, but to F1
  // (Step 10): -Zh(zeta)([t_lo]_1 + zeta^n [t_mid]_1 + zeta^2n
  // [t_hi]_1). Instead ([t_lo]_1 + zeta^n [t_mid]_1 + zeta^2n
  // [t_hi]_1) is added to F1 in Step 10 and the multiplication by
  // Zh(zeta) is accounted for by dividing by Zh(zeta) of r'(zeta) in
  // Step 8.
  template<typename ppT>
  void plonk_verifier<ppT>::step_nine(
				      const plonk_proof<ppT> proof,
				      const common_preprocessed_input<ppT> common_input
				      )
  {    
    // D1 is computed in 3 parts
    std::vector<libff::G1<ppT>> D1_part(3);

    Field alpha_power2 = libff::power(alpha, libff::bigint<1>(2));
    
    // compute D1_part[0]:    
    // (a_bar b_bar [q_M]_1 + a_bar [q_L]_1 + b_bar [q_R]_1 + c_bar [q_O]_1 + [q_C]_1) nu
    // Note: the paper omits the final multiplication by nu    
    std::vector<libff::G1<ppT>> curve_points_9
      {
       this->Q_polys_at_secret_g1[M],
       this->Q_polys_at_secret_g1[L],
       this->Q_polys_at_secret_g1[R],
       this->Q_polys_at_secret_g1[O],
       this->Q_polys_at_secret_g1[C]
      };
    std::vector<libff::Fr<ppT>> scalar_elements_9
      {
       proof.a_zeta * proof.b_zeta * this->nu,
       proof.a_zeta * this->nu,
       proof.b_zeta * this->nu,
       proof.c_zeta * this->nu,
       this->nu
      };
    D1_part[0] = plonk_multi_exp_G1<ppT>(curve_points_9, scalar_elements_9);

    // compute D1_part[1]:
    // ((a_bar + beta zeta + gamma)(b_bar + beta k1 zeta + gamma)(c_bar + beta k2 zeta + gamma) alpha + L1(zeta) alpha^2 + u) [z]_1
    Field D1_part1_scalar = 
      (proof.a_zeta + (this->beta * this->zeta) + this->gamma) *
      (proof.b_zeta + (this->beta * common_input.k1 * this->zeta) + this->gamma) *
      (proof.c_zeta + (this->beta * common_input.k2 * this->zeta) + this->gamma) * this->alpha * this->nu
      + this->L_0_zeta * alpha_power2 * this->nu + this->u;
    D1_part[1] = plonk_exp_G1<ppT>(proof.z_poly_at_secret_g1, D1_part1_scalar);

    // compute D1_part[2]:
    // (a_bar + beta s_sigma1_bar + gamma)(b_bar + beta s_sigma2_bar + gamma)alpha beta z_common_input.omega_roots_bar [s_sigma3]_1
    Field D1_part2_scalar = 
      ((proof.a_zeta + (this->beta * proof.S_0_zeta) + this->gamma) *
       (proof.b_zeta + (this->beta * proof.S_1_zeta) + this->gamma) *
       this->alpha * this->beta * proof.z_poly_xomega_zeta * this->nu) * Field(-1);
    D1_part[2] = plonk_exp_G1<ppT>(this->S_polys_at_secret_g1[2], D1_part2_scalar);

    // Compute D1 = D1_part[0] + D1_part[1] + D1_part[2]
    this->D1 = D1_part[0] + D1_part[1] + D1_part[2];
#ifdef DEBUG    
    printf("[%s:%d] D1_part[%d]\n", __FILE__, __LINE__, 0);
    D1_part[0].print();
    printf("[%s:%d] D1_part[%d]\n", __FILE__, __LINE__, 1);
    D1_part[1].print();
    printf("[%s:%d] D1_part[%d]\n", __FILE__, __LINE__, 2);
    D1_part[2].print();
    printf("[%s:%d] D1\n", __FILE__, __LINE__);
#endif // #ifdef DEBUG    
  }
  
  // Verifier Step 10: compute full batched polynomial commitment
  // [F]_1 = [D]_1 + v [a]_1 + v^2 [b]_1 + v^3 [c]_1 + v^4
  // [s_sigma_1]_1 + v^5 [s_sigma_2]_1 Note: to [F]_1 the erefernce
  // code also adds the term ([t_lo]_1 + zeta^n [t_mid]_1 + zeta^2n
  // [t_hi]_1) which is addedto [D]_1 in the paper (see commenst to
  // Steps 8,9)
  template<typename ppT>
  void plonk_verifier<ppT>::step_ten(
				     const plonk_proof<ppT> proof,
				     const common_preprocessed_input<ppT> common_input
				     )
  {
    Field zeta_power_n = libff::power(zeta, libff::bigint<1>(common_input.num_gates+2));
    Field zeta_power_2n = libff::power(zeta, libff::bigint<1>(2*(common_input.num_gates+2)));
    std::vector<Field> nu_power(7);
    for (size_t i = 0; i < nu_power.size(); ++i) {
      nu_power[i] = libff::power(this->nu, libff::bigint<1>(i));
    }    
    std::vector<libff::G1<ppT>> curve_points_10
      {
	proof.t_poly_at_secret_g1[lo], // nu^0
	  proof.t_poly_at_secret_g1[mid], // nu^0
	  proof.t_poly_at_secret_g1[hi], // nu^0
	  this->D1, // nu^1
	  proof.W_polys_blinded_at_secret_g1[a], // nu^2
	  proof.W_polys_blinded_at_secret_g1[b], // nu^3
	  proof.W_polys_blinded_at_secret_g1[c], // nu^4
	  this->S_polys_at_secret_g1[0], // nu^5
	  this->S_polys_at_secret_g1[1] // nu^6
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
    this->F1 = plonk_multi_exp_G1<ppT>(curve_points_10, scalar_elements_10);

  }
  
  // Verifier Step 11: compute group-encoded batch evaluation [E]_1
  template<typename ppT>
  void plonk_verifier<ppT>::step_eleven(
					const plonk_proof<ppT> proof
					)
  {
    std::vector<Field> nu_power(7);
    for (size_t i = 0; i < nu_power.size(); ++i) {
      nu_power[i] = libff::power(this->nu, libff::bigint<1>(i));
    }    
    std::vector<libff::G1<ppT>> curve_points_11
      {
	libff::G1<ppT>::one()
	  };
    std::vector<libff::Fr<ppT>> scalar_elements_11
      {
       this->r_prime_zeta +
       nu_power[1] * proof.r_zeta + // v^1
       nu_power[2] * proof.a_zeta + // v^2
       nu_power[3] * proof.b_zeta + // v^3
       nu_power[4] * proof.c_zeta + // v^4
       nu_power[5] * proof.S_0_zeta + // v^5
       nu_power[6] * proof.S_1_zeta + // v^6
       this->u * proof.z_poly_xomega_zeta
      };
    this->E1 = plonk_multi_exp_G1<ppT>(curve_points_11, scalar_elements_11);
  }
  
  // Verifier Step 12: batch validate all evaluations
  // 
  // Checks the following equality
  //
  // e( [W_zeta]_1 + u [W_{zeta common_input.omega_roots}]_1, [x]_2 ) *
  // e( -zeta [W_zeta ]_1 - u zeta common_input.omega_roots [W_{zeta common_input.omega_roots}]_1 - [F]_1 + [E]_1, [1]_2 )
  // = Field(1)
  //
  // Denoted as: 
  // e(first_lhs, second_lhs) * e(first_rhs, second_rhs) = 1
  //
  template<typename ppT>
  bool plonk_verifier<ppT>::step_twelve(
					const plonk_proof<ppT> proof,
					const srs<ppT> srs,
					const common_preprocessed_input<ppT> common_input
					)
  {
    // load test vectors for debug
#ifdef DEBUG
    plonk_example<ppT> example;
#endif // #ifdef DEBUG
    
    // add random element (noise) to the opening polynomials to check
    // that the pairing fails
    //    libff::G1<ppT> noise = libff::G1<ppT>::random_element();
    libff::G1<ppT> noise = libff::G1<ppT>::zero();
    
    std::vector<libff::G1<ppT>> curve_points_lhs
      {
       proof.W_zeta_at_secret + noise,
       proof.W_zeta_omega_at_secret
      };
    std::vector<libff::Fr<ppT>> scalar_elements_lhs
      {
       Field(1),
       this->u
      };
    libff::G1<ppT> pairing_first_lhs = plonk_multi_exp_G1<ppT>(curve_points_lhs, scalar_elements_lhs);
    libff::G2<ppT> pairing_second_lhs = srs.secret_powers_g2[1];
    
    std::vector<libff::G1<ppT>> curve_points_rhs
      {
       proof.W_zeta_at_secret,
       proof.W_zeta_omega_at_secret,
       this->F1,
       this->E1
      };
    std::vector<libff::Fr<ppT>> scalar_elements_rhs
      {
       // Warning! raise to the power of -1 to check e() * e()^-1 = 1
       Field(-1) * this->zeta,
       Field(-1) * this->u * this->zeta * common_input.omega_roots[base][1],
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
    bool b_accept = (result == libff::GT<ppT>::one());
    assert(b_accept);
    return b_accept;
  }
  
  // WARNING! The first three steps (as given in [GWC19] -- see
  // below) MUST be executed by the caller:
  //
  // Verifier Step 1: validate that elements belong to group G1
  // Verifier Step 2: validate that elements belong to scalar field Fr
  // Verifier Step 3: validate that the public input belongs to scalar field Fr        
  template<typename ppT>
  bool plonk_verifier<ppT>::verify_proof(
					 const plonk_proof<ppT> proof,
					 const srs<ppT> srs,
					 const common_preprocessed_input<ppT> common_input
					 )
  {

    // load test vector values form example for debug
#ifdef DEBUG
    plonk_example<ppT> example;
#endif // #ifdef DEBUG
    
    // compute verifier preprocessed input
    this->preprocessed_input(srs, common_input);
#ifdef DEBUG
    for (int i = 0; i < (int)common_input.Q_polys.size(); ++i) {
      printf("common_input.Q_polys_at_secret_G1[%d] \n", i);
      this->Q_polys_at_secret_g1[i].print();      
      libff::G1<ppT> Q_poly_at_secret_g1_i(this->Q_polys_at_secret_g1[i]);
      Q_poly_at_secret_g1_i.to_affine_coordinates();      
      assert(Q_poly_at_secret_g1_i.X == example.Q_polys_at_secret_g1[i][0]);
      assert(Q_poly_at_secret_g1_i.Y == example.Q_polys_at_secret_g1[i][1]);
    }
    for (int i = 0; i < (int)common_input.S_polys.size(); ++i) {
      printf("S_polys_at_secret_G1[%d] \n", i);
      this->S_polys_at_secret_g1[i].print();
      libff::G1<ppT> S_poly_at_secret_g1_i(this->S_polys_at_secret_g1[i]);
      S_poly_at_secret_g1_i.to_affine_coordinates();      
      assert(S_poly_at_secret_g1_i.X == example.S_polys_at_secret_g1[i][0]);
      assert(S_poly_at_secret_g1_i.Y == example.S_polys_at_secret_g1[i][1]);
    }
#endif // #ifdef DEBUG
    // Verifier Step 1: validate that elements belong to group G1
    // Executed by the caller
    
    // Verifier Step 2: validate that elements belong to scalar field Fr
    // Executed by the caller
    
    // Verifier Step 3: validate that the public input belongs to scalar field Fr
    // Executed by the caller

    // Verifier Step 4: compute challenges hashed transcript as in
    // prover description, from the common inputs, public input, and
    // elements of pi-SNARK (fixed to the test vectors for now)
    this->step_four();

    // Verifier Step 5: compute zero polynomial evaluation
    this->step_five(common_input);
#ifdef DEBUG
    printf("[%s:%d] zh_zeta ", __FILE__, __LINE__);
    this->zh_zeta.print();
    assert(this->zh_zeta == example.zh_zeta);
#endif // #ifdef DEBUG
    
    // Verifier Step 6: Compute Lagrange polynomial evaluation L1(zeta)
    // Note: the paper counts the L-polynomials from 1; we count from 0
    this->step_six(common_input);
#ifdef DEBUG
    printf("L_0_zeta ");
    this->L_0_zeta.print();
    assert(this->L_0_zeta == example.L_0_zeta);
#endif // #ifdef DEBUG    

    // Verifier Step 7: compute public input polynomial evaluation PI(zeta)
    this->step_seven(common_input);
#ifdef DEBUG
    printf("PI_zeta ");
    this->PI_zeta.print();
    assert(this->PI_zeta == example.PI_zeta);
#endif // #ifdef DEBUG
    
    // Verifier Step 8: compute quotient polynomial evaluation
    // r'(zeta) = r(zeta) - r0, where r0 is a constant term
    this->step_eight(proof);
#ifdef DEBUG
    assert(this->r_prime_zeta == example.r_prime_zeta);
#endif // #ifdef DEBUG    

    // Verifier Step 9: compute first part of batched polynomial
    // commitment [D]_1 
    this->step_nine(proof, common_input);
#ifdef DEBUG    
    this->D1.print();
    libff::G1<ppT> D1_aff(this->D1);
    D1_aff.to_affine_coordinates();      
    assert(D1_aff.X == example.D1[0]);
    assert(D1_aff.Y == example.D1[1]);
#endif // #ifdef DEBUG       

    // Verifier Step 10: compute full batched polynomial commitment
    // [F]_1
    this->step_ten(proof, common_input);
#ifdef DEBUG    
    printf("[%s:%d] F1\n", __FILE__, __LINE__);
    this->F1.print();
    libff::G1<ppT> F1_aff(this->F1);
    F1_aff.to_affine_coordinates();      
    assert(F1_aff.X == example.F1[0]);
    assert(F1_aff.Y == example.F1[1]);
#endif // #ifdef DEBUG    

    // Verifier Step 11: compute group-encoded batch evaluation [E]_1
    this->step_eleven(proof);
#ifdef DEBUG    
    printf("[%s:%d] E1\n", __FILE__, __LINE__);
    this->E1.print();
    libff::G1<ppT> E1_aff(this->E1);
    E1_aff.to_affine_coordinates();      
    assert(E1_aff.X == example.E1[0]);
    assert(E1_aff.Y == example.E1[1]);
#endif // #ifdef DEBUG

    // Verifier Step 12: batch validate all evaluations (check
    // pairing)
    bool b_accept = this->step_twelve(proof, srs, common_input);
    return b_accept;
  }
} // namespace libsnark

#endif // PLONK_PPZKSNARK_VERIFIER_TCC_
