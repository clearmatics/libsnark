/** @file
*****************************************************************************

Declaration of interfaces for ppzkSNARK proof system Plonk.

This includes:
- class srs
- class for proving key
- class for verification key
- class for key pair (proving key & verification key)
- class for proof
- generator algorithm / setup
- prover algorithm
- verifier algorithm 

The implementation instantiates the protocol of PlonK \[GWC19],

References:

\[GWC19]: 
"Plonk: Permutations over lagrange-bases for oecumenical noninteractive arguments of knowledge",
Ariel Gabizon, Zachary J. Williamson, and Oana Ciobotaru,
Cryptology ePrint Archive, Report 2019/953, 2019,
<https://eprint.iacr.org/2019/953>

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_HPP_
#define PLONK_PPZKSNARK_HPP_

#include <libff/algebra/curves/public_params.hpp>
#include <memory>

#define DEBUG 1

// number of generators for H, Hk1, Hk2
const size_t NUM_HGEN = 3;

namespace libsnark
{

  enum W_polys_id{a, b, c};
  enum Q_polys_id{L, R, M, O, C};
  enum t_polys_id{lo, mid, hi};
  enum omega_id{base, base_k1, base_k2};

template<typename FieldT> using polynomial = std::vector<FieldT>;

/************************ Common Preprocessed Input ***********************/

/**
 * Plonk common preprocessed input
 */
template<typename ppT> class common_preprocessed_input
{
public:
  using Field = libff::Fr<ppT>;
  // number of gates / constraints
  size_t num_gates;
  // number of selector polynomials (q-polynomials) (= 5 in the
  // vanilla Plonk proposal [GWC19])
  size_t num_qpolys;
  // Lagrange basis
  std::vector<polynomial<Field>> L_basis;
  // Public input polynomial
  polynomial<Field> PI_poly;
  // Circuit selector polynomials (Q-polynomials)
  std::vector<polynomial<Field>> Q_polys;
  // Permutation polynomials S_sigma_1, S_sigma_2, S_sigma_2 (see
  // [GWC19], Sect. 8.1)
  std::vector<polynomial<Field>> S_polys;
  // omega[0] are the n roots of unity, omega[1] are omega[0]*k1,
  // omega[2] are omega[0]*k2
  std::vector<std::vector<Field>> omega_roots;
  // H_gen contains the generators of H, k1 H and k2 H in one place
  // ie. omega, omega_k1 and omega_k2
  std::vector<Field> H_gen;
  // H_gen permuted according to the wire permutation
  std::vector<Field> H_gen_permute;
  // constants for H, k1 H, k2 H
  libff::Fr<ppT> k1;
  libff::Fr<ppT> k2;

  common_preprocessed_input(){};
  
};  

/********************************** SRS ***********************************/

/// The srs generated by the setup step. See kzg10.hpp .
template<typename ppT> class srs
{
public:
  /// Array of powers of secret \alpha, encoded in G1:
  /// [1]_1, [\alpha]_1, [\alpha^2]_1, ..., [\alpha^{n+2}]_1
  std::vector<libff::G1<ppT>> secret_powers_g1;
  
  /// Array of powers of secret \alpha, encoded in G2:
  /// [1]_2, [\alpha]_2
  std::vector<libff::G2<ppT>> secret_powers_g2;

  srs(){};
  
  srs(std::vector<libff::G1<ppT>> &&secret_powers_g1,
      std::vector<libff::G2<ppT>> &&secret_powers_g2);

  // for debug only
  void derive_from_secret(
			  const libff::Fr<ppT> secret,
			  size_t num_gates
			  );
};

  
/******************************** Proving key ********************************/
/**
 * A proving key for Plonk
 */
template<typename ppT> class plonk_proving_key
{
public:
  /// Array of powers of secret \alpha, encoded in G1:
  /// [1]_1, [\alpha]_1, [\alpha^2]_1, ..., [\alpha^{n+2}]_1
  std::vector<libff::G1<ppT>> secret_powers_g1;

  plonk_proving_key(){};
  plonk_proving_key(
		    std::vector<libff::G1<ppT>>&& secret_powers_g1)
    :secret_powers_g1(std::move(secret_powers_g1)){};
};

/******************************* Verification key ****************************/

/**
 * A verification key for Plonk
 */
template<typename ppT> class plonk_verification_key
{
public:
  /// Array of powers of secret \alpha, encoded in G2:
  /// [1]_2, [\alpha]_2
  std::vector<libff::G2<ppT>> secret_powers_g2;
  
  plonk_verification_key();
  plonk_verification_key(
		    std::vector<libff::G2<ppT>>&& secret_powers_g2)
    :secret_powers_g2(std::move(secret_powers_g2)){};
};

/********************************** Key pair *********************************/

/**
 * A key pair for Plonk, which consists of a proving key and a
 * verification key.
 */
template<typename ppT> class plonk_keypair
{
public:
    plonk_proving_key<ppT> pk;
    plonk_verification_key<ppT> vk;

    plonk_keypair() = default;
    plonk_keypair(const plonk_keypair<ppT> &other) =
        default;
    plonk_keypair(
        plonk_proving_key<ppT> &&pk,
        plonk_verification_key<ppT> &&vk)
        : pk(std::move(pk)), vk(std::move(vk))
    {
    }

    plonk_keypair(plonk_keypair<ppT> &&other) = default;
};
  
/*********************************** Proof ***********************************/

/**
 * A proof for the Plonk GG-ppzkSNARK.
 */
template<typename ppT> class plonk_proof
{
public:
  using Field = libff::Fr<ppT>;
  
  // Plonk proof Pi
  //
  // Pi ([a]_1, [b]_1, [c]_1, [z]_1,
  //     [t_lo]_1, [t_mi]_1, [t_hi]_1,
  //     \bar{a}, \bar{b}, \bar{c},
  //     \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
  //     [W_zeta]_1, [W_{zeta common_input.omega_roots}]_1
  //     r_zeta (*))
  
  // [a]_1, [b]_1, [c]_1
  std::vector<libff::G1<ppT>> W_polys_blinded_at_secret_g1;
  // [z]_1,
  libff::G1<ppT> z_poly_at_secret_g1;
  // [t_lo]_1, [t_mi]_1, [t_hi]_1
  std::vector<libff::G1<ppT>> t_poly_at_secret_g1;
  // \bar{a}, \bar{b}, \bar{c},
  Field a_zeta;
  Field b_zeta;
  Field c_zeta;
  // \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
  Field S_0_zeta;
  Field S_1_zeta;
  // \bar{z_w}
  Field z_poly_xomega_zeta;
  // [W_zeta]_1, [W_{zeta common_input.omega_roots}]_1  
  libff::G1<ppT> W_zeta_at_secret;
  libff::G1<ppT> W_zeta_omega_at_secret;
  // r_zeta (*)
  Field r_zeta;

  plonk_proof() {};
  plonk_proof(
	      std::vector<libff::G1<ppT>>& W_polys_blinded_at_secret_g1,
	      libff::G1<ppT>& z_poly_at_secret_g1,
	      std::vector<libff::G1<ppT>>& t_poly_at_secret_g1,
	      Field& a_zeta,
	      Field& b_zeta,
	      Field& c_zeta,
	      Field& S_0_zeta,
	      Field& S_1_zeta,
	      Field& z_poly_xomega_zeta,
	      libff::G1<ppT>& W_zeta_at_secret,
	      libff::G1<ppT>& W_zeta_omega_at_secret,
	      Field& r_zeta
	      ) : 
    W_polys_blinded_at_secret_g1(W_polys_blinded_at_secret_g1),
    z_poly_at_secret_g1(z_poly_at_secret_g1),
    t_poly_at_secret_g1(t_poly_at_secret_g1),
    a_zeta(a_zeta),
    b_zeta(b_zeta),
    c_zeta(c_zeta),
    S_0_zeta(S_0_zeta),
    S_1_zeta(S_1_zeta),
    z_poly_xomega_zeta(z_poly_xomega_zeta),
    W_zeta_at_secret(W_zeta_at_secret),
    W_zeta_omega_at_secret(W_zeta_omega_at_secret),
    r_zeta(r_zeta)
  {
  }
};

/**
 * Plonk prover. Computes object of class plonk_proof.
 */
template<typename ppT> class plonk_prover
{
private:  
  using Field = libff::Fr<ppT>;

  // round 0 (initialization)
  std::vector<libff::Fr<ppT>> zh_poly;
  //  libff::Fr<ppT> k1;
  //  libff::Fr<ppT> k2;
  polynomial<libff::Fr<ppT>> null_poly;
  polynomial<libff::Fr<ppT>> neg_one_poly; 
 
  // round 1
  std::vector<libff::Fr<ppT>> blind_scalars;
  std::vector<polynomial<libff::Fr<ppT>>> W_polys;
  std::vector<std::vector<libff::Fr<ppT>>> W_polys_blinded;
  std::vector<libff::G1<ppT>> W_polys_blinded_at_secret_g1;
  
  // round 2
  libff::Fr<ppT> beta;
  libff::Fr<ppT> gamma;
  polynomial<libff::Fr<ppT>> z_poly;
  libff::G1<ppT> z_poly_at_secret_g1;

  // round 3
  libff::Fr<ppT> alpha;
  std::vector<libff::Fr<ppT>> z_poly_xomega;
  std::vector<polynomial<libff::Fr<ppT>>> t_poly;
  polynomial<libff::Fr<ppT>> t_poly_long;
  std::vector<libff::G1<ppT>> t_poly_at_secret_g1;

  // round 4
  libff::Fr<ppT> zeta;
  libff::Fr<ppT> a_zeta;
  libff::Fr<ppT> b_zeta;
  libff::Fr<ppT> c_zeta;
  libff::Fr<ppT> S_0_zeta;
  libff::Fr<ppT> S_1_zeta;  
  libff::Fr<ppT> z_poly_xomega_zeta;
  libff::Fr<ppT> t_zeta;

  // round 5
  libff::Fr<ppT> nu;
  libff::Fr<ppT> r_zeta;
  libff::G1<ppT> W_zeta_at_secret;
  libff::G1<ppT> W_zeta_omega_at_secret;
  libff::Fr<ppT> u;

public:  
  // constructors: initialize round 0 variables
  plonk_prover() {};
  
  plonk_prover(
	       const common_preprocessed_input<ppT> common_input
	       );

  void compute_witness_polys(
			     const std::vector<libff::Fr<ppT>> witness,
			     const common_preprocessed_input<ppT> common_input
			     );    
  void round_one(
		 const std::vector<libff::Fr<ppT>> witness,
		 const common_preprocessed_input<ppT> common_input,
		 const srs<ppT> srs
		 );
  void round_two(
		 const std::vector<libff::Fr<ppT>> witness,
		 const common_preprocessed_input<ppT> common_input,
		 const srs<ppT> srs
		 );
  void round_three(
		   const common_preprocessed_input<ppT> common_input,
		   const srs<ppT> srs
		   );
  void round_four(
		  const common_preprocessed_input<ppT> common_input
		  );  
  void round_five(
		  const common_preprocessed_input<ppT> common_input,
		  const srs<ppT> srs
		  );
  plonk_proof<ppT> compute_proof(
				 const srs<ppT> srs,
				 const common_preprocessed_input<ppT> common_input
				 );
};

/***************************** Verifier *******************************/

/**
 * Plonk verifier. Verifies object of class plonk_proof.
 */
template<typename ppT> class plonk_verifier
{
  using Field = libff::Fr<ppT>;
  
private:
  // verifier preprocessed input
  std::vector<libff::G1<ppT>> Q_polys_at_secret_g1;
  std::vector<libff::G1<ppT>> S_polys_at_secret_g1;
  // secret * G2
  // evaluation of vanishing polynomial Zh at x=zeta i.e. Zh(zeta)
  Field zh_zeta;
  // Lagrange polynomial evaluation of polynomial L1 at x=zeta
  // i.e. L1(zeta)
  Field L_0_zeta;
  // Public input polynomial PI evaluated at x=zeta .e. PI(zeta)
  Field PI_zeta;
  // compute quotient polynomial evaluation r'(zeta) = r(zeta) - r0,
  // where r0 is a constant term Note:
  Field r_prime_zeta;
  // first part of batched polynomial commitment [D]_1
  libff::G1<ppT> D1;
  // full batched polynomial commitment [F]_1 = [D]_1 + v [a]_1 + v^2
  // [b]_1 + v^3 [c]_1 + v^4 [s_sigma_1]_1 + v^5 [s_sigma_2]_1
  libff::G1<ppT> F1;
  // group-encoded batch evaluation [E]_1
  libff::G1<ppT> E1;

  // challenges hashed transcript
  libff::Fr<ppT> beta;
  libff::Fr<ppT> gamma;
  libff::Fr<ppT> alpha;
  libff::Fr<ppT> zeta;
  libff::Fr<ppT> nu;
  libff::Fr<ppT> u;

public:
  plonk_verifier() {};
  
  void preprocessed_input(
			  const srs<ppT> srs,
			  const common_preprocessed_input<ppT> common_input
			  );

  void step_one(const plonk_proof<ppT> proof);
  void step_two(const plonk_proof<ppT> proof);
  void step_three(
		  const common_preprocessed_input<ppT> common_input
		  );
  void step_four();
  void step_five(
		 const common_preprocessed_input<ppT> common_input
		 );
  void step_six(
		const common_preprocessed_input<ppT> common_input
		);
  void step_seven(
		  const common_preprocessed_input<ppT> common_input
		  );
  void step_eight(
		  const plonk_proof<ppT> proof
		  );
  void step_nine(
		 const plonk_proof<ppT> proof,
		 const common_preprocessed_input<ppT> common_input
		 );
  void step_ten(
		const plonk_proof<ppT> proof,
		const common_preprocessed_input<ppT> common_input
		);
  void step_eleven(
		   const plonk_proof<ppT> proof
		   );
  bool step_twelve(
		   const plonk_proof<ppT> proof,
		   const srs<ppT> srs,
		   const common_preprocessed_input<ppT> common_input
		   );  
  bool verify_proof(
		    const plonk_proof<ppT> proof,
		    const srs<ppT> srs,
		    const common_preprocessed_input<ppT> common_input
		    );
};

/***************************** Main algorithms *******************************/


  
} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/plonk.tcc"

#endif // PLONK_PPZKSNARK_HPP_
