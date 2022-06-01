/** @file
*****************************************************************************

Declaration of Prover interfaces for ppzkSNARK proof system Plonk.

This includes:
- class for proof
- class for prover

References:

\[GWC19]:
"Plonk: Permutations over lagrange-bases for oecumenical noninteractive
arguments of knowledge", Ariel Gabizon, Zachary J. Williamson, and Oana
Ciobotaru, Cryptology ePrint Archive, Report 2019/953, 2019,
<https://eprint.iacr.org/2019/953>

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_PROVER_HPP_
#define PLONK_PPZKSNARK_PROVER_HPP_

namespace libsnark
{
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

    plonk_proof(){};
    plonk_proof(
        std::vector<libff::G1<ppT>> &W_polys_blinded_at_secret_g1,
        libff::G1<ppT> &z_poly_at_secret_g1,
        std::vector<libff::G1<ppT>> &t_poly_at_secret_g1,
        Field &a_zeta,
        Field &b_zeta,
        Field &c_zeta,
        Field &S_0_zeta,
        Field &S_1_zeta,
        Field &z_poly_xomega_zeta,
        libff::G1<ppT> &W_zeta_at_secret,
        libff::G1<ppT> &W_zeta_omega_at_secret,
        Field &r_zeta)
        : W_polys_blinded_at_secret_g1(W_polys_blinded_at_secret_g1)
        , z_poly_at_secret_g1(z_poly_at_secret_g1)
        , t_poly_at_secret_g1(t_poly_at_secret_g1)
        , a_zeta(a_zeta)
        , b_zeta(b_zeta)
        , c_zeta(c_zeta)
        , S_0_zeta(S_0_zeta)
        , S_1_zeta(S_1_zeta)
        , z_poly_xomega_zeta(z_poly_xomega_zeta)
        , W_zeta_at_secret(W_zeta_at_secret)
        , W_zeta_omega_at_secret(W_zeta_omega_at_secret)
        , r_zeta(r_zeta)
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
    plonk_prover(){};

    plonk_prover(const common_preprocessed_input<ppT> common_input);

    void compute_witness_polys(
        const std::vector<libff::Fr<ppT>> witness,
        const common_preprocessed_input<ppT> common_input);
    void round_one(
        const std::vector<libff::Fr<ppT>> witness,
        const common_preprocessed_input<ppT> common_input,
        const srs<ppT> srs);
    void round_two(
        const std::vector<libff::Fr<ppT>> witness,
        const common_preprocessed_input<ppT> common_input,
        const srs<ppT> srs);
    void round_three(
        const common_preprocessed_input<ppT> common_input, const srs<ppT> srs);
    void round_four(const common_preprocessed_input<ppT> common_input);
    void round_five(
        const common_preprocessed_input<ppT> common_input, const srs<ppT> srs);
    plonk_proof<ppT> compute_proof(
        const srs<ppT> srs, const common_preprocessed_input<ppT> common_input);

};

// Prover round 0 output
template<typename ppT>
struct round_zero_out_t {
  std::vector<libff::Fr<ppT>> zh_poly;
  polynomial<libff::Fr<ppT>> null_poly;
  polynomial<libff::Fr<ppT>> neg_one_poly;
};
  
// Prover round 1 output
template<typename ppT>
struct round_one_out_t {
  std::vector<libff::Fr<ppT>> blind_scalars;
  std::vector<polynomial<libff::Fr<ppT>>> W_polys;
  std::vector<std::vector<libff::Fr<ppT>>> W_polys_blinded;
  std::vector<libff::G1<ppT>> W_polys_blinded_at_secret_g1;
};

// Prover round 2 output
template<typename ppT>
struct round_two_out_t {
  libff::Fr<ppT> beta;
  libff::Fr<ppT> gamma;
  polynomial<libff::Fr<ppT>> z_poly;
  libff::G1<ppT> z_poly_at_secret_g1;
};

// Prover round 3 output
template<typename ppT>
struct round_three_out_t {
  libff::Fr<ppT> alpha;
  std::vector<libff::Fr<ppT>> z_poly_xomega;
  std::vector<polynomial<libff::Fr<ppT>>> t_poly;
  polynomial<libff::Fr<ppT>> t_poly_long;
  std::vector<libff::G1<ppT>> t_poly_at_secret_g1;
};
  
// Prover round 4 output
template<typename ppT>
struct round_four_out_t {
  libff::Fr<ppT> zeta;
  libff::Fr<ppT> a_zeta;
  libff::Fr<ppT> b_zeta;
  libff::Fr<ppT> c_zeta;
  libff::Fr<ppT> S_0_zeta;
  libff::Fr<ppT> S_1_zeta;
  libff::Fr<ppT> z_poly_xomega_zeta;
  libff::Fr<ppT> t_zeta;
};

/**
 * Plonk prover. Computes object of class plonk_proof.
 */
template<typename ppT> class plonk_prover_new
{
  using Field = libff::Fr<ppT>;
  
private:

  //  const round_zero_out_t<ppT> round_zero_out;
  
public:
  
  // constructors: initialize round 0 variables
  plonk_prover_new(){};

  // --- new ---

  static round_zero_out_t<ppT>
  round_zero(
	     const common_preprocessed_input<ppT> common_input);
  
  static round_one_out_t<ppT>
  round_one(
	    const round_zero_out_t<ppT> round_zero_out,
	    const std::vector<libff::Fr<ppT>> witness,
	    const common_preprocessed_input<ppT> common_input,
	    const srs<ppT> srs);
  
  static round_two_out_t<ppT>
  round_two(
	    const round_zero_out_t<ppT> round_zero_out,
	    const round_one_out_t<ppT> round_one_out,
	    const std::vector<libff::Fr<ppT>> witness,
	    const common_preprocessed_input<ppT> common_input,
	    const srs<ppT> srs);
  
  static round_three_out_t<ppT>
  round_three(
	      const round_zero_out_t<ppT> round_zero_out,
	      const round_one_out_t<ppT> round_one_out,
	      const round_two_out_t<ppT> round_two_out,  
	      const common_preprocessed_input<ppT> common_input,
	      const srs<ppT> srs);
  
  static round_four_out_t<ppT>
  round_four(
	     const round_one_out_t<ppT> round_one_out,
	     const round_three_out_t<ppT> round_three_out,
	     const common_preprocessed_input<ppT> common_input);
  
  static plonk_proof<ppT>
  compute_proof(
		const srs<ppT> srs,
		const common_preprocessed_input<ppT> common_input);
  
};
  
} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/prover.tcc"

#endif // PLONK_PPZKSNARK_PROVER_HPP_
