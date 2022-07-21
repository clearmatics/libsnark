/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_TCC_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_TCC_

/// Implementation of SRS interfaces for a ppzkSNARK for Plonk. See
/// srs.hpp .

namespace libsnark
{

// class usrs constructor
template<typename ppT>
usrs<ppT>::usrs(
    std::vector<libff::G1<ppT>> &&secret_powers_g1,
    std::vector<libff::G2<ppT>> &&secret_powers_g2)
    : secret_powers_g1(secret_powers_g1), secret_powers_g2(secret_powers_g2)
{
}

// class srs constructor
template<typename ppT>
srs<ppT>::srs(
    const size_t &num_gates,
    const size_t &num_qpolys,
    const polynomial<Field> &PI_poly,
    const std::vector<polynomial<Field>> &Q_polys,
    const std::vector<polynomial<Field>> &S_polys,
    const std::vector<std::vector<Field>> &omega_roots,
    const std::vector<Field> &H_gen,
    const std::vector<Field> &H_gen_permute,
    const libff::Fr<ppT> &k1,
    const libff::Fr<ppT> &k2,
    std::vector<libff::G1<ppT>> &&secret_powers_g1,
    std::vector<libff::G2<ppT>> &&secret_powers_g2,
    const polynomial<Field> &L_basis_zero,
    std::shared_ptr<libfqfft::evaluation_domain<libff::Fr<ppT>>> domain)
    : num_gates(num_gates)
    , num_qpolys(num_qpolys)
    , PI_poly(PI_poly)
    , Q_polys(Q_polys)
    , S_polys(S_polys)
    , omega_roots(omega_roots)
    , H_gen(H_gen)
    , H_gen_permute(H_gen_permute)
    , k1(k1)
    , k2(k2)
    , secret_powers_g1(secret_powers_g1)
    , secret_powers_g2(secret_powers_g2)
    , L_basis_zero(L_basis_zero)
    , domain(domain)
{
}

/// class plonk_verification_key
template<typename ppT>
plonk_verification_key<ppT>::plonk_verification_key(
    std::vector<libff::G2<ppT>> &&secret_powers_g2)
    : secret_powers_g2(std::move(secret_powers_g2)){};

/// class plonk_keypair constructor
template<typename ppT>
plonk_keypair<ppT>::plonk_keypair(
    plonk_proving_key<ppT> &&pk, plonk_verification_key<ppT> &&vk)
    : pk(std::move(pk)), vk(std::move(vk))
{
}

/// transcript_hasher constructor
template<typename ppT>
transcript_hasher<ppT>::transcript_hasher(size_t &istep) : istep(istep)
{
}

/// dummy implementation of get_hash that directly returns the
/// expected hard-coded hashes for the purposes of unit testing TODO
/// to be replaced by a call to a proper hash function e.g. SHA2,
/// BLAKE, etc.
template<typename ppT> libff::Fr<ppT> transcript_hasher<ppT>::get_hash()
{
    assert((this->istep >= 0) && (this->istep <= 5));
    using Field = libff::Fr<ppT>;

    Field beta = Field("3710899868510394644410941212967766116886736137326022751"
                       "891187938298987182388");
    Field gamma = Field("110379303840831945879077096653321168432672740458288022"
                        "49545114995763715746939");
    Field alpha = Field("379799789992747238930717819864848384921111623418803600"
                        "22719385400306128734648");
    Field zeta = Field("4327197228921839935583364394550235027071910395980312641"
                       "5018065799136107272465");
    Field nu = Field("275158598338697752421507265080923414294782807831923791651"
                     "55175653098691426347");
    Field u = Field("1781751143954696684632449211212056577828855388109883650570"
                    "6049265393896966778");
    if (this->istep == 0) {
        printf("[%s:%d] istep %d\n", __FILE__, __LINE__, (int)istep);
        this->istep++;
        return beta;
    }
    if (this->istep == 1) {
        printf("[%s:%d] istep %d\n", __FILE__, __LINE__, (int)istep);
        this->istep++;
        return gamma;
    }
    if (this->istep == 2) {
        printf("[%s:%d] istep %d\n", __FILE__, __LINE__, (int)istep);
        this->istep++;
        return alpha;
    }
    if (this->istep == 3) {
        printf("[%s:%d] istep %d\n", __FILE__, __LINE__, (int)istep);
        this->istep++;
        return zeta;
    }
    if (this->istep == 4) {
        printf("[%s:%d] istep %d\n", __FILE__, __LINE__, (int)istep);
        this->istep++;
        return nu;
    }
    if (this->istep == 5) {
        // reset step to 0
        printf("[%s:%d] istep %d\n", __FILE__, __LINE__, (int)istep);
        this->istep = 0;
        return u;
    }
    // error
    return 0;
}

/// Compute a universal srs (usrs). It is composed *only* of encoded
/// powers of the secret value in the group generator. Therefore a usrs
/// is independent of any particular circuit.
///
/// \note only for debug
template<typename ppT>
usrs<ppT> plonk_usrs_derive_from_secret(
    const libff::Fr<ppT> &secret, const size_t max_degree)
{
    // compute powers of secret times G1: 1*G1, secret^1*G1, secret^2*G1, ...
    const libff::bigint<libff::Fr<ppT>::num_limbs> secret_bigint =
        secret.as_bigint();
    const size_t window_size = std::max(
        libff::wnaf_opt_window_size<libff::G1<ppT>>(secret_bigint.num_bits()),
        1ul);
    const std::vector<long> naf =
        libff::find_wnaf<libff::Fr<ppT>::num_limbs>(window_size, secret_bigint);

    std::vector<libff::G1<ppT>> secret_powers_g1;
    secret_powers_g1.reserve(max_degree);
    libff::G1<ppT> secret_i_g1 = libff::G1<ppT>::one();
    secret_powers_g1.push_back(secret_i_g1);
    for (size_t i = 1; i < max_degree; ++i) {
        // secret^i * G1
        secret_i_g1 = libff::fixed_window_wnaf_exp<libff::G1<ppT>>(
            window_size, secret_i_g1, naf);
        secret_powers_g1.push_back(secret_i_g1);
    }

    // compute powers of secret times G2: 1*G2, secret^1*G2
    // Note: in Plonk we *always* have 2 encoded elemnts in G2
    std::vector<libff::G2<ppT>> secret_powers_g2;
    secret_powers_g2.reserve(2);
    // secret^0 * G2 = G2
    libff::G2<ppT> secret_0_g2 = libff::G2<ppT>::one();
    secret_powers_g2.push_back(secret_0_g2);
    // secret^1 * G2
    libff::G2<ppT> secret_1_g2 = secret * libff::G2<ppT>::one();
    secret_powers_g2.push_back(secret_1_g2);

    return usrs<ppT>(std::move(secret_powers_g1), std::move(secret_powers_g2));
}

/// Derive the (plain) SRS from the circuit description and the
/// USRS. The (plain) SRS is a specialization of the USRS for one
/// particular circuit i.e.
///
/// usrs = <encoded powers of secret>
/// srs = (proving_key, verificataion_key) = derive(usrs,
/// circuit_description)
template<typename ppT>
srs<ppT> plonk_srs_derive_from_usrs(
    const usrs<ppT> &usrs, const circuit_t<ppT> &circuit)
{
    // secret^i * G1
    std::vector<libff::G1<ppT>> secret_powers_g1;
    secret_powers_g1.reserve(circuit.num_gates + 3);
    for (size_t i = 0; i < (circuit.num_gates + 3); ++i) {
        secret_powers_g1.push_back(usrs.secret_powers_g1[i]);
    }
    // secret^i * G2
    std::vector<libff::G2<ppT>> secret_powers_g2;
    secret_powers_g2.reserve(2);
    for (size_t i = 0; i < 2; ++i) {
        secret_powers_g2.push_back(usrs.secret_powers_g2[i]);
    }

    std::shared_ptr<libfqfft::evaluation_domain<libff::Fr<ppT>>> domain =
        libfqfft::get_evaluation_domain<libff::Fr<ppT>>(circuit.num_gates);

    // compute 0-th Lagrange basis vector via inverse FFT
    polynomial<libff::Fr<ppT>> u(circuit.num_gates, libff::Fr<ppT>(0));
    u[0] = libff::Fr<ppT>(1);
    domain->iFFT(u);
    polynomial<libff::Fr<ppT>> L_basis_zero = u;

    srs<ppT> srs(
        circuit.num_gates,
        circuit.num_qpolys,
        circuit.PI_poly,
        circuit.Q_polys,
        circuit.S_polys,
        circuit.omega_roots,
        circuit.H_gen,
        circuit.H_gen_permute,
        circuit.k1,
        circuit.k2,
        std::move(secret_powers_g1),
        std::move(secret_powers_g2),
        std::move(L_basis_zero),
        domain);

    return srs;
}

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_TCC_
