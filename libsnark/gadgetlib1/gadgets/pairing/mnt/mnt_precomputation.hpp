/** @file
 *****************************************************************************

 Declaration of interfaces for pairing precomputation gadgets.

 The gadgets verify correct precomputation of values for the G1 and G2
 variables.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_PAIRING_MNT_MNT_PRECOMPUTATION_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_PAIRING_MNT_MNT_PRECOMPUTATION_HPP_

#include "libsnark/gadgetlib1/gadgets/curves/weierstrass_g1_gadget.hpp"
#include "libsnark/gadgetlib1/gadgets/curves/weierstrass_g2_gadget.hpp"
#include "libsnark/gadgetlib1/gadgets/pairing/pairing_params.hpp"

#include <libff/algebra/curves/mnt/mnt4/mnt4_init.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_init.hpp>
#include <memory>

namespace libsnark
{

// Forward declare some mnt-specific parameters
template<typename ppT> class mnt_pairing_params;

/**************************** G1 Precomputation ******************************/

/**
 * Not a gadget. It only holds values.
 */
template<typename ppT> class mnt_G1_precomputation
{
public:
    typedef libff::Fr<ppT> FieldT;
    typedef libff::Fqe<other_curve<ppT>> FqeT;
    typedef libff::Fqk<other_curve<ppT>> FqkT;

    std::shared_ptr<G1_variable<ppT>> P;
    std::shared_ptr<Fqe_variable<ppT>> PY_twist_squared;

    mnt_G1_precomputation();
    mnt_G1_precomputation(
        protoboard<FieldT> &pb,
        const libff::G1<other_curve<ppT>> &P,
        const std::string &annotation_prefix);
};

/**
 * Gadget that verifies correct precomputation of the G1 variable.
 */
template<typename ppT>
class mnt_precompute_G1_gadget : public gadget<libff::Fr<ppT>>
{
public:
    typedef libff::Fqe<other_curve<ppT>> FqeT;
    typedef libff::Fqk<other_curve<ppT>> FqkT;

    mnt_G1_precomputation<ppT> &precomp; // must be a reference.

    /* two possible pre-computations one for mnt4 and one for mnt6 */
    template<typename FieldT>
    mnt_precompute_G1_gadget(
        protoboard<FieldT> &pb,
        const G1_variable<ppT> &P,
        mnt_G1_precomputation<ppT> &precomp, // will allocate this inside
        const std::string &annotation_prefix,
        const typename std::enable_if<
            libff::Fqk<other_curve<ppT>>::extension_degree() == 4,
            FieldT>::type & = FieldT())
        : gadget<FieldT>(pb, annotation_prefix), precomp(precomp)
    {
        pb_linear_combination<FieldT> c0, c1;
        c0.assign(pb, P.Y * ((libff::mnt4_twist).squared().coeffs[0]));
        c1.assign(pb, P.Y * ((libff::mnt4_twist).squared().coeffs[1]));

        precomp.P.reset(new G1_variable<ppT>(P));
        precomp.PY_twist_squared.reset(new Fqe_variable<ppT>(
            pb, c0, c1, FMT(annotation_prefix, " PY_twist_squared")));
    }

    template<typename FieldT>
    mnt_precompute_G1_gadget(
        protoboard<FieldT> &pb,
        const G1_variable<ppT> &P,
        mnt_G1_precomputation<ppT> &precomp, // will allocate this inside
        const std::string &annotation_prefix,
        const typename std::enable_if<
            libff::Fqk<other_curve<ppT>>::extension_degree() == 6,
            FieldT>::type & = FieldT())
        : gadget<FieldT>(pb, annotation_prefix), precomp(precomp)
    {
        pb_linear_combination<FieldT> c0, c1, c2;
        c0.assign(pb, P.Y * ((libff::mnt6_twist).squared().coeffs[0]));
        c1.assign(pb, P.Y * ((libff::mnt6_twist).squared().coeffs[1]));
        c2.assign(pb, P.Y * ((libff::mnt6_twist).squared().coeffs[2]));

        precomp.P.reset(new G1_variable<ppT>(P));
        precomp.PY_twist_squared.reset(new Fqe_variable<ppT>(
            pb, c0, c1, c2, FMT(annotation_prefix, " PY_twist_squared")));
    }

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/**************************** G2 Precomputation ******************************/

/**
 * Not a gadget. It only holds values.
 */
template<typename ppT> class mnt_precompute_G2_gadget_coeffs
{
public:
    typedef libff::Fr<ppT> FieldT;
    typedef libff::Fqe<other_curve<ppT>> FqeT;
    typedef libff::Fqk<other_curve<ppT>> FqkT;

    std::shared_ptr<Fqe_variable<ppT>> RX;
    std::shared_ptr<Fqe_variable<ppT>> RY;
    std::shared_ptr<Fqe_variable<ppT>> gamma;
    std::shared_ptr<Fqe_variable<ppT>> gamma_X;

    mnt_precompute_G2_gadget_coeffs();
    mnt_precompute_G2_gadget_coeffs(
        protoboard<FieldT> &pb, const std::string &annotation_prefix);
    mnt_precompute_G2_gadget_coeffs(
        protoboard<FieldT> &pb,
        const G2_variable<ppT> &Q,
        const std::string &annotation_prefix);
};

/**
 * Not a gadget. It only holds values.
 */
template<typename ppT> class mnt_G2_precomputation
{
public:
    typedef libff::Fr<ppT> FieldT;
    typedef libff::Fqe<other_curve<ppT>> FqeT;
    typedef libff::Fqk<other_curve<ppT>> FqkT;

    std::shared_ptr<G2_variable<ppT>> Q;

    std::vector<std::shared_ptr<mnt_precompute_G2_gadget_coeffs<ppT>>> coeffs;

    mnt_G2_precomputation();
    mnt_G2_precomputation(
        protoboard<FieldT> &pb,
        const libff::G2<other_curve<ppT>> &Q_val,
        const std::string &annotation_prefix);
};

/**
 * Technical note:
 *
 * QX and QY -- X and Y coordinates of Q
 *
 * initialization:
 * coeffs[0].RX = QX
 * coeffs[0].RY = QY
 *
 * G2_precompute_doubling_step relates coeffs[i] and coeffs[i+1] as follows
 *
 * coeffs[i]
 * gamma = (3 * RX^2 + twist_coeff_a) * (2*RY).inverse()
 * gamma_X = gamma * RX
 *
 * coeffs[i+1]
 * RX = prev_gamma^2 - (2*prev_RX)
 * RY = prev_gamma * (prev_RX - RX) - prev_RY
 */
template<typename ppT>
class mnt_precompute_G2_gadget_doubling_step : public gadget<libff::Fr<ppT>>
{
public:
    typedef libff::Fr<ppT> FieldT;
    typedef libff::Fqe<other_curve<ppT>> FqeT;
    typedef libff::Fqk<other_curve<ppT>> FqkT;

    mnt_precompute_G2_gadget_coeffs<ppT> cur;
    mnt_precompute_G2_gadget_coeffs<ppT> next;

    std::shared_ptr<Fqe_variable<ppT>> RXsquared;
    std::shared_ptr<Fqe_sqr_gadget<ppT>> compute_RXsquared;
    std::shared_ptr<Fqe_variable<ppT>> three_RXsquared_plus_a;
    std::shared_ptr<Fqe_variable<ppT>> two_RY;
    std::shared_ptr<Fqe_mul_gadget<ppT>> compute_gamma;
    std::shared_ptr<Fqe_mul_gadget<ppT>> compute_gamma_X;

    std::shared_ptr<Fqe_variable<ppT>> next_RX_plus_two_RX;
    std::shared_ptr<Fqe_sqr_gadget<ppT>> compute_next_RX;

    std::shared_ptr<Fqe_variable<ppT>> RX_minus_next_RX;
    std::shared_ptr<Fqe_variable<ppT>> RY_plus_next_RY;
    std::shared_ptr<Fqe_mul_gadget<ppT>> compute_next_RY;

    mnt_precompute_G2_gadget_doubling_step(
        protoboard<FieldT> &pb,
        const mnt_precompute_G2_gadget_coeffs<ppT> &cur,
        const mnt_precompute_G2_gadget_coeffs<ppT> &next,
        const std::string &annotation_prefix);
    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/**
 * Technical note:
 *
 * G2_precompute_addition_step relates coeffs[i] and coeffs[i+1] as follows
 *
 * coeffs[i]
 * gamma = (RY - QY) * (RX - QX).inverse()
 * gamma_X = gamma * QX
 *
 * coeffs[i+1]
 * RX = prev_gamma^2 + (prev_RX + QX)
 * RY = prev_gamma * (prev_RX - RX) - prev_RY
 *
 * (where prev_ in [i+1] refer to things from [i])
 *
 * If invert_Q is set to true: use -QY in place of QY everywhere above.
 */
template<typename ppT>
class mnt_precompute_G2_gadget_addition_step : public gadget<libff::Fr<ppT>>
{
public:
    typedef libff::Fr<ppT> FieldT;
    typedef libff::Fqe<other_curve<ppT>> FqeT;
    typedef libff::Fqk<other_curve<ppT>> FqkT;

    bool invert_Q;
    mnt_precompute_G2_gadget_coeffs<ppT> cur;
    mnt_precompute_G2_gadget_coeffs<ppT> next;
    G2_variable<ppT> Q;

    std::shared_ptr<Fqe_variable<ppT>> RY_minus_QY;
    std::shared_ptr<Fqe_variable<ppT>> RX_minus_QX;
    std::shared_ptr<Fqe_mul_gadget<ppT>> compute_gamma;
    std::shared_ptr<Fqe_mul_gadget<ppT>> compute_gamma_X;

    std::shared_ptr<Fqe_variable<ppT>> next_RX_plus_RX_plus_QX;
    std::shared_ptr<Fqe_sqr_gadget<ppT>> compute_next_RX;

    std::shared_ptr<Fqe_variable<ppT>> RX_minus_next_RX;
    std::shared_ptr<Fqe_variable<ppT>> RY_plus_next_RY;
    std::shared_ptr<Fqe_mul_gadget<ppT>> compute_next_RY;

    mnt_precompute_G2_gadget_addition_step(
        protoboard<FieldT> &pb,
        const bool invert_Q,
        const mnt_precompute_G2_gadget_coeffs<ppT> &cur,
        const mnt_precompute_G2_gadget_coeffs<ppT> &next,
        const G2_variable<ppT> &Q,
        const std::string &annotation_prefix);
    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/**
 * Gadget that verifies correct precomputation of the G2 variable.
 */
template<typename ppT>
class mnt_precompute_G2_gadget : public gadget<libff::Fr<ppT>>
{
public:
    typedef libff::Fr<ppT> FieldT;
    typedef libff::Fqe<other_curve<ppT>> FqeT;
    typedef libff::Fqk<other_curve<ppT>> FqkT;

    std::vector<std::shared_ptr<mnt_precompute_G2_gadget_addition_step<ppT>>>
        addition_steps;
    std::vector<std::shared_ptr<mnt_precompute_G2_gadget_doubling_step<ppT>>>
        doubling_steps;

    size_t add_count;
    size_t dbl_count;

    mnt_G2_precomputation<ppT> &precomp; // important to have a reference here

    mnt_precompute_G2_gadget(
        protoboard<FieldT> &pb,
        const G2_variable<ppT> &Q,
        mnt_G2_precomputation<ppT> &precomp, // will allocate this inside
        const std::string &annotation_prefix);
    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

} // namespace libsnark

#include "libsnark/gadgetlib1/gadgets/pairing/mnt/mnt_precomputation.tcc"

#endif // LIBSNARK_GADGETLIB1_GADGETS_PAIRING_MNT_MNT_PRECOMPUTATION_HPP_
