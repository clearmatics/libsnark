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
//#include <libsnark/polynomial_commitments/kzg10.hpp>
//#include <libsnark/polynomial_commitments/tests/polynomial_commitment_test_utils.hpp>

const size_t MAX_DEGREE = 254;

/*

example circuit

P(x) = x**3 + x + 5 = 3

circuit has 6 gates + 2 dummy gates

gates

      L    R    O
1 mul   x * x  = v1
2 mul  v1 * x  = v2
3 add  v2 + x  = v3
4 con5  1 * 5  = 5
5 pin   1 * 35 = 35
6 add  v3 + 5  = 35
7 dum   /    /    /
8 dum   /    /    /

wire polynomials

w_l = [ x, v1, v2,  1,  1, v3,  /,  /] = [a1, a2, a3, a4, a5, a6, a7] = a
w_r = [ x,  x,  x,  5, 35,  5,  /,  /] = [b1, b2, b3, b4, b5, b6, b7] = b
w_o = [v1, v2, v3,  5, 35, 35,  /,  /] = [c1, c2, c3, c4, c5, c6, c7] = c

witness

x = 3 => v1 = 9, v2 = 27, v3 = 30

w_l = a = [ 3,  9, 27,  1,  1, 30,  0,  0]
w_r = b = [ 3,  3,  3,  5, 35,  5,  0,  0]
w_o = c = [ 9, 27, 30,  5, 35, 35,  0,  0]

W = w_l + w_r + w_o

encoding of plonk gates

add  = [1, 1,-1, 0, 0]
mul  = [0, 0,-1, 1, 0]
con5 = [0, 1, 0, 0, 5]
pi   = [0, 1, 0, 0, 0]
dum  = [0, 0, 0, 0, 0]

gate polynomials

(q_L * a) + (q_r * b) + (q_o * c) + (q_m * a * b) + (q_c) = 0

                       q_l q_r q_o q_m q_c
1 mul   x * x  = v1 : [  0,  0, -1,  1,  0]
2 mul  v1 * x  = v2 : [  0,  0, -1,  1,  0]
3 add  v2 + x  = v3 : [  1,  1, -1,  0,  0]
4 con5  1 * 5  = 5  : [  0,  1,  0,  0, -5]
5 pin   1 * 35 = 35 : [  0,  1,  0,  0,  0]
6 add  v2 + 5  = 35 : [  1,  1,  0, -1,  0]
7 dum   /    /    / : [  0,  0,  0,  0,  0]
8 dum   /    /    / : [  0,  0,  0,  0,  0]

q_l = [ 0,  0,  1,  0,  0,  1,  0,  0]
q_r = [ 0,  0,  1,  1,  1,  1,  0,  0]
q_o = [-1, -1, -1,  0,  0,  0,  0,  0]
q_m = [ 1,  1,  0,  0,  0, -1,  0,  0]
q_c = [ 0,  0,  0, -1,  0,  0,  0,  0]


  auto reference_string = SRS
  proving_key key = 
  program_witness witness = 

  polynomial w_l;
  polynomial w_r;
  polynomial w_o;
  polynomial q_l;
  polynomial q_r;
  polynomial q_o;
  polynomial q_c;
  polynomial q_m;


*/
  
namespace libsnark
{

template<typename ppT> void test_plonk()
{
  //    typedef libff::Fr<default_r1cs_ppzksnark_pp> FieldT;
  //    default_r1cs_ppzksnark_pp::init_public_params();
  //  Fp_model;

    // Execute all tests for the given curve.
    ppT::init_public_params();
    //    assert(0);
    printf("[%s:%d] Test OK\n", __FILE__, __LINE__);
}

TEST(TestPlonk, BLS12_381)
{
    test_plonk<libff::bls12_381_pp>();
}

} // namespace libsnark
