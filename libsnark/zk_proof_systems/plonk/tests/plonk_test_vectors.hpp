#ifndef __PLONK_TEST_VECTORS_HPP__
#define __PLONK_TEST_VECTORS_HPP__

namespace libsnark
{

template<typename ppT> class plonk_test_vectors
{
public:
  using Field = libff::Fr<ppT>;
  using BaseField = libff::Fq<ppT>;

  // output from round 3
  std::vector<std::vector<BaseField>> t_poly_at_secret_g1;
  
  void initialize();
};

} // namespace libsnark
  
#include <libsnark/zk_proof_systems/plonk/tests/plonk_test_vectors.cpp>

#endif // __PLONK_TEST_VECTORS_HPP__
