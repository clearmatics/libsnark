template<typename ppT> void test_plonk()
{
    // Execute all tests for the given curve.
    ppT::init_public_params();
    printf("[%s:%d] Test OK\n", __FILE__, __LINE__);    
}

TEST(TestPolynomialCommitments, ALT_BN128)
{
    test_plonk<libff::alt_bn128_pp>();
}
