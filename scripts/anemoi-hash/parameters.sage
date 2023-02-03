#!/usr/bin/sage
# -*- mode: python ; -*-

from sage.all import *
import hashlib
import itertools
import datetime

from constants import *

load('anemoi.sage')

def anemoi_instances_bls12_381(A):

    # - 256-bit security level instantiations
    # -- BLS12_381_SCALRFIELD
    # --- 1 col
    A_BLS_12_381_SCALARFIELD_1_COL_256_BITS = AnemoiPermutation(
        q=BLS12_381_SCALARFIELD,
        n_cols=1,
        security_level=256
    )
    A.append(
        ("A_BLS_12_381_SCALARFIELD_1_COL_256_BITS",
         A_BLS_12_381_SCALARFIELD_1_COL_256_BITS))
    # --- 2 col    
    A_BLS_12_381_SCALARFIELD_2_COL_256_BITS = AnemoiPermutation(
        q=BLS12_381_SCALARFIELD,
        n_cols=2,
        security_level=256
    )
    A.append(
        ("A_BLS_12_381_SCALARFIELD_2_COL_256_BITS",
         A_BLS_12_381_SCALARFIELD_2_COL_256_BITS))    
    # --- 3 col    
    A_BLS_12_381_SCALARFIELD_3_COL_256_BITS = AnemoiPermutation(
        q=BLS12_381_SCALARFIELD,
        n_cols=3,
        security_level=256
    )
    A.append(
        ("A_BLS_12_381_SCALARFIELD_3_COL_256_BITS",
         A_BLS_12_381_SCALARFIELD_3_COL_256_BITS))    
    # ---4 col    
    A_BLS_12_381_SCALARFIELD_4_COL_256_BITS = AnemoiPermutation(
        q=BLS12_381_SCALARFIELD,
        n_cols=4,
        security_level=256
    )
    A.append(
        ("A_BLS_12_381_SCALARFIELD_4_COL_256_BITS",
         A_BLS_12_381_SCALARFIELD_4_COL_256_BITS))
    
def anemoi_instances_bls12_377(A):

    # - 256-bit security level instantiations
    # -- BLS12_377_SCALRFIELD
    # --- 1 col
    A_BLS_12_377_SCALARFIELD_1_COL_256_BITS = AnemoiPermutation(
        q=BLS12_377_SCALARFIELD,
        n_cols=1,
        security_level=256
    )
    A.append(
        ("A_BLS_12_377_SCALARFIELD_1_COL_256_BITS",
         A_BLS_12_377_SCALARFIELD_1_COL_256_BITS))
    # --- 2 col    
    A_BLS_12_377_SCALARFIELD_2_COL_256_BITS = AnemoiPermutation(
        q=BLS12_377_SCALARFIELD,
        n_cols=2,
        security_level=256
    )
    A.append(
        ("A_BLS_12_377_SCALARFIELD_2_COL_256_BITS",
         A_BLS_12_377_SCALARFIELD_2_COL_256_BITS))    
    # --- 3 col    
    A_BLS_12_377_SCALARFIELD_3_COL_256_BITS = AnemoiPermutation(
        q=BLS12_377_SCALARFIELD,
        n_cols=3,
        security_level=256
    )
    A.append(
        ("A_BLS_12_377_SCALARFIELD_3_COL_256_BITS",
         A_BLS_12_377_SCALARFIELD_3_COL_256_BITS))    
    # ---4 col    
    A_BLS_12_377_SCALARFIELD_4_COL_256_BITS = AnemoiPermutation(
        q=BLS12_377_SCALARFIELD,
        n_cols=4,
        security_level=256
    )
    A.append(
        ("A_BLS_12_377_SCALARFIELD_4_COL_256_BITS",
         A_BLS_12_377_SCALARFIELD_4_COL_256_BITS))
    
def anemoi_instances_mnt4(A):

    # - 256-bit security level instantiations
    # -- MNT4_SCALRFIELD
    # --- 1 col
    A_MNT4_SCALARFIELD_1_COL_256_BITS = AnemoiPermutation(
        q=MNT4_SCALARFIELD,
        n_cols=1,
        security_level=256
    )
    A.append(
        ("A_MNT4_SCALARFIELD_1_COL_256_BITS",
         A_MNT4_SCALARFIELD_1_COL_256_BITS))
    # --- 2 col    
    A_MNT4_SCALARFIELD_2_COL_256_BITS = AnemoiPermutation(
        q=MNT4_SCALARFIELD,
        n_cols=2,
        security_level=256
    )
    A.append(
        ("A_MNT4_SCALARFIELD_2_COL_256_BITS",
         A_MNT4_SCALARFIELD_2_COL_256_BITS))    
    # --- 3 col    
    A_MNT4_SCALARFIELD_3_COL_256_BITS = AnemoiPermutation(
        q=MNT4_SCALARFIELD,
        n_cols=3,
        security_level=256
    )
    A.append(
        ("A_MNT4_SCALARFIELD_3_COL_256_BITS",
         A_MNT4_SCALARFIELD_3_COL_256_BITS))    
    # ---4 col    
    A_MNT4_SCALARFIELD_4_COL_256_BITS = AnemoiPermutation(
        q=MNT4_SCALARFIELD,
        n_cols=4,
        security_level=256
    )
    A.append(
        ("A_MNT4_SCALARFIELD_4_COL_256_BITS",
         A_MNT4_SCALARFIELD_4_COL_256_BITS))
    
def anemoi_instances_mnt6(A):

    # - 256-bit security level instantiations
    # -- MNT6_SCALRFIELD
    # --- 1 col
    A_MNT6_SCALARFIELD_1_COL_256_BITS = AnemoiPermutation(
        q=MNT6_SCALARFIELD,
        n_cols=1,
        security_level=256
    )
    A.append(
        ("A_MNT6_SCALARFIELD_1_COL_256_BITS",
         A_MNT6_SCALARFIELD_1_COL_256_BITS))
    # --- 2 col    
    A_MNT6_SCALARFIELD_2_COL_256_BITS = AnemoiPermutation(
        q=MNT6_SCALARFIELD,
        n_cols=2,
        security_level=256
    )
    A.append(
        ("A_MNT6_SCALARFIELD_2_COL_256_BITS",
         A_MNT6_SCALARFIELD_2_COL_256_BITS))    
    # --- 3 col    
    A_MNT6_SCALARFIELD_3_COL_256_BITS = AnemoiPermutation(
        q=MNT6_SCALARFIELD,
        n_cols=3,
        security_level=256
    )
    A.append(
        ("A_MNT6_SCALARFIELD_3_COL_256_BITS",
         A_MNT6_SCALARFIELD_3_COL_256_BITS))    
    # ---4 col    
    A_MNT6_SCALARFIELD_4_COL_256_BITS = AnemoiPermutation(
        q=MNT6_SCALARFIELD,
        n_cols=4,
        security_level=256
    )
    A.append(
        ("A_MNT6_SCALARFIELD_4_COL_256_BITS",
         A_MNT6_SCALARFIELD_4_COL_256_BITS))
    
def anemoi_instances_bw6_761(A):

    # - 256-bit security level instantiations
    # -- BW6_761_SCALRFIELD
    # --- 1 col
    A_BW6_761_SCALARFIELD_1_COL_256_BITS = AnemoiPermutation(
        q=BW6_761_SCALARFIELD,
        n_cols=1,
        security_level=256
    )
    A.append(
        ("A_BW6_761_SCALARFIELD_1_COL_256_BITS",
         A_BW6_761_SCALARFIELD_1_COL_256_BITS))
    # --- 2 col    
    A_BW6_761_SCALARFIELD_2_COL_256_BITS = AnemoiPermutation(
        q=BW6_761_SCALARFIELD,
        n_cols=2,
        security_level=256
    )
    A.append(
        ("A_BW6_761_SCALARFIELD_2_COL_256_BITS",
         A_BW6_761_SCALARFIELD_2_COL_256_BITS))    
    # --- 3 col    
    A_BW6_761_SCALARFIELD_3_COL_256_BITS = AnemoiPermutation(
        q=BW6_761_SCALARFIELD,
        n_cols=3,
        security_level=256
    )
    A.append(
        ("A_BW6_761_SCALARFIELD_3_COL_256_BITS",
         A_BW6_761_SCALARFIELD_3_COL_256_BITS))    
    # ---4 col    
    A_BW6_761_SCALARFIELD_4_COL_256_BITS = AnemoiPermutation(
        q=BW6_761_SCALARFIELD,
        n_cols=4,
        security_level=256
    )
    A.append(
        ("A_BW6_761_SCALARFIELD_4_COL_256_BITS",
         A_BW6_761_SCALARFIELD_4_COL_256_BITS))
    
def anemoi_instances_bn128(A):

    # - 256-bit security level instantiations
    # -- BN128_SCALRFIELD
    # --- 1 col
    A_BN128_SCALARFIELD_1_COL_256_BITS = AnemoiPermutation(
        q=BN128_SCALARFIELD,
        n_cols=1,
        security_level=256
    )
    A.append(
        ("A_BN128_SCALARFIELD_1_COL_256_BITS",
         A_BN128_SCALARFIELD_1_COL_256_BITS))
    # --- 2 col    
    A_BN128_SCALARFIELD_2_COL_256_BITS = AnemoiPermutation(
        q=BN128_SCALARFIELD,
        n_cols=2,
        security_level=256
    )
    A.append(
        ("A_BN128_SCALARFIELD_2_COL_256_BITS",
         A_BN128_SCALARFIELD_2_COL_256_BITS))    
    # --- 3 col    
    A_BN128_SCALARFIELD_3_COL_256_BITS = AnemoiPermutation(
        q=BN128_SCALARFIELD,
        n_cols=3,
        security_level=256
    )
    A.append(
        ("A_BN128_SCALARFIELD_3_COL_256_BITS",
         A_BN128_SCALARFIELD_3_COL_256_BITS))    
    # ---4 col    
    A_BN128_SCALARFIELD_4_COL_256_BITS = AnemoiPermutation(
        q=BN128_SCALARFIELD,
        n_cols=4,
        security_level=256
    )
    A.append(
        ("A_BN128_SCALARFIELD_4_COL_256_BITS",
         A_BN128_SCALARFIELD_4_COL_256_BITS))
    
def anemoi_instances_alt_bn128(A):

    # - 256-bit security level instantiations
    # -- ALT_BN128_SCALRFIELD
    # --- 1 col
    A_ALT_BN128_SCALARFIELD_1_COL_256_BITS = AnemoiPermutation(
        q=ALT_BN128_SCALARFIELD,
        n_cols=1,
        security_level=256
    )
    A.append(
        ("A_ALT_BN128_SCALARFIELD_1_COL_256_BITS",
         A_ALT_BN128_SCALARFIELD_1_COL_256_BITS))
    # --- 2 col    
    A_ALT_BN128_SCALARFIELD_2_COL_256_BITS = AnemoiPermutation(
        q=ALT_BN128_SCALARFIELD,
        n_cols=2,
        security_level=256
    )
    A.append(
        ("A_ALT_BN128_SCALARFIELD_2_COL_256_BITS",
         A_ALT_BN128_SCALARFIELD_2_COL_256_BITS))    
    # --- 3 col    
    A_ALT_BN128_SCALARFIELD_3_COL_256_BITS = AnemoiPermutation(
        q=ALT_BN128_SCALARFIELD,
        n_cols=3,
        security_level=256
    )
    A.append(
        ("A_ALT_BN128_SCALARFIELD_3_COL_256_BITS",
         A_ALT_BN128_SCALARFIELD_3_COL_256_BITS))    
    # ---4 col    
    A_ALT_BN128_SCALARFIELD_4_COL_256_BITS = AnemoiPermutation(
        q=ALT_BN128_SCALARFIELD,
        n_cols=4,
        security_level=256
    )
    A.append(
        ("A_ALT_BN128_SCALARFIELD_4_COL_256_BITS",
         A_ALT_BN128_SCALARFIELD_4_COL_256_BITS))
    
def anemoi_instances_stdout(instances):
    for i in range(len(instances)):
        # string name
        A_str = instances[i][0]
        # actual instance that can be called as A[i][1].*
        A = instances[i][1]
        zero = 0
        width = 100
        print("------------------------------------------------------")
        print("instance         : {}".format(A_str))
        print("prime field      : {}".format(A.prime_field))
        print("Fr modulus       : {}".format(A.q))
        print("n_cols           : {}".format(A.n_cols))
        print("n_rounds         : {}".format(A.n_rounds))
        print("security level   : {}".format(A.security_level))
        print("mult generator g : {}".format(A.g))
        print("Q power          : {}".format(A.QUAD))
        print("alpha            : {}".format(A.alpha))
        print("alpha_inv        : {}".format(A.alpha_inv))
        print("beta             : {}".format(A.beta))
        print("gamma            : {}".format(zero))
        print("delta            : {}".format(A.delta))
        print("matrix M         :\n{}".format(A.mat))
        print("constants C      :\n{}".format(A.C))
        print("constants D      :\n{}".format(A.D))

# same as output_parameters() but stores parameters to file
def anemoi_instances_to_file(instances):
    f = open("instances.txt", "w")
    e = datetime.datetime.now()
    f.write("This file was automatically generated with SAGE script parameters.sage on %s/%s/%s at %s:%s:%s\n" % (e.day, e.month, e.year, e.hour, e.minute, e.second))
    for i in range(len(instances)):
        A_str = instances[i][0]
        A = instances[i][1]
        zero = 0
        width = 100
        f.write("------------------------------------------------------")
        f.write("instance         : {}\n".format(A_str))
        f.write("prime field      : {}\n".format(A.prime_field))
        f.write("Fr modulus       : {}\n".format(A.q))
        f.write("n_cols           : {}\n".format(A.n_cols))
        f.write("n_rounds         : {}\n".format(A.n_rounds))
        f.write("security level   : {}\n".format(A.security_level))
        f.write("mult generator g : {}\n".format(A.g))
        f.write("Q power          : {}\n".format(A.QUAD))
        f.write("alpha            : {}\n".format(A.alpha))
        f.write("alpha_inv        : {}\n".format(A.alpha_inv))
        f.write("beta             : {}\n".format(A.beta))
        f.write("gamma            : {}\n".format(zero))
        f.write("delta            : {}\n".format(A.delta))
        f.write("matrix M         :\n{}\n".format(A.mat))
        f.write("constants C      :\n{}\n".format(A.C))
        f.write("constants D      :\n{}\n".format(A.D))

def anemoi_parameters_in_cpp_format_to_file(instances, filename, curve_ppT):
    f = open(filename, "w")
    e = datetime.datetime.now()
    f.write("// This file was automatically generated with SAGE script parameters.sage on %s/%s/%s at %s:%s:%s\n\n" % (e.day, e.month, e.year, e.hour, e.minute, e.second))

    f.write("// Anemoi parameters for curve {}\n".format(curve_ppT))
    # get just the first instance -- all instances from a given curve
    # share the same parameters
    A = instances[0][1]

    f.write("template<> class anemoi_parameters<libff::{}>\n".format(curve_ppT))
    f.write("{\npublic:\n")    
    f.write("using ppT = libff::{};\n".format(curve_ppT))
    f.write("using FieldT = libff::Fr<ppT>;\n")
    f.write("using BignumT = libff::bigint<FieldT::num_limbs>;\n")
    f.write("static const bool b_prime_field = false;\n")
    f.write("static constexpr size_t multiplicative_generator_g = {};\n".format(A.g))
    f.write("static constexpr size_t alpha = {};\n".format(A.alpha))
    f.write("static constexpr size_t beta = multiplicative_generator_g;\n")
    f.write("static constexpr size_t gamma = 0;\n")
    f.write("static constexpr size_t quad_exponent = {};\n".format(A.QUAD))
    f.write("static const BignumT alpha_inv;\n")
    f.write("static const BignumT delta;\n")
    f.write("static const std::vector<std::vector<BignumT>> C_constants_col_one;\n")
    f.write("static const std::vector<std::vector<BignumT>> D_constants_col_one;\n")
    f.write("static const std::vector<std::vector<BignumT>> C_constants_col_two;\n")
    f.write("static const std::vector<std::vector<BignumT>> D_constants_col_two;\n")
    f.write("static const std::vector<std::vector<BignumT>> C_constants_col_three;\n")
    f.write("static const std::vector<std::vector<BignumT>> D_constants_col_three;\n")
    f.write("static const std::vector<std::vector<BignumT>> C_constants_col_four;\n")
    f.write("static const std::vector<std::vector<BignumT>> D_constants_col_four;\n")
    f.write("};\n")
    
    f.write("\n")    
    f.write("const anemoi_parameters<libff::{}>::BignumT anemoi_parameters<libff::{}>::alpha_inv = anemoi_parameters<libff::{}>::BignumT(\"{}\");\n".format(curve_ppT, curve_ppT, curve_ppT, A.alpha_inv))
    
    f.write("\n")    
    f.write("const anemoi_parameters<libff::{}>::BignumT anemoi_parameters<libff::{}>::delta = anemoi_parameters<libff::{}>::BignumT(\"{}\");\n".format(curve_ppT, curve_ppT, curve_ppT, A.delta))
    
#    f.write("namespace libsnark \n{\n")
#    f.write("} // namespace libsnark")
    
def anemoi_constants_in_cpp_format_to_file(instances, filename, curve_ppT):
    f = open(filename, "a")
    f.write("\n")    
    i_str = ["one", "two", "three", "four", "five", "six"]
    for i in range(len(instances)):
        A_str = instances[i][0]
        A = instances[i][1]
        f.write("// C constants for L = {} columns\n".format(i+1))
        f.write("const std::vector<std::vector<anemoi_parameters<libff::{}>::BignumT>> anemoi_parameters<libff::{}>::C_constants_col_{} = ".format(curve_ppT, curve_ppT, i_str[i]))
        f.write("{\n")
        for iround in range(len(A.C)):
            f.write("{")
            for icol in range(len(A.C[iround])):
                f.write("anemoi_parameters<libff::{}>::BignumT(\"{}\")".format(curve_ppT, A.C[iround][icol]))
                if icol < (len(A.C[iround]) - 1):
                    f.write(", ")
            f.write("}")
            if iround < (len(A.C) - 1):
                f.write(",\n")
        f.write("\n};\n")
        f.write("// D constants for L = {} columns\n".format(i+1))
        f.write("const std::vector<std::vector<anemoi_parameters<libff::{}>::BignumT>> anemoi_parameters<libff::{}>::D_constants_col_{} = ".format(curve_ppT, curve_ppT, i_str[i]))
        f.write("{\n")
        for iround in range(len(A.D)):
            f.write("{")
            for icol in range(len(A.D[iround])):
                f.write("anemoi_parameters<libff::{}>::BignumT(\"{}\")".format(curve_ppT, A.D[iround][icol]))
                if icol < (len(A.D[iround]) - 1):
                    f.write(", ")
            f.write("}")
            if iround < (len(A.D) - 1):
                f.write(",\n")
        f.write("\n};\n")
        
                  
if __name__ == "__main__":
    # bls12_381
    if 0:
        A = []
        anemoi_instances_bls12_381(A)
        filename = "parameters_bls12_381.txt"
        curve_ppT = "bls12_381_pp"
        anemoi_parameters_in_cpp_format_to_file(A, filename, curve_ppT)
        anemoi_constants_in_cpp_format_to_file(A, filename, curve_ppT)
    # bls12_377
    if 0:
        A = []
        anemoi_instances_bls12_377(A)
        filename = "parameters_bls12_377.txt"
        curve_ppT = "bls12_377_pp"
        anemoi_parameters_in_cpp_format_to_file(A, filename, curve_ppT)
        anemoi_constants_in_cpp_format_to_file(A, filename, curve_ppT)
    # mnt4
    if 0:
        A = []
        anemoi_instances_mnt4(A)
        filename = "parameters_mnt4.txt"
        curve_ppT = "mnt4_pp"
        anemoi_parameters_in_cpp_format_to_file(A, filename, curve_ppT)
        anemoi_constants_in_cpp_format_to_file(A, filename, curve_ppT)
    # mnt6
    if 0:
        A = []
        anemoi_instances_mnt6(A)
        filename = "parameters_mnt6.txt"
        curve_ppT = "mnt6_pp"
        anemoi_parameters_in_cpp_format_to_file(A, filename, curve_ppT)
        anemoi_constants_in_cpp_format_to_file(A, filename, curve_ppT)
    # bw6_761 (WARNING! slow ~10 min.)
    if 0:
        A = []
        anemoi_instances_bw6_761(A)
        filename = "parameters_bw6_761.txt"
        curve_ppT = "bw6_761_pp"
        anemoi_parameters_in_cpp_format_to_file(A, filename, curve_ppT)
        anemoi_constants_in_cpp_format_to_file(A, filename, curve_ppT)
    # bn128
    if 1:
        A = []
        anemoi_instances_bn128(A)
        filename = "parameters_bn128.txt"
        curve_ppT = "bn128_pp"
        anemoi_parameters_in_cpp_format_to_file(A, filename, curve_ppT)
        anemoi_constants_in_cpp_format_to_file(A, filename, curve_ppT)
    # alt_bn128
    if 1:
        A = []
        anemoi_instances_alt_bn128(A)
        filename = "parameters_alt_bn128.txt"
        curve_ppT = "alt_bn128_pp"
        anemoi_parameters_in_cpp_format_to_file(A, filename, curve_ppT)
        anemoi_constants_in_cpp_format_to_file(A, filename, curve_ppT)