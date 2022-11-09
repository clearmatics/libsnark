#!/usr/bin/sage
# -*- mode: python ; -*-

from sage.all import *
import hashlib
import itertools

from constants import *

load('anemoi.sage')

def anemoi_selected_instances():

    # accumulating selected Anemoi instances
    A = []
    
    # 128-bit security level instantiations
    A_BLS_12_381_SCALARFIELD_1_COL_128_BITS = AnemoiPermutation(
        q=BLS12_381_SCALARFIELD,
        n_cols=1,
        security_level=128
    )
    A.append(
        ("A_BLS_12_381_SCALARFIELD_1_COL_128_BITS",
         A_BLS_12_381_SCALARFIELD_1_COL_128_BITS))
    A_BLS_12_381_SCALARFIELD_4_COL_128_BITS =AnemoiPermutation(
        q=BLS12_381_SCALARFIELD,
        n_cols=4,
        security_level=128
    )
    A.append(
        ("A_BLS_12_381_SCALARFIELD_4_COL_128_BITS",
         A_BLS_12_381_SCALARFIELD_4_COL_128_BITS))
    A_BLS_12_381_SCALARFIELD_6_COL_128_BITS =AnemoiPermutation(
        q=BLS12_381_SCALARFIELD,
        mat=matrix.circulant([1, 1, 3, 4, 5, 6]),
        n_cols=6,
        security_level=128)
    A.append(
        ("A_BLS_12_381_SCALARFIELD_6_COL_128_BITS",
         A_BLS_12_381_SCALARFIELD_6_COL_128_BITS))
    return A
    
def output_parameters():
    instances = anemoi_selected_instances()    
    for i in range(len(instances)):
#    for i in range(1):
        A_str = instances[i][0]
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
    return instances
                  
if __name__ == "__main__":
    A = output_parameters()


