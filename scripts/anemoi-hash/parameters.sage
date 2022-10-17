#!/usr/bin/sage
# -*- mode: python ; -*-

from sage.all import *
import hashlib
import itertools

from constants import *

load('anemoi.sage')
     
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


