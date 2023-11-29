"""
This is a module to store functions related to Polynomial arithmetic
where the coefficients are Galois Field Elements 

add_polys, mult_polys, poly_divide_r, remove_zeros
"""

from GF_and_Polynomials import *

def add_polys(GF, ax, bx):
    """
    Adds the polynomial ax and bx together

    GF (Galois_Field) -  Galois Field we are working in
    ax (Int[]) - first polynomial
    bx (Int[]) - second polynomial

    result (Int[]) - resultant polynomial
    """
    result = max([ax, bx], key=len)

    small_degree = min(len(ax), len(bx))
    for i in range(1, small_degree+1):
        result[-i] = GF.add(ax[-i], bx[-i])

    return result

def mult_polys(GF, ax, bx):
    """
    Multiplies the polynomial ax and bx together

    GF (Galois_Field) -  Galois Field we are working in
    ax (Int[]) - first polynomial
    bx (Int[]) - second polynomial

    prod (Int[]) - resultant product polynomial
    """

    len_ax = len(ax)
    len_bx = len(bx)
    prod = [0] * (len_ax + len_bx -1)

    for i in range(len_ax):
        for j in range(len_bx):
            to_add = GF.mult(ax[i], bx[j])
            prod[i+j] = GF.add(prod[i+j], to_add)

    return prod


def poly_divide_r(GF, mx, gx):
    """
    This takes the polynomial mx and divides it by polynomial gx to find the remainder rx
    
    GF (Galois_Field) -  Galois Field we are working in
    mx (Int[]) - Polynomial to be divided
    gx (Int[]) - Polynomial to divide by

    remainder (Int[]) - resultant remainder polynomial
    """
    remainder= mx
    gx_len = len(gx)

    # ie while we can keep dividing
    # mulitplies gx to make it same as most significant term as remainder in last step
    # then subtract to form new remainder
    while len(remainder) >= gx_len:
        # multiply gx by highest intermediate last term - we start with mx 
        to_subtract = [GF.mult(i, remainder[0]) for i in gx]

        for i in range(gx_len):
            remainder[i] = GF.add(remainder[i], to_subtract[i]) # as addition same as subtraction
        remainder = remove_zeros(remainder)
    
    return remainder

def remove_zeros(poly):
    """
    Take polynomial with leading 0 terms and removes so in simplist form
    poly (Int[]) - list of coefficients

    poly (Int[]) - resultant polynomial
    """
    non_zero = 0
    for i in poly:
        if i==0:
            non_zero +=1
        else:
            break
            
    return poly[non_zero:]