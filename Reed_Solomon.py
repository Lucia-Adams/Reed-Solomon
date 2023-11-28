import math
from GF_and_Polynomials import *


def gen_poly(GF, gen_degree):
    """
    GF (Galois_Field) -  Galois Field we are working in
    gen_degree (Int) - degree of the generator polynomial ie n-k
    
    result (Int[]) - list of generator polynomial with coefficnets as GF elements
    """
    gen_roots = GF.elem_table[1:gen_degree+1]

    # g(x) = (x-alpa^0)(x-alpha^2)...(x-alpha^n-k-1)
    # ie g(x) = (x+alpa^0)(x+alpha^2)...(x+alpha^n-k-1) as addition is subtratcion in GF
    # need to expand this but using the field multiplication

    # Take first root to get (x-alpha^0)
    result = [1, gen_roots.pop(0)]

    for r in gen_roots:
        mult_by_root = [GF.mult(i, r) for i in result] # multiplies result what we had by root
        result.append(0) # what we had multiplies eveything by x

        # go through from low to highest degree and add mult-by root to result 
        for i in range(1,len(mult_by_root)+1):
            result[-i] = GF.add(result[-i], mult_by_root[-i])

    return result

def remove_zeros(poly):
    """
    Take polynomial with leading 0 terms and removes so in simplist form
    poly (Int[]) - list of coefficients
    """
    non_zero = 0
    for i in poly:
        if i==0:
            non_zero +=1
        else:
            break
            
    return poly[non_zero:]

def add_polys(GF, ax, bx):
    """
    Adds the polynomial ax and bx together

    GF (Galois_Field) -  Galois Field we are working in
    ax (Int[]) - first polynomial
    bx (Int[]) - second polynomial
    """
    result = max([ax, bx], key=len)

    small_degree = min(len(ax), len(bx))
    for i in range(1, small_degree+1):
        result[-i] = GF.add(ax[-i], bx[-i])

    return result

def mult_polys(GF, ax, bx):

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
    


def message_to_int_list():
    """
    This might take say string of binary which should be mbits x message length k in length
    Then divides this up into mbits and makes a lits of ints 
    """
    None


# IMPLEMENTING THE ALGORITHM

def encode(RS_n, RS_k , RS_m, irreducible_p, tx):
    """
    Returns the message encoded by RS(n,k) with chosen irreducible polynomial

    RS_n (Int) - RS codeword length
    RS_k (Int) - RS original message length
    RS_m (Int) - Symbol length in bits 
    irreducible_p (Str) - Chosen irreducible polynomial for the Galois Field
    
    tx (Int[]) - List of integers of tx message (for now might convert from binary later)
    """

    irpol = Binary_Polynomial(irreducible_p)
    GF = Galois_Field(RS_m, irpol)
    #print(f"Element table: {GF.elem_table}")
    
    # Creates generator polnomial
    gx = gen_poly(GF, (RS_n - RS_k)) 
    print(f"Generator polynomial {gx}")

    # We have message 
    mx = tx.copy()
    mx.extend([0 for i in range(RS_n - RS_k)]) # this multiplies mx by x^(n-k)

    # Remainder of mx.x^(n-k) / gx and of degree (n-k)
    rx = poly_divide_r(GF, mx, gx)
    print(f"Remainder of mx.x^(n-k) / gx is {rx}")

    tx.extend(rx) # equivalent to adding remainder to mx.x^(n-k)

    print(f"Encoded message : {tx}")
    return tx


def GF_horner(GF, x, poly):
    """
    Evaluate polynomials using horners method within the Galois Field
    uses say 2x^3 - 6x^2 + 2x - 1 = ((2x - 6)x + 2)x - 1

    GF (Galois_Field) -  Galois Field we are working in
    x (Int) - Galois Field element we are evaluating at
    poly (Int[]) - The pollynomial we are evaluating
    """

    value = poly[0] # if degree 0 then this is the result

    for i in range(1, len(poly)):
        value = GF.mult(value, x) # multiply by what we evaluating at
        value = GF.add(value, poly[i])
    
    return value

    
def check_syndromes(GF, rx, gen_degree):
    """
    Check the syndromes of the transimitted message to find the errors

    GF (Galois_Field) -  Galois Field we are working in
    rx (Int[]) - the recieved message to decode
    gen_degree (Int) - degree of the generator polynomial ie n-k

    syndromes (Int[]) - list of integer syndromes for each root
    """

    # We evaluate rx at each alpha^i from the generator polynomial
    gen_roots = GF.elem_table[1:gen_degree+1]
    syndromes = []

    for root in gen_roots:
        syndromes.append(GF_horner(GF, root, rx))
 
    return syndromes

def berlekamp_massey(GF, syndromes, added_bits):
    """
    This finds the error locator polynomial - a solution to the linear equations of the syndromes 
    GF (Galois_Field) -  Galois Field we are working in
    syndromes (Int[]) - list of integer syndromes 
    added_bits (Int) - the extra bits added to make codeword ie 2t=n-k
    
    vx (Int[]) - the error locator polynomial
    """
    K = 1 # step parameter - used to track which syndome we are up to
    L = 0 # tracks order of equations
    vx = [1] # error locator polynomial - set to 1
    cx = [1,0] # correction polynomial - set to x

    while (K <= added_bits):
        #  print(f"\nStatus: K:{K} L:{L} cx:{cx} vx:{vx}")

        e = syndromes[K-1]
        for i in range(1, L+1):
            # note v(x) in literature indicies labelled 0 as lowest power element 
            # hence the -i-1 on the vx as 'going backwards'
            e = GF.add(e, GF.mult(syndromes[K-1-i], vx[-i-1]))

        # if e = 0, assumes no error and skips to increment K
        if e != 0:
            # new approximation for error locator polynomial - new_vx = vx + (e * cx)
            new_vx = [GF.mult(e,i) for i in cx]
            new_vx = add_polys(GF, vx, new_vx)

            if ((2*L) < K):
                L = K-L # limits L to number of syndromes we have
                cx = [GF.div(i, e) for i in vx] # divide cx by e 
            vx = new_vx

        cx.append(0) # equivalent to mulitplying by x
        K += 1
    
    # print(f"\nStep final : {K} {L} {cx} {vx}")
    return vx

def error_mag_poly(GF, vx, sx, added_bits):
    """
    This finds the error magnitude polynomial which is calcluated
    by using the error locator polynomial and the syndrome polynomial
    This is often denoted as capital omega in literature

    GF (Galois_Field) -  Galois Field we are working in
    vx (Int[]) - the error locator polynomial
    sx (Int[]) - the syndrome polynomial
    added_bits (Int) - the extra bits added to make codeword ie 2t=n-k

    """
    # calcualated via first multiplying syndrome polynomial and error locator
    qx = mult_polys(GF, sx, vx)

    # then the error magnitude polynomial is this modulo x^2t
    # ie remove anything of or past the power x^2t
    qx = qx[-added_bits:]
    qx = remove_zeros(qx)  
    return qx

def find_inv_roots(GF, vx):
    """
    The inverse of the roots of the berlekamp massey error locator polynomial
    provide the error locator numbers
    This goes through and finds these inverses of the roots 

    GF (Galois_Field) -  Galois Field we are working in
    vx (Int[]) - the error locator polynomial

    roots (Int[]) - List of roots of the polynomial
    """
    inverses = GF.inverse_table[1:]
    roots = []

    # loops through table of inverses to see what evaluates to 0 and is thus a root
    for i in range(len(inverses)):
        val = GF_horner(GF, inverses[i], vx)
        if val == 0:
            # the value i is the power of the alpha element to which the 
            # value stored at it is an inverse hence can just add i
            roots.append(i)

    return roots

def forney_algorithm(GF, error_locations, vx, qx):
    """
    This uses the error locator polynomial and the error value polynomial
    to calulate the error value
    The forney algorithm is actually based on lagrange interpolation!
    """

    # ecah error location is the power of alpha that is Xj in literarture
    for l in error_locations:
        l_elem = GF.elem_table[l+1]
        l_inv = GF.inverse_table[l+1]

        # evaluate qx at each error locators inverse
        qx_eval = GF_horner(GF, l_inv, qx)

        # evaluate derivative of vx at each error locators inverse
        # this is equivalent to setting even powers of vx to zero and 
        # then dividing by x = error locators inverse
        vx_deriv_eval = vx.copy()
        for i in range(1, len(vx)+1, 2):
            vx_deriv_eval[-i] = 0 # sets even powers to 0
        vx_deriv_eval.pop() # divides by x
        vx_deriv_eval = remove_zeros(vx_deriv_eval)
        vx_deriv_eval = GF_horner(GF, l_inv, vx_deriv_eval)
        
        error = GF.div(qx_eval, vx_deriv_eval)
        error = GF.mult(error, l_elem)
        
        print(error)



def decode(RS_n, RS_k , RS_m, irreducible_p, rx):
    """
    Returns the message decoded from RS(n,k) with chosen irreducible polynomial

    RS_n (Int) - RS codeword length
    RS_k (Int) - RS original message length
    RS_m (Int) - Symbol length in bits 
    irreducible_p (Str) - Chosen irreducible polynomial for the Galois Field
    
    rx (Int[]) - List of integers of rx message (for now might convert from binary later)
    """

    irpol = Binary_Polynomial(irreducible_p)
    GF = Galois_Field(RS_m, irpol)

    syndromes = check_syndromes(GF, rx, (RS_n-RS_k))
    print(f"Syndromes are: {syndromes}")

    # Calculate error-locator polynomial using berlekamp
    vx = berlekamp_massey(GF, syndromes, (RS_n-RS_k))
    print(f"Error locator polynomial: {vx}")

    #  Now solve the error locator polynomial 
    error_locations = find_inv_roots(GF, vx)
    print(f"Error locations are the powers: {error_locations}")

    sx = syndromes.copy()
    sx.reverse() # syndrome polynomial has coefficients of the syndromes
    qx = error_mag_poly(GF, vx, sx , (RS_n-RS_k))
    print(f"Error magnitude polynomial: {qx}")

    forney_algorithm(GF, error_locations, vx, qx)
    
    


# global definitions for RS(15, 11) with GF(16) and m is 4 bits
n = 15
k = 11
m = int(math.log(n+1,2))

message = [1,2,3,4,5,6,7,8,9,10,11]


n_DVB = 255
k_DVB = 239
m_DVB = int(math.log(n_DVB+1,2))



rx = encode(n, k , m , '10011', message)
# should perturb 
bleh = [1, 2, 3, 4, 5, 11, 7, 8, 9, 10, 11, 3, 1, 12, 12]
decode(n, k , m , '10011', bleh)



#encode(n_DVB, k_DVB ,m_DVB , '100011101')


