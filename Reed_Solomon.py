import math
import random
import textwrap
from GF_and_Polynomials import *
from GF_Polynomials import *

def gen_poly(GF, gen_degree):
    """
    GF (Galois_Field) -  Galois Field we are working in
    gen_degree (Int) - degree of the generator polynomial ie n-k
    
    result (Int[]) - list of generator polynomial with coefficients as GF elements
    """
    gen_roots = GF.elem_table[1:gen_degree+1]

    # g(x) = (x-alpa^0)(x-alpha^2)...(x-alpha^n-k-1)
    # ie g(x) = (x+alpa^0)(x+alpha^2)...(x+alpha^n-k-1) as addition is subtraction in GF
    # need to expand this but using the field multiplication
    # Take first root to get (x-alpha^0)
    result = [1, gen_roots.pop(0)]

    for r in gen_roots:
        mult_by_root = [GF.mult(i, r) for i in result] # multiplies result by root
        result.append(0) # multiplies result by x

        # go through from low to highest degree and add mult-by root to result 
        for i in range(1,len(mult_by_root)+1):
            result[-i] = GF.add(result[-i], mult_by_root[-i])

    return result

def GF_horner(GF, x, poly):
    """
    Evaluate polynomials using horners method within the Galois Field
    uses say 2x^3 - 6x^2 + 2x - 1 = ((2x - 6)x + 2)x - 1

    GF (Galois_Field) -  Galois Field we are working in
    x (Int) - Galois Field element we are evaluating at
    poly (Int[]) - The polynomial we are evaluating
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

        e = syndromes[K-1]
        for i in range(1, L+1):
            # note v(x) in literature indices labelled 0 as lowest power element 
            # hence the -i-1 on the vx as 'going backwards'
            e = GF.add(e, GF.mult(syndromes[K-1-i], vx[-i-1]))

        # if e = 0, assumes no error and skips to increment K
        if e != 0:
            # new approximation for error locator polynomial - new_vx = vx + (e * cx)
            new_vx = [GF.mult(e,i) for i in cx]
            new_vx = add_polys(GF, vx, new_vx)

            if ((2*L) < K):
                L = K-L # limits L to number of syndromes we have
                cx = [GF.div(i, e) for i in vx] # divide cx by error e
            vx = new_vx

        cx.append(0) # equivalent to mulitplying by x
        K += 1
    
    return vx

def error_mag_poly(GF, vx, sx, added_bits):
    """
    This finds the error magnitude polynomial which is calculated
    by using the error locator polynomial and the syndrome polynomial
    This is often denoted as capital omega in literature

    GF (Galois_Field) -  Galois Field we are working in
    vx (Int[]) - the error locator polynomial
    sx (Int[]) - the syndrome polynomial
    added_bits (Int) - the extra bits added to make codeword ie 2t=n-k

    qx (Int[]) - the error magnitude polynomial
    """
    # calculated via first multiplying syndrome polynomial and error locator
    qx = mult_polys(GF, sx, vx)

    # then the error magnitude polynomial is this modulo x^2t
    # ie remove anything of or past the power x^2t
    qx = qx[-added_bits:]
    qx = remove_zeros(qx)  
    return qx

def find_inv_roots(GF, vx):
    """
    The inverse of the roots of the Berlekamp Massey error locator polynomial
    provide the error locator numbers
    This goes through and finds these roots and their inverses

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

    GF (Galois_Field) -  Galois Field we are working in
    error_locations (Int[]) - powers in codeword polynomial at which the errors are located
    vx (Int[]) - the error locator polynomial
    qx (Int[]) - the error magnitude polynomial

    error_values (Int[]) - list of integer error values
    """
    error_values = []

    # each error location is the power of alpha that is Xj in literarture
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
        
        error_values.append(error)
    
    return error_values


def encode(RS_n, RS_k , RS_m, irreducible_p, tx, verbose):
    """
    Returns the message encoded by RS(n,k) with chosen irreducible polynomial

    RS_n (Int) - RS codeword length
    RS_k (Int) - RS original message length
    RS_m (Int) - Symbol length in bits 
    irreducible_p (Str) - Chosen irreducible polynomial for the Galois Field
    verbose (Boolean) - if set to True then will print out intermediate results
    
    tx (Int[]) - List of integers of tx message 
    """
    verbose_string = "\n -- ENCODING -- \n"

    irpol = Binary_Polynomial(irreducible_p)
    GF = Galois_Field(RS_m, irpol)
    
    # Creates generator polynomial
    gx = gen_poly(GF, (RS_n - RS_k)) 
    verbose_string += (f"Generator polynomial {gx}\n")

    mx = tx.copy()
    mx.extend([0 for i in range(RS_n - RS_k)]) # this multiplies mx by x^(n-k)

    # Remainder of mx.x^(n-k) / gx and of degree (n-k)
    rx = poly_divide_r(GF, mx, gx)
    verbose_string +=(f"Remainder of mx.x^(n-k) / gx is {rx}\n")

    tx.extend(rx) # equivalent to adding remainder to mx.x^(n-k)

    if verbose:
        print(verbose_string)

    return tx


def decode(RS_n, RS_k , RS_m, irreducible_p, rx, verbose):
    """
    Returns the message decoded from RS(n,k) with chosen irreducible polynomial

    RS_n (Int) - RS codeword length
    RS_k (Int) - RS original message length
    RS_m (Int) - Symbol length in bits 
    irreducible_p (Str) - Chosen irreducible polynomial for the Galois Field  
    rx (Int[]) - List of integers of rx message (for now might convert from binary later)
    verbose (Boolean) - if set to true then will print out intermediate results

    corrected (Int[]) - returns the original message with any errors corrected
    """
    verb_string = "\n -- DECODING -- \n"

    irpol = Binary_Polynomial(irreducible_p)
    GF = Galois_Field(RS_m, irpol)

    syndromes = check_syndromes(GF, rx, (RS_n-RS_k))
    verb_string += (f"Syndromes are: {syndromes}\n")

    # Calculate error-locator polynomial using Berlekamp
    vx = berlekamp_massey(GF, syndromes, (RS_n-RS_k))
    verb_string += (f"Error locator polynomial: {vx}\n")

    #  Now solve the error locator polynomial 
    #  this finds the powers in the codeword at which the errors occur
    error_locations = find_inv_roots(GF, vx)
    verb_string += (f"Error locations are at the message polynomial powers: {error_locations}\n")

    sx = syndromes.copy()
    sx.reverse() # syndrome polynomial has coefficients of the syndromes
    qx = error_mag_poly(GF, vx, sx , (RS_n-RS_k))
    verb_string += (f"Error magnitude polynomial: {qx}")

    errors = forney_algorithm(GF, error_locations, vx, qx)
    
    # Now we correct the errors!
    # Just add (=subtract) the errors from the recieved message
    for i in range(len(error_locations)):
        place = error_locations[i]
        rx[-place-1] = GF.add(rx[-place-1], errors[i])
    # rx is now the corrected polynomial!

    if verbose:
        print(verb_string)
     
    corrected = rx[:-(RS_n-RS_k)]
    return corrected


def bin_to_int_list(bin_file, sym_bits):
    """
    This takes a file of a binary number, then splits that binary into 
    symbols defined by sym_bits and puts the integer value of these chunks into a list
    This is used for demonstration purposes that the algorithm works for any binary data

    mes_str (String) - message string to convert
    sym_bits (Int) - size of symbol 

    mes_int (Int[]) - message as a list of ints of sym_bits bits
    """

    with open(bin_file, 'r') as file:
        bin_string = file.read()

    padding = (sym_bits - (len(bin_string) % sym_bits)) % sym_bits
    bin_string += '0'*padding
    
    symbol_str = textwrap.wrap(bin_string, 4)
    mes_int = [int(symb, base=2) for symb in symbol_str]

    return mes_int


def demo(verbose=True, rand=True):
    """
    This is to demonstrate using the functions with an example using 
    the file "binary.txt"

    verbose (Boolean) - If true, gives extra details on inner computations 
    rand (Boolean) - If set to True randomises the errors, if False uses set errors 
    """

    print("\nREED SOLOMON DEMONSTRATION")
    # Here we set our n,k and m parameters
    n = 15
    k = 11
    m = int(math.log(n+1,2))

    # demonstrate with the binary in binary.txt 
    message = bin_to_int_list("binary.txt", 4)
    # For the sake of this demo, we run RS once for one codeword/packet so
    # we make the length of the message to be the length of the RS codeword
    assert len(message) == k

    print(f"This demo uses the binary string of length k*m from binary.txt")
    print(f"Using RS(n={n},k={k}) encoding with m={m} bits\n")
    print(f"The message to send is: \n{message} as {m} bit integers")

    # encodes message using RS(n,k) = RS(15,11)
    # with irreducible polynomial x^4 + x + 1 ie '10011'
    rx = encode(n, k , m , '10011', message, verbose)
    print(f"The encoded message to send is: {rx}\n")

    print("-- TRANSMISSION -- \n")
    # Here we are simulating errors in transmission
    if rand:
        print("(Randomising errors)")
        max_errors = (n-k) // 2
        for i in range(max_errors):
            error_place = random.randint(1, n-1)
            error = random.randint(1, n)
            rx[error_place] = error    
    else:
        rx[5] = 11
        rx[12] = 1

    print(f"In 'transmission' this changes to : {rx}")
    corrected = decode(n, k , m , '10011', rx, verbose)
    print(f"\nThis is then corrected back to : {corrected}\n")


if __name__ == "__main__":
    # Here we run the demo! 
    demo()