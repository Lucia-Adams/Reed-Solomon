import math

# from numpy.polynomial import Polynomial as p


#  This is a polnomial but expressed in binary coefficients
class GF_Polynomial():
    # strc =  string coefficient
    def __init__(self, strc):
        self.strc = strc

    def __str__(self):
        return f"GF_Polynomial {self.strc}"
    
    def binary(self):
        # converts string in base 2 to integer base 10, then convert back to python binary
        return bin(int(self.strc, base=2))
    
    def denary(self):
        return int(self.strc, base=2)
    
    def string(self):
        return self.strc

    def degree(self):
        return len(self.strc) - self.strc.index('1') - 1



class Galois_Field():

    def __init__(self, gf_size_m, irpol):
        """
        gf_size_m (Int) - m for the Galois Field denoted GF(2^m) 
        irpol (Polynomial) -  irreducible polynomial of degree m
        """
        self.size_m = gf_size_m
        self.size = 2**gf_size_m
        self.irpol = irpol
        self.elem_table = self.generate_gf_elem_table()
        self.inverse_table = self.generate_gf_elem_inverse_table()

        assert irpol.degree() == gf_size_m


    def generate_gf_elem_table(self):
        """
        Generates table of all the elements with irreducible polynomial of degree m
        """
        # have first two elements of table standard, so need GF-2 more elements 
        table = [0, 1] 
        # This is the highest degree term of the irreducible polynomial is Polynomial('1' + ('0' * irpol.degree())) ie integer is 2^degree = size galois field
        high_deg = 2**(self.irpol.degree())
        # Holds integer version of rest of the irreducible polynomial ie eveything but it's highest degree
        rest = self.irpol.denary() - high_deg

        for i in range(high_deg-2):
            last_alpha = table[-1]
            new_alpha = table[-1] << 1  # multiply by 2

            if new_alpha & high_deg:
                # has overflown too high so subtract off 
                new_alpha -= high_deg
                # xor with rest of polynomial which is adding back what we took off as xor is addition
                new_alpha = new_alpha ^ rest 
        
            table.append(new_alpha)

        return table

    
    def generate_gf_elem_inverse_table(self):
        """
        One way to do GF division is by using inverses ie the element multiplied by its inverse produce 1 = alpha^0
        gf_elem_table (Int[]) - list of GF(2^m) elements length 2^m
        The inverse of alpha^m will be stored at index m+1 in the table (as first element is 0 for non exisitent logarithm)
        """

        gf_size = self.size
        units = self.elem_table[1:] 
        #  The element 0 has no inverse so define as 0 here
        inverse_table = [0]

        for i in range(0,gf_size-1):
            # go through each element alpha^i in table and calculate inverse
            inverse_index = (-i) % (gf_size-1) # this is the power of alpha that is the inverse element
            inverse_table.append(units[inverse_index])
        
        return inverse_table
    
    def add(self, a ,b):
        """
        Addition and subtarction and XOR are equivalent in the Galois Field 
        a,b (Int) - Galois Field elements to add
        """
        return a^b

    
    def mult(self, a,b):
        """
        This multiplies elements by adding their indicies modulo 2^m-1
        a,b (Int) - Galois Field elements to be multiplied together to return another GF element
        """
        # get mutliplicative elemnts ie remove 0 from begining so 1 is alpha^0 at index 0
        units = self.elem_table[1:] 

        if a==0 or b==0:
            return 0
        
        a_index = units.index(a)
        b_index = units.index(b)
        result_index = (a_index + b_index) % (self.size-1)

        return units[result_index]
    
    def div(self, a, b):
        """
        This divides elements by multiplying by inverses in the GF
        a,b (Int) - Galois Field elements to perform a / b
        """

        if a== 0 :
            return 0
        if b == 0:
            # Can't do this divide by 0!!
            return -1 
        
        units = self.elem_table[1:] 
        
        b_index = units.index(b)
        # b = alpha^m inverse element is stored at index m+1 
        b_inverse = self.inverse_table[b_index+1] 

        return self.mult(a, b_inverse)


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

def encode(RS_n, RS_k , RS_m, irreducible_p, message):
    """
    Returns the message encoded by RS(n,k) with chosen irreducible polynomial

    RS_n (Int) - RS codeword length
    RS_k (Int) - RS original message length
    RS_m (Int) - Symbol length in bits 
    irreducible_p (Str) - Chosen irreducible polynomial for the Galois Field
    
    message (Int[]) - List of integers of message (for now might convert from binary later)
    """

    irpol = GF_Polynomial(irreducible_p)
    GF = Galois_Field(RS_m, irpol)
    #print(f"Element table: {GF.elem_table}")
    
    # Creates generator polnomial
    gx = gen_poly(GF, (RS_n - RS_k)) 
    print(f"Generator polynomial {gx}")

    # We have message 
    mx = message.copy()
    mx.extend([0 for i in range(RS_n - RS_k)]) # this multiplies mx by x^(n-k)

    # Remainder of mx.x^(n-k) / gx and of degree (n-k)
    rx = poly_divide_r(GF, mx, gx)
    print(f"Remainder of mx.x^(n-k) / gx is {rx}")

    message.extend(rx) # equivalent to adding remainder to mx.x^(n-k)

    print(f"Encoded message : {message}")
    return message

    
def check_syndromes(rx, GF, gen_degree):
    """
    Check the syndromes of the transimitted message to find the errors

    rx (Int[]) - the recieved message to decode
    GF (Galois_Field) -  Galois Field we are working in
    gen_degree (Int) - degree of the generator polynomial ie n-k
    
    """

    # We evaluate rx at each alpha^i from the generator polynomial




# global definitions for RS(15, 11) with GF(16) and m is 4 bits
n = 15
k = 11
m = int(math.log(n+1,2))

message = [1,2,3,4,5,6,7,8,9,10,11]


n_DVB = 255
k_DVB = 239
m_DVB = int(math.log(n_DVB+1,2))



encode(n, k , m , '10011', message)

#encode(n_DVB, k_DVB ,m_DVB , '100011101')


