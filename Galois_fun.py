import math

# from numpy.polynomial import Polynomial as p


#  This is a polnomial but expressed in binary coefficients
class GF_Polynomial():
    # strc =  string coefficient
    def __init__(self, strc):
        self.strc = strc

    def __str__(self):
        return f"Polynomial {self.strc}"
    
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
        self.mult_table = self.generate_gf_elem_table()
        self.inverse_table = self.generate_gf_elem_inverse_table()

        assert irpol.degree() == m


    def generate_gf_elem_table(self):
        """
        Generates table of all the elements with irreducible polynomial of degree m
        """
        # have first two elements of table standard, so need GF-2 more elements 
        table = [0, 1] 
        # This is the highest degree term of the irreducible polynomial is Polynomial('1' + ('0' * irpol.degree())) ie integer is 2^degree = size galois field
        high_deg = 2**(irpol.degree())
        # Holds integer version of rest of the irreducible polynomial ie eveything but it's highest degree
        rest = irpol.denary() - high_deg

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
        units = self.mult_table[1:] 
        #  The element 0 has no inverse so define as 0 here
        inverse_table = [0]

        for i in range(0,gf_size-1):
            # go through each element alpha^i in table and calculate inverse
            inverse_index = (-i) % (gf_size-1) # this is the power of alpha that is the inverse element
            inverse_table.append(units[inverse_index])
        
        return inverse_table

    
    def mult(self, a,b):
        """
        This multiplies elements by adding their indicies modulo 2^m-1
        a,b (Int) - Galois Field elements to be multiplied together to return another GF element
        """
        # get mutliplicative elemnts ie remove 0 from begining so 1 is alpha^0 at index 0
        units = self.mult_table[1:] 

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
        
        units = self.mult_table[1:] 
        
        b_index = units.index(b)
        # b = alpha^m inverse element is stored at index m+1 
        b_inverse = self.inverse_table[b_index+1] 

        return self.mult(a, b_inverse)
    



def gen_poly(field_elem, gen_degree):
    """
    field_elem (Int[]) -  List of generated field elements
    gen_degree (Int) - degree of the generator polynomial ie n-k
    """
    
    gen_roots = field_elem[1:gen_degree+1]
    # g(x) = (x-alpa^0)(x-alpha^2)...(x-alpha^n-k-1)




# global definitions for RS(15, 11) with GF(16) and m is 4 bits
n = 15
k = 11
m = int(math.log(n+1,2))


# Irreducible polynomial is p(x) = x^4 + X + 1 and degree matches m 
#irpol = GF_Polynomial('100011101') for GF(256) ie GF(2^8) other example works!!
irpol = GF_Polynomial('10011')

GF = Galois_Field(m, irpol)
print(GF.mult_table)
print(GF.inverse_table)

print(GF.mult(10,13))

# maybe test another value for this
print(GF.div(11, 10))


