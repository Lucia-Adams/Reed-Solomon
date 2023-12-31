"""
This file stores the classes we use for Binary Polynomials, which are used for Galois
Field elements and a Galois Field Class and it's associated methods

Binary_Polynomial()
Galois_Field()
"""

class Binary_Polynomial():
    #  This is a polynomial but expressed in binary coefficients

    def __init__(self, strc):
        self.strc = strc # strc =  string coefficient

    def __str__(self):
        return f"Binary_Polynomial {self.strc}"
    
    def binary(self):
        # converts string in base 2 to integer base 10, then convert back to python binary
        return bin(int(self.strc, base=2))
    
    def denary(self):
        return int(self.strc, base=2)

    def degree(self):
        return len(self.strc) - self.strc.index('1') - 1


class Galois_Field():
    #  This is a class to represent a chosen Galois Field

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

        # checks this is true when defining to ensure it is valid
        assert irpol.degree() == gf_size_m


    def generate_gf_elem_table(self):
        """
        Generates table of all the elements with irreducible polynomial of degree m

        table (Int[]) - list of GF(2^m) elements length 2^m
        """
        # Have first two elements of table standard, so need GF-2 more elements 
        table = [0, 1] 
        # Get highest degree term of the irreducible polynomial 
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
        One way to do GF division is by using inverses ie the element multiplied by its 
        inverse produce 1 = alpha^0
        The inverse of alpha^m will be stored at index m+1 in the table 
        (as first element is 0 for non exisitent logarithm)

        inverse_table (Int[]) - list of GF(2^m) inverse elements
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
   
    def mult(self, a,b):
        """
        This multiplies elements by adding their indicies modulo 2^m-1
        a,b (Int) - Galois Field elements to be multiplied together to return another GF element

        returns (Int) result
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

        returns (Int) result
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

    def add(self, a ,b):
        """
        Addition and subtraction and XOR are equivalent in the Galois Field 
        a,b (Int) - Galois Field elements to add

        returns (Int) result
        """
        return a^b