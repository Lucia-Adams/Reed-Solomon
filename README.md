# Reed-Solomon
This is an implementation of the Reed Solomon algorithm in Python.

I have made an accompanying video to explain the algorithm here: https://youtu.be/UkDAMomEvxY

## Reed_Solomon.py
This is where the encoding and decoding functions are, as well as the sub-algorithms such as the
Berlekamp-Massey Algorithm, Forney algorithm etc.

> I have also included a demo function for RS(15,11). This reads in the file `binary.txt` which contains one codewords worth of binary data to exhibit the function working. You can set this to be verbose which will show you the substeps. It also allows you to choose whether the errors are random or not. The non-random errors follow the main example given within the guide by C.K.P Clarke*, should you want to compare my implementation with their explanations.

## GF_and_Polynomials.py
In this file I have implemented two classes
 - `Binary_Polynomial` is for polynomials with binary coefficients
 - `Galois_Field` is to represent the Galois field and I have implemented addition, multiplication and division. (Addition and Subtraction are equivalent in this case as we work modulo 2)

## GF_Polynomials.py
This file includes functions for polynomials with Galois Field element coefficients. I use
the Galois field operations I wrote in `GF_and_Polynomials.py`


*C.K.P Clarke paper on Reed-Solomon: https://downloads.bbc.co.uk/rd/pubs/whp/whp-pdf-files/WHP031.pdf