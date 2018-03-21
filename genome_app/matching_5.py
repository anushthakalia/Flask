from decimal import *
getcontext().prec = 8

#return mass of given amino acid
def mass_of_aminoacid(aa):
    mass_table = { 'A':71.03711,
                   'C':103.00919,
                   'D':115.02694,
                   'E':129.04259,
                   'F':147.06841,
                   'G':57.02146,
                   'H':137.05891,
                   'I':113.08406,
                   'K':128.09496,
                   'L':113.08406,
                   'M':131.04049,
                   'N':114.04293,
                   'P':97.05276,
                   'Q':128.05858,
                   'R':156.10111,
                   'S':87.03203,
                   'T':101.04768,
                   'U':150.95363,
                   'V':99.06841,
                   'W':186.07931,
                   'Y':163.06333 }

    aa = aa.upper()

    #FASTA notation : ambiguous amino acids
    if 'B' in aa:
        print('Ambiguity: B can be either Asparagine (N) or Aspartic acid (D)!')
        return None
    if 'Z' in aa:
        print('Ambiguity: Z can be either  Glutamine (Q) or Glutamic acid (E)!')
        return None

    mass = 0
    for i in aa:
        try:
            mass += mass_table[i]
        except KeyError:
            print('Error: Could not find a mass for an amino acid %s.' % i)
            return None

    # Return the sum of the monoisotopic masses.
    return mass


#Calculate mass for all possible cuts in peptide
def possible_masses(p):
    masses = []

    for i in range(len(p)):
        masses.append(Decimal(mass_of_aminoacid(p[:i])))
        masses.append(Decimal(mass_of_aminoacid(p[i:])))

    return masses


#Given 2 spectra, use their Minkowski difference to find shift value and maximum multiplicity(hence the shared peaks count)
def compare_spectra_convolution(S1,S2):  #assuming S1 and S2 are lists of floating point integers
    minkowski_difference = dict()
    for x in S1:
        for y in S2:
            z = x-y
            if z in minkowski_difference:
                minkowski_difference[z] += 1
            else:
                minkowski_difference[z] = 1
    max_multiplicity = max(minkowski_difference.values())
    #max_x = list(minkowski_difference.keys())[list(minkowski_difference.values()).index(max_multiplicity)]
    return max_multiplicity

def sec_1(filename1,filename2):
    spectrum_file_path = filename1
    db_file_path = filename2
    max_mul= 0
    max_peptide = ''

    # Read the integer, n, the peptides, and the complete spectrum.
    with open(spectrum_file_path, 'r') as infile:
        spectrum = [Decimal(i) for i in infile.readlines()]
    with open(db_file_path, 'r') as infile1:
        n = int(infile1.readline())
        peptides = [infile1.readline().strip() for i in range(n)]

    # Comapre each peptide to the given spectrum
    for i in peptides:
        mul = compare_spectra_convolution(possible_masses(i), spectrum)
        if mul >= max_mul:
            max_mul = mul
            max_peptide = i

    return max_peptide,max_mul

def main():
    # Variables to hold the largest multiplicity and associated peptide.
    max_peptide,max_mul=sec_1("input_spectrum.txt","database_peptide.txt")
    print('Most similar peptide to our spectrum=')
    print(max_peptide)
    print('No. of shared peaks=', max_mul)


if __name__ == '__main__':
    main()