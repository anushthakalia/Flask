import sys

def mass_to_aa(val, tolerance=0.001):
    ''' Returns the amino acid corresponding to a given mass. '''

    # The monoisotopic masses of each
    aa_table = { 71.03711:'A',
                 103.00919:'C',
                 115.02694:'D',
                 129.04259:'E',
                 147.06841:'F',
                 57.02146:'G',
                 137.05891:'H',
                 113.08406:'I',
                 128.09496:'K',
                 113.08406:'L',
                 131.04049:'M',
                 114.04293:'N',
                 97.05276:'P',
                 128.05858:'Q',
                 156.10111:'R',
                 87.03203:'S',
                 101.04768:'T',
                 150.95363:'U',
                 99.06841:'V',
                 186.07931:'W',
                 163.06333:'Y' }

    # Keep track of the closest match to the given mass. Admittedly this is
    # only useful in certain circumstances...
    closest = ['', 999]

    for mass, aa in aa_table.items():
        diff = abs(val - mass)
        if diff < closest[1]:
            closest = [aa, diff]

        # Return if a match is found.
        if diff < tolerance:
            return aa

    # Print a warning message if no match is found.
    #print('Note: Could not find an amino acid with monoisotopic mass %.5f.' % val)
    #print(' '*6 + 'Closest match is', closest[0], '(mass difference %5f).' % closest[1])


def build_peptide(l, peptide='', aa=0):
    ''' Given a dictionary of fragment masses, with the next highest fragment
    mass and an amino acid representing the gap between them, iterably build a
    peptide by starting with the smallest mass.
    '''
    if aa == 0:
        aa = min(l)

    if aa not in l:
        return peptide
    else:
        for i in l[aa]:
            return build_peptide(l, peptide+i[0], i[1])


def peptide_from_spectrum(l):
    # Create a directed graph of each mass and their possible associated amino
    # acids.
    pairs = {}
    for i in range(len(l)):
        for j in range(i, len(l)):
            aa = mass_to_aa(l[j]-l[i])
            if aa:
                if l[i] in pairs:
                    pairs[l[i]].append((aa, l[j]))
                else:
                    pairs[l[i]] = [(aa, l[j])]

    # Iterably build the peptide starting from the smallest mass.
    peptide = build_peptide(pairs)

    # Return the completed peptide of length n.
    return peptide


def main():
    
    with open('rosalind_sgra.txt', 'r') as infile:
        l = list(map(float, infile.readlines()))

    print(peptide_from_spectrum(l))


if __name__ == '__main__':
    main()