from flask import Flask, render_template,request
import FASTA as fasta
import os
import DNA_operations as RevComp
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import sys
from decimal import *
getcontext().prec = 8


app= Flask(__name__)

@app.route('/')
def index():
    return render_template("homepage.html")

@app.route('/strAssembler',methods=['GET','POST'])

def strAssembler():
	if request.method=='POST':
		if request.form.get('cb2'):
			filename=request.form['filename']
			superstring=LongestCommonSuperstring( [fasta[1] for fasta in fasta.ReadFASTA(filename)])
			return render_template("form.html",superstring=superstring)
		elif request.form.get('cb1'):
			text=request.form['reads']
			with open('input.txt','w') as f:
				f.write(text)
			superstring=LongestCommonSuperstring( [fasta[1] for fasta in fasta.ReadFASTA('input.txt')])
			os.remove('input.txt')
			return render_template("form.html",superstring=superstring)
		else:
			return render_template("form.html")

	return render_template("form.html")

@app.route('/graphAssembler',methods=['GET','POST'])
def graphAssembler():
    if request.method=='POST':
        if request.form.get('cb2'):
            filename=request.form['filename']
            superstring=debruijn(filename)
            make_graph()
            return render_template("formg.html",superstring=superstring)   
        elif request.form.get('cb1'):
            text=request.form['reads']
            with open('input.txt','w') as f:
                f.write(text)
            superstring=debruijn('input.txt')
            os.remove('input.txt')
            return render_template("formg.html",superstring=superstring)
        else:
            return render_template("formg.html")

    return render_template("formg.html")
    
@app.route('/spectrum',methods=['GET','POST'])
def spectrum():
    if request.method=='POST':
        if request.form.get('cb2'):
            filename=request.form['filename']
            with open(filename, 'r') as myfile:
                sequence=myfile.read().replace('\n', '')
            struct = final_spec(sequence)
            superstring = write_structure(struct,sequence)
            
            return render_template("forms.html",superstring=superstring)   
        elif request.form.get('cb1'):
            text=request.form['string']
            with open('input.txt','w') as f:
                f.write(text)
            with open('input.txt', 'r') as myfile:
                sequence=myfile.read().replace('\n', '')
            struct = final_spec(sequence)
            superstring = write_structure(struct,sequence)
            os.remove('input.txt')
            return render_template("forms.html",superstring=superstring)
        else:
            return render_template("forms.html")

    return render_template("forms.html")
    
@app.route('/sec1',methods=['GET','POST'])
def sec1():
    if request.method=='POST':
        if request.form.get('cb2'):
            filename1=request.form['filename1']
            filename2=request.form['filename2']

            max_peptide,max_mul=sec_1(filename1,filename2)
            superstring=[max_peptide,max_mul]
            
            
            return render_template("formr.html",superstring=superstring)   
        elif request.form.get('cb1'):
            text1=request.form['string1']
            text2=request.form['string2']
            with open('input1.txt','w') as f:
                f.write(text1)
            with open('input2.txt','w') as f:
                f.write(text2)
            max_peptide,max_mul=sec_1('input1.txt','input2.txt')
            superstring=[max_peptide,max_mul]
            os.remove('input1.txt')
            os.remove('input2.txt')
            return render_template("formr.html",superstring=superstring)
        else:
            return render_template("formr.html")

    return render_template("formr.html")

@app.route('/sec2',methods=['GET','POST'])
def sec2():
        if request.form.get('cb2'):
            filename=request.form['filename']
            with open(filename, 'r') as infile:
                l = list(map(float, infile.readlines()))
            
            superstring=peptide_from_spectrum(l)
            
            return render_template("formrr.html",superstring=superstring)   
        elif request.form.get('cb1'):
            text=request.form['string']
            with open('input.txt','w') as f:
                f.write(text)
            with open('input.txt', 'r') as infile:
                l = list(map(float, infile.readlines()))
            
            superstring=peptide_from_spectrum(l)
            os.remove('input.txt')
            return render_template("formrr.html",superstring=superstring)
        else:
            return render_template("formrr.html")
        return render_template("formrr.html")

def MergeMaxOverlap(str_list):
    '''Given a list of strings, returns the list of strings with the two strings having maximum overlap merged into a single string.'''
    max_length = -1

    for prefix_index in range(len(str_list)):
        # Don't compare a string with itself.
        for suffix_index in [num for num in range(len(str_list)) if num != prefix_index]:
            # Name the two selected strings for code readability.
            prefix, suffix = str_list[prefix_index], str_list[suffix_index]
            # Begin finding the maximum overlap between the prefix and suffix strings.
            i = 0
            while prefix[i:] != suffix[0:len(prefix[i:])]:
                i += 1
            # Store the overlap length and string indicies if they the longest thus far.
            if len(prefix) - i > max_length:
                max_length = len(prefix) - i
                max_indicies = [prefix_index, suffix_index]

    # Return all strings without maximum overlap, and the merged string with maximum overlap.
    return [str_list[j] for j in range(len(str_list)) if j not in max_indicies] + [str_list[max_indicies[0]] + str_list[max_indicies[1]][max_length:]]


def LongestCommonSuperstring(str_list):
    '''Return the longest common superstring from a list of strings.'''
    # Note: In general, not necessarily the longest but a good approximation at worst.
    # For the Rosalind problem we're given the condition that there exists a unique way to reconstruct the entire chromosome from these reads by gluing 
    # together pairs of reads that overlap by more than half their length, so this method should be optimal.

    # Repeatedly find the merge strings with the maximum overlap until we have one string.
    while len(str_list) > 1:
        str_list = MergeMaxOverlap(str_list)

    return str_list[0]

def debruijn(filename):
    with open(filename) as input_data:
        k_mers = [line.strip() for line in input_data.readlines()]

# Get the edge elements.
    DBG_edge_elmts = set()
    for kmer in k_mers:
        DBG_edge_elmts.add(kmer)
        DBG_edge_elmts.add(RevComp.ReverseComplementDNA(kmer))

# Create the edges.
    k = len(k_mers[0])
    edge = lambda elmt: '('+elmt[0:k-1]+', '+elmt[1:k]+')'
    DBG_edges = [edge(elmt) for elmt in DBG_edge_elmts]

# Write and save the adjacency list.
    l='\n'.join(DBG_edges)
    map = str.maketrans('', '', '(),')
    DBG_edges=[ element.translate(map) for element in DBG_edges]
    with open('DBRU.txt', 'w') as output_file:
        output_file.write('\n'.join(DBG_edges))
        return l


def make_graph():
    ins = open( "DBRU.txt", "r" )
    data = [tuple(str(n) for n in line.split()) for line in ins]
    G = nx.MultiDiGraph(data)
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos,  node_color = 'r',node_size=1000)
    nx.draw_networkx_edges(G, pos,  arrows=True)
    nx.draw_networkx_labels(G,pos)
    plt.savefig('/home/cic/Documents/Flask/genome_app/static/graph.png', format="PNG")



def pair(i,j):
    if ((i=='A') & (j=='U'))|((i=='U') & (j=='A')):
        return 2
    elif ((i=='G') & (j=='C'))|((i=='C') & (j=='G')):
        return 3
    elif ((i=='G') & (j=='U'))|((i=='U') & (j=='G')):
        return 1
    else:
        return 0


def Fill(DP, i, j,sequence):
    structure = []
    
    structure.append(DP[i+1][i+j])
    structure.append(DP[i][i+j-1])
    structure.append(DP[i+1][i+j-1]+pair(sequence[i],sequence[i+j]))
    
    if i+3<=i+j:
        tmp=[]
        for k in range(i+1,j+i):
            tmp.append(DP[i,k]+DP[k+1,i+j]);
        structure.append(max(tmp))

    return max(structure)


def traceback(DP,i,j,pair1,sequence):
  if i<j:
    if DP[i,j]==DP[i+1,j]:
      traceback(DP,i+1,j,pair1,sequence)
    elif DP[i,j]==DP[i,j-1]:
      traceback(DP,i,j-1,pair1,sequence)
    elif DP[i,j]==(DP[i+1,j-1]+pair(sequence[i],sequence[j])):
      pair1.append([i,j,str(sequence[i]),str(sequence[j])])
      traceback(DP, i+1,j-1,pair1,sequence);
    else:
      for k in range(i+1,j):
        if DP[i,j]==DP[i,k]+DP[k+1,j]:
          traceback(DP, i,k,pair1, sequence);
          traceback(DP, k+1,j,pair1, sequence);
          break;
  return pair1;

def write_structure(structure,sequence):
    dot_bracket = ["." for _ in range(len(sequence))]
    #print structure
    for s in structure:
        dot_bracket[min(s)] = "("
        dot_bracket[max(s)] = ")"
    return "".join(dot_bracket)


def final_spec(sequence):
  


  N = len(sequence)
  structure = []

  DP = np.zeros((N,N))
  count=0
  for j in range(1,N):
      for i in range(0,N-j):

          DP[i][i+j]=Fill(DP,i,j,sequence)

  struct=[]
  pair1=traceback(DP,0,N-1,[],sequence)
  #print "max # of folding pairs: ",len(pair1);
  for x in range(0,len(pair1)):
      
      pair1[x][1]=pair1[x][1]-pair1[x][0]-1
      if (pair1[x][1]>=3) :
          #print '%d %d %s==%s' % (pair1[x][0],pair1[x][1]+pair1[x][0]+1,pair1[x][2],pair1[x][3]);
          struct.append((pair1[x][0],pair1[x][1]+pair1[x][0]+1))
  return struct


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
if __name__=="__main__":
	app.run()
    