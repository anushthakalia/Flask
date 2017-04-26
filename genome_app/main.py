from flask import Flask, render_template,request
import FASTA as fasta
import os
import DNA_operations as RevComp
import networkx as nx
import matplotlib.pyplot as plt


app= Flask(__name__)

@app.route('/strAssembler',methods=['GET','POST'])

def strAssembler():
	if request.method=='POST':
		if request.form.get('cb2'):
			filename=request.form['filename']
			superstring=LongestCommonSuperstring( [fasta[1] for fasta in fasta.ReadFASTA(filename)])
			return render_template("result.html",superstring=superstring)
		elif request.form.get('cb1'):
			text=request.form['reads']
			with open('input.txt','w') as f:
				f.write(text)
			superstring=LongestCommonSuperstring( [fasta[1] for fasta in fasta.ReadFASTA('input.txt')])
			os.remove('input.txt')
			return render_template("result.html",superstring=superstring)
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
            return render_template("resultg.html",superstring=superstring)   
        elif request.form.get('cb1'):
            text=request.form['reads']
            with open('input.txt','w') as f:
                f.write(text)
            superstring=debruijn('input.txt')
            os.remove('input.txt')
            return render_template("resultg.html",superstring=superstring)
        else:
            return render_template("formg.html")

    return render_template("formg.html")
    

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

if __name__=="__main__":
	app.run()
    