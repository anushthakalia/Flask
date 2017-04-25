from flask import Flask, render_template,request
import FASTA as fasta
import os


app= Flask(__name__)

@app.route('/send',methods=['GET','POST'])

def send():
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

if __name__=="__main__":
	app.run()