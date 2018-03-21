import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys


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

def main():
  with open('spect.txt', 'r') as myfile:
    sequence=myfile.read().replace('\n', '')

  struct = final_spec(sequence)
  dot_brac_notation = write_structure(struct,sequence)
  print(dot_brac_notation)

if __name__=="__main__":
  main()