#! /usr/bin/env python
"""dynamic programming algorithm to count the number of compatible k-mutants
with a given structure and sequence."""

import sys
PRINT=0

def m1(x,y):
	"""For a nucleotide pair (x, y), let m1(x, y) be the number of single-point
	mutations in (x, y), that keep base pair compatibility."""
	if x=='A' and y=='U': m=1
	elif x=='C' and y=='G': m=1
	elif x=='G' and y=='C': m=1
	elif x=='G' and y=='U': m=2
	elif x=='U' and y=='A': m=1
	elif x=='U' and y=='G': m=2
	else: 
		print "sequence and structure are not compatible!"
		sys.exit(1)
	return m

def m2(x,y):
	"""For a nucleotide pair (x, y), let m2(x, y) be the number of two-point
	mutations in (x, y), that keep base pair compatibility."""
	if x=='A' and y=='U': m=4
	elif x=='C' and y=='G': m=4
	elif x=='G' and y=='C': m=4
	elif x=='G' and y=='U': m=3
	elif x=='U' and y=='A': m=4
	elif x=='U' and y=='G': m=3
	else: 
		print "sequence and structure are not compatible!"
		sys.exit(1)
	return m

def computeNumOfCompatibleStr_draft(seq,s,K):
	"""first implementation of the algorithm-- NOT USED"""
	n=len(s)
	baseD = getBpList(s)
	z = {}
	for i in range(n):
		for j in range(n):
			for r in range(-2,K+1):
				z[(i,j,r)] = 0.0
	for i in range(n):
		z[(i,i,0)] = 1.0
		z[(i,i,1)] = 3.0
	for d in range(1,n):
		for i in range(0,n-d):
			for r in range(K+1):
				j=i+d
				#if (r<=j-i+1):
				if baseD[i]==j:
					z[(i,j,r)] += z[(i+1,j-1,r)] + z[(i+1,j-1,r-1)] * m1(seq[i],seq[j]) + z[(i+1,j-1,r-2)] * m2(seq[i],seq[j])
				elif baseD[j]!=-1 and baseD[j]<j: #s[j]==')'
					z[(i,j,r)] += z[(i,j-1,r)] + z[(i,j-1,r-1)] * m1(seq[baseD[j]],seq[j]) +  z[(i,j-1,r-2)] * m2(seq[baseD[j]],seq[j])
				elif (baseD[i]==-1 and baseD[j]==-1):
					z[(i,j,r)] += 3*z[(i,j-1,r-1)] + z[(i,j-1,r)]
				elif baseD[i]>i:
					z[(i,j,r)] += z[(i+1,j,r)]
				elif baseD[j]!=-1 and baseD[j]>j: #s[j]=='('
					z[(i,j,r)] += z[(i,j-1,r)]
				elif baseD[i]!=-1 and baseD[i]<i: #s[i]==')'
					z[(i,j,r)] += z[(i+1,j,r)] + z[(i+1,j,r-1)] * m1(seq[baseD[i]],seq[i]) +  z[(i+1,j,r-2)] * m2(seq[baseD[i]],seq[i])
				if PRINT: print "%d\t%d\t%d\t%f\n" %(i,j,r,z[(i,j,r)])
	return z[(0,n-1,K)]


def computeNumOfCompatibleStr(seq,s,K):
	"""count the number of compatible K-mutants of seq that are compatible
		with a given structure s"""
	n=len(s)
	baseD = getBpList(s)
	z = {}
	for i in range(n):
		for j in range(n):
			for r in range(-2,K+1):
				z[(i,j,r)] = 0.0
	for i in range(n):
		if s[i] ==')' or s[i] =='(':
			z[(i,i,1)] = 0.0
			z[(i,i,0)] = 1.0
		else:
			z[(i,i,1)] = 3.0
			z[(i,i,0)] = 1.0
	for d in range(1,n):
		for i in range(0,n-d):
			for r in range(K+1):
				j=i+d
				#if (r<=j-i+1):
				#if baseD[i]==j:
					#z[(i,j,r)] += z[(i+1,j-1,r)] + z[(i+1,j-1,r-1)] * m1(seq[i],seq[j]) + z[(i+1,j-1,r-2)] * m2(seq[i],seq[j])
				if baseD[j]!=-1 and baseD[j]<j: #s[j]==')'
					z[(i,j,r)] += z[(i,j-1,r)] + z[(i,j-1,r-1)] * m1(seq[baseD[j]],seq[j]) +  z[(i,j-1,r-2)] * m2(seq[baseD[j]],seq[j])
				elif baseD[j]!=-1 and baseD[j]>j: #s[j]=='('
					z[(i,j,r)] += z[(i,j-1,r)]
				elif (baseD[j]==-1):
					z[(i,j,r)] += 3*z[(i,j-1,r-1)] + z[(i,j-1,r)]
				if PRINT: print "%d\t%d\t%d\t%f\n" %(i,j,r,z[(i,j,r)])
	return z[(0,n-1,K)]
	
def getBpList(s):
	"""given a structure in dot bracket notation return the
	list of its base pairs"""
	
  q = []; basedict={}
  for i in range(len(s)):
	  basedict[i]=-1
  for i in range(len(s)):
    if s[i]=="(":
      q.append(i)
    elif s[i]==")":
      p = q.pop()
      basedict[p]=i;
      basedict[i]=p;
  #print BpList[::-1], pBaseList[::-1]
  return basedict
				
if __name__ == '__main__':
	if len(sys.argv) != 4:
		print """dynamic programming algorithm to count the number of compatible k-mutants
				 with a given structure and sequence
				 
				 Usage: %s sequence structure k
				 
				 Example:./RNACompatible.py  CCCCCCCAAAAAGGGGGGG '(((((((.....)))))))' 5""" %sys.argv[0]
		sys.exit(1)
	seq = sys.argv[1]
	s = sys.argv[2]
	K = int(sys.argv[3])
	if(K>len(s)):
		print "Error: k>len(s)"
		sys.exit(1)
	print computeNumOfCompatibleStr(seq,s,K)
