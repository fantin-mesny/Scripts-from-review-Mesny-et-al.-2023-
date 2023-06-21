import sys
import os
import argparse


def get_params(argv):
	parser = argparse.ArgumentParser(description='rename sequences in fasta to match tree branches')
	parser.add_argument('-i', '--i', help="fasta from pal2nal", required=True)
	a = parser.parse_args()
	return a

if __name__ == '__main__':
	a = get_params(sys.argv[1:])

	newfile=''
	with open(a.i,'r') as inp:
		lines=inp.readlines()
	for line in lines:
		if '>' in line:
			newfile+=line.replace('jgi',line.split('|')[1]+'_jgi') #.replace('|','_').replace('.','_').replace('#','_')
		else:
			newfile+=line
	with open(a.i.replace('.pal2nal','.2.pal2nal'),'w+') as outp:
		outp.write(newfile)
		print('output:'+a.i.replace('.pal2nal','.2.pal2nal'))