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
		if 'Root' in line:
			newfile+=line.replace('>Root',line.split('_')[0]+'_Root')
		else:
			newfile+=line
	with open(a.i.replace('.pal2nal','.2.pal2nal'),'w+') as outp:
		outp.write(newfile)