
import os
import glob

os.chdir('/home/ubuntu/mytest')


store=[]
store={}
# store={'ShortRead':[]}

for f in glob.glob('/home/ubuntu/AAVengeR/*.R'):
	with open(f,'r') as infile:
		for l in infile:

			if l.startswith('library(') == True:
				l=l.replace('library(','')
				l=l.replace(')\n','')

				if l not in store.keys():
					store[l] = []

				store[l].append(f)

print(store)

del store['rtracklayer']


with open('getRversion.R','w') as outfile:
	for i in store:
		outfile.write(f'library({i})\n')
	os.system('print(sessionInfo())')

os.system('head -n 100 getRversion.R')

os.system('Rscript getRversion.R')

print(set(store))
