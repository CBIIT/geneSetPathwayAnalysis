import os
import sys 
import re

import pprint 
pp = pprint.PrettyPrinter(indent=4)	

def unique(seq): 
	# order preserving
	checked = []
	for e in seq:
		if e not in checked:
			checked.append(e)
	return checked

with open("pharmgkb_attr.txt") as f:
	lines = f.readlines()
	
pharmgkb_mapping = {}	

for line in lines:
	items = line.rstrip().split("\t")
	if items[0] == "Gene":
		pharmgkb_mapping[items[1]] = [items[2]]
		#DEBUG
		#print items[1] + "\t" + items[2]

#DEBUG
#pp.pprint(pharmgkb_mapping.keys())

gene_sets = []

os.chdir("./pharmgkb_owl")

for file in os.listdir("."):
	cur_gene_set = []
	if file.endswith(".sif"):
		#DEBUG
		#print file
		with open(file) as f:
			lines = f.readlines()
			for line in lines:
					items = line.rstrip().split("\t")
					#DEBUG
					#pp.pprint(items)
					intA = ""
					intB = ""
					match = re.search(r'Protein.(PA[0-9]+)', items[0])
					if match:
						intA = match.group(1)
					match = re.search(r'Protein.(PA[0-9]+)', items[2])
					if match:
						intB = match.group(1)
					#DEBUG
					#print intA + "\t" + intB
					if intA in pharmgkb_mapping.keys():
						# pharmgkb_mapping[intA] returns a list 
						# so grab the first element
						cur_gene_set.append(pharmgkb_mapping[intA][0])
					if intB in pharmgkb_mapping.keys(): 
						cur_gene_set.append(pharmgkb_mapping[intB][0])
		tmp = [file, "http://www.pharmgkb.org/"]
		cur_gene_set = tmp + unique(cur_gene_set)				
	gene_sets.append(cur_gene_set)

for cur_gene_set in gene_sets:
	if cur_gene_set: 
		for item in cur_gene_set: 
			#DEBUG
			#pp.pprint(item)
			if item: 
				sys.stdout.write(item + "\t")
		sys.stdout.write("\n")


	
