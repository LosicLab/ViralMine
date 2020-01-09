#!/bin/python

# Read through GT_Scores files in the BLAST window output and compile fractional levels of 
# HPV Genotypes for each patient. Can use a "bootstrapping" contig-size normalization method, or
# can use a 'rigorous' mode that uses raw bitscores. The rigorous mode is the default, as BLAST bitscore
# takes into account sequence length already. The bootstrapping should be useful if you want to treat 
# all contigs (no matter size) as equal in determining viral genotype.


import argparse as ap
import re

def get_arguments():
	parser = ap.ArgumentParser(description="Fractional HBV genotype qunatifier by contig bitscore matches")
	parser.add_argument("-f", help="HBV Window Blast output file (should be concatenated list of contig scores if multiple samples assessed)", required=True, type=str)
	parser.add_argument("-method", help="Mixed genotype calculation method, either 'rig' or 'boot'. This will either normalize bitscores before fractional scoring or not. Default is 'rig'.", required=False, type=str, default='rig', choices=['rig','boot'])
	parser.add_argument("-threshold", help="Fraction of total bitscore value that a genotype must reach to be qualified as a coinfection type. Default is 0.1 (10%).", required=False, type=float, default=0.1)
	return parser.parse_args()
args = get_arguments()

##### Bootstrapping option #####
def gt_normalization(cscores, contig_len):  #accepts a dict of gt scores per contig:

	gts = cscores.split('\n')[:-1]
	gt_scores = dict()
	for gt in gts:
		if gt != '':
			if re.findall(" [1-9]",gt) != []:
				gt_scores[gt[0:5]] = float(gt.split(' ')[1])
			else:
				gt_scores[gt[0:5]] = 0
		else:
			pass

	nscores = dict()
	for gt in gt_scores:
		nscores[gt] = gt_scores[gt]/contig_len

	return nscores


##### Rigorous option #####
def gt_frac_score(cscores):

	gts = cscores.split('\n')[:-1]
	gt_scores = dict()
	#print(gts)
	for gt in gts:
		if gt != '':
			#print(gt)
			if re.findall(" [1-9]",gt) != []:
				gt_scores[gt[0:5]] = float(gt.split(' ')[1])
			else:
				gt_scores[gt[0:5]] = 0
		else:
			pass

	#print(gt_scores)
	return gt_scores	

######################
### Main function: ###
######################

file = open(args.f, 'r') #This is a concatenated file of all GT_Scores files
method = args.method
threshold = args.threshold
line = file.readline()

contigs = []
while line:
	contig = []
	contig.append(line)
	line = file.readline()
	while line.startswith(">") != True and line != "":
		contig.append(line)
		line=file.readline()

	contig = ''.join(contig)
	contigs.append(contig)

file.close()

patients = set([x.split('_')[0][2:] for x in contigs])
#print(patients)

viral_gtness = dict()
for pat in patients:
	# pull out all contigs that match patient ID:
	pat_contig_set = [ctg for ctg in contigs if ctg.find(pat) != -1]
	# Now start contig by contig analysis:
	tot_nscores = {"HPV16":0,"HPV18":0,"HPV33":0,"HPV45":0,"HPV31":0,"HPV58":0,"HPV52":0,"HPV35":0,"HPV59":0,"HPV56":0,"HPV51":0,"HPV39":0,"HPV73":0,"HPV68":0,"HPV82":0}
	#print(pat_contig_set)

	if method is not 'rig':
		for ctg in pat_contig_set:
			clen = int(ctg.split('length: ')[1].split('\n')[0])
			cscores = '\n'.join(ctg.split('\n')[1:])
			nscores = gt_normalization(cscores,clen)
			for key in nscores:
				tot_nscores[key] = (tot_nscores[key]) + nscores[key]
			#print(nscores)

	elif method is 'rig':
		for ctg in pat_contig_set:
			cscores = '\n'.join(ctg.split('\n')[1:])
			nscores = gt_frac_score(cscores)
			for key in nscores:
				tot_nscores[key] = (tot_nscores[key]) + nscores[key]
			#print(nscores)

	total_nscore_all = sum(tot_nscores.values())
	#print(total_nscore_all)
	for key in tot_nscores:
		tot_nscores[key] = round((tot_nscores[key])/total_nscore_all, 4)

	viral_gtness[pat] = tot_nscores

#print(viral_gtness)

with open('./output_table.tsv','w') as out:
	out.write("Patient\tHPV16\tHPV18\tHPV33\tHPV45\tHPV31\tHPV58\tHPV52\tHPV35\tHPV59\tHPV56\tHPV51\tHPV39\tHPV73\tHPV68\tHPV82\n")
	for pat in patients:
		vals = []
		for key in viral_gtness[pat]:
			vals.append(str(viral_gtness[pat][key]))

		line = str(pat+'\t'+'\t'.join(vals)+'\n')
		out.write(line)

with open('./pat_coinfect.tsv','w') as out:
	for pat in patients:
		coinfs = []
		for key in viral_gtness[pat]:
			if viral_gtness[pat][key] > threshold:
				coinfs.append(key)

		line = str(pat+'\t'+','.join(coinfs)+'\n')
		out.write(line)







	



