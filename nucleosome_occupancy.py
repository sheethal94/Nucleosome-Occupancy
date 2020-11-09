"""
name: Sheethal Umesh Nagalakshmi
lang: python3
file: nucleosome_occupancy.py
desc: Finds nculeosome occupancy

"""

import pandas as pd
import numpy as np

def occupancy(infile1,infile2,outfile):
	
	"""
		:param infile1: Input file containing nucleosome dyad positions 
		:param infile2: Input file containing transcription factor peaks
		:param outfile:	Output file containing nucleosome ocuupancy within CHIP-loci range
		:return: Array of nucleosome occupancy
		
	"""
	# Reading nucleosome positions file
	data = pd.read_csv(infile1, sep = '\t', names = ['chrom','chrStart','chrEnd'],
	skiprows=1, engine='python') 
	data[['chrStart','chrEnd']].astype(int64)
	# End position of the last nucleosome
	min_value = data['chrStart'].min()
	# Start position of the first nucleosome
	max_value = data['chrEnd'].max()
	
	# Reading TF CHIP-seq file containing chromosome number and TF binding sites
	data1 = pd.read_csv(infile2, sep = '\t', names = ['chrom','first_pos','last_pos'])
	
	# Storing the number of peaks
	no_rows = data1.shape[0]
	
	# Intitializing an array to store nucleosomes 
	nucleosome_number = np.zeros(max_value+1,dtype = np.int)
	for index,row in data.iterrows():
		start,end = row['chrStart'],row['chrEnd']
		# Storing the number of nucleosomes found in each position 
		nucleosome_number[start:end] += 1
		
	nucleosome_number.astype(int)
	result_array = ([])
	final_nucleosome_count = ([])

	for index,row in data1.iterrows():
		start,stop = row['first_pos'],row['last_pos']
		# Extracting nucleosomes within the CHIP-loci range
		count = nucleosome_number[start:stop] 
		result_array = np.append(result_array,count)
	nucleosome_count = np.zeros(no_rows*2000 - len(result_array))
	
	# Storing the aligned nucleosomes to the Transcriptional start sites 
	final_nucleosome_count = np.concatenate((result_array,nucleosome_count))
	nucleosome_count = np.reshape(final_nucleosome_count,(-1,2000))
	data2 = pd.DataFrame(nucleosome_count)
	data2.astype(int)
	data3 = data2.sum(axis=0)
	
	# Writing the number of nucleosomes in each position to output file
	pd.DataFrame(data3).to_csv(outfile, sep='\t', index = False)

def main():
	infiles1 = ['chemical_pos/1_pos.txt','chemical_pos/2_pos.txt']
	infiles2 = ['TF_chipseq/TF/chr1.txt','TF_chipseq/TF/chr2.txt']
	outfiles = ['TF_chipseq/TF/Array/chr1.txt','TF_chipseq/TF/Array/chr2.txt']	
	l = len(infiles1)
	for i in range(l):		
		infile1 = infiles1[i]
		infile2 = infiles2[i]
		outfile = outfiles[i]
		occupancy(infile1,infile2,outfile)
			
		
if __name__ == '__main__':
	main()


























	





