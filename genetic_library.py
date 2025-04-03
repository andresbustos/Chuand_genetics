1#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 16:54:31 2020

Library with functions to analyze genetic sequences. Mostly Fastq files.

In developement, now guaranteed to be 100% correct, sorry.

@author: Andres de Bustos Molina
"""

import pandas as pd
import re
import collections
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import Bio.Seq
import Bio.SeqUtils
from Bio import SeqIO
import sys
import os
from itertools import chain
import json
import glob
import operator
import random
from itertools import product, permutations

#names and configuration lists
nucleotides      = ['A','C','G','T']
beautiful_colors = ['MediumTurquoise','MediumSlateBlue', 'MidnightBlue','MediumVioletRed','Orange','Olive','DarkSlateGray','Fuchsia','DarkMagenta','LightSkyBlue','LimeGreen','MediumBlue','MediumSeaGreen','MediumOrchid','Salmon','Red','Yellow','indigo','blue','darkorange','darkolivegreen','saddlebrown','slategrey','cyan','darkcyan','chocolate','darkkhaki','darkslategrey','blueviolet','royalblue','deeppink','gold']


def read_parameters(parname):
	"""
	Reads the parameters.json file and writes on the screen
	"""
	with open(parname,'rt') as f:
		p=json.load(f)
	for folder in [p['outfolder']]:
		if not os.path.exists(folder):
			os.mkdir(folder)
	os.system('cp ' + parname + ' ' + p['outfolder']+'/')

	print('\n - Parameters - \n')
	for k in p.keys():
		print(k,'\t', p[k])

	#make a list of the fastq to analyze
	p.update({'fastq_names':[a for a in sorted(glob.glob(os.path.join(p['infolder'], '*.fastq')))]})
	other_fastqs = glob.glob(os.path.join(os.path.join(p['infolder'], '*/*.fastq'))) #find all fastqs in the subfolders
	for f in other_fastqs:
		p['fastq_names'].append(f)


	if "dict_conf" in p.keys():
		for pp in p['dict_conf'].keys():
			p['dict_conf'][pp].update({"files": [a for a in sorted(glob.glob(os.path.join(p['dict_conf'][pp]['dir'], '*.fastq')))]})
			#print('PMuS files:',pp, '-->',p['dict_conf'][pp],os.path.join(p['dict_conf'][pp]['dir'], '*.fastq'), [a for a in sorted(glob.glob(os.path.join(p['dict_conf'][pp]['dir'], '*.fastq')))])
			#also append all fastq files in all subfolders of p['dict_conf'][pp]['dir']
			other_fastqs = glob.glob(os.path.join(p['dict_conf'][pp]['dir'], '*/*.fastq'))
			for f in other_fastqs:
				p['dict_conf'][pp]['files'].append(f)
	return p

def fastq_to_df(filename):
	"""
	Reads a file with fastq format and returns a dataframe with the data
	"""
	with open(filename,'rt') as f:
		lines = f.readlines()
		if len(lines)%4 != 0:
			print('File not in fasq format\n')
	df              = pd.DataFrame(columns=['id','sequence','comments','quality'])
	df['id']        = lines[0::4]
	df['sequence']  = [Bio.Seq.Seq(a) for a in lines[1::4]]
	df['comments']  = lines[2::4]
	df['quality']   = lines[3::4]
	return df

def analyze_fastq(parameters, filename):
	"""
	Takes a FASTQ filename and peerforms statistics.
	All output files are albeled with the FAST filename
	"""

	print('\nAnalyzing file',filename )

	df             = pd.DataFrame()
	df['sequence'] = list(SeqIO.parse(filename, "fastq"))
	name           = os.path.basename(filename)[:-6]
	outfolder      = os.path.join(parameters['outfolder'], 'DESCRIPTION_'+name)
	if not os.path.exists(outfolder):
		os.mkdir(outfolder)

	#histogram on the sequences length
	lens   = [len(a.seq) for a in df['sequence']]
	av_len = np.mean(lens)
	plt.figure(dpi=500)
	plt.hist(lens, bins=25, alpha=0.75,histtype='bar',rwidth=0.9)
	plt.yscale('log')
	plt.title('Length of the Sequences')
	plt.xlabel('Lenght of the Sequence')
	plt.ylabel('Number of Sequences')
	plt.savefig(os.path.join(outfolder,name + '_histogram_length.png'))
	plt.clf()

	#histogram on the abundancy of the sequences
	cont = collections.Counter([str(a.seq) for a in df['sequence']]).most_common()
	values = [int(c[1]) for c in cont]
	plt.figure(dpi=500)
	plt.hist(values, bins=50, alpha=0.75,histtype='bar',rwidth=0.9)
	plt.yscale('log')
	plt.title('Historgam of Abundancy of the Sequences')
	plt.xlabel('Abundancy')
	plt.ylabel('Number of sequences')
	plt.savefig(os.path.join(outfolder, name + '_histogram_sequence_abundancy.png'))
	plt.clf()
	values = None

	#write the data to a csv file
	with open( os.path.join(outfolder, name + '_abundancies.csv'), 'wt') as f:
		f.write( "sequence,abundancy\n" )
		for k,v in  cont:
			f.write( "{},{}\n".format(k,v) )

	#histogram on the average quality of each sequence
	q = [ np.mean(a.letter_annotations['phred_quality']) for a in df['sequence'] ]
	plt.figure(dpi=500)
	plt.hist(q, bins=25, alpha=0.75,histtype='bar',rwidth=0.9)
	plt.yscale('log')
	plt.title('Average Quality of the Sequences')
	plt.xlabel('Phred quality')
	plt.ylabel('Number of sequences')
	plt.savefig(os.path.join(outfolder, name + '_histogram_sequence_quality.png'))
	plt.clf()
	q = None

	#histogram on the average quality of each nucleotide
	q = list( chain.from_iterable( [ a.letter_annotations['phred_quality'] for a in df['sequence']] ))
	num_nucl = len(q)
	av_q     = np.mean(q)
	plt.figure(dpi=500)
	plt.hist(q, bins=max(q)-min(q), alpha=0.75,histtype='bar',rwidth=0.9)
	plt.yscale('log')
	plt.title('Quality of the Nucleotides')
	plt.xlabel('Phred quality')
	plt.ylabel('Number of nucleotides')
	plt.savefig(os.path.join(outfolder,name + '_histogram_nucleotide_quality.png'))
	plt.clf()
	q = None

	#average nucleotide quality of each position on the FASTQ file.
	Maxlen = max([ len(a.seq) for a in df['sequence']])
	position = [a for a in range(1,Maxlen+1)]
	p_total  = np.zeros(Maxlen)
	p_total2 = np.zeros(Maxlen)
	count    = np.zeros(Maxlen)

	for s in df['sequence']:
		l = len(s.seq)
		count[0:l]    += 1
		p_total[0:l]  += s.letter_annotations['phred_quality']
		p_total2[0:l] += [a**2 for a in s.letter_annotations['phred_quality']]

	errors  = p_total2/count - [(p_total[i]/count[i])**2 for i in range(len(p_total))]
	errors  = errors**0.5

	fig, ax1 = plt.subplots(dpi=500, figsize=[9,4])
	plt.title('Av. Quality of the Nucleotides at each Position')
#	ax1.plot(position, p_total/count,  'bo-', alpha=0.4, label='Phred quality')
	ax1.errorbar(position, p_total/count, yerr=errors, fmt= 'bo-', alpha=0.4, label='Phred quality')

	ax1.set_xlabel('Position')
	ax1.set_ylabel('Prhed quality')
	ax2 = ax1.twinx()
	ax2.plot(position, count,  'r-', alpha=0.4, label='Num. sequences')
	ax2.set_ylabel('Num. sequences')
	fig.legend()
	plt.savefig(os.path.join(outfolder,name + '_nucleotide_quality_VS_position.png'))
	plt.clf()

	#average nucleotide quality of each position in ths POSITION INTERVAL on the FASTQ file.
	Maxlen   = parameters['pos_interval'][1] - parameters['pos_interval'][0]
	position = [str(a) for a in range(parameters['pos_interval'][0],parameters['pos_interval'][1]+1)]
	p_total  = np.zeros(Maxlen+1)
	p_total2 = np.zeros(Maxlen+1)
	count    = 0

	for s in df['sequence']:
		if  parameters['pos_interval'][1]+1 > len(s.seq):
			continue
		count += 1
		p = s.letter_annotations['phred_quality'][ parameters['pos_interval'][0]: parameters['pos_interval'][1]+1]
		p_total[:]  += p
		p_total2[:] += [a**2 for a in p]

	errors = p_total2/count - [(a/count)**2 for a in p_total]
	errors = errors**0.5

	fig, ax1 = plt.subplots(dpi=500, figsize=[6,4])
	plt.title('Av. Quality of the Nucleotides at the selected Positions')
	ax1.errorbar(position, p_total/count,yerr=errors, fmt='bo-', alpha=0.4, label='Phred quality')
	ax1.set_xticklabels([str(a) if int(a)%10==0 else '' for a in position ])
	ax1.set_xlabel('Position')
	ax1.set_ylabel('Prhed quality')
	plt.savefig(os.path.join(outfolder,name + '_nucleotide_quality_VS_sel_position.png'))
	plt.clf()

	#text file with some basic information
	f = open(os.path.join(outfolder,name + '_MISC.log'), 'wt')
	s = 'Basic information for the file '+ name +'.fastq\n\n'
	s += "Number of sequences                                 = " + str(len(df))+'\n'
	s += "Number of nucleotides                               = " + str(num_nucl) +'\n'
	s += "Average sequence length                             = " + str(av_len)+ '\n'
	s += "Average nucleotide quality                          = " + str(av_q)+ '\n'
	s += "Sequences found with the selected position interval = " + str(count)+ '\n'


	c = collections.Counter( chain.from_iterable( [ a for a in df['sequence']] ) )
	c = collections.OrderedDict(sorted(c.items()))
	num_nucl =  sum(c.values())
	s += "\n Composition of the " + str(num_nucl) + " nucleotides:\n"
	for k in c.keys():
		s += '\t' + k + ' --> ' + str(c[k]) + '  (' + str(round(100.*c[k]/num_nucl,2)) + '%) \n'
	f.write(s)
	f.close()

def merge_fastq(fastq_list, new_name ='newfile.fastq'):
	"""
	Combines a list of fastq files into a single one with the new name profived.

	Parameters
	----------
	fastq_list : TYPE
		DESCRIPTION.
	new_name : TYPE, optional
		DESCRIPTION. The default is 'newfile.fastq'.

	Returns
	-------
	None.

	"""
	print('\nMerging FASTQ files ...')
	try:
		os.ystem('rm '+ new_name)
	except:
		pass

	lines = []
	for file in  fastq_list:
	    f = open(file, 'rt')
	    a = f.readlines()
	    lines.extend(a[1::2])
	    f.close()

	string   = new_name[:-5]
	f        = open(new_name, 'wt')
	for i in range(int(len(lines)/2)):
	    f.write('@' + string + str(i+1) + ' ' + str(i+1) +'/1\n')
	    f.write(lines[2*i])
	    f.write('+\n')
	    f.write(lines[2*i+1])
	f.close()
	print('Done\n')

def BCS(par):
	"""
	Searches the barcodes defined in the library string in the fastq files
	"""

	lib = par['lib_string']
	lib = lib.replace('N','[A,C,G,T]')
	lib = lib.replace('S','[C,G]')

	fff = open(os.path.join(par['outfolder'], 'ALL_SEARCHES.txt'), 'wt')

	for f in par['fastq_names']:
		name           = os.path.basename(f[:-6])
		outfolder      = os.path.join(par['outfolder'], 'SEARCHES_'+name)
		if not os.path.exists(outfolder):
			os.mkdir(outfolder)

		df             = pd.DataFrame({'sequence': list(SeqIO.parse(f, "fastq"))})
		df['in_lib']   = [ re.findall(lib, str(s.seq)) for s in df['sequence'] ]

		ff =  open(os.path.join(outfolder, 'SEARCH_DATA.txt'),'wt')
		s  = '** Some data about the search of ' + par['lib_string'] + ' in file ' + f + ' ** \n'
		s += '\nFound ' + str( len( [a[0] for a in df['in_lib'] if a!=[]]) ) + ' sequences in a FASTQ with '  + str(len(df)) + ' entries\n'
		ff.write(s)
		ff.close()

		fff.write(s)

		quality_plot(df, par, outfolder)
		s = repetitions_stuff(df['in_lib'], outfolder, par)
		fff.write(s+'\n\n')
		composition_stuff(df['in_lib'], outfolder)

	fff.close()

def quality_plot(df, par, outfolder):
	"""
	Plots the average quality of the lib sequences found in the FASTQ file,
	as a function of the nucleotide position in the lib_string.
	"""
#	from re import finditer

	from matplotlib.patches import Rectangle
	N      = len(par['lib_string'])
	q_sum  = np.zeros(N)
	q_sum2 = np.zeros(N)
	cont   = 0

	for f_seq, lib_seq in zip(df['sequence'],df['in_lib']):
		if len(lib_seq)==0:
			continue

		match = [a for a in re.finditer(lib_seq[0],str(f_seq.seq))]
		if (len(match)>1):
			print('More than one lib string found in an entry of the fastq. This should not happen. Exiting.')
			sys.exit()

		q      = f_seq.letter_annotations['phred_quality'][match[0].span()[0]: match[0].span()[1]]
		q_sum  += q
		q_sum2 += [a**2 for a in q]
		cont   += 1

		#figure

	N_pos = [a for a in re.finditer('N', par['lib_string'])]
	S_pos = [a for a in re.finditer('S', par['lib_string'])]
	errors = q_sum2/(cont) - [(a/(cont))**2 for a in q_sum]
	errors = errors**0.5
	Merror = max(max(errors), 0.5)
	ymin = min(q_sum/cont)-1-Merror  #for the rectangles
	R_h  = max(q_sum/cont) - ymin+1+Merror #for the rectangles

	fig = plt.figure(dpi=500, figsize=(9,4))
	ax  = fig.add_subplot(111)
	plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
	ax.errorbar([a for a in range(1, N+1)], q_sum/cont, yerr=[a for a in errors], fmt='bo', alpha=0.4 )
	plt.title('Quality of the sequences found in the FASTQ')
	ax.set_xlabel('Position and nucleotides (lib_string)')
	ax.set_ylabel('Phred quality')
	plt.xticks(ticks=[a for a in range(N)], labels=[str(i+1)+'\n'+a  if (i+1)%5==0 else '\n'+a for i,a in enumerate(par['lib_string'])])
	for p in N_pos:
		ax.add_patch(Rectangle( xy=(p.span()[0]-0.5, ymin ), width=1, height=R_h, linewidth=0, color='red', alpha=0.25, fill=True)
)
	for p in S_pos:
		ax.add_patch(Rectangle( xy=(p.span()[0]-0.5, ymin ), width=1, height=R_h, linewidth=0,color='blue', alpha=0.25, fill=True)
)
	plt.savefig(os.path.join(outfolder,'av_quality_found_sequences.png'))
	plt.clf()

def repetitions_stuff(sequences, outfolder, par):
	"""
	To detect repetitions of the sequences of the library in the FAST file
	"""
	sequences   = [a[0] for a in sequences if a!=[]]
	count       = collections.Counter( sequences )

	#divide the sequences in thow sets, maj and min
	df0   =  pd.DataFrame(data={'seq':[a for a in count.keys()], 'count':[a for a in count.values()]} )
	df0.sort_values(by=['count'], ascending=False, inplace=True)

	limit     = par['BCS_threshold']*sum(df0['count'])
	acc_total = 0
	sets      = ['min' for _ in range(len(df0))]

	for i,c in enumerate(df0['count']):
		if acc_total > limit:
			break
		acc_total += c
		sets[i]='maj'
	df0['set']=sets
	df0.to_csv(os.path.join(outfolder,'sequence_counts.csv'),	index=False)
	sets_counter = collections.Counter(df0['set'])

   #some data to a file
	f =  open(os.path.join(outfolder, 'SEARCH_DATA.txt'),'at')
	s =  '\n#Most common sequence ocurrences in the fastq file'
	s += '\n#sequence_in_the library, num_ocurrences\n'
	for c in count.most_common(10):
		s += c[0] + ','+ str(c[1])+'\n'
	s += '\nNumber of different sequences in the maj group = ' + str(sets_counter['maj']) +'\n'
	s += '\nNumber of different sequences in the min group = ' + str(sets_counter['min']) +'\n'

	f.write(s)
	f.close()

	#figure
	plt.figure(dpi=500)

	if len(count)==1:
		N = list(count.values())[0]
	elif len(count)==0:
		N = 1
		print('Not enough data to do repetitiions stuff.')
		return 0
	else:
		N = max(count.values())

	h = np.histogram([a for a in count.values()], bins=N)
	plt.plot([a for a in range(1,N+1)],h[0], 'bo', alpha=0.4 )
	plt.yscale('log')
	plt.title('Total:' + str(sum(count.values())) + ' sequences')
	plt.xlabel('Abundance')
	plt.ylabel('Number of unique bars')
	plt.savefig(os.path.join(outfolder,'histogram_repetition.png'))
	plt.clf()

	return s

def composition_stuff(sequences, outfolder):
	"""
	Composition [A,C,G or T] for every position in the sequence
	"""

	#analyze the data anf count
	sequences = [a[0] for a in sequences if a!=[]]
	a         = np.array( [list(s) for s in sequences])
	composition = []

	if len(a) == 0:
		print('No data available to do composition_stuff\n')
		return 0

	for loc in range(a.shape[1]):
		composition.append(collections.Counter(a[:,loc]))

	if len(composition) == 0:
		print('Not enough data to calculate the composition.')
		return 0

	histos = {}
	for nucl in nucleotides:
		d = {nucl : [ b[nucl] if nucl in b.keys() else 0 for b in composition ] }
		histos.update(d)

	#figure
	fig, ax = plt.subplots(dpi=500)
	width   = 0.65
	xvalues = [b+1 for b in range(a.shape[1])]
	ax.bar(xvalues, histos['A'], width, label='A', color='g', alpha=0.4)
	bot = [ histos['A'][i] for i in range(a.shape[1])]
	ax.bar(xvalues, histos['C'], width, bottom=bot,label='C', color='r', alpha=0.4)
	bot = [ bot[i] + histos['C'][i] for i in range(a.shape[1])]
	ax.bar(xvalues, histos['G'], width, bottom=bot,label='G', color='b', alpha=0.4)
	bot = [ bot[i] + histos['G'][i] for i in range(a.shape[1])]
	ax.bar(xvalues, histos['T'], width, bottom=bot,label='T', color='orange', alpha=0.4)
	plt.xlabel('Nucleotide position in the sequence')
	plt.ylabel('Appearances')
	ax.legend()
	plt.savefig(os.path.join(outfolder,'histogram_composition.png'))

def abundancy_filter(df0, col, threshold):
	"""
	* filter df0 according to column col. The criteria is the abundancy of col.
	* theshold is a number in (0,1)
	* we keep the entries that represent  a fraction 'threshold' of the total number of entries
	* All ties are kept.

	"""

	cont = collections.Counter(df0[col]).most_common()
	print('\tFound ' + str(len(cont)) + ' diferent sequences')

	limit = int(len(df0))*threshold
	isel  = 0
	acc_sum = 0


	for ic in range(len(cont)):
		acc_sum += cont[ic][1]
		if acc_sum >= limit:
			isel = ic+1
	#	print('\tLimit reached, cutting the number of sequences to ',isel)
			break

	for isel0 in range(isel,len(cont)):
		if cont[isel0][1]<cont[isel][1]:
	#				print('\tFinal number of sequences:', isel0)
			isel = isel0
			break

	sel_data = [c[0] for c in cont[:isel]]
# 	s += '\tUnique sequences to keep ' +str(len(sel_sequences)) + '\n'
	df0 = df0[df0[col].isin(sel_data)]
# 	s += '\tAfter filtering there are ' + str(len(df)) + ' sequences'  + '\n'

	return df0

def representativity_filter(df0, col, threshold):
	"""
	* filter df0 according to column col. The criteria is the abundancy of col.
	* theshold is a number in (0,1)
	* we keep the entries that represent  a fraction 'threshold' of the total number of entries
	* All ties are kept.

	"""

	cont = collections.Counter(df0[col])
	print('\tFound ' + str(len(cont)) + ' diferent sequences')
	total = len(df0)
	df0['representativity'] = [float(cont[c])/total*100 for c in df0[col]]

	return df0[df0['representativity']>=threshold]

def PMuS(par):
	"""
	Performs Point Mutation Search according to the input parameters
	"""
#	par = read_parameters('parameters_PMuS.json')

	for mut in par['dict_conf'].keys():

		for f in par['dict_conf'][mut]['files']:

			s = '\nPMuS on file ' + f  + '\n'
			name           = os.path.basename(f[:-6])
			outfolder      = os.path.join(par['outfolder'], 'PMuS_'+name)
			if not os.path.exists(outfolder):
				os.mkdir(outfolder)

			df             = pd.DataFrame({'sequence': list(SeqIO.parse(f, "fastq"))})
			df['sequence'] = [str(a.seq) for a in df['sequence']]
			#limit to the number of sequences that represent > 99% of the total sequences in the file
			s +=  '\tBefore filtering there are ' + str(len(df)) + ' sequences'  + '\n'
			df = representativity_filter(df, 'sequence', par['rep_filter_thr'])
			s += '\tUnique sequences to keep ' +str(len([a for a in set(df['sequence'])])) + '\n'

			s += '\tAfter filtering there are ' + str(len(df)) + ' sequences'  + '\n'


			expr = par['dict_conf'][mut]['bc'].upper()
			aux  =  [a for a in par['dict_conf'][mut]['mutation'].keys() ]
			aux = '('+'|'.join(aux)+')'

			expr = expr.replace('X',aux)

			#find sequencies and save the results in the string
			df['match']   = [ re.findall(expr, s) for s in df['sequence'] ]
			sequences     = [a[0] for a in df['match']  if a!=[]]
			s  +=  '\tFound a total of ' + str(len(sequences)) + ' sequences (with possible repetitions)\n'
			count         = collections.Counter( sequences )

			suma = 0
			plot_dict = {} #for plots
			s += '\tRegex = ' + expr +  '  \n'
			for k in  par['dict_conf'][mut]['mutation'].keys():
				if k in count.keys():
					s+= '\t    ' +str(k) + '    ' + str(count[k]) + ' <--> ' +str( count[k]/len(df)*100.)+  ' %' + '\n'
					suma += count[k]
					plot_dict.update({par['dict_conf'][mut]['mutation'][k]:count[k]})
				else:
					s += '\t    ' +str(k) + '    0 <--> ' +  '0.00 %'  + '\n'
					plot_dict.update({par['dict_conf'][mut]['mutation'][k]:0})

			s += '\t    OTHER  ' + str(len(df) - suma) +  '    <--> ' + str( (len(df) - suma)/len(df)*100. ) + ' %'  + '\n'
			plot_dict.update({'OTHER':len(df)-suma})

			#bar plot with the previous results
			fig, ax = plt.subplots(dpi=500, figsize=[2,4])
			plt.subplots_adjust(left=0.5)
			plt.title('PMuS')
			ax.set_xticks([])
			y_offset = 0
			colors = ['blue','green', 'orange','red' ] + beautiful_colors
			for i,k in enumerate(plot_dict.keys()):
				if  plot_dict[k] > 0:
					plt.bar(0, plot_dict[	k], 0.5, bottom=y_offset, color=colors[i], alpha=0.5)
					plt.text(-0.13, y_offset+plot_dict[k]*0.4, k)
					y_offset += plot_dict[k]
			fig.savefig(os.path.join(outfolder, os.path.basename(f[:-6])+'_barplot_PMuS.png'))

			#file with the PMuS results, for further plot
			file = open(os.path.join(outfolder, 'PMuS.csv'), 'wt')
			file.write('mutation,count,percentage\n')
			total = sum(plot_dict.values())
			for k in plot_dict.keys():
				file.write(k+','+ str(plot_dict[k])+',' + str(plot_dict[k]*100./total) + '\n' )
			file.close()


			#filtering
			mask = [len(a)==0 for a in df['match']]

			#now save in a file the most common found sequences
			count = collections.Counter(df[[not a for a in mask]]['sequence'])
			acc_sum = 0
			with open( os.path.join(outfolder, os.path.basename(f[:-6])+'_most_common_sequences.csv'), 'wt') as file:
				file.write( "sequence,abundancy\n" )
				for k,v in  count.most_common(6):
					file.write( "{},{}\n".format(k,v) )
					acc_sum += v
				file.write( "{},{}\n".format('Rest of detected sequences',len(df[[not a for a in mask]])-acc_sum) )

			#write a file with information with the seauences that appear in 'OTHERS'
			count_others = collections.Counter(df[mask]['sequence']).most_common()
			#write the data to a csv file
			with open( os.path.join(outfolder, os.path.basename(f[:-6])+'_OTHER_abundancies.csv'), 'wt') as file:
				file.write( "sequence,abundancy\n" )
				for k,v in  count_others:
					file.write( "{},{}\n".format(k,v) )

			#print the onfo on the screen and save it in a log file
			print(s)
			file = open(os.path.join(outfolder, os.path.basename(f[:-6])+'_PMuS.log'),'wt')
			file.write(s)
			file.close()

##############################################################################
def transversal_BCS_analysis(par, N=6):
	"""
	Creates plots with the data of all the fastqs
	"""
	#plot with the N  most common occurrences

 	# par = parameters

	colorlist = beautiful_colors + [c for c in mcolors.CSS4_COLORS.keys() if c not in beautiful_colors]

	#label for the output data (png and dat file)
	assert(N>0)
	label = str(N).zfill(4)

	icolor     = 0
	color_dict = {}
	df_values =pd.DataFrame(columns = ['fastq_file', 'sequence','percentage'])

	fig, ax = plt.subplots(dpi=600) #, figsize=[2,4])

	max_len = max([len(os.path.basename(a)[0:-6]) for a in par['fastq_names']])
# 	plt.subplots_adjust(left=0.1, bottom=0.23)
	plt.subplots_adjust(left=0.1, bottom=0.23*max_len/18.)

# 	plt.subplots_adjust(left=0.1, bottom=0.23)
	plt.title('The '+str(N)+' most common sequences in each file')
	ax.set_xticks([])
	ax.set_ylim(0,108)
	ax.set_ylabel('Representativity in %')
	ax = plt.gca()
	if len(par['fastq_names']) <=6:
		ax.set_xlim(-0.75,5.75)
		fsize=7
	else:
		fsize=5

	rep = [] #list of numbers of the representativity of the selected sequences in the file [%].

	for ifile, f in enumerate(par['fastq_names']):
		name           = os.path.basename(f[:-6])
		df             = pd.read_csv(os.path.join(par['outfolder'], 'SEARCHES_'+name,'sequence_counts.csv')	)
		norm           = df['count'].sum()
		y_offset       = 0
		rep.append(100.*df['count'][:N].sum() / df['count'].sum())

		df = df.head(N)

		for seq, val in zip(df['seq'], df['count']):
			if seq not in color_dict.keys():
				color_dict.update({seq:{'color':colorlist[icolor], 'label':'s-'+str(icolor)} })
				icolor += 1
			value     = val*100./norm
			plt.bar(ifile,value, 0.8, bottom=y_offset, color=color_dict[seq]['color'], alpha=0.5, linewidth=1, edgecolor='k')
			plt.text(ifile-0.1, y_offset+value*0.4, color_dict[seq]['label'], fontsize=6)
			y_offset += value

		#grey block with the rest of sequences
		plt.bar(ifile,100.-y_offset, 0.8, bottom=y_offset, color='#E0E0E0', alpha=0.5, linewidth=1, edgecolor='k') #,hatch='x')
		plt.text(ifile-0.3, 103, r'$R_T=$'+str(rep[ifile].round(1)) + '%', fontsize=fsize)

		#store all data to use when writing data file
		df             = pd.read_csv(os.path.join(par['outfolder'], 'SEARCHES_'+name,'sequence_counts.csv')	)
		df_temp = pd.DataFrame(columns = ['fastq_file', 'sequence','percentage'])
		df_temp['sequence']            = df['seq']
		df_temp['fastq_file']          = [f for _ in range(len(df_temp['sequence']))]
		df_values = pd.concat([df_values, df_temp])

	plt.xticks([a for a in range(len(par['fastq_names']))], [os.path.basename(a)[0:-6] for a in par['fastq_names']], rotation=-45	, fontsize=7)
	fig.savefig(os.path.join(par['outfolder'], 'BCS_barplots_'+label+'.png'))
	plt.clf()

	#make again the previous plot, but NAKED
	fig, ax = plt.subplots(dpi=600)
	plt.subplots_adjust(left=0.1, bottom=0.23*max_len/18.)
	ax.set_xticks([])
	ax.set_ylim(0,108)
	ax.set_ylabel('Representativity in %')
	if len(par['fastq_names']) <=6:
		ax.set_xlim(-0.75,5.75)
	plt.xticks([])

	rep = [] #list of numbers of the representativity of the selected sequences in the file [%].

	for ifile, f in enumerate(par['fastq_names']):
		print(ifile)
		name           = os.path.basename(f[:-6])
		df             = pd.read_csv(os.path.join(par['outfolder'], 'SEARCHES_'+name,'sequence_counts.csv')	)
		norm           = df['count'].sum()
		y_offset       = 0
		rep.append(100.*df['count'][:N].sum() / df['count'].sum())

		df = df.head(N)
		for seq, val in zip(df['seq'], df['count']):
			if seq not in color_dict.keys():
				color_dict.update({seq:{'color':colorlist[icolor], 'label':'s-'+str(icolor)} })
				icolor += 1
			value     = val*100./norm
			plt.bar(ifile, value, 0.8, bottom=y_offset, color=color_dict[seq]['color'], alpha=0.5, linewidth=1, edgecolor='k')
			y_offset += value

		#grey block with the rest of sequences
		plt.bar(ifile,100.-y_offset, 0.8, bottom=y_offset, color='#E0E0E0', alpha=0.5, linewidth=1, edgecolor='k') #,hatch='x')
	ax = plt.gca()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(True)
	plt.yticks([0,20,40,60,80,100])
	fig.savefig(os.path.join(par['outfolder'], 'BCS_barplots_'+label+'_naked.png'))
	plt.clf()



	#make again the previous plot, but CIRCULAR
	fig = plt.figure(dpi=600)
	ax  = plt.subplot(111,polar=True)
	ax.set_xticks([])

	angle_cont = 0.5
	d_angle    = 2.* np.pi / ( len(par['fastq_names']) + 1.)
	d_angle_sep = d_angle / len(par['fastq_names'])
	min_r      = 50.
	max_r      = 100.
	for ifile, f in enumerate(par['fastq_names']):
		name           = os.path.basename(f[:-6])
		df             = pd.read_csv(os.path.join(par['outfolder'], 'SEARCHES_'+name,'sequence_counts.csv')	)
		norm           = df['count'].sum()
		y_offset       = 0

		df = df.head(N)

		for seq, val in zip(df['seq'], df['count']):
			if seq not in color_dict.keys():
				color_dict.update({seq:{'color':colorlist[icolor], 'label':'s-'+str(icolor)} })
				icolor += 1
			value     = (val*100./norm)*(max_r-min_r)/max_r
			ax.bar(angle_cont*d_angle + ifile*d_angle_sep, value, d_angle, bottom=min_r+y_offset, color=color_dict[seq]['color'], alpha=0.5, linewidth=0.4, edgecolor='k')
			y_offset += value
		#grey block with the rest of sequences
		ax.bar(angle_cont*d_angle+ ifile*d_angle_sep, max_r-min_r-y_offset, d_angle, bottom=min_r+y_offset, color='#E0E0E0', alpha=0.5, linewidth=0.4, edgecolor='k') #,hatch='x')
		angle_cont += 1

	plt.axis('off')
	fig.savefig(os.path.join(par['outfolder'], 'BCS_barplots_'+label+'_naked_circular.png'))
	plt.clf()


	df_values= df_values.set_index(['fastq_file', 'sequence'])
	#write data into a file
	ff  = open( os.path.join(par['outfolder'], 'BCS_barplots_data_'+label+'.dat'),'wt' )
	ff.write('sequence,label,color,')
	for fastqname in par['fastq_names']:
		ff.write('Percentage_'+fastqname+',')
	ff.write('\n')
	for seq in color_dict.keys():
		ff.write(seq+','+color_dict[seq]['label']+','+ color_dict[seq]['color'])
		for fastq_file in par['fastq_names']:
			if (fastq_file,seq) in df_values.index:
				ff.write(',' + str( df_values.loc[fastq_file,seq]['percentage']))
			else:
				ff.write(',' + '0.0')
		ff.write('\n')
	ff.close()


	#now we track the sequences in par['Ref_fastq'] in the files par["List_fastq"], to produce 'Sequence_tracking' results
	results_dict = {}

	df_ref       = pd.read_csv(os.path.join(par['outfolder'], 'SEARCHES_'+par['Ref_fastq'][:-6],'sequence_counts.csv'), index_col='seq'	)
	norm         = df_ref['count'].sum()
	df_ref       = df_ref[:N]
	if norm == 0:
		norm = 1.

	results_dict[par['Ref_fastq']]  = [a / norm *100. for a in df_ref['count'].values]

	for f in par['List_fastqs']:
		results_dict[f] = []
		df0  = pd.read_csv(os.path.join(par['outfolder'], 'SEARCHES_'+ f[:-6],'sequence_counts.csv'), index_col='seq'	)
		norm = df0['count'].sum()
		for s in df_ref.index:
			if s in df0.index:
				a = df0.loc[s]['count'] /norm *100.
			else:
				a=0.0
			results_dict[f].append(a)

	#now we got the data, we plot it

	fig,ax = plt.subplots(dpi=600)
	max_len = max( [len(a) for a in results_dict.keys()]) - 6 #The 6 is for the .fastq extension that we remove
	plt.subplots_adjust(left=0.12, bottom=0.23*max_len/17.)
	ax.set_ylabel('Representativity in %')
	for i in range(N):
		y  = [results_dict[k][i] for k in results_dict.keys()]
		ax.plot(y, '-o', label='s'+str(i), alpha=0.6)
	plt.legend()
	plt.xticks([a for a in range(len(results_dict.keys()))], [a[:-6] for a in results_dict.keys()], rotation=-45, fontsize=7)

	plt.savefig(os.path.join(par['outfolder'],'Sequence_tracking.png'))

	df_out = pd.DataFrame(results_dict)
	df_out.index = df_ref.index
	df_out.to_csv(os.path.join(par['outfolder'], 'Sequence_tracking.csv'))

##############################################################################
def transversal_PMuS_analysis(par):
	"""
	Creates plots with the data of all the PMuS analysis previously done
	"""

	#plot with the N  most common occurrences
#	par  = gl.read_parameters('parameters_PMuS.json')

	colorlist = beautiful_colors + [c for c in mcolors.CSS4_COLORS.keys() if c not in beautiful_colors]


	for mut_folder in par['dict_conf'].keys():

		df_all = pd.DataFrame(columns=['file','mutation','count','percentage']) #df to save all PMus.csv data into a single file

		icolor     = 0
		color_dict = {}

		fig, ax = plt.subplots(dpi=500) #, figsize=[2,4])
		plt.subplots_adjust(left=0.1, bottom=0.28)
		plt.title('PMuS composition')
		ax.set_xticks([])
		ax.set_ylabel('Representativity in %')

		for ifile,f in enumerate(par['dict_conf'][mut_folder]['files']):
			name           = os.path.basename(f[:-6])
			datafolder     = os.path.join(par['outfolder'], 'PMuS_'+name)
			df             = pd.read_csv(os.path.join(datafolder, 'PMuS.csv'))
			df['file']     = [f for _ in range(len(df))]
			df_all         = pd.concat([df_all, df])
			y_offset       = 0

			for mut, val in zip(df['mutation'], df['percentage']):
				if mut not in color_dict.keys():
					color_dict.update({mut:{'color':colorlist[icolor]} })
					icolor += 1
				if val > 0:
					plt.bar(ifile,val, 0.45, bottom=y_offset, color=color_dict[mut]['color'], alpha=0.5, linewidth=1, edgecolor='k')
					plt.text(ifile-0.1, y_offset+val*0.4, mut, fontsize=6, rotation=90)
					y_offset += val

		plt.xticks([a for a in range(len(par['dict_conf'][mut_folder]['files']))], [os.path.basename(a)[0:-6] for a in par['dict_conf'][mut_folder]['files']], rotation=90	, fontsize=6)
		fig.savefig(os.path.join(par['outfolder'], str(mut_folder)+ '_PMuS_barplots.png'))

		df_all.to_csv(os.path.join(par['outfolder'], str(mut_folder)+'_all_PMuS.csv' ), index=False)

def find_coincidences(par):
	"""
	Creates an array with % coincidences: outputfolder/percentage_tracking.csv

	There the element [i,j] represents the fracion of seuqences in file i that appears in file j (in %)
	"""

	fastq_list = par['List_fastqs']
	if par['Ref_fastq'] not in par['List_fastqs']:
		fastq_list = [par['Ref_fastq']] + par['List_fastqs']
	fastq_list = sorted(fastq_list)

	Nfiles  = len(fastq_list)
	results           = np.zeros([Nfiles, Nfiles]) #for percetage_tracking: Sequence_tracking.csv
	results_bis       = np.zeros([Nfiles, Nfiles]) #for percetage_tracking_absolute: percentage_tracking.csv
	results_tris      = np.zeros([Nfiles, Nfiles]) #for tracking_absolute: unique_sequencies.csv
	unique_sequencies = np.zeros(Nfiles)
	df_list = []

	for f in fastq_list:
		df_list.append( pd.read_csv(os.path.join(par['outfolder'], 'SEARCHES_'+f[:-6],'sequence_counts.csv'), index_col='seq'))

	for i in range(Nfiles):
		unique_sequencies[i] = len(df_list[i])
		for j in range(i+1):
			df = df_list[i].join(df_list[j], how='inner',lsuffix='_i',rsuffix='_j')
			if i==j:
				results[i,j] = 100.
				results_bis[i,j] = 100.
			else:
				results[i,j]     = sum([ min(a,b) for a,b in zip(df['count_i'], df['count_j'])]) / df_list[i]['count'].sum()*100.
				results_bis[i,j] = len(df) / len(df_list[j])*100.
				results[j,i]     = sum([ min(a,b) for a,b in zip(df['count_i'], df['count_j'])]) / df_list[j]['count'].sum()*100.
				results_bis[j,i] = len(df) / len(df_list[i])*100.
			results_tris[i,j] = len(df)
			results_tris[j,i] = len(df)

	#save results
	df = pd.DataFrame(index=fastq_list, columns=fastq_list, data=results).rename_axis('Ref_fastq')
	df.to_csv(os.path.join(par['outfolder'],'percentage_tracking.csv'),sep=',')

	#plot the data
	fig, ax = plt.subplots(dpi=800)
	max_len = max( [len(a) for a in fastq_list]) - 6 #The 6 is for the .fastq extension that we remove
	plt.subplots_adjust(left=0.015*max_len + 0.06, bottom=0.015*max_len+0.01)
	cmap ='rainbow'
	results  = results[::-1,:]
	im = ax.pcolormesh(results, cmap=cmap, alpha=0.7)
	ax.set_xticklabels([a[:-6] for a in fastq_list], fontsize=7, rotation=90)
	ax.set_yticklabels([a[:-6] for a in fastq_list[::-1]], fontsize=7, rotation=0)
	ax.set_xticks([a + 0.5 for a in range(Nfiles)])
	ax.set_yticks([a + 0.5 for a in range(Nfiles)])
	ax.set_ylabel('Reference')
	fig.colorbar(im, ax=ax)
	plt.savefig(os.path.join(par['outfolder'],'percentage_tracking.png'))

	#save results bis
	df = pd.DataFrame(index=fastq_list, columns=fastq_list, data=results_bis).rename_axis('Ref_fastq')
	df.to_csv(os.path.join(par['outfolder'],'percentage_tracking_absolute.csv'),sep=',')
	#plot the data
	fig, ax = plt.subplots(dpi=800)
	max_len = max( [len(a) for a in fastq_list]) - 6 #The 6 is for the .fastq extension that we remove
	plt.subplots_adjust(left=0.015*max_len + 0.06, bottom=0.015*max_len+0.01)
	cmap ='rainbow'
	results_bis  = results_bis[::-1,:]
	im = ax.pcolormesh(results_bis, cmap=cmap, alpha=0.7)
	ax.set_xticklabels([a[:-6] for a in fastq_list], fontsize=7, rotation=90)
	ax.set_yticklabels([a[:-6] for a in fastq_list[::-1]], fontsize=7, rotation=0)
	ax.set_xticks([a + 0.5 for a in range(Nfiles)])
	ax.set_yticks([a + 0.5 for a in range(Nfiles)])
	ax.set_ylabel('Reference')
	fig.colorbar(im, ax=ax)
	plt.savefig(os.path.join(par['outfolder'],'percentage_tracking_absolute.png'))


	df = pd.DataFrame(index=fastq_list, columns=fastq_list, data=results_tris, dtype=int).rename_axis('Ref_fastq')
	df.to_csv(os.path.join(par['outfolder'],'common_sequencies.csv'),sep=',')

	df = pd.DataFrame(index=fastq_list, columns=['unique_sequencies'], data=unique_sequencies, dtype=int).rename_axis('fastq')
	df.to_csv(os.path.join(par['outfolder'],'unique_sequencies.csv'),sep=',')