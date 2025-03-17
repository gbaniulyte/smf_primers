from cmath import nan
import os
import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
from datetime import datetime
print('Start:', datetime.now().strftime("%d-%m-%Y %Hh:%Mm"))

#primer parameters
ampl_size, tm_range, primer_size, no_gcg_min_range = [300, 500], [62, 66], [18, 36], 12
min_nt_from_motif, tm_difference, y_count, y_3end = 20, 3, 3, 10
def rev_compl(sequence):
	sequence = sequence.upper()
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'Y':'R'}
	reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))
	return reverse_complement
def bisulf_convert_Y(seq):
	bases = str(seq).upper()
	bases = bases.replace("CG", "YG").replace("GC", "GY").replace("C", "T")
	return bases
def gcg_positions(seq):
	gcg_list = []
	for pos in range(0, len(seq)+1):
		di_nt = seq[pos:pos+2]
		if di_nt == 'GC':
			gcg_list.append(pos+1)
		elif di_nt == 'CG':
			gcg_list.append(pos)
	return list(set(gcg_list)) #remove duplicates that arise due to GCG motifs
#try finding a region of x nt with no GC/CG; add rev primer ranges, diff gcg list
def gcg_free_regions(gcg_pos_list, primer_region_list):
	ranges = []
	for i, nt in enumerate(gcg_pos_list[:-1]):
		if gcg_pos_list[i+1] - nt >= no_gcg_min_range:
			if nt and gcg_pos_list[i+1] in primer_region_list:
				ranges.append([nt + 1, gcg_pos_list[i+1] -1])
	return ranges
def all_possible_primers(gcg_free_ranges, seq, motif_start_end):
	ups_primers = {}
	ds_primers = {} #plus orientation first
	if gcg_free_ranges != []: #move out of function
		for r in gcg_free_ranges: # r = [31, 57]
			for end in range(r[1], r[0] + 1, -1): #nt = 31, 32 etc
				for lenght in range(primer_size[0], primer_size[1] + 1): #lenght = 18, 19 etc..
					start = end - lenght
					primer = seq[start : end]
					try:
						tm = mt.Tm_NN(primer, nn_table = mt.DNA_NN4, Na = 0, K = 50, Tris = 25, Mg = 2, dNTPs = 0.2, saltcorr = 6)
					except IndexError:
						break
					if tm_range[0] < tm <= tm_range[1]:
						if start < motif_start_end[0] - primer_size[1] and primer.count('Y') <= 4: #add selection against Y in the 3' end???
							ups_primers.update({primer:[tm, start]})
						elif start > motif_start_end[1] + primer_size[1] and primer.count('Y') <= 4:
							ds_primers.update({primer:[tm, start]})
					elif tm_range[1] < tm:
						break
	return ups_primers, ds_primers
#filter out primers that have Y within x nt from 3' end
def filter_y_3end(primer_list, nt_from_3end):
	reduced_list = {}
	for primer_pair, tm_size in primer_list.items():
		if 'Y' not in primer_pair[0][-nt_from_3end:] and 'Y' not in primer_pair[1][-nt_from_3end:]:
			reduced_list.update({primer_pair:tm_size})
	return reduced_list
#filter primers that have at least 1 G at the 3end
def gc_3end(primer_list):
	reduced_list = {}
	for primer_pair, tm_size in primer_list.items():
		if 'G' in primer_pair[0][-1] and 'G' in primer_pair[1][-1]:
			reduced_list.update({primer_pair:tm_size})
		if reduced_list == {}:
			if 'G' in primer_pair[0][-1] or 'G' in primer_pair[1][-1]:
				reduced_list.update({primer_pair:tm_size})
	return reduced_list
#then try and pair them for the closest Tm, then by the best distance
def pair_primers(us_primer_dict, ds_primer_dict):
	primer_pairs = {}
	for us_primer, us_tm_start in us_primer_dict.items():
		for ds_primer, ds_tm_start in ds_primer_dict.items():
			tm_diff = abs(us_tm_start[0] - ds_tm_start[0])
			if tm_diff <= tm_difference:
				amplicon_size = ds_tm_start[1] + len(ds_primer) - us_tm_start[1]
				if ampl_size[0] <= amplicon_size <= ampl_size[1]:
					ups_primer_end = us_tm_start[1] + len(us_primer)
					primer_pairs.update({(us_primer, ds_primer) : [us_tm_start[0], ds_tm_start[0], amplicon_size, ups_primer_end, ds_tm_start[1]]})
	return primer_pairs
#convert genomic coordinates into relative coordinates
def rel_coord_convert(dhs_coord, ups_seq, motif_coord, motif_strand):
	dhs_len = dhs_coord[1] - dhs_coord[0]
	if ups_seq == 'NA':
		rel_dhs = [0, dhs_len - 1]
		rel_motif = [motif_coord[0] - dhs_coord[0], dhs_coord[1] - motif_coord[1]]
		full_seq = rel_dhs
	else:
		extra_len = len(ups_seq) #sequence added before/after DHS lenght to make rgion xbp
		full_seq = [0, extra_len + dhs_len + extra_len]
		rel_dhs = [extra_len, extra_len + dhs_len - 1]
		motif_start = motif_coord[0] - dhs_coord[0] + extra_len
		rel_motif = [motif_start, motif_start + motif_coord[1] - motif_coord[0]]
	if motif_strand == '-':
		rel_dhs = [full_seq[1] - rel_dhs[1], full_seq[1] - rel_dhs[0]] #check math
		rel_motif = [full_seq[1] - rel_motif[1], full_seq[1] - rel_motif[0]]
	return full_seq, rel_dhs, rel_motif
#filter out primers by best distance from the DHS
def filter_by_dist_from_dhs(primer_list, dhs_coord):
	reduced_list = {}
	for pair, param in primer_list.items():
		if param[3] + len(pair[0]) < dhs_coord[0] and param[4] + len(pair[1]) > dhs_coord[1]:
			reduced_list.update({pair:param})
	#if can't find both primers outside DHS get at least one outside DHS
	if reduced_list == {}:
		for pair, param in primer_list.items():
			if param[3] + len(pair[0]) < dhs_coord[0] or param[4] + len(pair[1]) > dhs_coord[1]:
				reduced_list.update({pair:param})
	return reduced_list
#final filter, select a set from by the longest amplicon size
def long_ampl(primer_list):
	ampl = 0
	primer_set = {}
	for pair, param in primer_list.items():
		if param[2] > ampl:
			ampl = param[2]
			primer_set = {pair:param}
	#try:
	return primer_set
#recalculate distance from motif based on primer start position
def calc_motif_dist(primer_list, rel_motif):
	if 'NA' not in primer_list:
		new_coord = [rel_motif[0] - primer_list[5], primer_list[6] - rel_motif[1]]
		primer_list[5] = new_coord[0]
		primer_list[6] = new_coord[1]
		primer_list[1] = rev_compl(primer_list[1])
	return primer_list

if __name__ == '__main__': #needed for modin to avoid freeze_support() error

	import multiprocessing.popen_spawn_win32
	import modin.pandas as pd

	def seq_for_primer(row):
		ups_dhs, ds_dhs, dhs, dhs_coord, motif_strand = str(row.Extend_seq_ups), str(row.Extend_seq_ds), row.DHSseq, [row.Start, row.End], row.strand
		try:
			motif_coord = [row.Start_p53, row.End_p53]
		except AttributeError:
			motif_coord = [row.Start_ATF4_motif, row.End_ATF4_motif]
		if ups_dhs == 'NA':
			seq = dhs
		else:
			seq = ups_dhs + dhs + ds_dhs
		if motif_strand == '-':  # design primer for the strand that motif is on
			seq = rev_compl(seq)
		gcg = gcg_positions(seq.upper())
		converted_seq = bisulf_convert_Y(seq.upper())
		rel_seq_coord, rel_dhs_coord, rel_motif_coord = rel_coord_convert(dhs_coord, ups_dhs, motif_coord, motif_strand)
		#postitions in the sequence that are outside of motif aka where primers can be designed
		primer_pos_list = list(range(
			0, rel_motif_coord[0])) + list(range(rel_motif_coord[1], rel_seq_coord[1] + 1))
		ranges = gcg_free_regions(gcg, primer_pos_list)
		all_primers = all_possible_primers(ranges, converted_seq, rel_motif_coord)
		paired_primers = pair_primers(all_primers[0], all_primers[1])
		if paired_primers != {}:
			gc_end_primers = gc_3end(paired_primers)
			if gc_end_primers != {}:
				#trying different filters:
				gc_distDHS_primers = filter_by_dist_from_dhs(
					gc_end_primers, rel_dhs_coord)
				if gc_distDHS_primers != {}:
					gc_distDHS_noY_primers = filter_y_3end(
						gc_distDHS_primers, y_3end)
					if gc_distDHS_noY_primers != {}:
						final_primer_set = long_ampl(gc_distDHS_noY_primers)
					else:
						final_primer_set = long_ampl(gc_distDHS_primers)
				else:
					final_primer_set = long_ampl(gc_end_primers)
			else:
				final_primer_set = long_ampl(paired_primers)
		#add something to select best primer
		else:
			# check how many is needed
			final_primer_set = {('NA', 'NA'): ['NA', 'NA', 'NA', 'NA', 'NA']}
		primer_info_list = []
		for k, v in final_primer_set.items():
			primer_info_list.extend(k)  # unlists keys from tuple
			primer_info_list.extend(v)  # unlists keys from list
		primer_info_list = calc_motif_dist(primer_info_list, rel_motif_coord)
		return primer_info_list

	# stops converting 'NA' as Nan
	tbl = pd.read_csv('input.csv', keep_default_na=False)
	#make new table with col for primer seq, tm, dist from motif, dist from DHS, amplicon size, use numpy vectorized
	tbl[['Fwd_primer', 'Rev_primer', 'Tm_fwd', 'Tm_rev', 'Amplicon_size', 'Distance_from_motif_ups', 'Distance_from_motif_ds']] = tbl.apply(seq_for_primer, axis=1, result_type='expand')
	tbl['GC_in_motif'] = tbl['motif_seq_p53'].apply(lambda x: 'Yes'  if 'GC' in x.upper() else 'No')
	tbl['No_of_GC_in_DHS'] = tbl['DHSseq'].apply(lambda x: x.upper().count('GC')) #orientation doesn't matter because GC on + strand is GC on - strand
	tbl.to_csv('output.csv')
	print('End:', datetime.now().strftime("%d-%m-%Y %Hh:%Mm"))



