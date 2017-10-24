# output:
# -volcano plot of significance by distance of means
# -histogram of distributions of 'interesting' mutations
# input:
# -dataframe with mean distances and pvalues
# -mutation to growth rate dictionary
#
# parameters:
# filename (a .csv pandas dataframe)
# version (will be appended to the end of the plot file name)
# OPTIONAL PARAMS (used to select "interesting" mutations)
# optional: the maximum p-value for mutations with a higher growth rate than wt
# optional: the minimum distance between growth rates (for mutations with higher growth rate than wt)
# optional: the maximum p-value for mutations with a lower growth rate than wt
# optional: the minimum distance between growth rates (for mutations with lower growth rate than wt)



from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import numpy as np
import pickle as pic
import pandas as pd
import sys


def distance_plot(df, df_neg, df_pos, file_prefix, version, parameters):
	'''
	plotting p-values against distance from WT mean
	p-values are displayed on a negative log scale
	'''
	plt.figure(0)

	# plot full data frame
	plt.scatter(df['mean_dist'], df['p_value'], s=4, color='grey')

	# highlight the interesting mutations
	plt.scatter(df_neg['mean_dist'], df_neg['p_value'], s=4, color='orange')
	plt.scatter(df_pos['mean_dist'], df_pos['p_value'], s=4, color='orange')

	plt.axvline(0, color='black', lw=1)  # display vertical axis
	plt.title('Difference from wt')
	plt.xlabel('distance from wt mean')
	plt.ylabel('p-value (on negative log scale)')

	current_ax = plt.gca()
	current_ax.set_yscale('log')
	current_ax.set_position((.1, .3, .8, .6))  # make room for parameters in caption
	plt.figtext(.08, .08, "parameters: " + ', '.join(str(i) for i in parameters))

	# plot rectangle around the interesting positive values (not set up for neg log scale)
	# pos_d_max = df['mean_dist'].max()
	# neg_d_max = df['mean_dist'].min()
	# current_ax.add_patch(Rectangle((pos_d_min, 0), pos_d_max - pos_d_min, pos_p_max, fill='blue'))
	# plot rectanble around the interesting negative values
	# current_ax.add_patch(Rectangle((neg_d_max, 0), neg_d_min - neg_d_max, neg_p_max, fill=False))

	current_ax.set_ylim(bottom=10**-5)  # set min for p-value range (this will appear at the top)
	current_ax.set_ylim(top=1)
	current_ax.set_ylim(current_ax.get_ylim()[::-1])  # flip vertical axis

	plt.savefig(file_prefix + '_volcano_' + version + '.png', dpi=300)


def distribution_plot(fit_dict, mutations, file_prefix, version, parameters):
	'''
	plotting distributions of interesting mutations against WT distribution
	will likely need to adjust number of bins, alpha (opacity), ranges, etc.
	'''
	plt.figure(1)
	bins = np.linspace(-.3, .2, 250)  # should be based on min and max slope
	
	axes = plt.gca()
	axes.set_yscale('log')

	#plot the wildtype
	plt.hist(fit_dict['WT'], bins, label='WT', color='grey')
	colors = cm.rainbow(np.linspace(0, 1, len(mutations)))  # generate a colormap for the mutations

	for i, mutant in enumerate(mutations):
		plt.hist(fit_dict[mutant], bins, label=mutant, color=colors[i], alpha=.5)

	plt.legend(loc='upper right')
	plt.title('Growth rates of selected mutants')
	plt.xlabel('growth rate distribution')
	axes.set_position((.1, .3, .8, .6))  # make room at the bottom of the graph for the parameters
	plt.figtext(.08, .08, "parameters: " + ', '.join(str(i) for i in parameters))

	plt.savefig(file_prefix + '_distribution_' + version + '.png', dpi=300)


def fdr(df):
	'''
	determine the cut-off p-value based on a false_rate
	returns dataframe with only significant values
	'''
	false_rate = .05
	hyp_count = len(df.index)
	fdr_valid = [False]*hyp_count
	df = df.sort_values('p_value', na_position='last')
	for i, pv in enumerate(df['p_value']):
		if pv <= false_rate * (1+i) / hyp_count:  # +1 to correct for rank
			fdr_valid[i] = True
		else:
			break

	df['FDR_valid'] = fdr_valid
	return df[df['FDR_valid'] == True]


def reference_df(df, mut_dict):
	'''determine how many barcodes cover each mutant'''
	barcod_cov = []
	for mut in df["Mutation"]:
		barcod_cov.append((len(mut_dict[mut[1:-1]])))  # handle the weird extra quotes
	df['Barcode_cov'] = barcod_cov

	# if you want to return a filtered version (say only mutations with 4+ barcodes):
	# return df[df['Barcode_cov'] >= 4]

	return df



def main():
	args = sys.argv[1:]
	file_pre= args[0]  # sample to be analyzed (eg 'fitness_ctrl_1')
	ver = args[1]  # version will be appended to output files to prevent overwriting
	params = args[2:]

	# read in dataframe of mutant, distance from WT, and p-values
	filename = file_pre + '_Ttest.csv'
	df = pd.read_csv(filename, sep=',')
	df = df.rename(columns = {'Mean Distance': "mean_dist", 'p-value': "p_value"})

	# uses false discovery rate to determine max significant p-value
	# if you want to restrict the rest of the script to only significant values:
	#  replace df with df_sig
	df_sig = fdr(df)

	# read in the dictionary of growth rates by mutation
	dict_file = 'mutDict_' + file_pre + '.p'
	fit_dict = pic.load(open(dict_file, "rb"))

	# if no parameters for mutant selection are passed in, make no selections
	# you can also extract the max p-value from df_sig and use that: df_sig['p_value'].max()
	# to pick the distance, look at df_sig['mean_dist'] or check the output volcano plot
	if params == []:
		params = [0, df['mean_dist'].max(), 0, df['mean_dist'].min()]
	pos_p = float(params[0])
	pos_d = float(params[1])
	neg_p = float(params[2])
	neg_d = float(params[3])

	# identify the interesting mutations, retrieve their slope distributions, plot them
	neg = df[(df["mean_dist"] < neg_d) & (df["p_value"] < neg_p)]
	pos = df[(df["mean_dist"] > pos_d) & (df["p_value"] < pos_p)]
	cool_pos = pos['Mutation'].tolist()
	cool_neg = neg['Mutation'].tolist()

	# add barcode coverage to the dataframe (how many barcodes per mutation)
	# this is for reference - it seems that many 'interesting & significant' mutations
	#  are covered by only two barcodes. These may not be significant at all.
	# these are saved along with their version number
	# The parameters are not written to the file.
	pos = reference_df(pos, fit_dict)
	neg = reference_df(neg, fit_dict)

	pos.to_csv('pos' + '_' + file_pre + ver + '.csv', sep=',')
	neg.to_csv('neg' + '_' + file_pre + ver + '.csv', sep=',')


	# generate a volcano and distribution plot
	# you will need to edit these functions for a good display
	# these are saved as .png with their version file prefix, plot type, and version
	# the parameters are written on the bottom of each plot
	cool_mut = [i[1:-1] for i in cool_pos] + [j[1:-1] for j in cool_neg]

	distance_plot(df, neg, pos, file_pre, ver, params)  # volcano plot
	distribution_plot(fit_dict, cool_mut, file_pre, ver, params) 


if __name__ == '__main__':
    main()
