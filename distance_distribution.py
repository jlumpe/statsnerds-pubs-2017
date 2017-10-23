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
	# adjust to take negative log of pvalue
	plt.figure(0)
	plt.scatter(df['mean_dist'], df['p_value'], s=4, color='grey')
	plt.scatter(df_neg['mean_dist'], df_neg['p_value'], s=4, color='orange')
	plt.scatter(df_pos['mean_dist'], df_pos['p_value'], s=4, color='orange')
	plt.axvline(0, color='black', lw=1)
	plt.title('Difference from wt')
	plt.xlabel('distance from wt mean')
	plt.ylabel('p-value (on negative log scale)')

	
	pos_d_max = df['mean_dist'].max()
	neg_d_max = df['mean_dist'].min()

	current_ax = plt.gca()
	current_ax.set_yscale('log')
	current_ax.set_position((.1, .3, .8, .6))
	plt.figtext(.08, .08, "parameters: " + ', '.join(str(i) for i in parameters))

	# plot rectangle around the interesting positive values
	# current_ax.add_patch(Rectangle((pos_d_min, 0), pos_d_max - pos_d_min, pos_p_max, fill='blue'))
	# plot rectanble around the interesting negative values
	# current_ax.add_patch(Rectangle((neg_d_max, 0), neg_d_min - neg_d_max, neg_p_max, fill=False))
	current_ax.set_ylim(bottom=10**-5)
	current_ax.set_ylim(top=1)
	current_ax.set_ylim(current_ax.get_ylim()[::-1])

	plt.savefig(file_prefix + '_volcano_' + version + '.png', dpi=300)


def distribution_plot(mutations, file_prefix, version, parameters):

	# read in the dictionary of slopes
	dict_file = 'mutDict_' + file_prefix + '.p'
	fit_dict = pic.load(open(dict_file, "rb"))

	plt.figure(1)
	bins = np.linspace(-.3, .2, 250)  # should be based on min and max slope
	colors = cm.rainbow(np.linspace(0, 1, len(mutations)))
	axes = plt.gca()
	axes.set_yscale('log')
	#plot the wildtype
	plt.hist(fit_dict['WT'], bins, label='WT', color='grey')
	for i, mutant in enumerate(mutations):
		plt.hist(fit_dict[mutant], bins, label=mutant, color=colors[i])
	plt.legend(loc='upper right')
	plt.title('Growth rates of selected mutants')
	plt.xlabel('growth rate distribution')
	axes.set_position((.1, .3, .8, .6))
	plt.figtext(.08, .08, "parameters: " + ', '.join(str(i) for i in parameters))

	plt.savefig(file_prefix + '_distribution_' + version + '.png', dpi=300)


def main():
	args = sys.argv[1:]
	file_pre= args[0]
	filename = file_pre + '_Ttest.csv'
	ver = args[1]
	params = args[2:]

	df = pd.read_csv(filename, sep=',', index_col=1)
	df = df.rename(columns = {'Mean Distance': "mean_dist", 'p-value': "p_value"})
	
	# if no parameters for mutant selection are passed in, make no selections
	if params == []:
		params = [0, df['mean_dist'].max(), 0, df['mean_dist'].min()]
	pos_p = float(params[0])
	pos_d = float(params[1])
	neg_p = float(params[2])
	neg_d = float(params[3])

	# identify the interesting mutations, retrieve their slope distributions, plot them
	pos = df[(df["mean_dist"] > pos_d) & (df["p_value"] < pos_p)]
	cool_pos = pos.index.tolist()
	neg = df[(df["mean_dist"] < neg_d) & (df["p_value"] < neg_p)]
	cool_neg = neg.index.values.tolist()

	cool_mut = [i[1:-1] for i in cool_pos] + [j[1:-1] for j in cool_neg]

	distance_plot(df, neg, pos, file_pre, ver, params)  # volcano plot
	distribution_plot(cool_mut, file_pre, ver, params) 


if __name__ == '__main__':
    main()
