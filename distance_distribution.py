# output:
# -volcano plot of significance by distance of means
# -histogram of distributions of 'interesting' mutations
# input:
# -dataframe with mean distances and pvalues
# -mutation to growth rate dictionary


from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import numpy as np
import pickle as pic
import pandas as pd


def distance_plot(df, pos_p_min, pos_d_min, neg_p_min, neg_d_min):
	plt.figure(0)
	plt.scatter(df['mean_dist'], df['p_value'])
	plt.axvline(0, color='black', lw=1)
	plt.title('Difference from wt')
	plt.xlabel('distance from wt mean')
	plt.ylabel('p-value')
	
	pos_d_max = df['mean_dist'].max()
	neg_d_max = df['mean_dist'].min()

	current_ax = plt.gca()
	# plot rectangle around the interesting positive values
	current_ax.add_patch(Rectangle((pos_d_min, pos_p_min), pos_d_max - pos_d_min, 1 - pos_p_min, fill=False))
	# plot rectanble around the interesting negative values
	current_ax.add_patch(Rectangle((neg_d_min, neg_p_min), neg_d_max - neg_d_min, 1 - neg_p_min, fill=False))
	plt.savefig('volcano.png')

def distribution_plot(mutations):
	# TODO use real dictionary
	# improve colors
	# set alpha less than one in histo plot to handle overlap

	# read in the dictionary of slopes
	# fit_dict = pic.load("fitness_dict.pkl", "rb")
	fit_dict = {'6_T': [1, 2, 2.05, 3], 
	'8_S': [5, 5, 7, 4 ,5.5, 6.7], 
	'1_M': [-1, -1.1, -1.1, -1, -.9],
	'WT': [.57, .56, .5, .5, .47, .44, .44, .46, .4, .6, .63, .63, .68, .7]
	}
	plt.figure(1)
	bins = np.linspace(-10, 10, 100)  # should be based on min and max slope
	colors = cm.rainbow(np.linspace(0, 1, len(mutations) + 1))
	#plot the wildtype
	plt.hist(fit_dict['WT'], bins, label='WT', color=colors[0])
	for i, mutant in enumerate(mutations):
		plt.hist(fit_dict[mutant], bins, label=mutant, color=colors[i + 1])
	plt.legend(loc='upper right')
	plt.savefig('mutant_dist.png')

def main():
	df = pd.read_csv('dummy_Ttest.csv', sep=',', index_col=1)
	df = df.rename(columns = {'Mean Distance': "mean_dist", 'p-value': "p_value"})
	
	# pick coordinates for the interesting mutations
	# use dummy values for now
	# TODO where should we set these coordinates
	pos_p = .7
	pos_d = 2
	neg_p = .5
	neg_d = -4

	distance_plot(df, pos_p, pos_d, neg_p, neg_d)

	# identify the interesting mutations, retrieve their slope distributions, plot them
	pos = df[(df["mean_dist"] > pos_d) & (df["p_value"] > pos_p)]
	cool_pos = pos.index.tolist()
	neg = df[(df["mean_dist"] < neg_d) & (df["p_value"] > neg_p)]
	cool_neg = neg.index.values.tolist()

	cool_mut = [i[1:-1] for i in cool_pos] + [j[1:-1] for j in cool_neg]


	distribution_plot(cool_mut)


if __name__ == '__main__':
    main()
