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


def distance_plot(df, df_neg, df_pos):
	# adjust to take negative log of pvalue
	plt.figure(0)
	print(df_neg)
	print(df_pos)
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

	# plot rectangle around the interesting positive values
	# current_ax.add_patch(Rectangle((pos_d_min, 0), pos_d_max - pos_d_min, pos_p_max, fill='blue'))
	# plot rectanble around the interesting negative values
	# current_ax.add_patch(Rectangle((neg_d_max, 0), neg_d_min - neg_d_max, neg_p_max, fill=False))
	# current_ax.set_ylim(bottom=10**-28)
	# current_ax.set_ylim(top=1)
	current_ax.set_ylim(current_ax.get_ylim()[::-1])

	plt.savefig('volcano.png', dpi=300)


def distribution_plot(mutations):

	# read in the dictionary of slopes
	fit_dict = pic.load(open("fitness_ctrl_1.pkl", "rb"))

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
	plt.savefig('mutant_dist.png', dpi=300)


def main():
	df = pd.read_csv('fitness_ctrl_1_Ttest.csv', sep=',', index_col=1)
	df = df.rename(columns = {'Mean Distance': "mean_dist", 'p-value': "p_value"})
	
	# pick coordinates for the interesting mutations
	# use dummy values for now
	# TODO where should we set these coordinates
	pos_p = .05
	pos_d = .06
	neg_p = .05
	neg_d = -.1

	# identify the interesting mutations, retrieve their slope distributions, plot them
	pos = df[(df["mean_dist"] > pos_d) & (df["p_value"] < pos_p)]
	cool_pos = pos.index.tolist()
	neg = df[(df["mean_dist"] < neg_d) & (df["p_value"] < neg_p)]
	cool_neg = neg.index.values.tolist()

	cool_mut = [i[1:-1] for i in cool_pos] + [j[1:-1] for j in cool_neg]

	distance_plot(df, neg, pos)
	distribution_plot(cool_mut)


if __name__ == '__main__':
    main()
