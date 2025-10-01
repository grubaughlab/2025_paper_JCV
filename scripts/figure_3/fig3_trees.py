## import necessary packages
import baltic as bt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from io import StringIO
import pandas as pd
import numpy as np
from matplotlib.dates import YearLocator
from matplotlib.ticker import MultipleLocator, AutoMinorLocator 

# check default style
plt.style.use("default")

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white' 

# set general formatting
font = {'family' : 'Helvetica',
'weight' : 'bold',
'size' : 15}
mpl.rcParams.update({"svg.fonttype": 'none', 'text.usetex': False})

# load trees by segment
a_gengenus = bt.loadNewick('/Users/elliebourgikos/Desktop/figures/fig3/a_gengenus_fig_filtered.newick')
b_gengenus = bt.loadNewick('/Users/elliebourgikos/Desktop/figures/fig3/b_gengenus_fig_filtered.newick')

# set absolute time by defining most recent tip date
a_gengenus.setAbsoluteTime(2022.545205479452)
b_gengenus.setAbsoluteTime(2022.4520547945206)

# define mosquito colors
gengenus_colors = {
    "Multi_Aedes": "#000000",
    "Multi_Anopheles": "#41242b",
    "Multi_Other": "#b9a3a3",
    "Uni_Aedes": "#86cacc",
    "Uni_Coquillettidia": "#4f7f84"
}

# define the full state names
gengenus_names = {
    "Multi_Aedes": "Multivoltine Aedes",
    "Multi_Anopheles": "Multivoltine Anopheles",
    "Multi_Other": "Multivoltine Other",
    "Uni_Aedes": "Univoltine Aedes", 
    "Uni_Coquillettidia": "Univoltine Coquillettidia"
}

# create a legend for each lineage
patch_list = []
for abbr, color in gengenus_colors.items():
    full_name = gengenus_names[abbr]
    patch_list.append(mpatches.Patch(color=color, label=full_name))



# create a figure and axis for the tree
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 15

# define two side-by-side plots
fig, ax = plt.subplots(1, 2, figsize=(12, 6), facecolor='white', sharex=True)

# make sure you can iterate based on ax
if isinstance(ax, plt.Axes):
    ax = [ax]

# set figure background to white
fig.patch.set_facecolor('white')

# define the color function
def c_func(k):
    return gengenus_colors.get(k.traits.get('gen_genus', "unknown"), "#c5c6c7")

# define aesthetic attributes of the plots
x_attr=lambda k: k.absoluteTime 
c_func = lambda k: gengenus_colors.get(k.traits.get('gen_genus', "unknown"), "#c5c6c7")
b_func=lambda k: 1 ## set the branch width

# define fixed manual branch width and tip size
fixed_branch_width = 1
fixed_tip_size = 20

# plot trees
for ax_, tree, title in zip(
    ax,
    [a_gengenus, b_gengenus],
    ["Lineage A", "Lineage B"]):

    tree.plotTree(ax_, x_attr=x_attr, width=fixed_branch_width, colour=c_func)
    tree.plotPoints(ax_, x_attr=x_attr, size=fixed_tip_size, zorder=100, colour=c_func)

    # axes
    ax_.xaxis.tick_bottom()
    ax_.xaxis.set_major_locator(MultipleLocator(40))
    ax_.xaxis.set_minor_locator(AutoMinorLocator()) 
    ax_.tick_params(axis='x', which='major', length=4)
    ax_.yaxis.tick_left()
    [ax_.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
    ax_.tick_params(axis='y', size=0)
    ax_.set_yticklabels([])
    ax_.set_title(title)
    ax_.set_facecolor('white')


# add lineage legends to lineage a (ax[0])
ax[0].legend(handles=patch_list, loc="lower left", fontsize=15, frameon=False)

# plot
plt.tight_layout()
plt.savefig("figure4.png", facecolor="white", bbox_inches="tight", dpi=300, transparent=False)
plt.show()