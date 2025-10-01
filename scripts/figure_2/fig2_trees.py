## import necessary packages
import baltic as bt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from io import StringIO
import pandas as pd
import numpy as np

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
a_state = bt.loadNexus('/Users/elliebourgikos/Desktop/figures/fig2/a_state_fig.nex')
b_state = bt.loadNexus('/Users/elliebourgikos/Desktop/figures/fig2/b_state_fig.nex')

# set absolute time by defining most recent tip date
a_state.setAbsoluteTime(2023.7397260273972)
b_state.setAbsoluteTime(2023.4109589041095)

# define state colors
state_colors = {
    "CT": "#5a2790",
    "NY": "#c59fd4",
    "MA": "#b35a9e",
    "PA": "#7887c0",
    "95% HPD": "#f8dcaa"  
}

state_colors_b = {
    "CT": "#5a2790",
    "NY": "#c59fd4",
    "MA": "#b35a9e"
}

# define the full state names
state_names = {
    "CT": "Connecticut",
    "NY": "New York",
    "MA": "Massachusetts",
    "PA": "Pennsylvania", 
    "95% HPD": "Connecticut Intro Time"
}

state_names_b = {
    "CT": "Connecticut",
    "NY": "New York",
    "MA": "Massachusetts"
}

# create a legend for each lineage
patch_list = []
for abbr, color in state_colors.items():
    full_name = state_names[abbr]
    patch_list.append(mpatches.Patch(color=color, label=full_name))

patch_list_b = []
for abbr, color in state_colors_b.items():
    full_name = state_names_b[abbr]
    patch_list_b.append(mpatches.Patch(color=color, label=full_name))

# check all node heights to ensure that the node actually exists
# heights = sorted(set([node.height for node in a_state.Objects]))
# for h in heights:
#    print(f"{h:.4f}")


# identify introduction nodes to plot
a_target_height = 229.89368843051
tolerance = 0.01

matches = []
for node in a_state.Objects:
    if abs(node.height - a_target_height) < tolerance:
        matches.append((node.height, node))

if not matches:
    print("No node found near that height.")
    # show nearby heights for sanity check
    for node in a_state.Objects:
        if abs(node.height - a_target_height) < 0.1:
            print(f"Nearby node: {node.height:.10f}")
    exit()

# check that the node was correctly identified
a_intro_node = None
for node in a_state.Objects:
    if abs(node.height - a_target_height) < tolerance:
        a_intro_node = node
        break

if a_intro_node is None:
    print("No node found near that height.")

# identify the hpd area
a_hpd_lower, a_hpd_upper = a_intro_node.traits['height_95%_HPD']


# identify introduction nodes to plot on b tree
b_target_height = 0.0000
tolerance = 0.01

b_intro_node = None
for node in b_state.Objects:
    if abs(node.height - b_target_height) < tolerance:
        b_intro_node = node
        break

if b_intro_node is None:
    print("No node found near that height.")

# identify the hpd area
b_hpd_lower, b_hpd_upper = b_intro_node.traits['height_95%_HPD']


# create a figure and axis for the tree
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 15

# define the color function to return color based on state
def c_func(k):
    return state_colors.get(k.traits.get('state', "unknown"), "#c5c6c7")


# define two side-by-side plots
fig, ax = plt.subplots(1, 2, figsize=(12, 6), facecolor='white', sharex=True)

# make sure you can iterate based on ax
if isinstance(ax, plt.Axes):
    ax = [ax]

# set figure background to white
fig.patch.set_facecolor('white')

# define the color function
def c_func(k):
    return state_colors.get(k.traits.get('state', "unknown"), "#c5c6c7")

# define aesthetic attributes of the plots
x_attr=lambda k: k.absoluteTime ## x coordinate of branches will be time
c_func = lambda k: state_colors.get(k.traits.get('state', "unknown"), "#c5c6c7") ## color tips by state
b_func=lambda k: 1 ## set the branch width

# define fixed manual branch width and tip size
fixed_branch_width = 1
fixed_tip_size = 20

# Convert HPD bounds from node height to absolute time
a_hpd_abs = (
    a_state.mostRecent - a_hpd_upper,
    a_state.mostRecent - a_hpd_lower
)

b_hpd_abs = (
    b_state.mostRecent - b_hpd_upper,
    b_state.mostRecent - b_hpd_lower
)


# plot trees with hpd bars
for ax_, tree, title, hpd_abs in zip(
    ax,
    [a_state, b_state],
    ["Lineage A", "Lineage B"],
    [a_hpd_abs, b_hpd_abs]):

    tree.plotTree(ax_, x_attr=x_attr, width=fixed_branch_width, colour=c_func)
    tree.plotPoints(ax_, x_attr=x_attr, size=fixed_tip_size, zorder=100, colour=c_func)

    # axes
    ax_.xaxis.tick_bottom()
    ax_.yaxis.tick_left()
    [ax_.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
    ax_.tick_params(axis='y', size=0)
    ax_.set_yticklabels([])
    ax_.set_title(title)
    ax_.set_facecolor('white')

    # add the hpd region
    ax_.axvspan(hpd_abs[0], hpd_abs[1], color='#f8dcaa', alpha=0.3, zorder=0)


# add lineage legends to lineage a (ax[0])
ax[0].legend(handles=patch_list, loc="lower left", fontsize=15, frameon=False)

# plot
plt.tight_layout()
plt.savefig("figure2.png", facecolor="white", bbox_inches="tight", dpi=600, transparent=False)
plt.show()
