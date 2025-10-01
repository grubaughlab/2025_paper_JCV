import baltic as bt
import re
import io
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde
import numpy as np
import random
import itertools

# identify outlier tips to exclude
excluded_ct_tips = {"318", "228", "354", "348"}

# translate tip names to numbers
def parse_translate_block(filepath):
    translate_map = {}
    with open(filepath, 'r') as f:
        in_translate = False
        for line in f:
            line = line.strip()
            if line.startswith("Translate"):
                in_translate = True
                continue
            if in_translate:
                if line == ";":
                    break
                line = line.rstrip(',')
                m = re.match(r'(\d+)\s+(.+)', line)
                if m:
                    tip_id, full_label = m.groups()
                    translate_map[tip_id] = full_label
    return translate_map

# extract all ct tips (without the outliers)
def get_ct_tip_ids(translate_map, exclude_ids):
    ct_tip_ids = []
    for tip_id, label in translate_map.items():
        state = label.split('|')[0]
        if state == "CT" and tip_id not in exclude_ids:
            ct_tip_ids.append(tip_id)
    return ct_tip_ids

# clean newick strings from beast all trees
def extract_clean_newick(raw_line):
    if "=" in raw_line:
        raw_line = raw_line.split("=", 1)[1].strip()
    raw_line = re.sub(r"\[&[^\]]*\]", "", raw_line)
    start = raw_line.find("(")
    if start == -1:
        return None
    return raw_line[start:]

# function for leaf names under node of interest
def get_leaf_names(node):
    if hasattr(node, 'getExternal') and callable(node.getExternal):
        leaves = node.getExternal()
        if not leaves:
            return set()
        if isinstance(leaves[0], str):
            return set(leaves)
        else:
            return set(leaf.name for leaf in leaves)

    def recursive_leaves(n):
        if not hasattr(n, 'children') or not n.children:
            if isinstance(n, str):
                return [n]
            else:
                return [getattr(n, 'name', None)]
        else:
            leaves_list = []
            for child in n.children:
                leaves_list.extend(recursive_leaves(child))
            return leaves_list

    leaves = recursive_leaves(node)
    return set([l for l in leaves if l is not None])

# identify mrca node height of monophyletic ct clade
def get_CT_mrca_node_height(tree, ct_tip_codes):
    try:
        ct_tips = [n for n in tree.Objects if n.is_leaf() and n.name in ct_tip_codes]
        if len(ct_tips) < 2:
            print("Not enough CT tips to determine MRCA")
            return None

        # find mrca even if not monophyletic
        mrca_node = tree.commonAncestor(ct_tips)
        return mrca_node.height

    except Exception as e:
        print(f"Error finding MRCA of CT tips: {e}")
        return None


# load in files and run function
file_path = "a_state_dta_ALL.trees"
sample_size = 1000 # only sampling 1000 of 7000 posterior trees, randomly
random_seed = 123
most_recent_sampling_year = 2023.7397260273972  # most recent a lineage sample date

random.seed(random_seed)

translate_map = parse_translate_block(file_path)
print(f"Parsed {len(translate_map)} tips from Translate block")

ct_tip_codes = get_ct_tip_ids(translate_map, excluded_ct_tips)
print(f"Number of CT tips (excluding outliers): {len(ct_tip_codes)}")

with open(file_path, "r") as f:
    all_raw_trees = [line.strip() for line in f if line.startswith("tree")]

total_trees = len(all_raw_trees)
print(f"Total trees in file: {total_trees}")

if sample_size > total_trees:
    raise ValueError("Sample size exceeds number of trees in the file.")

sampled_trees = random.sample(all_raw_trees, sample_size)
print(f"Randomly sampled {len(sampled_trees)} trees with seed {random_seed}\n")

converted_times = []

for i, raw in enumerate(sampled_trees):
    cleaned = extract_clean_newick(raw)
    if not cleaned:
        print(f"Skipping tree {i}: malformed Newick")
        continue

    try:
        t = bt.loadNewick(io.StringIO(cleaned))

        # most recent tip date
        tip_heights = [leaf.height for leaf in t.getExternal() if hasattr(leaf, 'height')]
        if not tip_heights:
            print(f"Tree {i}: No tip heights found")
            continue
        max_tip_height = max(tip_heights)

        # find mrca height of the ct clade
        mrca_height = get_CT_mrca_node_height(t, ct_tip_codes)
        if mrca_height is not None:
            year = most_recent_sampling_year - (max_tip_height - mrca_height)
            converted_times.append(year)
        else:
            print(f"Tree {i}: No valid CT MRCA found")
    except Exception as e:
        print(f"Tree {i} failed: {e}")

print(f"\nCollected {len(converted_times)} adjusted MRCA node heights")
if converted_times:
    print(f"First 5 adjusted heights (years): {converted_times[:5]}") # check heights to make sure they look correct

# plot kde of the node heights (adjust to time)
converted_times = np.random.normal(1950, 10, 300)

if converted_times.any():
    x = np.linspace(min(converted_times) - 10, max(converted_times) + 10, 500)
    kde = gaussian_kde(converted_times)
    y = kde(x)

    # Build a soft gradient colormap from muted red to transparent
    cmap = LinearSegmentedColormap.from_list(
        "soft_red", [(1, 1, 1, 0), ("#c5454e")], N=256
    )

    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    ax.patch.set_alpha(0)  # transparent background

    # Fill with gradient by plotting many thin horizontal strips
    for i in range(len(x) - 1):
        plt.fill_between(
            [x[i], x[i + 1]],
            0,
            [y[i], y[i + 1]],
            color=cmap(0.6),  # adjust intensity of gradient
            alpha=0.8
        )

    # KDE line on top
    plt.plot(x, y, color="#5a2790", linewidth=3)

    # aesthetics
    plt.xlim(1900, 2000)
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    plt.savefig("a_kde_gradient.png", dpi=600, transparent=True)
    plt.show()



