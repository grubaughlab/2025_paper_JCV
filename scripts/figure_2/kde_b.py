import baltic as bt
import re
import io
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde
import numpy as np
import random
import itertools

# clean newick strings from beast all trees
def extract_clean_newick(raw_line):
    if "=" in raw_line:
        raw_line = raw_line.split("=", 1)[1].strip()
    raw_line = re.sub(r"\[&[^\]]*\]", "", raw_line)
    start = raw_line.find("(")
    if start == -1:
        return None
    return raw_line[start:]

# load in files and set parameters
file_path = "b_state_dta.trees"
sample_size = 1000
random_seed = 123 # keep seed across analysis
most_recent_sampling_date = 2023.4109589041095

random.seed(random_seed)

# load in trees and identify root node to target
with open(file_path, "r") as f:
    all_raw_trees = [line.strip() for line in f if line.startswith("tree")]

total_trees = len(all_raw_trees)
print(f"Total trees in file: {total_trees}")

if sample_size > total_trees:
    raise ValueError("Sample size exceeds number of trees in the file.")

sampled_trees = random.sample(all_raw_trees, sample_size)
print(f"Randomly sampled {len(sampled_trees)} trees\n")

converted_root_times = []

for i, raw in enumerate(sampled_trees):
    cleaned = extract_clean_newick(raw)
    if not cleaned:
        print(f"Skipping tree {i}: malformed Newick")
        continue

    try:
        t = bt.loadNewick(io.StringIO(cleaned))
        root_height = t.root.height
        tip_heights = [leaf.height for leaf in t.getExternal()]
        max_tip_height = max(tip_heights)
        calendar_root_time = most_recent_sampling_date - (max_tip_height - root_height)
        converted_root_times.append(calendar_root_time)

    except Exception as e:
        print(f"Tree {i} failed: {e}")

print(f"\nCollected {len(converted_root_times)} root times")

# plot kde of the node heights (adjust to time)
if converted_root_times:
    x = np.linspace(min(converted_root_times) - 1, max(converted_root_times) + 1, 500)
    kde = gaussian_kde(converted_root_times)
    y = kde(x)

    lower, upper = np.percentile(converted_root_times, [2.5, 97.5])
    median = np.median(converted_root_times)

    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    ax.patch.set_alpha(0)  # transparent background

    # Gradient fill under curve (muted red → transparent)
    cmap = LinearSegmentedColormap.from_list(
        "soft_red", [(1, 1, 1, 0), "#c5454e"], N=256
    )

    for i in range(len(x) - 1):
        plt.fill_between(
            [x[i], x[i + 1]],
            0,
            [y[i], y[i + 1]],
            color=cmap(0.7),  # darker near curve
            alpha=0.8
        )

    # KDE line
    plt.plot(x, y, color="#5a2790", linewidth=3)

    # aesthetics
    plt.xlim(1750, 1950)
    plt.xticks([])
    plt.title("")
    plt.xlabel("")
    plt.yticks([])
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig("b_kde_plot.png", dpi=600, transparent=True)
    plt.show()
else:
    print("No valid root times collected, skipping KDE plot.")
