# load necessary packages
from ete3 import Tree
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from collections import defaultdict, Counter
import numpy as np
from scipy.stats import entropy
from datetime import datetime

# set font and sizes
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

# date extraction function
def extract_date(label):
    """Extract date tuple (YYYY, MM, DD) from tip label in format ...|YYYY-MM-DD"""
    parts = label.split("|")
    if len(parts) >= 2:
        match = re.match(r"(\d{4})-(\d{2})-(\d{2})", parts[-1].strip())
        if match:
            return tuple(map(int, match.groups()))
    return None

# convert to decimal year function
def decimal_year(date_tuple):
    """Convert (YYYY, MM, DD) tuple to decimal year float"""
    year, month, day = date_tuple
    dt = datetime(year, month, day)
    year_start = datetime(year, 1, 1)
    year_end = datetime(year + 1, 1, 1)
    return year + (dt - year_start).total_seconds() / (year_end - year_start).total_seconds()

# identify node ID function
def assign_node_ids(tree):
    for i, node in enumerate(tree.traverse("postorder")):
        node.add_feature("node_id", i)

# reward parsing function
def parse_translate_block(filepath):
    translate = {}
    inside_translate = False
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.lower().startswith("translate"):
                inside_translate = True
                continue
            if inside_translate:
                if line.endswith(";"):
                    inside_translate = False
                else:
                    match = re.match(r"(\d+)\s+(['\"]?)([^,'\";]+)\2,?", line)
                    if match:
                        number, _, name = match.groups()
                        translate[number] = name
    return translate

# assign absolute times based on heights
def set_absolute_times(tree, latest_decimal, latest_tip):
    """
    Assign absolute times to each node in tree.
    Time zero is latest_decimal (most recent tip).
    Node absoluteTime = latest_decimal - distance to latest_tip
    """
    # distance from root to latest tip
    dist_latest_tip = tree.get_distance(latest_tip)
    
    # root absolute time is latest_decimal - dist_latest_tip
    root_abs_time = latest_decimal - dist_latest_tip
    tree.add_feature("absoluteTime", root_abs_time)
    
    # traversal to assign absoluteTime to every node
    for node in tree.traverse("preorder"):
        if node.is_root():
            node.absoluteTime = root_abs_time
        else:
            node.absoluteTime = node.up.absoluteTime + node.dist


# load tree
tree_file = "a_reward_MCC.tree"

# parse tree
with open(tree_file, "r") as f:
    for line in f:
        if line.strip().lower().startswith("tree"):
            tree_line = line.strip()
            break

tree_newick = re.sub(r"^tree\s+\w+\s+=\s+\[&R\]\s+", "", tree_line, flags=re.IGNORECASE)
tree_stripped = re.sub(r"\[&[^\]]*\]", "", tree_newick)
t = Tree(tree_stripped, format=1)

assign_node_ids(t)

translate_map = parse_translate_block(tree_file)
for leaf in t.iter_leaves():
    if leaf.name in translate_map:
        leaf.name = translate_map[leaf.name]


# parse the reward annotations
branch_annotations = re.findall(r"\[&([^\]]+)\]:([0-9.eE+-]+)", tree_newick)
postorder_nodes = list(t.traverse("postorder"))[1:]  # skip root
branch_rewards = {}

if len(postorder_nodes) != len(branch_annotations):
    print(f"Warning: branches ({len(postorder_nodes)}) != annotations ({len(branch_annotations)})")

for node, (ann_str, length_str) in zip(postorder_nodes, branch_annotations):
    branch_length = float(length_str)
    rewards = {}
    if ann_str:
        matches = re.findall(r"gen_genus\.reward_([A-Za-z0-9_]+)=(\d+\.\d+(?:[eE][+-]?\d+)?)", ann_str)
        for loc, val in matches:
            if not loc.endswith("_median"):
                rewards[loc] = float(val)
    parent_id = node.up.node_id if node.up else None
    child_id = node.node_id
    if parent_id is not None:
        branch_rewards[(parent_id, child_id)] = (rewards, branch_length)


# find most recent tip date
tip_dates = [(leaf, extract_date(leaf.name)) for leaf in t.get_leaves()]
tip_dates = [(leaf, d) for leaf, d in tip_dates if d is not None]

if not tip_dates:
    raise ValueError("No valid dates found in tip labels.")

latest_tip, latest_date = max(tip_dates, key=lambda x: x[1])
latest_decimal = decimal_year(latest_date)

print(f"Latest tip: {latest_tip.name}, date: {latest_date}, decimal: {latest_decimal:.3f}")

# assign absolute time
set_absolute_times(t, latest_decimal, latest_tip)

# cumulative rewards by decimal year
reward_by_year = defaultdict(lambda: defaultdict(float))

for (parent_id, child_id), (rewards, branch_length) in branch_rewards.items():
    child_node = next(n for n in t.traverse("postorder") if n.node_id == child_id)
    if not hasattr(child_node, "absoluteTime") or not hasattr(child_node.up, "absoluteTime"):
        continue
    start_time = child_node.up.absoluteTime
    end_time = child_node.absoluteTime
    duration = abs(end_time - start_time)
    if duration == 0:
        continue

    # split branches into smaller decimal steps
    steps = max(1, int(duration * 10))  # 0.1 year resolution
    for i in range(steps):
        frac = i / steps
        time_point = start_time + frac * (end_time - start_time)
        year_bin = round(time_point, 1)
        for k, v in rewards.items():
            reward_by_year[year_bin][k] += v / steps


# define manual colors and labels
manual_colors = {
    "Multi_Aedes": "#000000",
    "Multi_Anopheles": "#593b40",
    "Multi_Other": "#b9a3a3",
    "Uni_Aedes": "#86cacc",
    "Uni_Coquillettidia": "#4f7f84"
}

manual_labels = {
    "Multi_Aedes": "Multivoltine Aedes",
    "Multi_Anopheles": "Multivoltine Anopheles",
    "Multi_Other": "Multivoltine Other",
    "Uni_Aedes": "Univoltine Aedes",
    "Uni_Coquillettidia": "Univoltine Coquillettidia"
}

sorted_years = sorted(reward_by_year.keys())
reward_categories = sorted({k for v in reward_by_year.values() for k in v})

# set colors
colors = [manual_colors.get(cat, "#cccccc") for cat in reward_categories]
labels = [manual_labels.get(cat, cat) for cat in reward_categories]

# create matrix for reward proportions
matrix = []
for y in sorted_years:
    year_data = reward_by_year[y]
    total = sum(year_data.values())
    if total > 0:
        proportions = [year_data.get(cat, 0) / total for cat in reward_categories]
    else:
        proportions = [0] * len(reward_categories)
    matrix.append(proportions)

matrix = np.array(matrix)

# plot rewards across the full time period of the tree
if not any(matrix.flatten()):
    print("All proportions are zero — cannot plot stackplot.")
else:
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.stackplot(sorted_years, matrix.T, labels=labels, colors=colors)
    ax.legend(loc="upper left", fontsize=10, title="Reward category")

    ax.set_title("Lineage A: Mosquito Reward Proportions Over Time")

    ax.set_xlabel("Year", fontweight="bold")
    ax.set_ylabel("Proportion", fontweight="bold")

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False) # remove outer spines

    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.set_xlim(min(sorted_years), max(sorted_years))
    ax.set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig("a_rewards_full.png", facecolor="white", bbox_inches="tight", dpi=600, transparent=False)
    plt.close()

    # plot only rewards from 1997 to present
    filtered_years = [year for year in sorted_years if year >= 1997]
    if filtered_years:
        filtered_matrix = matrix[[sorted_years.index(y) for y in filtered_years]]
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.stackplot(filtered_years, filtered_matrix.T, labels=labels, colors=colors)

        # ax.set_xlabel("Year", fontweight="bold")
        # ax.set_ylabel("Markov Reward Proportion", fontweight="bold")

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.set_xlim(min(filtered_years), max(filtered_years))
        ax.set_ylim(0, 1)
        plt.tight_layout()
        plt.savefig("a_rewards_since1997.png", facecolor="white", bbox_inches="tight", dpi=600, transparent=False)
        plt.close()
    else:
        print(" No data available for years ≥ 1997.")


# statistics to report
# collapse into years (calculating reward proportion by years)
annual_rewards = defaultdict(lambda: defaultdict(float))

for year_bin, categories in reward_by_year.items():
    year_int = int(year_bin)  # calendar year
    for cat, val in categories.items():
        annual_rewards[year_int][cat] += val

# convert to data frame
annual_df = pd.DataFrame.from_dict(annual_rewards, orient="index").fillna(0)
annual_df = annual_df.sort_index()

# calculate annual reward proportions
annual_prop = annual_df.div(annual_df.sum(axis=1), axis=0).fillna(0)

# save files as csvs
annual_df.to_csv("a_rewards_absolute_time.csv")
annual_prop.to_csv("a_rewards_proportions.csv") 

# summary statistics across full time period
total_time_per_category = annual_df.sum()
mean_annual_proportion = annual_prop.mean()

# write summary tables
summary_df = pd.DataFrame({
    "Total_time": total_time_per_category,
    "Mean_annual_proportion": mean_annual_proportion
})
summary_df.to_csv("a_rewards_summary.csv")

# print("\n=== Reportable Statistics ===")
# print("Absolute reward time per year (head):")
# print(annual_df.head())
# print("\nProportions per year (head):")
# print(annual_prop.head())
# print("\nSummary totals:")
# print(summary_df)


