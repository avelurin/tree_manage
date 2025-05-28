# Tree Manager (tree_manajer)

A Python toolkit for phylogenetic tree manipulation and analysis using DendroPy.

## Features

### 1. Tree Loading & Saving
- Load trees from Newick and Nexus formats
- Save trees to Newick format
- Automatic handling of taxon name formatting

### 2. Tree Manipulation
- Filter trees by keeping specific taxa
- Extract taxa labels from trees
- Find Most Recent Common Ancestor (MRCA)
- Add new species to existing trees
- Balance branch lengths
- Scale branch lengths

### 3. Advanced Features
- Adjust branch lengths based on reference trees
- ASCII tree visualization
- Support for both rooted and unrooted trees

## Installation

1. Create a virtual environment (recommended):
```bash
python -m venv env
source env/bin/activate  # On Linux/Mac
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Quick Start

```python
from tree_manajer import load_newick_tree, filter_tree, save_tree

# Load a tree
tree = load_newick_tree("your_tree.nwk")

# Filter the tree to keep only specific species
target_species = ["Species1", "Species2", "Species3"]
filtered_tree = filter_tree(tree, target_species)

# Save the modified tree
save_tree(filtered_tree, "filtered_tree.nwk")
```

## Function Documentation

### Loading and Saving
- `load_newick_tree(file_path)`: Load a tree from Newick format
- `load_nexus_tree(nexus_file)`: Load a tree from Nexus format
- `save_tree(tree, output_file)`: Save a tree to Newick format

### Tree Operations
- `filter_tree(tree, target_species)`: Keep only specified taxa
- `extract_taxa(tree)`: Get all taxa labels
- `find_mrca_loose(tree, taxon_labels)`: Find MRCA of specified taxa
- `add_species(tree, new_species, sister_species)`: Add a new species
- `balance_leaf_depths(tree, length_base_sp)`: Balance leaf depths
- `adjust_branch_lengths(old_work_tree, ref_small_tree)`: Adjust branch lengths
- `scale_branch_lengths(tree, factor, in_place=False)`: Scale all branch lengths

## Notes
- All functions handle errors gracefully with appropriate error messages
- Branch length manipulations preserve tree topology
- Taxa names are automatically cleaned (dots replaced with spaces)
