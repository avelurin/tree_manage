from dendropy import Tree

# ==============================
# 1️ LOADING AND SAVING TREES
# ==============================

def load_newick_tree(file_path):
    """Loads a phylogenetic tree from a Newick or Nexus file."""
    tree = Tree.get(path=file_path, schema="newick")
        # Fix taxon names
    for taxon in tree.taxon_namespace:
        taxon.label = taxon.label.replace(".", " ")  # Replace dots with spaces
    return tree

def load_nexus_tree(nexus_file):
    """
    Loads a phylogenetic tree from a Nexus file.

    Args:
        nexus_file (str): Path to the Nexus file.

    Returns:
        Tree: A DendroPy Tree object.
    """
    try:
        tree = Tree.get(path=nexus_file, schema="nexus")
            # Fix taxon names
        for taxon in tree.taxon_namespace:
            taxon.label = taxon.label.replace(".", " ")  # Replace dots with spaces
        print(f"Done: Nexus tree loaded from {nexus_file}")
        return tree
    except Exception as e:
        print(f"Fail: Could not load the Nexus tree. Error: {e}")
        return None

def save_tree(tree, output_file="updated_tree.nwk"):
    """Saves the tree to a Newick file."""
    tree.write(path=output_file, schema="newick")
    print(f"Done: Tree saved to {output_file}")

# ==============================
# 2️ TREE FILTERING
# ==============================

def filter_tree(tree, target_species):
    """
    Filters the tree by keeping only taxa in `target_species`.
    It preserves the original case of taxon names.
    """

    # Create a dictionary to match taxon names without case sensitivity
    target_species_dict = {taxon.lower(): taxon for taxon in target_species}

    # Identify taxa to remove (checking without case sensitivity)
    species_to_remove = [
        taxon.label for taxon in tree.taxon_namespace if taxon.label.lower() not in target_species_dict
    ]

    print(f"Removing {len(species_to_remove)} taxa from the tree.")

    for label in species_to_remove:
        node = tree.find_node_with_taxon_label(label)
        
        if node:
            if node.parent_node:  # Ensure the node is not the root
                tree.prune_subtree(node)
            else:
                print(f"Warning: Skipping root node '{label}', cannot be pruned.")

    return tree


def extract_taxa(tree):
    """
    Extracts all taxa labels from a phylogenetic tree.

    Args:
        tree (Tree): A DendroPy tree object.

    Returns:
        set: A set of unique taxa labels.
    """
    taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}  # Ensure no None values
    print(f"Done: Extracted {len(taxa)} taxa from the tree.")
    return taxa  # Returns as a set for efficient operations


# ==============================
# 3️ FINDING MRCA
# ==============================

def find_mrca_loose(tree, taxon_labels):
    """
    Finds the MRCA (Most Recent Common Ancestor) for a list of taxa, even if extra species are present.
    """
    tree.is_rooted = True

    found_species = [species for species in taxon_labels if tree.find_node_with_taxon_label(species)]

    if len(found_species) < len(taxon_labels):
        print(f"Warning: Some taxa are missing from the tree: {set(taxon_labels) - set(found_species)}")
        return None

    best_node = None
    for node in tree.postorder_node_iter():
        descendants = {t.taxon.label for t in node.leaf_nodes() if t.taxon}
        if set(taxon_labels).issubset(descendants):
            best_node = node
            break

    if best_node:
        print(f"Done: MRCA found for {taxon_labels}.")
        return best_node

    print(f"Fail: Could not find MRCA for {taxon_labels}.")
    return None

# ==============================
# 4️ ADDING A NEW SPECIES
# ==============================

def add_species(tree, new_species, sister_species):
    """
    Adds a new species to the tree as a sister to a given taxon or clade.
    If the species already exists, it does nothing.
    """
    if tree.find_node_with_taxon_label(new_species):
        print(f"Fail: {new_species} already exists in the tree.")
        return

    tree.is_rooted = True

    if isinstance(sister_species, str):
        sister_node = tree.find_node_with_taxon_label(sister_species)
    elif isinstance(sister_species, list):
        sister_node = find_mrca_loose(tree, sister_species)
    else:
        print("Fail: `sister_species` should be a string (single taxon) or a list (clade).")
        return

    if not sister_node:
        print("Fail: Could not find a node for insertion.")
        return

    new_internal_node = sister_node.parent_node.new_child()
    new_internal_node.edge_length = sister_node.edge_length / 2
    sister_node.edge_length = new_internal_node.edge_length
    sister_node.parent_node = new_internal_node

    taxon_new = tree.taxon_namespace.require_taxon(label=new_species)
    new_node = new_internal_node.new_child(taxon=taxon_new)
    new_node.edge_length = new_internal_node.edge_length

    print(f"Done: {new_species} added as sister to {sister_species}.")

# ==============================
# 5️ BALANCING BRANCH LENGTHS
# ==============================

def distance_from_root(node):
    """Computes the distance from the root to a given node."""
    distance = 0.0
    while node is not None and node.parent_node is not None:
        if node.edge_length is not None:
            distance += node.edge_length
        node = node.parent_node
    return distance

def balance_leaf_depths(tree, length_base_sp):
    """
    Adjusts branch lengths so that all leaves are at the same total distance from the root.
    """
    node = tree.find_node_with_taxon_label(length_base_sp)
    if not node:
        print(f"Fail: {length_base_sp} not found in the tree.")
        return

    target_depth = distance_from_root(node)

    for leaf in tree.leaf_node_iter():
        current_depth = distance_from_root(leaf)
        diff = target_depth - current_depth

        if diff > 0:
            leaf.edge_length += diff
        elif diff < 0:
            leaf.edge_length = max(0.1, leaf.edge_length + diff)

    print("Done: Branch lengths adjusted.")


# ==============================
# 6️ ADJUSTING BRANCH LENGTHS BASED ON REFERENCE TREE
# ==============================



def adjust_branch_lengths(old_work_tree, ref_small_tree):
    """
    Replaces a matching subtree in work_tree with ref_small_tree, ensuring proper scaling.

    Parameters:
    - work_tree: DendroPy Tree (the larger tree to modify)
    - ref_small_tree: DendroPy Tree (the reference small tree to insert)

    Returns:
    - Modified work_tree with the adjusted subtree.
    """
    work_tree = old_work_tree.clone()

    def distance_from_root(node):
        """Computes the distance from the root to a given node."""
        distance = 0.0
        while node is not None and node.parent_node is not None:
            if node.edge_length is not None:
                distance += node.edge_length
            node = node.parent_node
        return distance

    def get_subtree_height(node):
        """Returns the maximum distance from the given node to any leaf in its subtree."""
        return max(distance_from_root(leaf) - distance_from_root(node) for leaf in node.leaf_nodes())

    def find_mrca_in_work_tree(work_tree, ref_small_tree):
        """Finds the MRCA (Most Recent Common Ancestor) in work_tree that corresponds to ref_small_tree."""
        real_small_taxa = {node.taxon.label for node in ref_small_tree.leaf_node_iter() if node.taxon}
        large_taxa = {taxon.label for taxon in work_tree.taxon_namespace}
        common_taxa = real_small_taxa.intersection(large_taxa)

        if len(common_taxa) < 2:
            print("Error: Not enough common taxa between small and large trees.")
            return None

        return work_tree.mrca(taxon_labels=common_taxa)

    def scale_branch_lengths(ref_tree, scale_factor):
        """Scales the branch lengths of ref_tree by a given factor."""
        for node in ref_tree.preorder_node_iter():
            if node.edge_length:
                node.edge_length *= scale_factor

    def replace_subtree(work_tree, ref_small_tree, mrca_large):
        """Replaces the subtree in work_tree with ref_small_tree while preserving overall height."""
        parent_large = mrca_large.parent_node
        if not parent_large:
            print("Error: Cannot replace the root of the large tree.")
            return

        # Compute height before replacement
        original_height = get_subtree_height(mrca_large)

        # Get the taxa before replacement
        original_taxa = {node.taxon.label for node in mrca_large.leaf_nodes() if node.taxon}

        # Remove the existing subtree
        parent_large.remove_child(mrca_large)

        # Clone ref_small_tree and scale it
        new_subtree = ref_small_tree.clone()
        new_height = get_subtree_height(new_subtree.seed_node)

        if new_height == 0:
            print("Warning: New subtree height is 0, skipping scaling.")
        else:
            scale_factor = original_height / new_height
            scale_branch_lengths(new_subtree, scale_factor)

        # Attach the new subtree
        new_subtree_root = new_subtree.seed_node
        parent_large.add_child(new_subtree_root)
        new_subtree_root.edge_length = mrca_large.edge_length  # Preserve the connection length

        # Get the taxa after replacement
        updated_taxa = {node.taxon.label for node in new_subtree.leaf_nodes() if node.taxon}

        # Find missing taxa
        removed_taxa = original_taxa - updated_taxa

        print("Subtree replacement completed successfully.")
        return removed_taxa


    # Step 1: Find MRCA in the large tree
    mrca_large = find_mrca_in_work_tree(work_tree, ref_small_tree)
    if not mrca_large:
        print("Error: No matching MRCA found in work_tree.")
        return None

    # Step 2: Replace the subtree
    removed_taxa = replace_subtree(work_tree, ref_small_tree, mrca_large)
    print(f"Removed taxa: {removed_taxa}" if removed_taxa else "No taxa were removed.")


    return work_tree

# ==============================
# 7 VISUALIZATION
# ==============================

def display_ascii_tree(tree):
    """Prints an ASCII representation of the tree."""
    print(tree.as_ascii_plot())
    print(tree)


# ==============================
# 8️ SCALING BRANCH LENGTHS
# ==============================

def scale_branch_lengths(tree, factor, in_place=False):
    """
    Divides each branch length in the tree by the given factor.

    Args:
        tree (Tree): A DendroPy Tree object.
        factor (float): The divisor for branch lengths (must be > 0).
        in_place (bool): If True, modifies the original tree;
                         otherwise, works on a cloned copy.

    Returns:
        Tree: The tree with all branch lengths scaled.
    """
    if factor <= 0:
        raise ValueError("factor must be > 0")

    # work on a clone if not modifying in place
    if not in_place:
        tree = tree.clone()

    # iterate over all nodes in preorder and scale their edge lengths
    for nd in tree.preorder_node_iter():
        if nd.edge_length is not None:
            nd.edge_length = nd.edge_length / factor

    return tree
