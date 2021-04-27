import sys
import operator
import ete3

from scripts_of import trees2ologs_of
from scripts_of import parallel_task_manager

class Finalise(object):
    def __enter__(self):
        pass
    def __exit__(self, type, value, traceback):
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()

def GetRoot(tree, species_tree_rooted, GeneToSpecies):
        roots = trees2ologs_of.GetRoots(tree, species_tree_rooted, GeneToSpecies)
        m = roots[0][0]
        print(roots)
        roots = [r for score, r in roots if score==m]
        print("%d qual-scoring roots" % len(roots))
        if len(roots) > 0:
            root_dists = [r.get_closest_leaf()[1] for r in roots]
            i, _ = max(enumerate(root_dists), key=operator.itemgetter(1))
            print("Using: " + "|".join(roots[i].get_leaf_names()))
            return roots[i]
        else:
            return None # single species tree

def process_trees(sp_tree_fn, gene_trees, GeneToSpecies = trees2ologs_of.GeneToSpecies_dash):
    species_tree_rooted = ete3.Tree(sp_tree_fn, format=3)
    for t_fn in gene_trees:
        tree = ete3.Tree(t_fn)
        # Root tree
        root = GetRoot(tree, species_tree_rooted, GeneToSpecies)
        if root == None: 
            raise Exception("No root")
        if root != tree:
            tree.set_outgroup(root)
        print("Tree to analyse:")
        print(tree)
        # get levels in species tree
        speciesObserved = set([GeneToSpecies(g) for g in tree.get_leaf_names()])
        if len(speciesObserved) == 1:
            print("Only one species, exiting")
            return
        n = species_tree_rooted
        children = n.get_children()
        leaves = [set(ch.get_leaf_names()) for ch in children]
        have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]
        while sum(have) < 2:
            n = children[have.index(True)]
            children = n.get_children()
            leaves = [set(ch.get_leaf_names()) for ch in children]
            have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]
        # print(leaves) # leaves are the two splits we're looking for
        # There might be more than two groups, consider each in turn as A and the rest as B
        n_ops = len(leaves)
        if n_ops != 2:
            raise Exception("Non-binary tree")
        if n_ops == 2:
            # the two options are identical
            n_ops = 1
            i = 0
        # for i in xrange(n_ops):
        t1 = leaves[i]
        t1 = t1.intersection(speciesObserved)
        n_t1_req = 2 if len(t1) >= 2 else 1
        t2 = set.union(*[l for j,l in enumerate(leaves) if j!=i])
        t2 = t2.intersection(speciesObserved)
        n_t2_req = 2 if len(t2) >= 2 else 1
        trees2ologs_of.StoreSpeciesSets(tree, GeneToSpecies) # ingroup/outgroup identification
        print("Splitting requirements:")
        print((t1, t2))
        print((n_t1_req, n_t2_req))
        # find all possible locations in the gene tree at which the root should be
        T = {True,}
        F = {False,}
        TF = set([True, False])
        print("Partitions:")
        partitions = set()
        for m in tree.traverse('postorder'):
            if m.is_leaf():
                m.add_feature("done", False)
                continue
            if all(ch.done for ch in m.get_children()):
                m.add_feature("done", True)
                continue
            print([len(ch) for ch in m.get_children()])
            print([ch.done for ch in m.get_children()])
            if len(m.sp_down.intersection(t1))>=n_t1_req and len(m.sp_down.intersection(t2))>=n_t2_req:
                if any(ch.done for ch in m.get_children()):
                    print("*** Overriding ***")
                    overridden = [ch for ch in m.get_children() if ch.done]
                    print("Adding", len(m), m.get_leaf_names())
                    to_remove = []
                    for ch in overridden:
                        remove_leaves = set(ch.get_leaf_names())
                        # traverse down to the previous partitions that were stored and remove them
                        # or just remove anyhing that has these genes
                        for p in partitions:
                            if len(remove_leaves.intersection(p.get_leaf_names())) > 0:
                                print("Removing", len(p), "|".join(p.get_leaf_names()))
                                to_remove.append(p)
                    for p in to_remove:
                        partitions.remove(p)
                print(len(m))
                print(m)
                partitions.add(m)
                m.add_feature("done", True)
            else:
                m.add_feature("done", False)
        L = [len(p) for p in partitions]
        print(L)
        for p in partitions:
            print("|".join(p.get_leaf_names()))
        print((len(tree), sum(L)))

if __name__ == "__main__":
    with Finalise():
        sp_tree = sys.argv[1]
        gene_trees = sys.argv[2]
        try:
            t = ete3.Tree(gene_trees)
            gene_trees = [gene_trees, ]
        except IOError:
            gene_trees = list(glob.glob(gene_trees + "/*"))
        except Exception as e:
            print("Unexpected exception")
            print(type(e))
            print(e)
        process_trees(sp_tree, gene_trees)
