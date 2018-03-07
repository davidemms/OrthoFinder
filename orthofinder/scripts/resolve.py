# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 14:15:17 2017

@author: david
"""  
import parallel_task_manager

import glob
import argparse
import tree as tree_lib

import trees2ologs_of as om1
    
def DetachAndCleanup(top, n):
    """
    top - the node that all the nodes in the current analysis will be below
    nUp - the node above n originally, from here to nTop needs the cached data updated
    """
    nUp = n.up
    n = n.detach()
    ch = nUp.get_children()
    if len(ch) == 1:
        if nUp.is_root():
            # can't delete the root, need to detach the child nodes and reattach them to the root
            for c in ch[0].get_children():
                nUp.add_child(c.detach())
            ch[0].delete()
            return n, nUp, top
        else:
            if top == nUp: top = nUp.up 
            nUpNew = nUp.up
            nUp.delete()
            return n, nUpNew, top
    return n, nUp, top

def GraftAndUpdate(top, n, s):
    """
    Remove n from its location in the tree and make it sister to s. Update only the sp_down data. NOT the sp_up data.
    Args:
        nTop - node being reconcilled
        n - node to move
        s - node to make n sister to.
    Returns:
        The root of the tree
    Implementation:
        n and s could be anywhere in relation to one another, need to change all the species sets that could be affected
    """
    n, nUp, top = DetachAndCleanup(top, n)
    parent = s.up
    if parent == None:
        new = tree_lib.TreeNode()
        new.add_feature("sp_up", set())
        parent = new
        top = new
        new.add_child(n)
        s = s.detach()
        new.add_child(s)
        s.sp_up = n.sp_down
    else:
        new = parent.add_child()
        new.add_child(n)
        s = s.detach()
        new.add_child(s)
    
    # sort out distances, new node will have a default distance of 1.0, correct this
    new.dist = 0.1*s.dist
    s.dist   = 0.9*s.dist
    
    # sort out sp_down data - can do a more efficient routine later if necessary (updating only what needs updating). Otherwise:
    # all nodes on path from nTop to s need sp_down updating
    new.add_feature("sp_down", s.sp_down.union(n.sp_down))
    parent.sp_down = parent.sp_down.union(n.sp_down)
    if not parent == top:
        r = parent.up
        while r != top:
            r.sp_down = set.union(*[ch.sp_down for ch in r.get_children()])
            r = r.up
    # all nodes on path from nTop to n need sp_down updating
    r = nUp
    if r != None:
        while r != top:
            r.sp_down = set.union(*[ch.sp_down for ch in r.get_children()])
            r = r.up
    return new.get_tree_root()
    
def GraftTripartAndUpdate(nTop, s1, s2, p):
    """
    Make s1 and s2 sister clades and hang the from parent node p
    Args:
        s1, s2 - nodes to be made sisters of each other
        p - node to hang (s1,s2) from
    """
    d1 = s1.up.dist
    s1, s1Parent , nTop = DetachAndCleanup(nTop, s1)   
    d2 = s2.up.dist
    s2, s2Parent, nTop = DetachAndCleanup(nTop, s2)   
    new = p.add_child()
    new.dist = 0.5*(d1+d2)
    new.add_child(s1)
    new.add_child(s2)
    new.add_feature("sp_down", s1.sp_down.union(s2.sp_down))
    n = s1Parent
    while n != nTop:
        n.sp_down = set.union(*[ch.sp_down for ch in n.get_children()])
        n = n.up
    n = s2Parent
    while n != nTop:
        n.sp_down = set.union(*[ch.sp_down for ch in n.get_children()])
        n = n.up
    nTop.sp_down = set.union(*[ch.sp_down for ch in nTop.get_children()])
    return new.get_tree_root()
        
def ContainsMonophyletic(e, O, i, I):
    """
    Check if O is a monophyletic clade in e. Returns (isMono, n) where n is the node if O is monophletic otherwise n is None
    Args:
        e - the node to check
        O - the species set to look for
        i - the current iteration 
        I - the one-past-the-last iteration. Number of checks will be M-m
    Returns:
        qSuccess, node, depth
        where:
            qSuccess - True if found species set O to be monophyletic
            node - node of the monophyletic clade if qSuccess else the node it made it down to
            depth - depth of the node, 0=(e is O), 1=(it's child is O), etc.
    """
    E = e.sp_down
    if E == O:
        return True, e, i     
    i+=1
    if i == I or e.is_leaf(): return False, e, i-1
    ch = e.get_children()
    inChild = [bool(O & c.sp_down) for c in ch]
    if sum(inChild) > 1: return False, e, i-1         # bipart : cheange to if sum(inChild) == len(ch): return False, None 
    return ContainsMonophyletic(ch[inChild.index(True)], O, i, I)  # bipart need to check more than one potentially

def check_monophyly(node, taxa):
    """ This should be used as a wrapper of the ete method since that method returns false for a single node. It should return true
    Args:
        node - the node under which to search
        taxa - the list of taxa
    Returns
        bool - are the taxa monophyletic
    """
    if len(taxa) == 1:
        return (list(taxa)[0] in node)
    else: 
        return node.check_monophyly(target_attr='name', values=taxa)[0]

t = True
f = False

# dA = 1, dB = 0
# Subcases (S&O, S&U)
c10_a = {(f,f), (t,f), (t,t)}      
c10_b = {(f,t),}

# dA = 1, dB = 1
# Subcases (S&O, S&W, S&U)
c11_a = {(f,f,f), (t,t,t), (t,t,f), (t,f,t)}
c11_b = {(f,t,f), (f,f,t)}
c11_c = {(f,t,t),}
 
# dA = 2, dB = 0
# Subcases (S&V, S&Y, S&O, V&Y)
#c20_a = {(f,f,f,f),}
#c20_b = {(f,f,f,t), (f,f,t,f), (f,f,t,t), (f,t,f,t), (f,t,t,t), (t,f,t,t), (t,t,f,f), (t,t,f,t), (t,t,t,t)}
c20_a = {(f,f,f,f),(f,f,t,f),}
c20_b = {(f,f,f,t), (f,f,t,t), (f,t,f,t), (f,t,t,t), (t,f,t,t), (t,t,f,f), (t,t,f,t), (t,t,t,t)}
c20_c = {(f,t,f,f), (f,t,t,f), (t,f,f,f), (t,f,f,t), (t,f,t,f)}

def resolve(n, M):
    """This is an improved but also more staightforward implementation of resolve. It also considers far more sub-cases.
    It comes from writing up and generalising the first version.  
    
    Implementation:
        lower case letters are nodes, uppercase letters are species sets
    """
    if "recon" in n.features:
        sis = n.recon
        n.del_feature('recon')
        if check_monophyly(n, sis[0]) and check_monophyly(n, sis[1]):
            s0 = n.get_common_ancestor(sis[0]) if len(sis[0]) > 1 else (n&sis[0][0])
            s1 = n.get_common_ancestor(sis[1]) if len(sis[1]) > 1 else (n&sis[1][0])
            return GraftAndUpdate(n, s0, s1)
    elif "recon_2" in n.features:
        """ (destination, target, move_out)"""
        moves = n.recon_2
        n.del_feature('recon_2')
        if check_monophyly(n, moves[0]) and check_monophyly(n, moves[1]) and check_monophyly(n, moves[2]):
            s0 = n.get_common_ancestor(moves[0]) if len(moves[0]) > 1 else (n&moves[0][0])
            s1 = n.get_common_ancestor(moves[1]) if len(moves[1]) > 1 else (n&moves[1][0])
            d = s1.get_sisters()[0].dist
            tree = GraftAndUpdate(n, s1, s0)
            tree.dist = d
            s3 = tree.get_common_ancestor(moves[0] + moves[1])
            sister = s3.up.up
            top = sister.up if sister.up != None else sister
            s2 = n.get_common_ancestor(moves[2]) if len(moves[2]) > 1 else (n&moves[2][0])
            return GraftAndUpdate(top, s2, sister)
        
    ch = n.get_children()
    if len(ch) != 2: return n.get_tree_root()    #bipart
    X,Y = [c.sp_down for c in ch]
    O = X&Y
    if not O: return n.get_tree_root()# no overlap, nothing to do
    if X == Y: return n.get_tree_root()  # identical species sets, they're paralogues
    sis = n.get_sisters()
    s = sis[0] if len(sis) == 1 else None
    S = set() if (n.is_root() or len(n.get_sisters()) > 1) else set.union(*[nother.sp_down for nother in n.up.get_children() if nother != n]) # works for non-binary too
    
    successA, nA, dA = ContainsMonophyletic(ch[0], O, 0, 3)   
    successB, nB, dB = ContainsMonophyletic(ch[1], O, 0, 3)    
    if dA+dB > 2:
        if dA == 2:
            successA, nA, dA = ContainsMonophyletic(ch[0], O, 0, 2)   
        if dB == 2:
            successB, nB, dB = ContainsMonophyletic(ch[1], O, 0, 2)                  
        
    if (dA==0 and dB==1) or (dA==1 and dB==0): 
        # lowest letter is the one that branches
        if dA == 1:
            a, b = ch
        else:
            b, a = ch
        ch = a.get_children()
        u, v = ch
        U,V = [c.sp_down for c in ch]
        if O & U: 
            v,u = ch
            U = u.sp_down
        case = (bool(S&O), bool(S&U))  
        if case in c10_a:
            return GraftAndUpdate(n, b, v)
        elif case in c10_b:
            n.up.add_feature("recon", (b.get_leaf_names(),s.get_leaf_names())) 
            return n.get_tree_root()
        return n.get_tree_root()
    elif (dA==1 and dB==1): 
        cha = ch[0].get_children()
        U,V = [c.sp_down for c in cha]
        if O & U: 
            v, u = cha
            TEMP = U
            U = V
            V = TEMP
        else:
            u, v = cha
        chb = ch[1].get_children()
        W,X = [c.sp_down for c in chb]
        if O & W: 
            x, w = chb
            TEMP = X
            X = W
            W = TEMP
        else:
            w, x = chb
        case = (bool(S&O), bool(S&W), bool(S&U))
        if case in c11_a:
            return GraftAndUpdate(n, x, v)
        elif case in c11_b:
            if W&S:
                w = u
                x = v
            n.up.add_feature("recon_2", (x.get_leaf_names(),s.get_leaf_names(),w.get_leaf_names())) 
            return n.get_tree_root()
        elif case in c11_c:
            n.up.add_feature("recon", (s.get_leaf_names(), w.get_leaf_names()))     
            return n.get_tree_root()
        return n.get_tree_root()
    elif (dA==2 and dB==0) or (dA==0 and dB==2): 
        if dA == 2:
            a,b = ch
        else:
            b,a = ch
        cha = a.get_children()
        u, v = cha
        U,V = [c.sp_down for c in cha]
        if O & V: 
            v, u = cha
            U =  u.sp_down
            V =  v.sp_down
        chu = u.get_children()[:2]  # fix
        y,z = chu
        Y,Z = [c.sp_down for c in chu]
        if O & Y:
            z,y = chu
            Y = y.sp_down
        case = (bool(S&V), bool(S&Y), bool(S&O), bool(V&Y))
        if case in c20_a:
            return GraftAndUpdate(n, b, z)
        elif case in c20_b:
            return GraftAndUpdate(n, v, b)
        elif case in c20_c:
            n.up.add_feature("recon", (b.get_leaf_names(), s.get_leaf_names()))     
            return n.get_tree_root()
        return n.get_tree_root()
    return n.get_tree_root()

def NumberOfOrthologues(tree, GeneToSpecies):
    import numpy as np
    import itertools
    species = list(set(map(GeneToSpecies, tree.get_leaf_names())))
    genes = tree.get_leaf_names()
    gDict = {gene:i for i,gene in enumerate(genes)}
    sDict = {sp:i for i,sp in enumerate(species)}
    orthologues = np.zeros((len(genes), len(species)))
    nOrtho = 0
    for n in tree.traverse('postorder'):
        if n.is_leaf(): continue
        ch = n.get_children()        
        if len(ch) == 2: 
            l0 = ch[0].get_leaf_names()
            l1 = ch[1].get_leaf_names()
            s0 = {GeneToSpecies(l) for l in l0}
            s1 = {GeneToSpecies(l) for l in l1}
            if s0&s1: continue
            nOrtho += len(l0)*len(l1)
            for g0 in l0:
                for s in s1:
                    orthologues[gDict[g0], sDict[s]] = 1
            for g1 in l1:
                for s in s0:
                    orthologues[gDict[g1], sDict[s]] = 1
        elif len(ch) > 2:
            for ch0, ch1 in itertools.combinations(ch, 2):
                l0 = ch0.get_leaf_names()
                l1 = ch1.get_leaf_names()
                s0 = {GeneToSpecies(l) for l in l0}
                s1 = {GeneToSpecies(l) for l in l1}
                if s0&s1: continue
                nOrtho += len(l0)*len(l1)
                for g0 in l0:
                    for s in s1:
                        orthologues[gDict[g0], sDict[s]] = 1
                for g1 in l1:
                    for s in s0:
                        orthologues[gDict[g1], sDict[s]] = 1
                
#    print(nOrtho)
#    N = (len(species)-1)*len(genes)
#    print((sum(sum(orthologues)),float(N)))
#    print(sum(sum(orthologues)))

class Finalise(object):
    def __enter__(self):
        pass
    def __exit__(self, type, value, traceback):
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()
        
def Resolve_Main(trees_fn_or_dir, species_tree_rooted_fn, GeneToSpecies, qTest):   
    """
    Resolves the single tree or trees in the directory trees_fn. If no species tree is provided then the gene tree is assumed to 
    be rooted
    Args:
        trees_fn - tree filename or directory containing trees
        species_tree_rooted_fn - species tree used to root the tree. If None then the gene tree is assumed to be rooted already
    """
    if species_tree_rooted_fn != None:
        species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
        
    qDir = True
    try:
        tree = tree_lib.Tree(trees_fn_or_dir)
        qDir = False
    except:
        try:
            tree = tree_lib.Tree(trees_fn_or_dir, format=3)
            qDir = False
        except:
            pass
    trees = glob.glob(trees_fn_or_dir + "/*") if qDir else [trees_fn_or_dir] 
    for trees_fn in trees:
        # Root using species tree if provided, otherwise tree should have been rooted already
        tree = tree_lib.Tree(trees_fn)
        if len(tree) == 1: continue
        if species_tree_rooted_fn != None:
            tree.prune(tree.get_leaf_names())
            roots = om1.GetRoots(tree, species_tree_rooted, GeneToSpecies)
            if len(roots) == 0: continue
            # Pick the first root for now
            root = roots[0]
            if root != tree:
                tree.set_outgroup(root)
            tree.write(outfile=trees_fn+"_rooted.tre")
        om1.StoreSpeciesSets(tree, GeneToSpecies)  
        if qTest:
            # Just perform a single reconciliation (including making move that is flagged but not initially applied in the reconciliation_
            for iNode, n in enumerate(tree.traverse()):
                n.dist = float(iNode * 0.01)
            ch = tree.get_children()
            if len(ch) != 2:
                print("ERROR: Not binary")
                return tree
            print(tree)
            n = ch[0] if len(ch[0]) > len(ch[1]) else ch[1]
            tree = resolve(n, GeneToSpecies)
            print(tree)
            if "recon" in tree.features or "recon_2" in tree.features:
                tree = resolve(tree, GeneToSpecies)
                print(tree)
        else:
            # Perform full reconciliation
            for n in tree.traverse("postorder"):
                tree = resolve(n, GeneToSpecies)
            NumberOfOrthologues(tree, GeneToSpecies)    
        tree.write(outfile=(trees_fn + ".rec.tre"))           
           
if __name__ == "__main__":
    with Finalise():
        parser = argparse.ArgumentParser()
        parser.add_argument("gene_tree")
        parser.add_argument("-r", "--rooted_species_tree")
        parser.add_argument("-s", "--separator", choices=("dot", "dash", "second_dash", "3rd_dash", "hyphen"), help="Separator been species name and gene name in gene tree taxa")
        parser.add_argument("-t", "--test", action="store_true", help="Perform a single operation on the largest node one step down--allows testing of the method")
        args = parser.parse_args()
        Resolve_Main(args.gene_tree, args.rooted_species_tree, om1.GetGeneToSpeciesMap(args), args.test)
