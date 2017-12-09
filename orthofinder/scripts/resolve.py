# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 14:15:17 2017

@author: david
"""  
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
    Implementation:
        n and s could be anywhere in relation to one another, need to change all the species sets that could be affected
    """
#    c1, c2 = top.get_children()
#    print((len(c1), len(c2)))
#    print(top)
#    print(n)
#    print(s)
    n, nUp, top = DetachAndCleanup(top, n)
    parent = s.up
    try:
        new = parent.add_child()
    except:
        print(top.get_tree_root())
        print(top)
        print(n)
        print(s)
        raise
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
    return True
    
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
    return True
        
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

depth = 2 # check down to A, B, C & D. Can't just change this without checking rest of code

def resolve(n, M):
    """
    Rearrange node in more parsimonious configuration. 
    Cost of Moves+Dupliciations+Losses after is guaranteed to be lower than cost of Dups+Losses before.
    M=1, D=1, L=1
    Args:
        n - tree_lib node with overlapping species sets on it's children
        M - function taking a gene name and returning the species name
        
    Implementation:
        X, Y - the sets of species either side of the node of interest
        S - the set of species sister to this node
        (A, B) = X 
        (C, D) = Y
        Lowercase letters are the corresponding nodes
        #bipart - cases which will need to be generalised to deal with non-binary nodes
        
    Overall plan:
        First deal with all cases where the tree is binary
            - Must also set the spsecies sets correctly (this must be tested in the unit tests)
        Then deal with non-binary cases
    """
    if "recon" in n.features:
        sis = n.recon
        n.del_feature('recon')
        if n.check_monophyly(target_attr='name', values=sis[0])[0] and n.check_monophyly(target_attr='name', values=sis[1])[0]:
            s0 = n.get_common_ancestor(sis[0])
            s1 = n.get_common_ancestor(sis[1])
            return GraftAndUpdate(n, s0, s1)
    ch = n.get_children()
    if len(ch) != 2: return     #bipart
    X,Y = [c.sp_down for c in ch]
    if not (X & Y): return # no overlap, nothing to do
    if X == Y: return   # identical species sets, they're paralogues
    sis = n.get_sisters()
    s = sis[0] if len(sis) == 1 else None
    S = None if n.is_root() else set.union(*[nother.sp_down for nother in n.up.get_children() if nother != n]) # works for non-binary too
    
    if len(X) == 1 or len(Y) == 1:
        """
        Part 1
        
           /-S
        --|
          |   /-Y
           \-|
             |   /-B
              \-|
                 \-A
                 
        i.e. Y is a single species
        Question, does a moving y decrease the overall reconciliation cost?
        Know:
            Y is a single species 
            |X| > 1 (since o/w X=Y given there is an overlap and |Y| = 1. Have returned if X=Y) 
        """
        if len(Y) == 1:
            # then know len(X) != 1
            x, y = ch
        elif len(X) == 1:  
            # then know len(Y) != 1
            TEMP = X
            X = Y
            Y = TEMP
            y, x = ch
        TEMP = x.get_children()
        if len(TEMP) > 2: return False #bipart
        A,B = [t.sp_down for t in TEMP]
        a,b = TEMP
        
        # Cases
        # 1. WLOG Y \in A & Y \notin B
        # If instead Y \in A and Y \in B - can't do anything, would have been resolved by now if it could have been
        inA = bool(Y & A)
        inB = bool(Y & B)
        if inA + inB == 1:
            # Case 1
            if inB:
                TEMP = A
                A = B
                B = TEMP
                TEMP = a
                a = b
                b = TEMP
            # Case 1a. Y = A
            # Case 1b. Y < A
            if len(A) == 1:
                # Case 1a, Y == A
                # Possibilities: 
                #   i.  B == S. Should move y up
                #   ii. B != S. y becomes sister of a
                if S!=None and B & S:
                    # The move should be done later, remove the exception once that's been done
                    # Mark node for later rearrangement
                    if s!= None: n.up.add_feature("recon", (y.get_leaf_names(),s.get_leaf_names()))    # bipart: if polytomy for S species, they should all be put together
                    return False
                else:
                    return GraftAndUpdate(n, y, a) 
            else:
                # Case 1b
                O = Y # Y is a single species
                successX, nX, dX = ContainsMonophyletic(x, O, 0, depth + 1) # can go one deeper if len(Y) == 1 since requires one fewer move on that half
                if not successX: return False
#                print("Do I need to consider doing something if ¬successX?")
                assert(dX != 0) # Can't have dX == 0
                if dX == 1:
                    return GraftAndUpdate(n, nX, y) 
                elif dX == 2:
                    # If there the sister and sister once removed overlap then it's a duplication before (sister, Y) 
#                    outerN = next(nn for nn in x.get_children() if not (nn.sp_down & Y))
                     # This code may not be used any more
#                    print(y.get_leaf_names())
#                    print(nX.get_leaf_names())
#                    print(outerN.get_leaf_names())
#                    GraftAndUpdate(n, y, outerN) 
                    return GraftAndUpdate(n, y, nX) 
                else:
                    raise NotImplemented # shouldn't get here
                    
                    
#        else:
#            # Case 2
#            # Y must only be on one of the children of A, otherwise too many moves
#            inCh = [Y & ch_.sp_down for ch_ in a.get_children()]  #bipart
#            if sum(inCh) > 1: return False
#            # check the node is entirely Y
#            nThis = a.get_children()[0 if inCh[0] else 1] # bipart
#            if len(nThis.sp_down) != 1: return False
#            # Doesn't matter at this point if we bring nThis up or move y down, same orthologues will be recovered
#            return GraftAndUpdate(n, y, nThis)
    else:
        """
        Part 2
        
           /-S
          |
        --|      /-C
          |   /-|
          |  |   \-D
           \-|
             |   /-A
              \-|
                 \-B
        """
        # Need to deal with the case that there is a monophyletic clade of species in both X and Y and this is the only overlap
        # Otherwise there is nothing to do (current plan)
        #       What are the more complicated rearrangements that could be done. Would any of them lower the recon cost?
        #       if there were duplicated species across A & B, they would have already been resolved at the lower node if it were possible
        # Implementation: Check if the overlapping species are monophyletic in X and Y, if so then rearrange if it lowers recon 
        # cost (i.e. the clades are not too deep)
        
        O = X.intersection(Y)
        x, y = ch
        successX, nX, dX = ContainsMonophyletic(x, O, 0, depth)
        successY, nY, dY = ContainsMonophyletic(y, O, 0, depth)
        # regardless if we get to the monophyletic clade, it still lowers the recon cost to put them together
        if depth != 2:
            raise NotImplemented
        if dX == 0 and dY == 0:
            return False
        elif dX == 1 and dY == 1:
            # transform to tripartition ((nX,nY), XOther, yOther) if dX==1 and dY == 1
            return GraftTripartAndUpdate(n, nX, nY, n) 
        elif dX==0:
            sis = nY.get_sisters()
            if len(sis) != 1: return False # bipart
            if S != None and sis[0].sp_down & S: 
                # mark nX and s to be put together as sister species later in the tree traversal
                if s!= None: n.up.add_feature("recon", (nX.get_leaf_names(),s.get_leaf_names()))    # bipart: if polytomy for S species, they should all be put together
                return False # different rearrangement with next node up is better
            return GraftAndUpdate(n, nX, nY)
        elif dY==0:
            sis = nX.get_sisters()
            if len(sis) != 1: return False # bipart
            if S != None and sis[0].sp_down & S: 
                # mark nY and s to be put together as sister species later in the tree traversal
                if s!= None: n.up.add_feature("recon", (nY.get_leaf_names(),s.get_leaf_names()))    # bipart: if polytomy for S species, they should all be put together
                return False # different rearrangement with next node up is better
            return GraftAndUpdate(n, nY, nX)
        else:
            raise NotImplemented

def Resolve_Main(trees_fn, species_tree_rooted_fn, GeneToSpecies):
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    try:
        tree = tree_lib.Tree(trees_fn)
    except:
        tree = tree_lib.Tree(trees_fn, format=3)
    tree.prune(tree.get_leaf_names())
    if len(tree) == 1: return tree
#        print(tree)
    roots, j, GeneMap = om1.GetRoots(tree, species_tree_rooted, GeneToSpecies)
    if len(roots) == 0: return tree
    # Pick the first root for now
    root = roots[0]
    if root != tree:
        tree.set_outgroup(root)
    om1.StoreSpeciesSets(tree, GeneToSpecies)
    for n in tree.traverse("postorder"):
        resolve(n, GeneToSpecies)
    return tree
           
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gene_tree")
    parser.add_argument("rooted_species_tree")
    parser.add_argument("-s", "--separator", choices=("dot", "dash", "second_dash", "3rd_dash", "hyphen"), help="Separator been species name and gene name in gene tree taxa")
    args = parser.parse_args()
    tree = Resolve_Main(args.gene_tree, args.rooted_species_tree, om1.GetGeneToSpeciesMap(args))
    tree.write(outfile=args.gene_tree + ".rec.tre")
                

