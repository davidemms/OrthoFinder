# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 11:00:11 2016

@author: david
"""
import numpy as np
from scipy.special import beta
 
"""
T - topology of the tree
B - set of branches in the tree
D - set of all duplications, {d_i | i \in B}
d_i - the duplications on branch i = (x_i, y_i)
O_i(j), where i, j \in B if the function O_i:B -> {->, r, <-} mapping a branch, j, to the orientation it would have
        if i where the root.
 
The probability of a given branch being the root (given the tree topology and the set of duplications) is
the product of each of the implied/required orientations of all the branches in the tree given only the duplications
on each branch individual and the topology of the tree, P_i(o | d_i, T). i ia the branch in question, o is its orientation, 
d_i are the numbers of duplications on the branch. 

Pclade = P(clade | D, T)
Pc = P_j(O_i(j) | d_j,T)

A.
==
       clade1
     /
root
     \
       clade2
       
p =       Pclade1 * Proot * Pclade2
     -----------------------------------
     \Pi_j (Pclade1 * Proot * Pclade2)_j


B.
==
          clade_i
        /
       i
      /   
clade
      \
       j
        \
          clade_j

Pclade = Pclade_i*Pc(O(i)) * Pclade_j*Pc(O(j))

C.
==

Pc = P_i(o | d_i, T) = P_i(d_i | o, T) P_i(o | T)
                       -------------------------- ,  X =  \sum_{o \in {<-, r, ->} P_i(d_i | o) P_i(o | T)
                                 X
(calculated by Ps_o_G_d, stored in Pc)                                 
                                 
                                 
where P_i(d_i | o, T) = P(d_i | o)  - independent of the branch, independent of the tree (only depends on orientation)
   (could potentially alter it so it is different for, say, terminal branches)

P_i(o | T) - is simply counted on the tree (what proportion of roots would give this orientation)     
      

================================ Branch Models ======================================

Poisson Model
-------------
           

"""


spTreeFormat = 1      

_lf = [0.000000000000000, 0.000000000000000, 0.693147180559945, 1.791759469228055, 3.178053830347946, 4.787491742782046, 6.579251212010101, 8.525161361065415, 10.604602902745251, 12.801827480081469, 15.104412573075516, 17.502307845873887, 19.987214495661885, 22.552163853123421, 25.191221182738683, 27.899271383840894, 30.671860106080675, 33.505073450136891, 36.395445208033053, 39.339884187199495, 42.335616460753485, 45.380138898476908, 48.471181351835227, 51.606675567764377, 54.784729398112319, 58.003605222980518, 61.261701761002001, 64.557538627006323, 67.889743137181526, 71.257038967168000, 74.658236348830158, 78.092223553315307, 81.557959456115029, 85.054467017581516, 88.580827542197682, 92.136175603687079, 95.719694542143202, 99.330612454787428, 102.968198614513810, 106.631760260643450, 110.320639714757390, 114.034211781461690, 117.771881399745060, 121.533081515438640, 125.317271149356880, 129.123933639127240, 132.952575035616290, 136.802722637326350, 140.673923648234250, 144.565743946344900, 148.477766951773020, 152.409592584497350, 156.360836303078800, 160.331128216630930, 164.320112263195170, 168.327445448427650, 172.352797139162820, 176.395848406997370, 180.456291417543780, 184.533828861449510, 188.628173423671600, 192.739047287844900, 196.866181672889980, 201.009316399281570, 205.168199482641200, 209.342586752536820, 213.532241494563270, 217.736934113954250, 221.956441819130360, 226.190548323727570, 230.439043565776930, 234.701723442818260, 238.978389561834350, 243.268849002982730, 247.572914096186910, 251.890402209723190, 256.221135550009480, 260.564940971863220, 264.921649798552780, 269.291097651019810, 273.673124285693690, 278.067573440366120, 282.474292687630400, 286.893133295426990, 291.323950094270290, 295.766601350760600, 300.220948647014100, 304.686856765668720, 309.164193580146900, 313.652829949878990, 318.152639620209300, 322.663499126726210, 327.185287703775200, 331.717887196928470, 336.261181979198450, 340.815058870798960, 345.379407062266860, 349.954118040770250, 354.539085519440790, 359.134205369575340, 363.739375555563470, 368.354496072404690, 372.979468885689020, 377.614197873918670, 382.258588773060010, 386.912549123217560, 391.575988217329610, 396.248817051791490, 400.930948278915760, 405.622296161144900, 410.322776526937280, 415.032306728249580, 419.750805599544780, 424.478193418257090, 429.214391866651570, 433.959323995014870, 438.712914186121170, 443.475088120918940, 448.245772745384610, 453.024896238496130, 457.812387981278110, 462.608178526874890, 467.412199571608080, 472.224383926980520, 477.044665492585580, 481.872979229887900, 486.709261136839360, 491.553448223298010, 496.405478487217580, 501.265290891579240, 506.132825342034830, 511.008022665236070, 515.890824587822520, 520.781173716044240, 525.679013515995050, 530.584288294433580, 535.496943180169520, 540.416924105997740, 545.344177791154950, 550.278651724285620, 555.220294146894960, 560.169054037273100, 565.124881094874350, 570.087725725134190, 575.057539024710200, 580.034272767130800, 585.017879388839220, 590.008311975617860, 595.005524249382010, 600.009470555327430, 605.020105849423770, 610.037385686238740, 615.061266207084940, 620.091704128477430, 625.128656730891070, 630.172081847810200, 635.221937855059760, 640.278183660408100, 645.340778693435030, 650.409682895655240, 655.484856710889060, 660.566261075873510, 665.653857411105950, 670.747607611912710, 675.847474039736880, 680.953419513637530, 686.065407301994010, 691.183401114410800, 696.307365093814040, 701.437263808737160, 706.573062245787470, 711.714725802289990, 716.862220279103440, 722.015511873601330, 727.174567172815840, 732.339353146739310, 737.509837141777440, 742.685986874351220, 747.867770424643370, 753.055156230484160, 758.248113081374300, 763.446610112640200, 768.650616799717000, 773.860102952558460, 779.075038710167410, 784.295394535245690, 789.521141208958970, 794.752249825813460, 799.988691788643450, 805.230438803703120, 810.477462875863580, 815.729736303910160, 820.987231675937890, 826.249921864842800, 831.517780023906310, 836.790779582469900, 842.068894241700490, 847.352097970438420, 852.640365001133090, 857.933669825857460, 863.231987192405430, 868.535292100464630, 873.843559797865740, 879.156765776907600, 884.474885770751830, 889.797895749890240, 895.125771918679900, 900.458490711945270, 905.796028791646340, 911.138363043611210, 916.485470574328820, 921.837328707804890, 927.193914982476710, 932.555207148186240, 937.921183163208070, 943.291821191335660, 948.667099599019820, 954.046996952560450, 959.431492015349480, 964.820563745165940, 970.214191291518320, 975.612353993036210, 981.015031374908400, 986.422203146368590, 991.833849198223450, 997.249949600427840, 1002.670484599700300, 1008.095434617181700, 1013.524780246136200, 1018.958502249690200, 1024.396581558613400, 1029.838999269135500, 1035.285736640801600, 1040.736775094367400, 1046.192096209724900, 1051.651681723869200, 1057.115513528895000, 1062.583573670030100, 1068.055844343701400, 1073.532307895632800, 1079.012946818975000, 1084.497743752465600, 1089.986681478622400, 1095.479742921962700, 1100.976911147256000, 1106.478169357800900, 1111.983500893733000, 1117.492889230361000, 1123.006317976526100, 1128.523770872990800, 1134.045231790853000, 1139.570684729984800, 1145.100113817496100, 1150.633503306223700, 1156.170837573242400]

def LogFactorial(n):
    if n < 0: raise NotImplementedError()
    elif n > 254:
        x = float(n + 1)
        return (x - 0.5)*np.log(x) - x + 0.5*np.log(2*np.pi) + 1.0/(12.0*x);
    else:
        return _lf[n]    
        
# With thanks to "Numerically Stable Hidden Markov Model Implementation", Tobias P. Mann 2006 for extended logarithm approach
# define None := log(0) 

def eexp(x):
    if x == None:
        return 0.
    else:
        return np.exp(x)

def eln(x):
    if x == 0.:
        return None
    else:
        return np.log(x)

def elnsum(elnx, elny):
    if (elnx == None) or (elny == None):
        return elny if elnx == None else elnx
    elif elnx > elny:
        return elnx + eln(1 + np.exp(elny-elnx))  
    else:
        return elny + eln(1 + np.exp(elnx-elny))
        
def elnprod(elnx, elny):
    if elnx == None or elny == None: 
        return None
    else:
        return elnx + elny
        
def elndiv(elnx, elny):
    if elnx == None:
        return None
    elif elny == None:
        print((elnx, elny))
        raise NotImplemented()
    else:
        return elnx - elny

def GetSpeciesName(name):
    if name.count("_") == 2:
        names = name.split("_")[1:]
    elif name.count("_") == 1:
        names = name.split("_")
    elif name.count("_") == 0:
        names = [name[0], name[1:]]
    else:
        names = (name[0], name[1:].lower())
    return names[0][0].upper() + ". " + names[1] 

def lnpoisson(k, r):
    """ k - events, r - rate"""
    return  elndiv(k * eln(r) - r, LogFactorial(k))

def get_bipartitions(tree):
     """ 
     It returns the set of all possible partitions under a
     node. Note that current implementation is quite inefficient
     when used in very large trees.

     t = Tree("((a, b), e);")
     partitions = t.get_partitions()

     # Will return: 
     # (a,b,e), ()
     # (a,e), (b)
     # etc
     """
     all_leaves = frozenset(tree.get_leaf_names())
     all_partitions = set()
     all_partitions.add(frozenset([all_leaves, frozenset()]))
     for n in tree.iter_descendants():
        p1 = frozenset(n.get_leaf_names())
        p2 = frozenset(all_leaves - p1)
        all_partitions.add(frozenset([p1, p2]))
     return all_partitions

# ===============================================================================================================================         

class BranchProbModel_corrected(object):
    """
    Probability of ->, <-, r shouldn't depend on tree 
    """
    @staticmethod
    def P_o_G_T(toA, A, B):
        """
        P(orientation | Tree)
        
        Inputs:
        toA - True: o=ToA, False: o=FromA, None: o=root
        A - cladeA
        B - cladeB
        P(o|T)
        
        Note, will be 0 for terminal branch and time moving away from it
        """
        nA = 2.*len(A) - 2   # number of branches in A part of tree
        nB = 2.*len(B) - 2   # number of branches in B part of tree
        n = nA + nB + 1     # number of branches in tree
        if toA == None:
            return 1./n
        elif (toA and len(B) == 1) or (not toA and len(A) == 1):
            return 0.
        else:
            return 0.5*(n-1)/n
            
# ===============================================================================================================================

class PoissonModel_IntergrateBranchLenthsSumFP(BranchProbModel_corrected):
    """Integrate over the possible position of the root along the length of the root branch and sum over the distributions of
    false positives
    """
    def __init__(self, rel_FP_rate_nonterm=0.001, rel_FP_rate_term=0.001, qExcludeTerminals=False):
        """
        rel_FP_rate - the rate of false positives relative to the rate of true positives, \lambda_{FP} = a \lambda_{TP}
        """
        self.a = rel_FP_rate_nonterm    
        self.a_term = rel_FP_rate_term    
        self.qExcludeTerminals = qExcludeTerminals        
        
    def P_d_G_o(self, m, n, toM, qTerminal=False):
        return eexp(self.lnP_d_G_o(m, n, toM, qTerminal))
           
    def lnP_d_G_o(self, m, n, toM, qTerminal=False):
        """Pnt(d|o)
        m - duplications to M
        n - duplications away from N
        toM - True: o=ToM, False: o=FromM, None: o=root
        """
        N = m+n             # total number of events
        if N == 0:
            return 1.
        alpha = self.a_term if qTerminal else self.a
        ltp = N/(1+alpha)  # lambda_TP
        lfp = alpha * ltp  # lambda_FP
        lltp = eln(ltp)
        la = eln(alpha)
        if toM == None:
            ln_pTot = eln(0.)
            for s in range(m+1):
                for t in range(n+1):
                    X = m - s + t
                    Y = n - t + s
                    integral = beta(X+1, Y+1)
                    ln_pTot = elnsum(ln_pTot, elnprod(eln(integral), (m+n)*lltp - ltp - lfp + (s+t)*la - LogFactorial(m-s) - LogFactorial(s) - LogFactorial(n-t) - LogFactorial(t)))
            return ln_pTot
        elif toM:
            return elnprod(lnpoisson(m ,ltp), lnpoisson(n, lfp))
        else:
            return elnprod(lnpoisson(m ,lfp), lnpoisson(n, ltp))
    
    def Ps_o_G_d(self, A, B, countsA, countsB):
        """ 
        *** Resolve Numerical Overflows: 3 Use logarithms *** 
        See Base class decription for remaining method documentation
        """
        if self.qExcludeTerminals:
            if len(A) == 1:
                countsA = 0
            if len(B) == 1:
                countsB = 0
        isTerminal = len(A) == 1 or len(B) == 1
        lnx = elnprod(self.lnP_d_G_o(countsA, countsB, True, isTerminal), eln(self.P_o_G_T(True, A, B)))  # P(d|->A, T)
        lny = elnprod(self.lnP_d_G_o(countsA, countsB, False, isTerminal), eln(self.P_o_G_T(False, A, B)))  # P(d|<-A, T) 
        lnz = elnprod(self.lnP_d_G_o(countsA, countsB, None, isTerminal), eln(self.P_o_G_T(None, A, B)))  # P(d|r, T)
        lntot = elnsum(elnsum(lnx, lnz), lny)  # both x or y could be small but not x & z (or y & z)
        x = eexp(elndiv(lnx, lntot))
        y = eexp(elndiv(lny, lntot))
        z = eexp(elndiv(lnz, lntot))
        return (x, y, z)  
        
         
# ===============================================================================================================================

class PoissonModel_WithTeminalModel(BranchProbModel_corrected):
    """Integrate over the possible position of the root along the length of the root branch and sum over the distributions of
    false positives
    """
    def __init__(self, rel_FP_rate_nonterm, TP_rate_term, FP_rate_term):
        """
        rel_FP_rate - the rate of false positives relative to the rate of true positives, \lambda_{FP} = a \lambda_{TP}
        """
        self.a = rel_FP_rate_nonterm    
        self.TP_rate_term = TP_rate_term         
        self.FP_rate_term = FP_rate_term         
        
    def P_d_G_o(self, m, n, toM, qTerminal=False):
        if qTerminal:
            return eexp(self.lnP_d_G_o_term(m, n))
        else:
            return eexp(self.lnP_d_G_o_nonterm(m, n))
            
    def lnP_d_G_o_term(self, m_inward, qRoot):
        """
        if root: duplications away from species arrive at absolute rate self.TP_rate_term 
        or
        if not root: duplications away from species arrive at absolute rate self.FP_rate_term 
        """
        if qRoot:
            return lnpoisson(m_inward, self.TP_rate_term)
        else:
            return lnpoisson(m_inward, self.FP_rate_term)
            
    def lnP_d_G_o(self, m, n, toM):
        """Pnt(d|o)
        m - duplications to M
        n - duplications away from N
        toM - True: o=ToM, False: o=FromM, None: o=root
        """
        N = m+n             # total number of events
        if N == 0:
            return 1.
        alpha = self.a
        ltp = N/(1+alpha)  # lambda_TP
        lfp = alpha * ltp  # lambda_FP
        lltp = eln(ltp)
        la = eln(alpha)
        if toM == None:
            ln_pTot = eln(0.)
            for s in range(m+1):
                for t in range(n+1):
                    X = m - s + t
                    Y = n - t + s
                    integral = beta(X+1, Y+1)
                    ln_pTot = elnsum(ln_pTot, elnprod(eln(integral), (m+n)*lltp - ltp - lfp + (s+t)*la - LogFactorial(m-s) - LogFactorial(s) - LogFactorial(n-t) - LogFactorial(t)))
            return ln_pTot
        elif toM:
            return elnprod(lnpoisson(m ,ltp), lnpoisson(n, lfp))
        else:
            return elnprod(lnpoisson(m ,lfp), lnpoisson(n, ltp))
    
    def Ps_o_G_d(self, A, B, countsA, countsB):
        if len(A) == 1 or len(B) == 1:
            # Terminal duplications
            if len(A) == 1:
                lnx = elnprod(self.lnP_d_G_o_term(countsB, False), eln(self.P_o_G_T(True, A, B)))
                lny = None
                lnz = elnprod(self.lnP_d_G_o_term(countsB, True), eln(self.P_o_G_T(None, A, B)))
            else:
                lnx = None
                lny = elnprod(self.lnP_d_G_o_term(countsA, False), eln(self.P_o_G_T(False, A, B)))
                lnz = elnprod(self.lnP_d_G_o_term(countsA, True), eln(self.P_o_G_T(None, A, B)))
        else:
            lnx = elnprod(self.lnP_d_G_o(countsA, countsB, True), eln(self.P_o_G_T(True, A, B)))  # P(d|->A, T)
            lny = elnprod(self.lnP_d_G_o(countsA, countsB, False), eln(self.P_o_G_T(False, A, B)))  # P(d|<-A, T) 
            lnz = elnprod(self.lnP_d_G_o(countsA, countsB, None), eln(self.P_o_G_T(None, A, B)))  # P(d|r, T)
        lntot = elnsum(elnsum(lnx, lnz), lny)  # both x or y could be small but not x & z (or y & z)
        x = eexp(elndiv(lnx, lntot))
        y = eexp(elndiv(lny, lntot))
        z = eexp(elndiv(lnz, lntot))
        return (x, y, z) 
           
# ===============================================================================================================================
# Tree-level model/calculation

def Pc(node, p_cond):
    x = frozenset(node.get_leaf_names())
#    print((p_cond[x][0], x, 'edge'))
    return p_cond[x][0]

def Pc_root(node, p_cond):
    return p_cond[frozenset(node.get_leaf_names())][2]
    
def P_clade(node, probs, p_cond):
    cluster = frozenset(node.get_leaf_names())
    if len(cluster) == 1: return 1.
    if cluster in probs:
        return probs[cluster] 
    else:
        l, r = node.children
        p = P_clade(l, probs, p_cond)*Pc(l, p_cond) * P_clade(r, probs, p_cond)*Pc(r, p_cond)
        probs[cluster] = p
#        print((p, cluster)) 
        return p
 
# ===============================================================================================================================


def GetBranchProbs(supported_clusters_counter, biparts, probModel): 
    p_cond = dict()
    for A, B in biparts:
        if len(A) == 0 or len(B) == 0: continue
        x,y,z = probModel.Ps_o_G_d(A, B, supported_clusters_counter[A], supported_clusters_counter[B])
        p_cond[A] = (x,y,z)
        p_cond[B] = (y,x,z)
    return p_cond
    
def GetFinalProbs(biparts, p_cond, tree):
    probs_clades = dict()
    p_final = dict() # unnormalised until the end
    total = 0.
    for A, B in biparts:
        if len(A) == 0 or len(B) == 0: continue
#        print(A)
        # tree for calculation
        if len(A) == 1:
            n2 = tree & list(A)[0]
        else:
            n2 = tree.get_common_ancestor(A)
        # cluster A could straddle root
        if n2 == tree:
            if len(B) == 1:
                n2 = tree & list(B)[0]
            else:
                n2 = tree.get_common_ancestor(B)
        if n2 != tree: tree.set_outgroup(n2)
        l, r = tree.children
        p = P_clade(l, probs_clades, p_cond) * Pc_root(l, p_cond) * P_clade(r, probs_clades, p_cond)
        total += p
        p_final[A] = p
#        print(p, A, 'final')
    p_final = {k:(v/total) for k, v in p_final.items()}
    return p_final
    
def GetTerminalRates(allSpecies, clades, supported_clusters_counter):
    nAntiTerm = len(allSpecies) - 1
    termDups = sorted([v for k,v in supported_clusters_counter.items() if len(k) == nAntiTerm], reverse=True)
    nSp = float(len(allSpecies))
    ntp = sum(termDups[:2]) / 2. 
    nfp = sum(termDups)
    rtp = (ntp + 1)
    rfp = (nfp + 1)/nSp
#    print("Terminal inward rate = %f" % rtp)
#    print("Terminal inward false-positive rate = %f" % rfp)
    return rtp, rfp, 
    
def GetAlpha(allSpecies, clades, supported_clusters_counter):
    nSp = len(allSpecies)
    nTerminal = sum([v for k,v in supported_clusters_counter.items() if (len(k) == 1 or len(k) == nSp-1)])
    nNonTerminal = sum([v for k,v in supported_clusters_counter.items() if not (len(k) == 1 or len(k) == nSp-1)])
    contradictions_term = dict()
    contradictions_nonterm = dict()
    for clade in clades:
        clade_p = allSpecies.difference(clade)
        against_term = 0
        against_nonterm = 0
        for observed, n in supported_clusters_counter.items():
            if (not observed.issubset(clade)) and (not observed.issubset(clade_p)):
                if len(observed) == 1 or len(allSpecies.difference(observed)) == 1:
                    against_term += n
                else:
                    against_nonterm += n
        contradictions_term[clade] = against_term
        contradictions_nonterm[clade] = against_nonterm
    m = min([contradictions_nonterm[k] + contradictions_term[k] for k in clades])
    n = sum(supported_clusters_counter.values())
#    nSupport = n-m
    nFP_term = []
    nFP_nonterm = []
    for clade in clades:
        if contradictions_nonterm[clade] + contradictions_term[clade] == m:
            nFP_nonterm.append(contradictions_nonterm[clade])
            nFP_term.append(contradictions_term[clade])
    nFP_mean = np.mean(nFP_nonterm)
    alpha_nonTerm = float(nFP_mean + 1.) / float(nNonTerminal - nFP_mean + 1.) 
    nFP_mean = np.mean(nFP_term)
    alpha_term = float(np.mean(nFP_term) + 1.) / float(nTerminal - nFP_mean) 
#    print("Observed %d duplications. %d support the best root(s) and %d contradict them." % (nTerminal + nNonTerminal, nSupport, nTerminal + nNonTerminal - nSupport))
#    nNonTrivial = sum([v for k,v in supported_clusters_counter.items() if len(k) != 1])
#    print("%d non-trivial duplications" % nNonTrivial)
#    print("alpha (non-terminal) = %f" % alpha_nonTerm)
#    print("alpha (terminal) = %f" % alpha_term)
    return alpha_nonTerm, alpha_term


def GetProbabilities(species_tree, allSpecies, clades, supported_clusters_counter):
    biparts = get_bipartitions(species_tree)
    alpha_nonTerm, alpha_term = GetAlpha(allSpecies, clades, supported_clusters_counter)
    tp_term, fp_term = GetTerminalRates(allSpecies, clades, supported_clusters_counter)
    probModel = PoissonModel_WithTeminalModel(alpha_nonTerm/10, tp_term, fp_term)
#    probModel = PoissonModel_IntergrateBranchLenthsSumFP(alpha_nonTerm/10., alpha_term/10., qExcludeTerminals = True)
    p_cond = GetBranchProbs(supported_clusters_counter, biparts, probModel)
    return GetFinalProbs(biparts, p_cond, species_tree)    
           
