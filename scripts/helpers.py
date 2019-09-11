import networkx as nx

def homopoly(kmer):
    """Check if kmer is a homopolymer"""
    last=None
    rv=1
    h=1
    for j in range(len(kmer)):
        if kmer[j] == last:
            h +=1
        else:
            rv=max(rv, h)
            h=1
            last = kmer[j]
    return(max(rv, h))

def findTiling(hsp):
    """Find a maximum scoring subset of local hits
    
    This is done by generating a directed graph. Every match is represented by a
    Start (S...) and End (E...) node which are connected by an edge with the 
    score as weight in both directions. Compatible (in the right direction, and 
    non-overlapping) edges are connected with zero weight edges. A maximum score
    path trough the graph than represents the set of maximum scoring 
    non-overlapping matches.
    """
    G=nx.DiGraph()
    nodes = []
    b_edges = []
    for h, tHsp in enumerate(hsp):
        nodes.extend(["S%i" % h, "E%i" % h])
        if tHsp[0] < tHsp[1]:
            #forward
            b_edges.append(("S%i" % h, "E%i" % h, -tHsp[2])) #use negative bit score as edge weight to be able to use minimum path algorithm to find maximum path
        else:
            #reverse
            b_edges.append(("E%i" % h, "S%i" % h, -tHsp[2]))
    G.add_nodes_from(nodes)
    G.add_weighted_edges_from(b_edges)
    #create edges between all non-overlapping matches that are in the same direction
    c_edges = []
    for i in range(len(hsp)):
        for j in range(len(hsp)):
            if i==j:
                continue
            if (hsp[i][0] < hsp[i][1]) and (hsp[j][0] < hsp[j][1]) and (hsp[i][1] < hsp[j][0]):
                #forward compatibility edge
                c_edges.append(("E%i" % i, "S%i" % j, 0))
            elif (hsp[i][0] > hsp[i][1]) and (hsp[j][0] > hsp[j][1]) and (hsp[i][1] > hsp[j][0]):
                #reverse compatibility edge
                c_edges.append(("S%i" % j, "E%i" % i, 0))
    G.add_weighted_edges_from(c_edges)

    #find sources and sinks and add START and END node
    s_edges = []
    for node in G.nodes():
        if G.in_degree(node) == 0:
            s_edges.append(("START", node, 0))
        if G.out_degree(node) == 0:
            s_edges.append((node, "END", 0))
    G.add_weighted_edges_from(s_edges)
        
    #find shortest (ie. maximum score) path from start to end
    pre, dist = nx.floyd_warshall_predecessor_and_distance(G)
    if dist["END"]["START"] < dist["START"]["END"]:
        start = "END"
        end = "START"
    else:
        start = "START"
        end = "END"
    #"record" shortest reverse path by
    # traversing backword from end to start along the shortest path
    # by going always to the predecessor in the shortet path
    current = end
    path = []
    while pre["START"][current] != start:
        path.append(pre["START"][current])
        current = pre["START"][current]
    
    #return the reverse of the path, exluding the artificial start and end nodes
    return [int(p[1:]) for p in path[::-2]]

def lca(lineageStrings, stringency=1.0, 
        unidentified=["unidentified", "unclassified", "unknown"],
        ignoreIncertaeSedis=True, sizes=None):
    """Find lowest common ancestor
    
    Input is a set of strings representing classifications as path from the root
    node to the most specific available classification. Stringency define what 
    proportion of classifications need to concur at a certain level to be 
    accepted. unidetnfied gives names that should be ignored. ignoreIncertaSedis
    flag can be set to not accept incerta sedis classification. sizes can give
    wights to each classification.
    """
    lineage = []
    mLineages = []
    #remove bootstrap values ("(100)", "(75)", etc.) if any
    for mLin in [l.strip(";").split(";") for l in lineageStrings]:
        mLineages.append([])
        for entry in mLin:
             mLineages[-1].append(entry.split("(")[0])
    maxLinLen = max([len(m) for m in mLineages])
    active = [True]*len(mLineages)
    for i in range(maxLinLen):
        total = 0.0
        counts = {}
        for m, memberLin in enumerate(mLineages):
            if not active[m]:
                continue #ignore lineages that were deactivated further up in the tree
            if len(memberLin) <= i:
                active[m] = False
                continue #ignore lineages that are not this long
            name = memberLin[i].split("__")[-1]
            if name in unidentified:
                continue # ignoring unidentified entrys
            if ignoreIncertaeSedis and name.startswith("Incertae"):
                continue # ignoring Incertae sedis entries.
                         # NOTE: this will mean lineages end at the first Incerta sedis
            if sizes is None:
                tSize = 1
            else:
                tSize = sizes[m]
            total += tSize
            try:
                counts[memberLin[i]] += tSize
            except KeyError:
                counts[memberLin[i]] = tSize
        if not counts:
            #no valid lineage entrys found in this level
            break
        most=sorted(counts.items(), key=lambda x: x[1], reverse=True)[0]
        different = total - most[1]
        #accept the lineage entry if its proportion of all (valid) classifications higher than stringency setting
        if different/total <= (1.0-stringency):
            lineage.append(most[0]) #add the most apearing entry to the new lineage
            #deactivate all lineages that were different at this level
            for m, memberLin in enumerate(mLineages):
                if active[m] and memberLin[i] != most[0]:
                    active[m] = False
        else:
            break
    if len(lineage) == 0:
        lineage = ["unknown"]
    return ";".join(lineage)

