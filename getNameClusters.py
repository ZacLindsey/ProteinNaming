# -*- coding: utf-8 -*-
"""
Created on Thu Jun 02 08:39:13 2016

@author: zacli

Methods to look at a list of proteins and find other proteins which are named
differently but might be the same, using usearch clustering.

Information about the format of the .tab file can be found at:
http://www.drive5.com/usearch/manual/opt_uc.html
"""
#==============================================================================
# R VERSION (1.0)
#   > Re-designed for simplicity (6/14 - 6/16)
#   > Intended to generate output for R
#==============================================================================
def gather_clusters_R_ZL(proteinFile = "Backbones_4.csv", 
                      clusterFile = "Plasmids20-200kb-6-9-2016_Clusters.tab",
                      outfile = "ClusterGroups.csv",
                      checkProteins = False):
    import csv

    cRows = [] # Rows which contain a cluster summary
    sRows = [] # Rows which contain a centroid or hit 
    
    with open(clusterFile,"r") as clusterfile:
        readTab = csv.reader(clusterfile, delimiter="\t")
        for clustRow in readTab:
            if clustRow[0] == "C":
                cRows.append(clustRow)
            else:
                sRows.append(clustRow)
    
    titles = []
    Names = []
    ClusterIDs = []
    ProtInfo = {}
    index = -1
    
    with open(proteinFile,"r") as backbones:
        backboneFile = csv.reader(backbones)
        for backboneRow in backboneFile:
            titles.append(backboneRow[0])
            Names.append({})
            ClusterIDs.append(set())
            index += 1
            for protein in backboneRow: 
                for srow in sRows:
                    currProt = srow[8].lower()
                    cid = srow[1]
                    if protein.lower() in currProt:
                        if checkProteins==False:
                            ClusterIDs[index].add(cid)
                            ProtInfo[currProt] = [srow[2],srow[3]]                            
                    if cid in ClusterIDs[index]:
                        ProtInfo[currProt] = [srow[2],srow[3]]
                        try:
                            Names[index][cid].add(currProt)
                        except KeyError:
                            Names[index][cid] = set()
                            Names[index][cid].add(currProt) 
                        sRows.remove(srow)
                        
    with open(outfile, "w") as out:
        Writer = csv.writer(out)
        for i in range(len(Names)):
            for j in Names[i]:
                rows = []
                for k in Names[i][j]:
                    row = []
                    row.append(titles[i])
                    row.append(j)
                    row.append(k)
                    row.append(ProtInfo[k][0])
                    row.append(ProtInfo[k][1])
                    rows.append(row)
                rows.sort(key=lambda x: int(x[3]))
                for row in reversed(rows):
                    Writer.writerow(row)
                     
#==============================================================================
#   Method to gather groups of clusters (June 2016) v5.0 (6/14)
#   Inputs: 
#   -Backbones.csv (a manually created csv where each row is a diff.
#       protein, and each column contains alternate names for that protein)
#   -Plasmids2016Clusters2.tab: a .uc or .tab file in usearch cluster format
#
#   Outputs:
#   -MostCommonName_ClusterGroups.csv: a csv file containing the 5 most
#   common names for each protein as well as its count and seq length
#   -ClusterGroups.csv: a csv file containing columns for each cluster group.
#       The first row is the title of the group, the second is a count
#       of how many clusters have at least one of the backbone genes.
#   -byRow_ClusterGroups.csv: same as above, but with each row containing the
#       protein names for each group
#   -Summary.csv: contains counts for how many proteins are in each group,
#    how many clusters, etc.
#
#==============================================================================
def gather_clusters_ZL(proteinFile = "Backbones_3.csv", 
                      clusterFile = "Plasmids20-200kb-6-9-2016_Clusters.tab",
                      outfile = "ClusterGroups.csv",
                      checkCounts = True,
                      writeMostSimilar = True,
                      writeMainFile = True):
    import csv

    nameGroups = [] # One for each row in Backbones.csv
    clusterGroups = []
    clusterCount = [] # The number of unique clusters found for each row
    protCount = []
    cRows = [] # Rows which contain a cluster summary
    sRows = [] # Rows which contain a centroid or hit 
    titles = []
    protSeqLengthsGroup = []
    
    with open(clusterFile,"r") as clusterfile:
        readTab = csv.reader(clusterfile, delimiter="\t")
        for clustRow in readTab:
            if clustRow[0] == "C":
                cRows.append(clustRow)
            else:
                sRows.append(clustRow)
                  
    with open(proteinFile,"r") as backbones:
        read = csv.reader(backbones)
        for backboneRow in read:
            titles.append(backboneRow[0])
            protGroup = [] 
            protSeqLengths = {}
            clustNumbers = set() # Redundancy not necessary for ID numbers
            
            for protein in backboneRow: 
                if len(str(protein))>0: # Ignore blanks in the csv file
                
                    # If the protein is found in any cluster, record the
                    # cluster ID number.
                    for clustRow in sRows:
                        if protein.lower() in clustRow[8].lower():
                            clustNumbers.add(clustRow[1])
                            
                            if clustRow[8] not in protSeqLengths:
                                protSeqLengths[clustRow[8]] = clustRow[2]
                                    
                    # Go through the clusters again. Add any names which
                    # correspond to any recorded cluster.
                    for clustRow in sRows:
                        if clustRow[1] in clustNumbers:
                            protGroup.append(clustRow[8])
                            
                            if clustRow[8] not in protSeqLengths:
                                protSeqLengths[clustRow[8]] = clustRow[2]
                      
            nameGroups.append(list(protGroup))
            clusterGroups.append(list(clustNumbers))
            
            protSeqLengthsGroup.append(protSeqLengths)
            
            clusterCount.append(len(clustNumbers))             
            protCount.append(len(set(protGroup)))
            
    # Optional code to check stats on number of clusters, etc.
    if checkCounts:
        Counter = 0
        counters = []
        for cluster in clusterGroups:
            for row in cRows:
                if row[1] in cluster:
                    Counter += int(row[2])
            counters.append(Counter)
            Counter = 0
            
        with open("Summary_"+outfile,"w") as LengthFile:
            Writer = csv.writer(LengthFile, delimiter=",")
            Writer.writerow(["Protein Name:"]+titles)
            Writer.writerow(["num_proteins_countedFromClusterFile:"]+counters)
            Writer.writerow(["num_distinct_names:"]+protCount)
            Writer.writerow(["CountClustersPerGroup:"]+clusterCount)
    
    mostCommonProts = []
    mostCommonProteinCounts = [] # Number of occurrences of most common prot
    if writeMostSimilar:
        temp = nameGroups
        for nG in temp:
            mCP = []
            mCPc = []
            for i in range(5): # Get the five most common proteins
                try:
                    mCP.append(max(set(nG), key=nG.count))
                    mCPc.append(nG.count(mCP[i]))
                    nG[:] = [x for x in nG if x != mCP[i]]
                except ValueError:
                    mCP.append("")
                    mCPc.append(0)
            mostCommonProts.append(mCP)
            mostCommonProteinCounts.append(mCPc)
            
        with open("MostCommonNames_"+outfile,"w") as nameFile:
            Writer = csv.writer(nameFile, delimiter=",")
            titlesRow = []
            for t in titles: #This just makes blank spaces for other info
                titlesRow.append(t)
                titlesRow.append("Num_appearances_in_clusterfile")
                titlesRow.append("Length_of_sequence")
            Writer.writerow(titlesRow)
            for i in range(5):
                currentRow=[]
                for j in range(len(titles)):
                    currProt = mostCommonProts[j][i]
                    currentRow.append(currProt)
                    currentRow.append(mostCommonProteinCounts[j][i])
                    if currProt in protSeqLengthsGroup[j]:
                        seqLength = protSeqLengthsGroup[j][currProt]
                    else:
                        seqLength = "NA"
                    currentRow.append(seqLength)
                Writer.writerow(currentRow)
        
        # For use in R
        with open("R_names_"+outfile,"w") as nameFile:
            Writer = csv.writer(nameFile, delimiter=",")
            for i in range(len(titles)):
                for j in range(5):
                    prot = mostCommonProts[i][j]
                    if len(prot)>0:
                        wr = [titles[i],prot,protSeqLengthsGroup[i][prot]]
                    Writer.writerow(wr)
            
    # Write a file containing all of the proteins in each backbone group
    if writeMainFile:
        pslg = protSeqLengthsGroup #Shorthand name...
        # The following line simply casts the dict values to all ints
        pslg_2 = [dict([a, int(x)] for a, x in b.iteritems()) for b in pslg]
        newoutput = []
        
        with open("byRow_"+outfile,"w") as NewFile:
            Writer = csv.writer(NewFile,delimiter=",")
            for i in range(len(pslg_2)):
                pslg_s = sorted(pslg_2[i],key=pslg_2[i].get)
                op = [n+" (SeqLength: "+str(pslg_2[i][n])+")" for n in pslg_s]
                rOutput = list(reversed(op))
                Writer.writerow(rOutput)
                newoutput.append(rOutput)
                
        # NOTE: The following creates a transpose matrix.
        # I would recommend ignoring it and using Excel instead.
        Rows = []
        Rows.append(titles)
        Rows.append(clusterCount)
        for i in range(len(max(newoutput,key=len))):
            Rows.append([])
            for j in range(len(newoutput)):
                try:
                    if len(newoutput[j][i])>0:
                        Rows[i+2].append(newoutput[j][i])
                except IndexError:
                    Rows[i+2].append("")
                    
        with open(outfile,"w") as NewFile:
            Writer = csv.writer(NewFile,delimiter=",")
            for row in Rows:
                Writer.writerow(row)
                
    return nameGroups
    
#==============================================================================
    #OLD VERSION        
#==============================================================================
def gather_clusters_ZL_old(proteinFile = "Backbones_2.csv", 
                      clusterFile = "Plasmids5-9-2016Clusters.tab",
                      outfile = "ClusterGroups.csv",
                      checkCounts = True,
                      writeFile = True):
    import csv

    nameGroups = [] # One for each row in Backbones.csv
    clusterGroups = []
    clusterCount = [] # The number of unique clusters found for each row
    protCount = []
    cRows = [] # Rows which contain a cluster summary
    sRows = [] # Rows which contain a centroid or hit 
    titles = []
    
    with open(clusterFile,"r") as clusterfile:
        readTab = csv.reader(clusterfile, delimiter="\t")
        for clustRow in readTab:
            if clustRow[0] == "C":
                cRows.append(clustRow)
            else:
                sRows.append(clustRow)
                  
    with open(proteinFile,"r") as backbones:
        read = csv.reader(backbones)
        for backboneRow in read:
            titles.append(backboneRow[0])
            protGroup = set() # One new group for each row (because the diff.
            # names in each row are supposed to be the same protein)
            clustNumbers = set() # Sets used to avoid redundancy
            
            for protein in backboneRow: 
                if len(str(protein))>0: # Ignore blanks in the csv file
                
                    # If the protein is found in any cluster, record the
                    # cluster number.
                    for clustRow in sRows:
                        if protein.lower() in clustRow[8].lower():
                            clustNumbers.add(clustRow[1])
                                    
                    # Go through the clusters again. Add any names which
                    # correspond to any recorded cluster.
                    for clustRow in sRows:
                        if clustRow[1] in clustNumbers:
                            protGroup.add(clustRow[8].split('_:_')[0])
                                                   
            nameGroups.append(list(protGroup))
            clusterGroups.append(list(clustNumbers))
            
            clusterCount.append(len(clustNumbers))             
            protCount.append(len(protGroup)) 
    
    # Optional code to check stats on number of clusters, etc.
    if checkCounts:
        Counter = 0
        counters = []
        for cluster in clusterGroups:
            for row in cRows:
                if row[1] in cluster:
                    Counter += int(row[2])
            counters.append(Counter)
            Counter = 0
                
        with open("Summary_"+outfile,"w") as LengthFile:
            Writer = csv.writer(LengthFile, delimiter=",")
            Writer.writerow(["Protein Name:"]+titles)
            Writer.writerow(["num_proteins_countedFromClusterFile:"]+counters)
            Writer.writerow(["num_distinct_names:"]+protCount)
            Writer.writerow(["CountClustersPerGroup:"]+clusterCount)

    # Write a file containing all of the proteins in each backbone group
    if writeFile:
        with open("byRow_"+outfile,"w") as NewFile:
            Writer = csv.writer(NewFile,delimiter=",")
            for ng in nameGroups:
                Writer.writerow(ng)
                
        # NOTE: The following creates a transpose matrix.
        # I would recommend ignoring it and using Excel instead.
        Rows = []
        Rows.append(titles)
        Rows.append(clusterCount)
        for i in range(len(max(nameGroups,key=len))):
            Rows.append([])
            for j in range(len(nameGroups)):
                try:
                    if len(nameGroups[j][i])>0:
                        Rows[i+2].append(nameGroups[j][i])
                except IndexError:
                    Rows[i+2].append("")
                    
        with open(outfile,"w") as NewFile:
            Writer = csv.writer(NewFile,delimiter=",")
            for row in Rows:
                Writer.writerow(row)
                
    return nameGroups