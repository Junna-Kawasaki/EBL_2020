#!/usr/bin/env python

import sys, re, os
import pandas as pd
import numpy as np
import math
argvs = sys.argv
from ete3 import Tree, PhyloTree, NodeStyle, TreeStyle, add_face_to_node,TextFace, AttrFace, CircleFace, random_color, PieChartFace, COLOR_SCHEMES

# 0. input files and parameters
tree_path=argvs[1]
age_path = argvs[2]
borna_species = argvs[3]
output_path = argvs[4]

## add node attributions: e.g., EBL integration ages
def add_feature_f(s,tree,dic):
	for leaf in tree:
		try:
			leaf.add_feature(str(s), dic.get(str(leaf.name)))
		except:
			leaf.add_feature(str(s), "NA")
	return(tree)

# 1. read tree file
t = Tree(tree_path, format = 1)
leaf_l = [i.name for i in t]

# 2. read the file about EBL orthologous groups
age_d = {}
for line in open(age_path):
	group, name = line.strip().split("\t")
	if group not in age_d.keys():
		age_d[group] = []
	
	name = name.replace(";","_")
	name = name.replace("+","_")
	age_d[group].append(name)

# 3. Append attributes to each node
## name
name_d = {}
borna_l = ["Culter", "Carbo", "Orthoborna"]
virus_l = []
for node in t:
	name = node.name
	## ウイルスの分類
	if "Culter" in name:
		name_d[name] = "Cultervirus"
		virus_l.append(name)
	elif "Carbo" in name:
		name_d[name] = "Carbovirus"
		virus_l.append(name)
	elif "Orthoborna" in name:
		name_d[name] = "Orthobornavirus"
		virus_l.append(name)
	elif "Nyamiviridae" in name :
		name_d[name] = "Nyamiviridae"
		virus_l.append(name)
	else:
		species = re.sub(r"([A-Z][a-z_]+)_[A-Z]+.+", r"\1",name)
		name_d[name] = species

## orthologous groups
internal_nm = 0
no_leaf = 0
no_match = 0
ortholog_n = 0
ancestor_d = {}
no_group_d = {}
for g in age_d.keys():
	group_l = age_d[g]
	node_l = set(group_l) & set(leaf_l)
	if len(node_l) > 1:
		ancestor = t.get_common_ancestor(node_l)
		ancestor.add_features(nodetype="internal")
		leaf_a = [i.name for i in ancestor]
		# leaf_a = ancestor.get_leaf_names
		if len(set(leaf_a) - set(node_l)) == 0:
			ancestor.add_feature("internal_group",str(g))
			species_l = [re.sub(r'(.+)_[A-Z]{2}_.+',r'\1',n) for n in node_l]
			ancestor_d[g] = leaf_a
			internal_nm = internal_nm+1
		else:
			no_match +=1
			no_group_d[g] = leaf_a
			ancestor.add_feature("internal_group",str(g))
	else:
		no_leaf +=1

# 4. Calculate genetic distances among nodes
spe_u_l = []
spe_d_l = []
spe_l = []

for node in t.traverse():
    if "internal_group" in node.features:
        g = node.internal_group
        if g in ancestor_d.keys():
            node.add_feature("name", str(g))
            spe_l.append(node)
            #print (node.name)
    else:
        name = node.name
        for key in borna_l:
            if key in name:
                spe_l.append(node)

##
header = []
k = 0
for s1 in spe_l:
    h = s1.name
    if "Cultervirus" in s1.name:
        h = h.split("_")[-1]
    if "Carbovirus" in s1.name:
        k += 1
        h = h.split("_")[-1] + str(k)
    header.append(h)
    l_d = []
    for s2 in spe_l:
        if s1 == s2:
            d = 0
            l_d.append(d)
        else:
            d_d = t.get_distance(s1.name,s2.name)
            l_d.append(d_d)
    
    # extract minmum distance
    min_d = sorted(l_d)[1]
    spe_d_l.append(l_d)

df_d = pd.DataFrame(spe_d_l, index = [str(i) for i in header], columns=[str(i) for i in header])


## 5. Define genetic distance to separate different viral species
df_species = pd.read_table(borna_species,sep='\t', names=("node","name","species"))
borna_species_l = set([i.strip() for i  in df_species.species.to_list()])
rename_d = dict(zip(map(str, [i.strip() for i  in df_species.node]), map(str,[ i.strip() for i  in df_species.species])))
df_d = df_d.rename(columns=rename_d, index=rename_d)

def species_distance(df):
    criteria_l = []
    distance_d = {}
    for i in range(0,df.shape[0]):
        l = df.iloc[i,:]
        min_dis = sorted(l.tolist())[1]
        distance_d[l.name] = min_dis
        if l.name in borna_species_l or "Carbovirus" in l.name:
            criteria_l.append(min_dis)
    criteria_min = min(criteria_l)
    criteria_max = max(criteria_l)
    diverse_l = [name for name in distance_d.keys() 
                 if distance_d[name] > criteria_min and "EBL" in name]
    return (diverse_l, criteria_min, criteria_max)

diverse_d_l, criteria_d_min, criteria_d_max = species_distance(df_d)
print ("Genetic distance between different extant bornavirus species: " + str(round(criteria_d_min,3)))

## Output: genetic distances between ancient and modern bornaviruses
reorder = ["carbovirus", "clade1", "clade2", "carbo_EBL1", "carbo_EBL2",
           "orthobornavirus", "clade3", "clade4","clade5",
           "ortho_EBL1", "ortho_EBL2",
           "cultervirus", "clade6", "culter_EBL1", "culter_EBL2"]

rename_d = {'Carbovirus1':[0,"carbovirus"],
            'Carbovirus2':[0,"carbovirus"],
            'Cultervirus':[0,"cultervirus"],
            'EBLN0':["EBLN18","clade1"],
            'EBLN11': ["EBLN3", "clade5"],
            'EBLN12': ["EBLN11","clade2"],
            'EBLN13': ["EBLN21","clade6"],
            'EBLN14': ["EBLN20","clade6"],
            'EBLN15': ["EBLN17","clade2"],
            'EBLN16': ["EBLN6","clade6"],
            'EBLN17': ["EBLN12","clade2"],
            'EBLN18': ["EBLN22","clade2"],
            'EBLN2': ["EBLN19","clade1"],
            'EBLN21': ["EBLN24","clade3"],
            'EBLN22': ["EBLN5","clade4"],
            'EBLN23': ["EBLN25","ortho_EBL2"],
            'EBLN25': ["EBLN27","ortho_EBL2"],
            'EBLN27': ["EBLN29","ortho_EBL2"],
            'EBLN28': ["EBLN30", "clade5"],
            'EBLN3': ["EBLN14", "clade1"],
            'EBLN30': ["EBLN32", "clade5"],
            'EBLN31': ["EBLN33", "clade5"],
            'EBLN34': ["EBLN36", "culter_EBL2"],
            'EBLN35': ["EBLN37","culter_EBL2"],
            'EBLN36': ["EBLN38","culter_EBL1"],
            'EBLN38': ["EBLN40","carbo_EBL2"],
            'EBLN39': ["EBLN41","carbo_EBL1"],
            'EBLN4': ["EBLN4", "clade5"],
            'EBLN41': ["EBLN43","culter_EBL2"],
            'EBLN42': ["EBLN44","culter_EBL2"],
            'EBLN43': ["EBLN45","carbo_EBL2"],
            'EBLN45': ["EBLN47","ortho_EBL1"],
            'EBLN46': ["EBLN48","carbo_EBL2"],
            'EBLN47': ["EBLN49","carbo_EBL1"],
            'EBLN5': ["EBLN15","clade2"],
            'EBLN52': ["EBLN54","ortho_EBL1"],
            'EBLN57': ["EBLN59","carbo_EBL2"],
            'EBLN6': ["EBLN2","clade5"],
            'EBLN61': ["EBLN63","carbo_EBL1"],
            'EBLN62': ["EBLN64","ortho_EBL1"],
            'EBLN64': ["EBLN66","carbo_EBL1"],
            'EBLN7': ["EBLN1","clade5"],
            'EBLN8': ["EBLN10","clade2"],
            'EBLN9':["EBLN16","clade2"],
            'Elapid 1 orthobornavirus':[0,"orthobornavirus"],
            'Mammalian 1 orthobornavirus':[0,"orthobornavirus"],
            'Mammalian 2 orthobornavirus':[0,"orthobornavirus"],
            'Passeriform 1 orthobornavirus':[0,"orthobornavirus"],
            'Passeriform 2 orthobornavirus':[0,"orthobornavirus"],
            'Psittaciform 1 orthobornavirus':[0,"orthobornavirus"],
            'Psittaciform 2 orthobornavirus':[0,"orthobornavirus"],
            'Waterbird 1 orthobornavirus':[0,"orthobornavirus"]}

sort_l = []
rename_d1 = {}
for i in df_d.columns.tolist():
    info_l = rename_d[i]
    name,order = info_l[:2]
    
    ## rename
    if name == 0:
        rename = i
    else:
        rename = name
    
    rename_d1[i] = rename
    
    ## order
    o = reorder.index(order)
    sort_l.append(o)
    #print (i, name, rename, o)

df_d["sort"] = sort_l
df_d_s = df_d.sort_values("sort")
df_d_ss = df_d_s.drop("sort", axis=1)
df_d_st = df_d_ss.T
df_d_st["sort"] = sort_l
df_d_sts = df_d_st.sort_values("sort")
df_d_sss = df_d_sts.drop("sort", axis=1)

df_out = df_d_sss.rename(columns=rename_d1, index=rename_d1)
df_out.to_csv(output_path,sep="\t")

print ("Genetic distances between ancient and modern viruses were listed in " + str(output_path))

