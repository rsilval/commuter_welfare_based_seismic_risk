import networkx as nx
import csv
import pickle
import ita_cost_counting_paths
import mahmod as mah
import util
import numpy as np
import math
import sys
import codecs
import copy
import bd_test
from collections import defaultdict

#This is a function to compute distance between UTC coordinates
def distance_coordinates(coord1,coord2):
    #Input: Coordinates of point 1 and 2.
    #Output: Approximate distance for the Bay Area between both points
    x_dif = abs(coord1[0]-coord2[0])
    y_dif = abs(coord1[1]-coord2[1])
    lat_factor = 110
    long_factor = 85
    distance = ((x_dif*lat_factor)**2+(y_dif*long_factor)**2)**0.5
    return distance

def top_n(ranking,N):
    #Funtion that selects top N bridges out a dictionary that is a ranking
    norm_val = ranking.values()
    sorted_val = sorted(norm_val,reverse=True)
    rank_n = sorted_val[N-1]
    list_top = []
    for key in ranking:
        if ranking[key] >= rank_n:
            list_top.append(key)
    return list_top

#Function to assign centroids of census tracts into nodes of the transportation network
def centroidcvs_to_nodes(cvs_filename):
    #Input: CVS file from GIS
    #Output: Python dictionary. Keys is census tract (string), item is assigned node in network
    G = nx.read_gpickle("input/graphMTC_CentroidsLength3int.gpickle")
    pos = {node: (G.position[node][0], G.position[node][1]) for node in G.position}
    dict_block_node = {}
    node_ban = ['9886','9047','9398','9564','9282','20020','8977','8991','12236','1654']

    with open(cvs_filename) as csvfile:
        read_data = csv.reader(csvfile,delimiter=',')

        idx = 0
        for row in read_data:
            if idx > 0: #This skips the header
                block = float(row[3])
                pos_x = float(row[0])
                pos_y = float(row[1])
                coords = (pos_y,pos_x)

                dis_min = 1000

                for node in pos:

                    if node in node_ban:
                        continue

                    dist = distance_coordinates(pos[node],coords)

                    if dist < dis_min:
                        node_block = node
                        if node not in G.nodes():
                            print 'Error in position key'
                        dis_min = dist
                dict_block_node[block] = node_block
                pos.pop(node_block)
            idx += 1
    return dict_block_node

#Function that creates a dictionary of demand using LODES data. Goes from census tract to demand in nodes of the network
def create_demand(cvs_filename,tracts):
    #Input: 1- CVS information from LODES Data. 2- List of OD (tracts)
    #Output: 1- Dictionary with demand. Keys: Tuple with nodes. Item: Number for demand. 2-Dictionary with deaggregated data from the table. Key: (oridnates pair). Item (list of 6 items with fractions of item from columns)
    dict_demand = {}
    dict_deag_data = {}
    for tract1 in tracts:
        for tract2 in tracts:
            pair = (float(tract1),float(tract2))
            dict_demand[pair] = 0
            dict_deag_data[pair] =[0,0,0,0,0,0]
    count = len(open(cvs_filename).readlines())


    read_data = csv.reader(codecs.open(cvs_filename, 'rb'))
    demand_sum = 0
    for i in range(0,count):
        row = next(read_data)
        if i > 0: #Skip first line
            #print i
            block_d = str(row[0])
            tract_d = float(block_d[0:11])
            block_o = str(row[1])
            tract_o = float(block_o[0:11])
            if tract_o in tracts and tract_d in tracts:
                dict_demand[(tract_o,tract_d)] =  dict_demand[(tract_o,tract_d)] + float(row[2])
                features_list = dict_deag_data[(tract_o,tract_d)]
                for i in range(0,6):
                    features_list[i] = features_list[i] + float(row[6+i])
                dict_deag_data[(tract_o, tract_d)] = features_list
    return dict_demand,dict_deag_data


#Code to generate dictionary of nodes assigned to census tract nodes.

cvs_filename = 'CensusTractCentroids.csv'
dict_block_node = centroidcvs_to_nodes(cvs_filename)
pickle.dump(dict_block_node, open( "dict_block_node.pkl", "wb" ) )

#print 'Done Assigning centroid blocks to nodes'

tracts = dict_block_node.keys()
for i in range(0,len(tracts)):
    tracts[i] = float(tracts[i])

#Code to generate dictionary of OD pairs and demand

cvs_filename2 = 'input/ca_od_main_JT00_2017.csv'
dict_demand,dict_deag_data = create_demand(cvs_filename2,tracts)
pickle.dump(dict_demand, open( "dict_demand.pkl", "wb" ) )
pickle.dump(dict_deag_data, open( "dict_deag_data.pkl", "wb" ) )
print 'Done assigning dictionary of demand'

#This part of the code load precomputed information from previous lines of code

dict_demand = pickle.load( open( "dict_demand.pkl", "rb" ) )
dict_block_node = pickle.load( open( "dict_block_node.pkl", "rb" ) )
dict_deag_data = pickle.load( open( "dict_deag_data.pkl", "rb" ) )

#Input to paralelize in sherlock THIS IS THE ONLY INPUT REQUIRED TO RUN THE CODE
index_scenario = float(sys.argv[1])
index_scenario = int(index_scenario)
scenario = index_scenario #Scenario index in the file with all of the scenarios
scenario = int(scenario)
st_scenario = str(scenario)

#This segment of the code transforms the demand information into a format that is used by the ITA code. This only needs to be run once.

dict_demand_ita = {}
dict_demand_ita_deag = {}

for pair in dict_demand:
    or_node = dict_block_node[float(pair[0])]
    dict_demand_ita[or_node] = {}
    dict_demand_ita_deag[or_node] = {}

for pair in dict_demand:
    or_node = dict_block_node[float(pair[0])]
    de_node = dict_block_node[float(pair[1])]
    demand_pair = dict_demand[pair]
    demand_deag = dict_deag_data[pair]
    dict_local = dict_demand_ita[or_node]
    dict_local.update({de_node:demand_pair})
    dict_local_deag = dict_demand_ita_deag[or_node]
    dict_local_deag.update({de_node:demand_deag})
    dict_demand_ita[or_node] = dict_local
    dict_demand_ita_deag[or_node] = dict_local_deag

pickle.dump(dict_demand_ita, open( "dict_demand_ita.pkl", "wb" ) )
pickle.dump(dict_demand_ita_deag, open( "dict_demand_ita_deag.pkl", "wb" ) )


#Load precomputed information of demand in ITA format


dict_demand_ita = pickle.load( open( "dict_demand_ita.pkl", "rb" ) )
dict_demand_ita_deag = pickle.load( open( "dict_demand_ita_deag.pkl", "rb" ) )
counting = copy.deepcopy(dict_demand_ita)

for origin in dict_demand_ita:
    for destination in dict_demand_ita[origin]:
        counting[origin][destination] = [0,0]

idx_analysis = 2 #This is an index that chooses between: 1: Computing travel times by using the congested path of the first iteration. 2: Using the congested path of each iteration

#Load information of the grpah that represents the transportation network in SF bay area
G = nx.read_gpickle("input/graphMTC_CentroidsLength3int.gpickle")
#G = mah.add_superdistrict_centroids(G)
t_a_dict = nx.get_edge_attributes(G, 't_a')

#Fix floating edges

for u,v,d in G.edges(data=True):
    if d['t_0'] == 0:
        d['t_0'] = float('inf')

for u,v,d in G.edges(data=True):
    if d['t_a'] == 0:
        d['t_a'] = float('inf')

#This part of the code assigns traffic given a graph and LODES demand for the undamaged version of the graph
#We compute the value of trips between census tracts for the undamaged situation. This doesn't need to be run on each simulation

it = ita_cost_counting_paths.ITA(G, dict_demand_ita,counting,idx_analysis)
newG,counting_new = it.assign()
counting_new = it.time_paths()

nameb = 'G_nd.pkl'
pathb = 'output/'+nameb

with open(pathb ,'wb') as f:
    pickle.dump(G, f)

namec = 'counting_nd.pkl'
pathc = 'output/'+namec
with open(pathc ,'wb') as f:
    pickle.dump(counting_new, f)
print 'Done Undamaged'

with open(pathb, 'rb') as f:
    G_nd = pickle.load(f)

'''Section of the code that fixes nodes that do not have capacity to satisfy demand'''

#Check for node in demand:
capacity_dict = nx.get_edge_attributes(G,'capacity_0')
time_0 = nx.get_edge_attributes(G,'t_0')
time_a = nx.get_edge_attributes(G,'t_a')
for origin in dict_demand_ita:
    demand_origin = 0
    for destination in dict_demand_ita[origin]:
        demand_origin += dict_demand_ita[origin][destination]
    #Call neighbors
    neighbors = G_nd.neighbors(origin)
    if len(neighbors) == 1:
        edge_1 = (origin,neighbors[0])
        edge_2 = (neighbors[0],origin)
        if edge_1 in time_a:
            if time_a[edge_1] == float('inf'):
                time_a[edge_1] = 0
        if edge_2 in time_a:
            if time_a[edge_2] == float('inf'):
                time_a[edge_2] = 0
    if len(neighbors) == 2:
        for node_ne in neighbors:
            edge_1 = (origin,neighbors[0])
            edge_2 = (neighbors[0],origin)
            if edge_1 in time_a:
                if time_a[edge_1] == float('inf'):
                      time_a[edge_1] = 0
            if edge_2 in time_a:
                if time_a[edge_2] == float('inf'):
                    time_a[edge_2] = 0
    if len(neighbors) == 3:
        for node_ne in neighbors:
            edge_1 = (origin, neighbors[0])
            edge_2 = (neighbors[0], origin)
            flag = 0
            if edge_1 in time_a:
                if time_a[edge_1] == float('inf'):
                    flag +=1
            if edge_2 in time_a:
                if time_a[edge_2] == float('inf'):
                    flag += 1
        if flag == 3:
            for node_ne in neighbors:
                edge_1 = (origin, neighbors[0])
                edge_2 = (neighbors[0], origin)
                flag = 0
                if edge_1 in time_a:
                    if time_a[edge_1] == float('inf'):
                        time_a[edge_1] = 0
                if edge_2 in time_a:
                    if time_a[edge_2] == float('inf'):
                        time_a[edge_2] = 0
    if len(neighbors) == 4:
        for node_ne in neighbors:
            edge_1 = (origin, neighbors[0])
            edge_2 = (neighbors[0], origin)
            flag = 0
            if edge_1 in time_a:
                if time_a[edge_1] == float('inf'):
                    flag +=1
            if edge_2 in time_a:
                if time_a[edge_2] == float('inf'):
                    flag += 1
        if flag == 4:
            for node_ne in neighbors:
                edge_1 = (origin, neighbors[0])
                edge_2 = (neighbors[0], origin)
                flag = 0
                if edge_1 in time_a:
                    if time_a[edge_1] == float('inf'):
                        time_a[edge_1] = 0
                if edge_2 in time_a:
                    if time_a[edge_2] == float('inf'):
                        time_a[edge_2] = 0

    #Check that capacity of edges is bigger than demand
    capacity_local = 0
    for node_nei in neighbors:
        edge_1 = (origin,node_nei)
        edge_2 = (node_nei,origin)
        if edge_1 in capacity_dict:
            capacity_local += capacity_dict[edge_1]
        elif edge_2 in capacity_dict:
            capacity_local += capacity_dict[edge_2]
    if capacity_local < demand_origin:
        for node_nei in neighbors:
            edge_1 = (origin, node_nei)
            edge_2 = (node_nei, origin)
            if edge_1 in capacity_dict:
                capacity_dict[edge_1] = float('inf')
                time_0[edge_1] = 0

            elif edge_2 in capacity_dict:
                capacity_dict[edge_2] = float('inf')
                time_0[edge_2] = 0

nx.set_edge_attributes(G,'capacity_0',capacity_dict)
nx.set_edge_attributes(G,'t_0',time_0)
nx.set_edge_attributes(G,'t_a',time_a)

''' Generate Damaged Version of the outputs '''

#Load scenarios

sa_matrix = util.read_2dlist('input/sample_ground_motion_intensity_maps_road_only_filtered.txt',delimiter='\t')#Using Mahalia info
lnsas = []

with open('input/20140114_master_bridge_dict.pkl', 'rb') as f:
    master_dict = pickle.load(f)

#Generate damaged bridges
runs = {'loFrag': 'mod_lnSa', 'hiFrag':'com_lnSa', 'loCoef':.75, 'hiCoef':1.25} #Simplified retrofitting information

idx_retro = 0
n_retro = 100
if idx_retro == 0:
    retrofitted = {site:False for site in master_dict} #In this case no retroffitting action is taken.
elif idx_retro == 1:
    retrofitted = {}
    with open('ranking_bridges_welfare.pkl', 'rb') as f:
        ranking_retro = pickle.load(f)
    top_bridges = top_n(ranking_retro,n_retro)
    for site in master_dict:
        if site not in top_bridges:
            retrofitted[site] = False
        else:
            retrofitted[site] = True

elif idx_retro == 2:
    with open('NNsensitivity2.pkl', 'rb') as f:
        ranking_retro = pickle.load(f)
    top_bridges = top_n(ranking_retro, n_retro)
    retrofitted = {}
    for site in master_dict:
        if site not in top_bridges:
            retrofitted[site] = False
        else:
            retrofitted[site] = True

elif idx_retro == 3:
    with open('ranking_bridges_low_income.pkl', 'rb') as f:
        ranking_retro = pickle.load(f)
    top_bridges = top_n(ranking_retro, n_retro)
    retrofitted = {}
    for site in master_dict:
        if site not in top_bridges:
            retrofitted[site] = False
        else:
            retrofitted[site] = True

last_idx = scenario # This is the index defined previously for the seismic scenario

#Predefinition of counting dictionary for the damaged version

counting_damage = copy.deepcopy(dict_demand_ita)
for origin in dict_demand_ita:
    for destination in dict_demand_ita[origin]:
        counting_damage[origin][destination] = [0, 0]
idx = 0
for row in sa_matrix:
    if idx == last_idx:
        lnsas.append([math.log(float(sa)) for sa in row[3:1746]])
    idx +=1

#Obtain a Damage Realization from Mahalia's Model

last_idx, damaged_bridges_internal, damaged_bridges_new, num_damaged_bridges = mah.compute_damage(lnsas[0], master_dict, last_idx, retrofitted, runs)
path_bridge_internal = 'damaged_bridges_internal.pkl'
with open(path_bridge_internal, 'rb') as f:
     damaged_bridges_internal = pickle.load(f)
damaged_bridges_internal = map(str, damaged_bridges_internal)

#Updates graph with damaged bridges

G_dam, roads_out = mah.damage_highway_network(damaged_bridges_internal, G_nd, master_dict, last_idx)
namebridges = 'bridges_out' + st_scenario + '.pkl'
path_bridges = 'output/'+namebridges
with open(path_bridges ,'wb') as f:
     pickle.dump(damaged_bridges_internal, f)

#Assing traffic for damaged version of the network

it_dam = ita_cost_counting_paths.ITA(G_dam,dict_demand_ita,counting_damage,idx_analysis)
newG2, counting_damage2 = it_dam.assign()
counting_damage2 = it_dam.time_paths()

namedam = 'counting_damage_'+st_scenario+'.pkl'
pathdam = 'output/'+namedam
with open(pathdam ,'wb') as f:
    pickle.dump(counting_damage2, f)

name1 = 'G_dam_'+st_scenario+'.pkl'
path1 = 'output/'+name1
with open(path1 ,'wb') as f:
    pickle.dump(it_dam.G, f)

