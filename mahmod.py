#Author: Mahalia Miller
#Date: February 10, 2014

import pickle, random, pp, pdb, time, networkx
import scipy
import scipy.stats
import math
import util, ita, bd

random.seed(1)
def compute_flow(damaged_graph):
	'''compute max flow between a start and end'''
	s=  'sf' #another option is '1000001'
	t = 'oak' #other options are '1000002' and 'sfo'
	try:
		flow = networkx.maximum_flow(damaged_graph, s, t, capacity='capacity') #not supported by multigraph
	except networkx.exception.NetworkXError as e:
		print ('found an ERROR: ', e)
		pdb.set_trace()
	return flow


def compute_shortest_paths(damaged_graph, demand):
	return -1


def compute_tt_vmt(damaged_graph, demand):
	#start = time.time()
	it = ita.ITA(damaged_graph,demand)
	newG = it.assign()
	#print('time to assign: ', time.time()-start)
	travel_time, probe = util.find_travel_time(damaged_graph) 
	vmt = util.find_vmt(damaged_graph)  
	''' in the undamaged case, this should be around 172 million (http://www.mtc.ca.gov/maps_and_data/datamart/stats/vmt.htm) over the course of a day, so divide by 0.053 (see demand note in main). BUT our trip table has only around 11 million trips (instead of the 22 million mentioned here: http://www.mtc.ca.gov/maps_and_data/datamart/stats/baydemo.htm because we are looking at vehicle-driver only and not transit, walking, biking, being a passenger in a car, etc. So, that's **8-9 million vehicle-miles divided by 2, which equals around 4 million vehicle-miles!**
	'''
	return travel_time, vmt, probe

def compute_tt_vmtRSL(damaged_graph, demand):
	#start = time.time()
	it = ita.ITA(damaged_graph,demand)
	newG = it.assign()
	#print('time to assign: ', time.time()-start)
	travel_time, probe = util.find_travel_time_sfRSL(damaged_graph)
	vmt = util.find_vmt(damaged_graph)
	''' in the undamaged case, this should be around 172 million (http://www.mtc.ca.gov/maps_and_data/datamart/stats/vmt.htm) over the course of a day, so divide by 0.053 (see demand note in main). BUT our trip table has only around 11 million trips (instead of the 22 million mentioned here: http://www.mtc.ca.gov/maps_and_data/datamart/stats/baydemo.htm because we are looking at vehicle-driver only and not transit, walking, biking, being a passenger in a car, etc. So, that's **8-9 million vehicle-miles divided by 2, which equals around 4 million vehicle-miles!**
	'''
	return travel_time, vmt, probe

def add_superdistrict_centroids(G):
	'''adds 34 dummy nodes for superdistricts'''
	sd_table = util.read_2dlist('input/superdistricts_clean.csv', ',', False)
	#for each superdistrict, create a dummy node. Make 2 directed edges from the dummy node to real nodes. Make 2 directed edges from real edges to dummy nodes.
	for row in sd_table:
		i = int(row[0])
		G.add_node(str(1000000 + i))
		G.add_edge(str(1000000 + i), str(row[1]), capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
		G.add_edge(str(1000000 + i), str(row[2]), capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
		G.add_edge(str(row[3]), str(1000000 + i), capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
		G.add_edge(str(row[4]), str(1000000 + i), capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 

	#add a sf dummy node, an oakland dummy node, and a SFO dummy node for max flow
	G.add_node('sf')
	G.add_node('oak')
	G.add_node('sfo')
	G.add_edge('sf', '1000001', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('sf', '1000002', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('sf', '1000003', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('sf', '1000004', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('sf', '1000005', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('1000018', 'oak', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('1000019', 'oak', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('1000020', 'oak', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('6564', 'sfo', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('6563', 'sfo', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('6555', 'sfo', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('9591', 'sfo', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('6550', 'sfo', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	G.add_edge('6599', 'sfo', capacity_0 = 100000,  capacity = 100000, lanes =1 , bridges=[], distance_0=1, distance = 1, t_a=1, t_0=1, flow=0, dailyvolume=1) #capacity in vehicles over all lanes, travel time in seconds, length in miles, flow in 
	return G


def damage_bridges(scenario, master_dict, retrofitted, runs):
	'''This function damages bridges based on the ground shaking values (demand) and the structural capacity (capacity). It returns two lists (could be empty) with damaged bridges (same thing, just different bridge numbering'''
	from scipy.stats import norm
	damaged_bridges_new = []
	damaged_bridges_internal = []

	#first, highway bridges and overpasses
	beta = 0.6 #you may want to change this by removing this line and making it a dictionary lookup value 3 lines below
	for site in master_dict.keys(): #1-1889 in Matlab indices (start at 1)

		if retrofitted[site]:			
			daSa = runs['hiCoef']*master_dict[site][runs['hiFrag']]
			daSa = 3*master_dict[site][runs['hiFrag']]
		else:			
			daSa = runs['loCoef']*master_dict[site][runs['loFrag']]

		lnSa = scenario[master_dict[site]['new_id'] - 1]
		prob_at_least_ext = norm.cdf((1/float(beta)) * (lnSa - math.log(daSa)), 0, 1) #want to do moderate damage state instead of extensive damage state as we did here, then just change the key name here (see master_dict description)
		U = random.uniform(0, 1)
		if U <= prob_at_least_ext:
			damaged_bridges_new.append(master_dict[site]['new_id']) #1-1743
			damaged_bridges_internal.append(site) #1-1889

	num_damaged_bridges = sum([1 for i in damaged_bridges_new if i <= len(master_dict.keys())])
	return damaged_bridges_internal, damaged_bridges_new, num_damaged_bridges


def compute_damage(scenario, master_dict, index, retrofitted, runs):
	'''goes from ground-motion intensity map to damage map '''
	#figure out component damage for each ground-motion intensity map
	damaged_bridges_internal, damaged_bridges_new, num_damaged_bridges = damage_bridges(scenario, master_dict, retrofitted, runs) #e.g., [1, 89, 598] #num_bridges_out is highway bridges only
 	return index, damaged_bridges_internal, damaged_bridges_new, num_damaged_bridges


def damage_highway_network(damaged_bridges_internal, G, master_dict, index):
	'''damaged bridges is a list of the original ids (1-1889, not the new ids 1-1743!!!!!!!) '''
	biggest_id_of_interest = max([int(k) for k in master_dict.keys()])
	road_bridges_out = sum([1 for i in damaged_bridges_internal if int(i) <= biggest_id_of_interest])
	try:
		if len(damaged_bridges_internal) > 0:
			b = damaged_bridges_internal[0].lower()
	except AttributeError:
		raise('Sorry. You must use the original ids, which are strings')
	list_of_u_v = []
	counter = 0
	for site in damaged_bridges_internal:
		if int(site) <= 1889: #in original ids, not new ones since that is what is in the damaged bridges list
			affected_edges = master_dict[site]['a_b_pairs_direct'] + master_dict[site]['a_b_pairs_indirect']
			list_of_u_v += affected_edges
			for [u,v] in affected_edges:

				G[str(u)][str(v)]['t_a'] = float('inf')
				G[str(u)][str(v)]['capacity'] = 0 
				G[str(u)][str(v)]['distance'] = 20*G[str(u)][str(v)]['distance_0']
				counter += 1
	if counter < road_bridges_out:
		print damaged_bridges_internal
		print road_bridges_out
		print counter
	# assert counter >= road_bridges_out, 'we should impact an edge per bridge minimum'
	return G, road_bridges_out


def measure_performance(damaged_graph, demand):
	# returns flow, shortest_paths, travel_time, vmt
	flow = compute_flow(damaged_graph)
	shortest_paths = compute_shortest_paths(damaged_graph, demand)
	travel_time, vmt, probe = compute_tt_vmt(damaged_graph, demand)
	return flow, shortest_paths, travel_time, vmt, probe

def measure_performanceRSL(damaged_graph, demand):
	# returns flow, shortest_paths, travel_time, vmt
	flow = compute_flow(damaged_graph)
	shortest_paths = compute_shortest_paths(damaged_graph, demand)
	travel_time, vmt, probe = compute_tt_vmtRSL(damaged_graph, demand)
	return flow, shortest_paths, travel_time, vmt, probe


def compute_road_performance(G, damaged_bridges_internal, demand, no_damage_travel_time, no_damage_vmt, no_damage_flow, no_damage_shortest_path, no_damage_probe, master_dict, index):
	'''computes network performance after damaging the network based on which bridges are damaged'''
	#start_time = time.time()
	if G == None:
		G = get_graph()
	#figure out road network damage
	if len(damaged_bridges_internal) > 0:
		G, road_bridges_out = damage_highway_network(damaged_bridges_internal, G, master_dict, index)
		#figure out impact (performance metrics)
		flow, shortest_paths, travel_time, vmt, probe = measure_performance(G, demand)
		G = util.clean_up_graph(G) #absolutely critical. otherwise, damage from scenario gets added to damage from previous scenarios!
	else: #no bridges are damaged, so no need to do all the calculations
		flow = no_damage_flow
		shortest_paths = no_damage_shortest_path
		travel_time = no_damage_travel_time
		vmt = no_damage_vmt
		road_bridges_out = 0
		probe = no_damage_probe
	#print 'total network performance calculation time: ', time.time() - start_time
	return index, road_bridges_out, flow, shortest_paths, travel_time, vmt, probe

def compute_road_performanceRSL(G, damaged_bridges_internal, demand, no_damage_travel_time, no_damage_vmt, no_damage_flow, no_damage_shortest_path, no_damage_probe, master_dict, index):
	'''computes network performance after damaging the network based on which bridges are damaged'''
	#start_time = time.time()
	if G == None:
		G = get_graph()
	#figure out road network damage
	if len(damaged_bridges_internal) > 0:
		G, road_bridges_out = damage_highway_network(damaged_bridges_internal, G, master_dict, index)
		#figure out impact (performance metrics)
		flow, shortest_paths, travel_time, vmt, probe = measure_performanceRSL(G, demand)
		G = util.clean_up_graph(G) #absolutely critical. otherwise, damage from scenario gets added to damage from previous scenarios!
	else: #no bridges are damaged, so no need to do all the calculations
		flow = no_damage_flow
		shortest_paths = no_damage_shortest_path
		travel_time = no_damage_travel_time
		vmt = no_damage_vmt
		road_bridges_out = 0
		probe = no_damage_probe
	#print 'total network performance calculation time: ', time.time() - start_time
	return index, road_bridges_out, flow, shortest_paths, travel_time, vmt, probe


def get_graph():
	import networkx
	'''loads full mtc highway graph with dummy links and then adds a few fake centroidal nodes for max flow and traffic assignment'''

	G = networkx.read_gpickle("input/graphMTC_CentroidsLength3int.gpickle")
	G = add_superdistrict_centroids(G)

	assert not G.is_multigraph() # Directed! only one edge between nodes
	G = networkx.freeze(G) # prevents edges or nodes to be added or deleted
	return G


def main(retrofitted, runs, alias):
	'''this is the main file that runs from ground-motion intensity map to network performance measure. You will  need to adjust various things below, such as the ground motion files, performance measure info and more. you should not need to change, however, the functions that they call'''
	seed_num = 0 #USER ADJUSTS THIS! other value examples: 1,2, 11, 14, ...
	random.seed(seed_num) #set random number generator seed so we can repeat this process

	'''
		Ground-motion intensity map data
	'''
	#load the earthquake info	
	sa_matrix = util.read_2dlist('input/sample_ground_motion_intensity_maps_road_only_filtered.txt',delimiter='\t')

	# Import disaster events
	# this does approx. 2000 ground-motion intensity maps (reduced to 249). These are hazard consistent.
	# sa_matrix = util.read_2dlist('input/sample_ground_motion_intensity_maps_road_only_filtered.txt',delimiter='\t')

	lnsas = []
	magnitudes = []
	rates = []
	indices = []

	for row in sa_matrix:
		lnsas.append([log(float(sa)) for sa in row[4:]])
		indices.append(float(row[0]))
		magnitudes.append(row[2])
		rates.append(float(row[3]))

	# component (bridge) damage map data
	sets = 1 		# number of bridge damage maps per ground-motion intensity map. USER ADJUSTS THIS! other value examples: 3,9,18
	targets = range(0, len(lnsas)*sets) #define the damage map IDs you want to consider. Note: this currently does not require modification. Just change the number of sets above.
	print '\tGot '+str(len(lnsas))+' GM maps, and '+str(len(lnsas[0]))+" site locations"
	print '\tConsidering %d damage maps (%d per ground-motion intensity map).' % (int(len(targets)), int(sets))
	
	targets = range(0,60)  # temporary trial
	#first load the all-purpose dictionary linking info about the bridges
	with open('input/20140114_master_bridge_dict.pkl','rb') as f:
		master_dict = pickle.load(f) #has 1743 keys. One per highway bridge. (NOT BART)	
		'''
		dict where the keyranges from 1 to 1889 and then the value is another dictionary with the following keys: 
		loren_row_number: the row number in Loren Turner's table that has info on all CA bridges (where the header line is row 0)
		original_id: the original id (1-1889) 
		new_id: the new id that excludes filtered out bridges (1-1743). Bridges are filtered out if a.) no seismic capacity data AND non-transbay bridge or b.) not located by Jessica (no edge list). this id is the new value that is the column number for the lnsa simulations.
		jessica_id: the id number jessica used. it's also the number in arcgis.
		a_b_pairs_direct: list of (a,b) tuples that would be directly impacted by bridge damage (bridge is carrying these roads)
		a_b_pairs_indirect: ditto but roads under the indirectly impacted bridges
		edge_ids_direct: edge object IDS for edges that would be directly impacted by bridge damage
		edge_ids_indirect: ditto but roads under the indirectly impacted bridges
		mod_lnSa: median Sa for the moderate damage state. the dispersion (beta) for the lognormal distribution is 0.6. (See hazus/mceer method)
		ext_lnSa: median Sa for the extensive damage state. the dispersion (beta) for the lognormal distribution is 0.6. (See hazus/mceer method)
		com_lnSa: median Sa for the complete damage state. the dispersion (beta) for the lognormal distribution is 0.6. (See hazus/mceer method)
		'''	

	num_of_interest_bridges = len(master_dict)
	num_of_total_bridges = len(master_dict)

	# network damage map data 
	G = get_graph()
	t_a = networkx.get_edge_attributes(G, 't_a')
	for element in t_a:
		if t_a[element] < 5:
			t_a[element] = 300
	networkx.set_edge_attributes(G, 't_a', t_a)

	t_0 = networkx.get_edge_attributes(G, 't_0')
	for element in t_a:
		if t_0[element] < 5:
			t_0[element] = 300
	networkx.set_edge_attributes(G, 't_0', t_0)

	cap = networkx.get_edge_attributes(G, 'capacity')
	for element in cap:
		if cap[element] < 5:
			cap[element] = 600
	networkx.set_edge_attributes(G, 'capacity', cap)


	# network performance map data
	# compute what the travel time and vehicle-miles-traveled values are without any damage
	demand = bd.build_demand('input/BATS2000_34SuperD_TripTableData.csv', 'input/superdistricts_centroids_dummies.csv') #we just take a percentage in ita.py, namely  #to get morning flows, take 5.3% of daily driver values. 11.5/(4.5*6+11.5*10+14*4+4.5*4) from Figure S10 of http://www.nature.com/srep/2012/121220/srep01001/extref/srep01001-s1.pdf. Note: these are vehicle-driver trips only (not transit, biking, walking, etc.)
	# pre-compute the network performance measures when there is no damage to save time later
	no_damage_travel_time, no_damage_vmt, no_damage_probe = compute_tt_vmtRSL(G, demand)
	no_damage_flow = compute_flow(G)
	no_damage_shortest_path = -1
	G = util.clean_up_graph(G) #so the trips assigned don't hang around

	# This is an intervention for SP4SF // 'compute_tt_vmt' modified to report arc times in probe
	#with open('output/probe_'+alias+'_NoDam.pkl','wb') as f:
	#	pickle.dump(no_damage_probe, f) #has 1743 keys. One per highway bridge. (NOT BART)

	# actually run damage map creation
	ppservers = ()    #starting a super cool parallelization 
	# Creates jobserver with automatically detected number of workers
	job_server = pp.Server(ppservers=ppservers)
	#print "\t\tStarting pp with", job_server.get_ncpus(), "workers"
	# set up jobs
	jobs = []
	for i in targets:
		jobs.append(job_server.submit(compute_damage, (lnsas[i%len(lnsas)], master_dict, targets[i], retrofitted, runs, ), modules = ('random', 'math', ), depfuncs = (damage_bridges, ))) 

	# get the results that have already run
	bridge_array_new = []
	bridge_array_internal = []
	indices_array = []
	bridge_array_hwy_num = []

	for job in jobs:
		(index, damaged_bridges_internal, damaged_bridges_new, num_damaged_bridges_road) = job()
		bridge_array_internal.append(damaged_bridges_internal)
		bridge_array_new.append(damaged_bridges_new)
		indices_array.append(index)
		bridge_array_hwy_num.append(num_damaged_bridges_road)
	print '\tDamage maps made; computing performance...'

	# #################################################################
	# ################## actually run performance measure realization creation #######################
	ppservers = ()    
	# Creates jobserver with automatically detected number of workers
	job_server = pp.Server(ppservers=ppservers)
	#print "\t\tStarting pp with", job_server.get_ncpus(), "workers"
	# set up jobs
	jobs = []

	for i in targets:
		jobs.append(job_server.submit(compute_road_performanceRSL, (None, bridge_array_internal[i], demand, no_damage_travel_time, no_damage_vmt, no_damage_flow, no_damage_shortest_path, no_damage_probe, master_dict, targets[i], ), modules = ('networkx', 'time', 'pickle', 'pdb', 'util', 'random', 'math', 'ita', ), depfuncs = (get_graph, add_superdistrict_centroids, damage_bridges, damage_highway_network, measure_performanceRSL, compute_flow, compute_shortest_paths, compute_tt_vmtRSL, ))) # functions, modules

	# get the results that have already run and save them	X
	LesProbes = {}
	i = 0
	for job in jobs:
		(index,  road_bridges_out, flow, shortest_paths, travel_time, vmt, probe) = job()
		assert indices_array[i] == index, 'the damage maps should correspond to the performance measure realizations'
		assert bridge_array_hwy_num[i] == road_bridges_out, 'we should also have the same number of hwy bridges out'		
		LesProbes[i] = (indices[i], rates[i], probe) 	# first item lets you know what lnsa it is (regardless of set)
		#with open('output/probe_'+alias+'_scen_'+str(i)+'_.pkl', 'wb') as f:
		#	pickle.dump(LesProbes, f)
		i += 1	
	LesProbes[0] = (0,1,no_damage_probe)
	#LesProbes[0] = (0,1-sum(rates),no_damage_probe)

	print '\tNetwork performance calculated.'
	return LesProbes
