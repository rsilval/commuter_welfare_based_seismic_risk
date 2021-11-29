# GB NOTE -- last modified 21 January 2021 to correct welfare coefficients and welfare loss due to lost trips.

import pickle
import sys
import time

# 1- Dictionary using origin census tract as keys and the value as the mean travel time - damaged condition
# 2- Dictionary using destination census tract as keys and the value as the mean travel time - damaged condition
# 3- Dictionary using origin census tract as keys and the value as the mean travel time - undamaged condition (this is not scenario dependent)
# 4- Dictionary using destination census tract as keys and the value as the mean travel time - undamaged condition (this is not scenario dependent)
# 5- Dictionary using origin census tract as keys and the value as lost trips.
# 6- Dictionary using destination census tract as keys and the value as lost trips.
# 7- Dictionary using origin census tract as keys and the value as well-being loss.
# 8- Dictionary using destination census tract as keys and the value as well-being loss.
# 9- Dictionary using origin census tract as keys and the value as welfare loss.
# 10- Dictionary using destination census tract as keys and the value as welfare loss.
# 11 - Dictionary using origin census tract as keys and the value as cost of road network disruption using Hackl et al. (2018) formula.
# 12 - Dictionary using destination census tract as keys and the value as cost of road network disruption using Hackl et al. (2018) formula.

	#Define number of scenario to be used in naming of output files and loading the input files
	# scenario = float(sys.argv[1]) # GB MODIFICATION -- to run automatically
	# scenario = int(scenario) # GB MODIFICATION
	#scenario = 1000 #Change to reflect scenario that is wanted to be postprocessed

scenario = float(sys.argv[1])
scenario = int(scenario)
st_scenario = str(scenario)

upper_bound = 4
fl_upper_bound = float(upper_bound)

sf_specific = False


fl_upper_bound = float(upper_bound)
if fl_upper_bound % 1 == 0.5:
	st_upper_bound = str(int(fl_upper_bound * 10))
else:
	st_upper_bound = str(upper_bound)

if sf_specific == 'True':
	output_directory = 'output_sf_' + st_upper_bound + 'h/' # stores outputs based on regional income groups and associated social value of travel time
else:
	output_directory = 'output/'# + st_upper_bound + 'h/' # stores outputs based on original (LODES) income groups and associated social value of travel time

output_directory = 'output/'
#Load information for demand
dict_demand_ita = pickle.load( open( "dict_demand_ita.pkl", "rb" ) )

dict_demand_ita_deag_file = 'dict_demand_ita_deag.pkl'

with open(dict_demand_ita_deag_file, 'rb') as f:
	dict_demand_ita_deag = pickle.load(f)

dict_block_node = pickle.load(open( "dict_block_node.pkl", "rb" ))
dict_node_block = {v: k for k, v in dict_block_node.iteritems()}

value_demand_ita = 0
for origin in dict_demand_ita:
	for destination in dict_demand_ita[origin]:
		value_demand_ita += dict_demand_ita[origin][destination]

#Load information from traffic model in Module_T
namec = 'counting_nd.pkl'
pathc = 'output/'+namec
counting = pickle.load(open(pathc,"rb"))

trips = []
for origin in counting:
	for destination in counting[origin]:
		for i in range(0,int(counting[origin][destination][0])):
			trips.append(counting[origin][destination][1])



name_dam = 'counting_damage_'+st_scenario+'.pkl'
path_dam = 'output/'+name_dam
counting_damage = pickle.load(open(path_dam,"rb"))
'''
### Generating Travel dictionaries ###

census_tract_info_origin = {}
census_tract_info_destination = {}

for block in dict_block_node:
	census_tract_info_origin[block] = 0
	census_tract_info_destination[block] = 0

for origin in counting:
	block_origin = dict_node_block[origin]
	for destination in counting[origin]:
		block_des = dict_node_block[destination]
		demand_local = counting[origin][destination][0]
		census_tract_info_origin[block_origin] += demand_local
		census_tract_info_destination[block_des] += demand_local

pickle.dump(census_tract_info_origin, open( "census_tract_info_origin.pkl", "wb" ) )
pickle.dump(census_tract_info_destination, open( "census_tract_info_destination.pkl", "wb" ) )
'''
###### Average travel time for census tract #### (Variables 1,2,3 and 4 in list)

### Damaged version
dict_origin_avtime_dam = {}
dict_destination_avtime_dam_cumtime = {}
dict_destination_avtime_dam_demand = {}
dict_destination_avtime_dam = {}

for block in dict_block_node:
	dict_origin_avtime_dam[block] = 0
	dict_destination_avtime_dam[block] = 0
	dict_destination_avtime_dam_cumtime[block] = 0
	dict_destination_avtime_dam_demand[block] =0

# For origin census tract
total_demand_value = 0
for origin in counting_damage:
	total_demand_origin = 0
	time_cumulated = 0
	for destination in counting_damage[origin]:
		time_local = counting_damage[origin][destination][1]
		demand_local = counting_damage[origin][destination][0]
		total_demand_value += demand_local
		dest_block = dict_node_block[destination]
		if time_local < float('inf'):
			time_cumulated += time_local*demand_local
			total_demand_origin += demand_local
			dict_destination_avtime_dam_demand[dest_block] += demand_local
			dict_destination_avtime_dam_cumtime[dest_block] += time_local*demand_local

	if total_demand_origin >0:
		average_time = float(time_cumulated)/float(total_demand_origin)
	else:
		average_time = 0

	block_node = dict_node_block[origin]
	dict_origin_avtime_dam[block_node] = average_time


#List item Number 1
path_origin_ave = output_directory + 'origin_average_time'+st_scenario+'.pkl'
with open(path_origin_ave,'wb') as f:
	pickle.dump(dict_origin_avtime_dam, f)

# For destination census tract

for key in dict_destination_avtime_dam:
	if dict_destination_avtime_dam_demand[key] >0:
		dict_destination_avtime_dam[key] = float(dict_destination_avtime_dam_cumtime[key])/float(dict_destination_avtime_dam_demand[key])

#List item Number 2
path_destination_ave = output_directory+'destination_average_time'+st_scenario+'.pkl'
with open(path_destination_ave,'wb') as f:
	pickle.dump(dict_destination_avtime_dam, f)

### Undamaged Version

flag_ud = 1 #Flag for undamaged version

if flag_ud == 1:
	dict_origin_avtime_nd = {}
	dict_destination_avtime_nd_cumtime = {}
	dict_destination_avtime_nd_demand = {}
	dict_destination_avtime_nd = {}

	for block in dict_block_node:
		dict_origin_avtime_nd[block] = 0
		dict_destination_avtime_nd_cumtime[block] = 0
		dict_destination_avtime_nd_demand[block] = 0
		dict_destination_avtime_nd[block] = 0

	# For origin census tract

	for origin in counting:
		total_demand_origin = 0
		time_cumulated = 0
		for destination in counting[origin]:
			time_local = counting[origin][destination][1]
			demand_local = counting[origin][destination][0]
			dest_block = dict_node_block[destination]
			block_node = dict_node_block[origin]


			if time_local < float('inf'):
				time_cumulated += time_local*demand_local
				total_demand_origin += demand_local
				dict_destination_avtime_nd_demand[dest_block] += demand_local
				dict_destination_avtime_nd_cumtime[dest_block] += time_local

		if total_demand_origin >0:
			average_time = time_cumulated/total_demand_origin
		else:
			average_time = 0

		block_node = dict_node_block[origin]
		dict_origin_avtime_nd[block_node] = average_time

	#Saving variable number 3 in the list

	path_origin_ave_nd = output_directory+'origin_average_time_undamaged.pkl'
	with open(path_origin_ave_nd, 'wb') as f:
		pickle.dump(dict_origin_avtime_nd, f)

	# For destination census tract

	for key in dict_destination_avtime_nd:
		if dict_destination_avtime_nd_demand[key] > 0:
			dict_destination_avtime_nd[key] = dict_destination_avtime_nd_cumtime[key] /dict_destination_avtime_nd_demand[key]

	#Saving variable 4 in the list

	path_destination_ave_nd = output_directory+'destination_average_time_undamaged.pkl'
	with open(path_destination_ave_nd, 'wb') as f:
		pickle.dump(dict_destination_avtime_nd, f)

###### Lost trips #####

#Definition of threshold values
factors = 3600*float(upper_bound) # seconds
low_bound = 60*10 # seconds

dict_origin_block = {}
dict_destination_block = {}

if sf_specific == 'True':
	# GB MODIFICATION -- using regional values for social value of travel time, based on monthly incomes, then dividing
	# by 160 hours of work per month * 12 months per year to get hourly income
	#coef_low = 1 / (2 * 6.05 ** 0.26) # MIDPOINT OF INCOME RANGE
	#coef_medium = 1 / (2* 19.90 ** 0.26) # MIDPOINT OF INCOME RANGE
	#coef_high = 1 / (2* 157.71 ** 0.26) # MIDPOINT OF INCOME RANGE

	# GB MODIFICATION -- using regional values (medians) for social value of travel time, based on annual incomes,
	# then dividing by 160 hours of work per month * 12 months per year to get hourly income
	coef_low = 1 / (2 * 5.21 ** 0.26) # MEDIAN PERCENTILE OF RANGE
	coef_medium = 1 / (2* 19.27 ** 0.26) # MEDIAN PERCENTILE OF RANGE
	coef_high = 1 / (2* 57.3 ** 0.26) # MEDIAN PERCENTILE OF RANGE
else:
	#### Welfare loss ##### GB Modification to get units correct --  hourly income
	coef_low = 1 / (2 * 3.9 ** 0.26)
	coef_medium = 1 / (2 * 14.3 ** 0.26)
	coef_high = 1 / (2 * 52.1 ** 0.26)

dict_origin_welfare = {}
dict_destination_welfare = {}
dict_origin_welfare_low ={}
dict_origin_welfare_mid ={}
dict_origin_welfare_high = {}
dict_dest_welfare_low ={}
dict_dest_welfare_mid ={}
dict_dest_welfare_hig = {}


for block in dict_block_node:
	dict_origin_block[block] = 0
	dict_destination_block[block] = 0
	dict_origin_welfare[block] = 0
	dict_destination_welfare[block] = 0
	dict_origin_welfare_low[block] = 0
	dict_origin_welfare_mid[block] = 0
	dict_origin_welfare_high[block] = 0
	dict_dest_welfare_low[block] = 0
	dict_dest_welfare_mid[block] = 0
	dict_dest_welfare_hig[block] = 0

start = time.time()
factor = factors
for origin in counting:
	for destination in counting[origin]:
		original_time = counting[origin][destination][1]
		trips = counting[origin][destination][0] # trips demanded
		new_time = counting_damage[origin][destination][1]
		delta_time = new_time - original_time
		demand_deag = dict_demand_ita_deag[origin][destination]
		low_demand = demand_deag[0]
		mid_demand = demand_deag[1]
		hig_demand = demand_deag[2]
		welfare_loss_od = low_demand*coef_low*delta_time + mid_demand*coef_medium*delta_time + hig_demand*coef_high *delta_time
		welfare_loss_od_low = low_demand*coef_low*delta_time
		welfare_loss_od_mid = mid_demand*coef_medium*delta_time
		welfare_loss_od_high = hig_demand*coef_high *delta_time
		block_or = dict_node_block[origin]
		block_de = dict_node_block[destination]
		if new_time >= low_bound:  # Threshold for acceptable lower bound
			if new_time >= factor:  # Trips are lost
				if original_time < float('inf'):  # GB ADDITION -- otherwise we get nan values for welfare loss
					welfare_loss_od_low = (fl_upper_bound * 3600 - original_time) * low_demand * coef_low  # TODO -- change this
					# to reflect max. acceptable commute time - original commute time
					welfare_loss_od_mid = (fl_upper_bound * 3600 - original_time) * mid_demand * coef_medium
					welfare_loss_od_high = (fl_upper_bound * 3600 - original_time) * hig_demand * coef_high
				else:  # if trip could not be made even on undamaged network,
					# no welfare loss when trip cannot be made on undamaged network
					welfare_loss_od_low = 0
					welfare_loss_od_mid = 0
					welfare_loss_od_high = 0
				welfare_loss_od = welfare_loss_od_low + welfare_loss_od_mid + welfare_loss_od_high

				loss_before_or = dict_origin_block[block_or]
				loss_before_des = dict_destination_block[block_de]
				loss_new_or = loss_before_or + trips
				loss_new_de = loss_before_des + trips
				dict_origin_block[block_or] = loss_new_or
				dict_destination_block[block_de] = loss_new_de
		dict_dest_welfare_hig[block_de] += welfare_loss_od_high
		dict_dest_welfare_mid[block_de] += welfare_loss_od_mid
		dict_dest_welfare_low[block_de] += welfare_loss_od_low
		dict_origin_welfare_high[block_or] += welfare_loss_od_high
		dict_origin_welfare_mid[block_or] += welfare_loss_od_mid
		dict_origin_welfare_low[block_or] += welfare_loss_od_low
		dict_origin_welfare[block_or] += welfare_loss_od
		dict_destination_welfare[block_de] += welfare_loss_od

#Saving Variables 5 and 6 of the list

path_origin = output_directory+'dict_origin_loss_'+st_scenario+'.pkl'
path_destination = output_directory+'dict_destination_loss_'+st_scenario+'.pkl'

with open(path_origin,'wb') as f:
	pickle.dump(dict_origin_block, f)

with open(path_destination,'wb') as f:
	pickle.dump(dict_destination_block, f)

# Saving Variables 9 and 10

path_origin_wel = output_directory+'origin_welfare_'+st_scenario+'.pkl'
path_destination_wel = output_directory+'destination_welfare_'+st_scenario+'.pkl'
path_origin_wel_low = output_directory+'origin_welfare_low_'+st_scenario+'.pkl'
path_origin_wel_mid = output_directory+'origin_welfare_mid_'+st_scenario+'.pkl'
path_origin_wel_hig = output_directory+'origin_welfare_hig_'+st_scenario+'.pkl'

path_destination_wel_low = 'output_welfare_destination/'+'destination_welfare_low_'+st_scenario+'.pkl'
path_destination_wel_mid = 'output_welfare_destination/'+'destination_welfare_mid_'+st_scenario+'.pkl'
path_destination_wel_hig = 'output_welfare_destination/'+'destination_welfare_hig_'+st_scenario+'.pkl'

with open(path_origin_wel,'wb') as f:
	pickle.dump(dict_origin_welfare, f)

with open(path_origin_wel_low,'wb') as f:
	pickle.dump(dict_origin_welfare_low, f)

with open(path_origin_wel_mid,'wb') as f:
	pickle.dump(dict_origin_welfare_mid, f)

with open(path_origin_wel_hig,'wb') as f:
	pickle.dump(dict_origin_welfare_high, f)

with open(path_destination_wel,'wb') as f:
	pickle.dump(dict_destination_welfare, f)

with open(path_destination_wel_low,'wb') as f:
	pickle.dump(dict_dest_welfare_low, f)

with open(path_destination_wel_mid,'wb') as f:
	pickle.dump(dict_dest_welfare_mid, f)

with open(path_destination_wel_hig,'wb') as f:
	pickle.dump(dict_dest_welfare_hig, f)

##### Cost of road network disruption relative to undamaged network ##### GB ADDITION
dict_origin_cost = {}
dict_destination_cost = {} # cost of road network disruption is sum of cost of lost trips and cost of travel delays

for block in dict_block_node:
	dict_origin_cost[block] = 0
	dict_destination_cost[block] = 0

alpha = 48 # dollars per hour
beta = 78*8 # dollars per hour times hours

for origin in counting:
	for destination in counting[origin]:
		original_time = counting[origin][destination][1] # time from origin to destination on undamaged network
		trips = counting[origin][destination][0] # trips demanded on undamaged network
		new_time = counting_damage[origin][destination][1] # time from origin to destination on damaged n etwork
		new_trips = counting_damage[origin][destination][0] # trips demanded on damaged network
		if new_time < float('inf'):
			delta_time = max(0,(new_time - original_time)/3600)
		else:
			delta_time = 0 # all trips will be lost, so don't assign a cost for travel delays
			#print('new_time = ', new_time, ' and trips made on damaged network = ', new_trips)

		block_or = dict_node_block[origin]
		block_de = dict_node_block[destination]

		lost_trips_origin = dict_origin_block[block_or]
		lost_trips_destination = dict_destination_block[block_de]

		if lost_trips_origin < 0:
			#print('lost_trips_origin = ', lost_trips_origin)
			lost_trips_origin = 0

		cost_origin = alpha*delta_time*new_trips + beta*lost_trips_origin # hourly cost of road network disruption to origin
		cost_destination = alpha*delta_time*new_trips + beta*lost_trips_destination # hourly cost of road network disruption to destination
		if cost_origin < float('inf') and cost_destination < float('inf'):
			pass
		else:
			print(delta_time, lost_trips_origin, lost_trips_destination, cost_origin, cost_destination)
		dict_origin_cost[block_or] += cost_origin
		dict_destination_cost[block_de] += cost_destination

print('Done computing costs, and they are not infinite!')

path_origin = output_directory+'dict_origin_cost_'+st_scenario+'.pkl' # Cost of road network disruption
path_destination = output_directory+'dict_destination_cost_'+st_scenario+'.pkl' # Cost of road network disruption

with open(path_origin,'wb') as f:
	pickle.dump(dict_origin_cost, f)

with open(path_destination,'wb') as f:
	pickle.dump(dict_destination_cost, f)