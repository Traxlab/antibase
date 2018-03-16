import json
from pyteomics import mass
def createNodes(antiBase_map):
	nodes_list = []
	group = 1
	for key in antiBase_map.keys():
		temp_node = {}
		temp_node["id"] = key + " , " + str(mass.calculate_mass(formula=key))
		temp_node["group"] = group
		nodes_list.append(temp_node)
		for i in range(0, len(antiBase_map[key][1])):
			temp_node = {}
			temp_node["id"] = str(antiBase_map[key][1][i]) + " , " + str(antiBase_map[key][3][i]) + " , " +str(antiBase_map[key][0][i])
			temp_node["group"] = group
			nodes_list.append(temp_node)
		group += 1
	return nodes_list

def createLinks(antiBase_map):
	links_list = []
	for key in antiBase_map.keys():
		for j in range(0, len(antiBase_map[key][1])):
			temp_link = {}
			temp_link["source"] = key + " , " + str(mass.calculate_mass(formula=key))
			temp_link["target"] = str(antiBase_map[key][1][j]) + " , " + str(antiBase_map[key][3][j]) + " , " +str(antiBase_map[key][0][j])
			temp_link["value"] = 10
			links_list.append(temp_link)
	return links_list
def makeJson(antiBase_map):
	nodes_list = createNodes(antiBase_map)
	links_list = createLinks(antiBase_map)
	json_list = {}
	json_list["nodes"] = nodes_list
	json_list["links"] = links_list

	with open('data.json', 'w') as fp:
		json.dump(json_list,fp)
