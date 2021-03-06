#!/bin/python
from list_elements import *
from pymatgen import MPRester
from fractions import gcd

# General declarations
my_api_key = 'MpHqeT5fxhS8y4uF'
metals = list_metals#['V','Ru','Ir','Ti']
m = MPRester(my_api_key)
outfile = 'allelasticity.txt'
basic_properties = ['task_id','pretty_formula','reduced_cell_formula','unit_cell_formula','spacegroup.number','energy','formation_energy_per_atom','density']
elasticity_properties_prefix = ['elasticity']
all_properties = basic_properties + elasticity_properties_prefix
elasticity_properties = ['G_Reuss','G_VRH','G_Voigt','K_Reuss','K_VRH','K_Voigt','poisson_ratio','universal_anisotropy']
basic_properties_noucf = basic_properties
basic_properties_noucf.remove('reduced_cell_formula') # Do not include reduced_cell_formula in the text file (need it for composition only)
basic_properties_noucf.remove('unit_cell_formula') # Do not include unit_cell_formula in the text file (need it for energy scaling only)

# All pretty_formula and spacegroup.number dict
all_pf = [] # Used to check repeated entries (same pretty_formula and spacegroup.number)
all_sg = [] # Used to check repeated entries (same pretty_formula and spacegroup.number)

# Generate data file
with open(outfile,'w') as ofile:
	for prop in basic_properties_noucf:
		ofile.write(prop + '\t')
	for prop in elasticity_properties:
		ofile.write(prop + '\t')
	for element in list_allelements:
		ofile.write(element + '\t')
	ofile.write('\n')

	# Query options
	criteria = {"elasticity": {"$exists": True}, "spacegroup.number": {"$exists": True}} # All compounds with elasticity properties

	# All properties: query response
	results = m.query(criteria=criteria, properties=all_properties)

	# Write all properties values to file
	resultInd = 1
	for result in results:
		if result["spacegroup.number"] > 0:
			# Check if current entry exists in psfg_dict
			isDiffEntry = True
			for i in range(len(all_pf)):
				if ((result["pretty_formula"] == all_pf[i]) and (result["spacegroup.number"] == all_sg[i])):
					isDiffEntry = False
					break
			if isDiffEntry:
				# Basic properties
				unit_cell_formula = result['unit_cell_formula']
				elementsgcd = reduce(gcd, unit_cell_formula.values())
				for prop in basic_properties_noucf:
					if prop != "energy":
						ofile.write(str(result[prop]) + '\t')
					else:
						ofile.write(str(result[prop]/elementsgcd) + '\t')
				
				# Elasticity properties
				for prop in elasticity_properties:
					if (result['elasticity'] != None) and (prop in result['elasticity']):
						if result['elasticity'][prop] == None:
							ofile.write('NA\t')
						else:
							ofile.write(str(result['elasticity'][prop]) + '\t')
					else:
						ofile.write('NA\t')
			
				# Composition
				reduced_cell_formula = result['reduced_cell_formula']
				for element in list_allelements:
					if (element in reduced_cell_formula):
						ofile.write('{:.0f}'.format(reduced_cell_formula[element]) + '\t')
					else:
						ofile.write(str('0') + '\t')
				ofile.write('\n')

				# Update all pretty_formula and spacegroup.number entries
				all_pf.append(result["pretty_formula"])
				all_sg.append(result["spacegroup.number"])
				
				print "Finished writing " + result['pretty_formula'] + "(" + str(resultInd) + "/" + str(len(results)) + ")"
				resultInd += 1