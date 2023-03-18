#!/usr/bin/env python3

"""
The program is used to calculate the number of hydrogen bonded rings in the pure water system.\n

USAGE: ./nucleation_tracker.py -f filename.gro [-r 10 -c 0 -r 10 -m vertex -d]\n
       -f = input a gro file with correct box dimension
       -r = size of the rings; 9 or 10 membered closed rings; by default to 8 membered polygons
       -c = '0' for minimal ring counting and '1' for total non-self-intersecting rings; by default to minimal rings
       -m = algorithm for ring pruning [vertex, hbond, hbondAngle, torsion]; by default to hbondAngle
       -d = generates only directional rings tracking proton acceptor waters
       -e = turns on energy definition of hydrogen bonding; by default -2.0 kcal/mol hbond energy
       -x = Output an .xyz trajectory file containing the ring center locations and type (using atomic number of elements for size). (default=off)
"""

import sys
import getopt
import os
import os.path
from os import path
import shutil
import math
import numpy as np
import linecache

def usage():
	print(__doc__)

def beginCalc(inputFileName, max_ring, ring_closure, algorithm):
	distance_tol = 3.5    # hbond distance tolerance (in Ångström)
	angle_tol = 30        # hbond angle tolerance (in degrees)
	angle_tol_rad = angle_tol * 3.1415926536 / 180.0  # hbond angle tolerance in radians
	tolerance = 10e-5
	
	num = 0
	hbonds = 0
	isHBond = 0
	hbondAngle = 0
	hbondDist = 0
	boxLength = [None]*3
	posVec = [None]*3		# three elements 1D vector list
	scaledVec = [None]*3
	vec1 = [None]*3
	vec2 = [None]*3
	hmat = [ [0]*3 for i in range(3)]
	invhmat = [ [0]*3 for i in range(3)]
	line_count = 3
	ref_ndx = 0
	is_H1 = 0
	is_H2 = 0
	oPosX = []
	oPosY = []
	oPosZ = []
	hbondIndex = []
	threeAtomWater = False
	fourAtomWater = False
	
	if os.path.isfile("ringsCount.dat"):
		# shutil.copy("ringsCount.dat", "ringsCount2.dat")
		os.remove("ringsCount.dat")
	if os.path.isfile("rings_location.xyz"):
		# shutil.copy("rings_location.xyz", "rings_location2.xyz")	
		os.remove("rings_location.xyz")
	
	print("\nEnumerating the closed rings in the system . . . \n")
	outFile = open("ringsCount.dat", 'w')
	outFile.write("{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}".format("Frames", "hbonds", "rings_3", "rings_4", "rings_5", "rings_6", "rings_7", "rings_8"))
	if max_ring > 8:
		outFile.write("{:>10s}".format("rings_9"))
		if max_ring > 9:
			outFile.write("{:>10s}".format("rings_10"))
	outFile.write("\n")
	
	def energy_defn(pos1, pos2, i, j, isHBond, etol = -2.0):
		# NOT FOR DIRECTIONALITY!!!!
		
		# a function to identify if there is an HBond between 3 point water models
		# following an energy evaluation
		isHBond = 0
		
		posVec = [None]*3
		charge_vec1 = [None]*9
		charge_vec2 = [None]*9
		q_vec1 = [None]*3
		q_vec2 = [None]*3
		tmp_vec = [None]*3
		scaledVec = [None]*3
		
		mag = 0
		lj_pot = 0
		el_pot = 0
		m_offset = 0.15
		# 'etol': -2 kcal/mol(more range) show greater ring population than -4 kcal/mol; non-ideal less stable hbond has high energy(+)
		energy_tol = etol   # kcal/mol
		coulomb_const = 332.0636  # kcal Å/(mol e^2)
		q_val = 1.040
		e_val = 0.1550
		s_val = 3.15365
		e_val *= 4.0
		
		q_vec1[0] = -q_val
		q_vec1[1] = q_val*0.5
		q_vec1[2] = q_val*0.5
		q_vec2[0] = -q_val
		q_vec2[1] = q_val*0.5
		q_vec2[2] = q_val*0.5
		
		# TIP4P default: determine the massless site location
		# to build the 3 site charge vectors
		# for water 1;  pos1 = [Ox, Oy, Oz, H1x, H1y, H1z, H2x, H2y, H2z]
		for k in range(3):
			tmp_vec[k] = pos1[3+k] + pos1[6+k] - 2*pos1[k]
			mag += tmp_vec[k]*tmp_vec[k]

		mag = math.sqrt(mag)
		for k in range(3):
			tmp_vec[k] /= mag
			tmp_vec[k] *= m_offset
			charge_vec1[k] = pos1[k] + tmp_vec[k]
			charge_vec1[3+k] = pos1[3+k]
			charge_vec1[6+k] = pos1[6+k]

		mag = 0
		# for water 2
		for k in range(3):
			tmp_vec[k] = pos2[3+k] + pos2[6+k] - 2*pos2[k]
			mag += tmp_vec[k]*tmp_vec[k]

		mag = math.sqrt(mag)
		for k in range(3):
			tmp_vec[k] /= mag
			tmp_vec[k] *= m_offset
			charge_vec2[k] = pos2[k] + tmp_vec[k]
			charge_vec2[3+k] = pos2[3+k]
			charge_vec2[6+k] = pos2[6+k]

		# first, we do the LJ potential with the O-O vector
		posVec[0] = pos2[0]-pos1[0]
		posVec[1] = pos2[1]-pos1[1]
		posVec[2] = pos2[2]-pos1[2]
		
		# do vector wrapping of periodic boundary conditions
		for k in range(3):
			scaledVec[k] = posVec[k] * invhmat[k][k]
			scaledVec[k] -= round(scaledVec[k])
			posVec[k] = scaledVec[k] * hmat[k][k]

		
		r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2]
		ri = 1.0/math.sqrt(r2)
		ri6 = 1.0/(r2*r2*r2)
		s_val = pow(s_val,6)
		lj_pot = e_val*s_val*((s_val*ri6*ri6) - ri6)		# V_LJ(r) = 4*epsilon[(sigma/r)^12 - (sigma/r)^6]
		
		# now do the Coulombic potential double loop [water1 vs water2]
		for m in range(3):
			for l in range(3):
		        # do vector wrapping of periodic boundary conditions in the loop
				for k in range(3):
					posVec[k] = charge_vec2[3*m+k] - charge_vec1[3*l+k]
					scaledVec[k] = posVec[k] * invhmat[k][k]
					scaledVec[k] -= round(scaledVec[k])
					posVec[k] = scaledVec[k] * hmat[k][k]

				r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2]
				r = math.sqrt(r2)
				ri = 1.0 / r
				el_pot += coulomb_const*q_vec1[m]*q_vec2[l]*ri

		pot = lj_pot + el_pot

		if (pot < energy_tol):
			if sorted([i,j]) not in hbondIndex:
				hbondIndex.append(sorted([i,j]))
				isHBond = 1
			hbond_ndx[i].append(j)
			#bondvec[0] = pot
			#bondvec[1] = pot
			#bondvec[2] = pot

		return isHBond

	# Speedy et al. 'Network Topology in Simulated Water' criteria to take only 4 hbonds and prune additional hbonds if any for water hbonded to more than 4 neighbors
	def hbond_tetra(pos1, pos2, i, j, isHBond, etol = -2.0):
		# NOT FOR DIRECTIONALITY!!!!
		
		# a function to identify if there is an HBond between 3 point water models
		# following an energy evaluation
		isHBond = 0
		
		posVec = [None]*3
		charge_vec1 = [None]*9
		charge_vec2 = [None]*9
		q_vec1 = [None]*3
		q_vec2 = [None]*3
		tmp_vec = [None]*3
		scaledVec = [None]*3
		
		mag = 0
		lj_pot = 0
		el_pot = 0
		m_offset = 0.15
		energy_tol = etol   # kcal/mol
		coulomb_const = 332.0636  # kcal Å/(mol e^2)
		q_val = 1.040
		e_val = 0.1550
		s_val = 3.15365
		e_val *= 4.0
		
		q_vec1[0] = -q_val
		q_vec1[1] = q_val*0.5
		q_vec1[2] = q_val*0.5
		q_vec2[0] = -q_val
		q_vec2[1] = q_val*0.5
		q_vec2[2] = q_val*0.5
		
		# TIP4P default: determine the massless site location
		# to build the 3 site charge vectors
		# for water 1
		for k in range(3):
			tmp_vec[k] = pos1[3+k] + pos1[6+k] - 2*pos1[k]
			mag += tmp_vec[k]*tmp_vec[k]

		mag = math.sqrt(mag)
		for k in range(3):
			tmp_vec[k] /= mag
			tmp_vec[k] *= m_offset
			charge_vec1[k] = pos1[k] + tmp_vec[k]
			charge_vec1[3+k] = pos1[3+k]
			charge_vec1[6+k] = pos1[6+k]

		mag = 0
		# for water 2
		for k in range(3):
			tmp_vec[k] = pos2[3+k] + pos2[6+k] - 2*pos2[k]
			mag += tmp_vec[k]*tmp_vec[k]

		mag = math.sqrt(mag)
		for k in range(3):
			tmp_vec[k] /= mag
			tmp_vec[k] *= m_offset
			charge_vec2[k] = pos2[k] + tmp_vec[k]
			charge_vec2[3+k] = pos2[3+k]
			charge_vec2[6+k] = pos2[6+k]

		# first, we do the LJ potential with the O-O vector
		posVec[0] = pos2[0]-pos1[0]
		posVec[1] = pos2[1]-pos1[1]
		posVec[2] = pos2[2]-pos1[2]
		
		# do vector wrapping of periodic boundary conditions
		for k in range(3):
			scaledVec[k] = posVec[k] * invhmat[k][k]
			scaledVec[k] -= round(scaledVec[k])
			posVec[k] = scaledVec[k] * hmat[k][k]

		
		r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2]
		ri = 1.0/math.sqrt(r2)
		ri6 = 1.0/(r2*r2*r2)
		s_val = pow(s_val,6)
		s_val = pow(s_val,6)
		lj_pot = e_val*s_val*((s_val*ri6*ri6) - ri6)		# V_LJ(r) = 4*epsilon[(sigma/r)^12 - (sigma/r)^6]
		
		# now do the Coulombic potential double loop [water1 vs water2]
		for m in range(3):
			for l in range(3):
		        # do vector wrapping of periodic boundary conditions in the loop
				for k in range(3):
					posVec[k] = charge_vec2[3*m+k] - charge_vec1[3*l+k]
					scaledVec[k] = posVec[k] * invhmat[k][k]
					scaledVec[k] -= round(scaledVec[k])
					posVec[k] = scaledVec[k] * hmat[k][k]

				r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2]
				r = math.sqrt(r2)
				ri = 1.0 / r
				el_pot += coulomb_const*q_vec1[m]*q_vec2[l]*ri

		pot = lj_pot + el_pot

		# for comparing hbond energies
		hbond_pot.append(pot)

		return hbond_pot

	# Luzar and Chandler geometric definition of hydrogen bonding
	def IsHBonding(pos1, pos2, i, j, isHBond):
		# first, we do the O-O vector
		posVec[0] = pos2[0]-pos1[0]
		posVec[1] = pos2[1]-pos1[1]
		posVec[2] = pos2[2]-pos1[2]
	
		# do vector wrapping of periodic boundary conditions
		for k in range(3):
			scaledVec[k] = posVec[k] * invhmat[k][k]
			scaledVec[k] -= round(scaledVec[k])
			posVec[k] = scaledVec[k] * hmat[k][k]
	
		# here, we take care of normalization
		# r2 = xVal*xVal + yVal*yVal + zVal*zVal
		r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2]
		r = math.sqrt(r2)
		dist_OO = r
		if dist_OO < distance_tol:		# Luzar & Chandler distance criteria
			# convert to unit vector for finding dot product
			ri = 1.0 / r
			vec1[0] = posVec[0] * ri
			vec1[1] = posVec[1] * ri
			vec1[2] = posVec[2] * ri
			
			# the O-H1
			posVec[0] = pos1[3]-pos1[0]
			posVec[1] = pos1[4]-pos1[1]
			posVec[2] = pos1[5]-pos1[2]
			
			# do vector wrapping of periodic boundary conditions
			for k in range(3):
				scaledVec[k] = posVec[k] * invhmat[k][k]
				scaledVec[k] -= round(scaledVec[k])
				posVec[k] = scaledVec[k] * hmat[k][k]
	
			# here, we take care of normalization
			# r2 = xVal*xVal + yVal*yVal + zVal*zVal
			r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2]
			r = math.sqrt(r2)
			ri = 1.0 / r
			
			vec2[0] = posVec[0] * ri
			vec2[1] = posVec[1] * ri
			vec2[2] = posVec[2] * ri
			
			dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]))
			if (dotProduct > 1.0):
				tempAngle = 0.0
			else:
				tempAngle = math.acos(dotProduct)
			
			# Luzar & Chandler hbond cone angle criteria
			if (tempAngle <= angle_tol_rad):
				if sorted([i,j]) not in hbondIndex:			# to remove hbond repeatation
					hbondIndex.append(sorted([i,j]))
					isHBond = 1
				if not directionality:
					hbond_ndx[i].append(j)
			else:
				# continue down the rabbit hole
				# the O-H2
				posVec[0] = pos1[6]-pos1[0]
				posVec[1] = pos1[7]-pos1[1]
				posVec[2] = pos1[8]-pos1[2]
				
				# do vector wrapping of periodic boundary conditions
				for k in range(3):
					scaledVec[k] = posVec[k] * invhmat[k][k]
					scaledVec[k] -= round(scaledVec[k])
					posVec[k] = scaledVec[k] * hmat[k][k]
	
				# here, we take care of normalization
				# r2 = xVal*xVal + yVal*yVal + zVal*zVal
				r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2]
				r = math.sqrt(r2)
				ri = 1.0 / r
				
				vec2[0] = posVec[0] * ri
				vec2[1] = posVec[1] * ri
				vec2[2] = posVec[2] * ri
			
				dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]))
				if (dotProduct > 1.0):
					tempAngle = 0.0
				else:
					tempAngle = math.acos(dotProduct)
			
				if (tempAngle <= angle_tol_rad):
					if sorted([i,j]) not in hbondIndex:
						hbondIndex.append(sorted([i,j]))
						isHBond = 1
					if not directionality:
						hbond_ndx[i].append(j)
				else:
					# continue further down the rabbit hole
					# considering only acceptors
					
					# flip the O-O vector
					vec1[0] *= -1
					vec1[1] *= -1
					vec1[2] *= -1
					
					# the O-H1
					posVec[0] = pos2[3]-pos2[0]
					posVec[1] = pos2[4]-pos2[1]
					posVec[2] = pos2[5]-pos2[2]
					
					# do vector wrapping of periodic boundary conditions
					for k in range(3):
						scaledVec[k] = posVec[k] * invhmat[k][k]
						scaledVec[k] -= round(scaledVec[k])
						posVec[k] = scaledVec[k] * hmat[k][k]
	
					# here, we take care of normalization
					# r2 = xVal*xVal + yVal*yVal + zVal*zVal
					r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2]
					r = math.sqrt(r2)
					ri = 1.0 / r
					
					vec2[0] = posVec[0] * ri
					vec2[1] = posVec[1] * ri
					vec2[2] = posVec[2] * ri
					
					dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]))
					if (dotProduct > 1.0):
						tempAngle = 0.0
					else:
						tempAngle = math.acos(dotProduct)
					
					if (tempAngle <= angle_tol_rad):
						if sorted([i,j]) not in hbondIndex:
							hbondIndex.append(sorted([i,j]))
							isHBond = 1
						hbond_ndx[i].append(j)
					else:
						# last level... I promise...
						
						# the O-H2
						posVec[0] = pos2[6]-pos2[0]
						posVec[1] = pos2[7]-pos2[1]
						posVec[2] = pos2[8]-pos2[2]
						
						# do vector wrapping of periodic boundary conditions
						for k in range(3):
							scaledVec[k] = posVec[k] * invhmat[k][k]
							scaledVec[k] -= round(scaledVec[k])
							posVec[k] = scaledVec[k] * hmat[k][k]
	
						# here, we take care of normalization
						# r2 = xVal*xVal + yVal*yVal + zVal*zVal
						r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2]
						r = math.sqrt(r2)
						ri = 1.0 / r
						
						vec2[0] = posVec[0] * ri
						vec2[1] = posVec[1] * ri
						vec2[2] = posVec[2] * ri
						
						dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]))
						if (dotProduct > 1.0):
							tempAngle = 0.0
						else:
							tempAngle = math.acos(dotProduct)
						
						if (tempAngle <= angle_tol_rad):
							if sorted([i,j]) not in hbondIndex:
								hbondIndex.append(sorted([i,j]))
								isHBond = 1
							hbond_ndx[i].append(j)
		return isHBond
	
	pos = 0
	frame_count = 0
	with open (inputFileName, 'r') as inFile:
		nAtoms = int(linecache.getline(inputFileName, 2))
		dummy_atom = linecache.getline(inputFileName, 6)
		if str(dummy_atom[10:15]).strip() == "MW":
			fourAtomWater = True
			nMols = int(nAtoms/4)
		else:
			threeAtomWater = True
			nMols = nAtoms/3
	
		# the 9 positio:ns are sorted by O(x,y,z), H1(x,y,z), H2(x,y,z)
		water = [[0 for _ in range(9)] for _ in range(int(nMols))]	# '-1' because '0' is also part of index
		# assume a water can for hbonds with 10 neighbors in worst case senario in liquid phase
		# '-1' because '0' is also part of index
		hbond_ndx = [[] for _ in range(int(nMols))]
	
		next(inFile)
		next(inFile)
		# grab the water atom positions
		for i, line in enumerate(inFile):
			if str(line[10:15]).strip() in ["OW", "O"]:
				water[pos][0] = float(line[20:28])*10 
				water[pos][1] = float(line[28:36])*10
				water[pos][2] = float(line[36:44])*10
				oPosX.append(float(line[20:28])*10)
				oPosY.append(float(line[28:36])*10)
				oPosZ.append(float(line[36:44])*10)
				is_H1 = 1
		
			elif bool(is_H1):
				water[pos][3] = float(line[20:28])*10 
				water[pos][4] = float(line[28:36])*10
				water[pos][5] = float(line[36:44])*10
				is_H1 = 0
				is_H2 = 1

			elif bool(is_H2):
				water[pos][6] = float(line[20:28])*10 
				water[pos][7] = float(line[28:36])*10
				water[pos][8] = float(line[36:44])*10
				is_H2 = 0
				pos += 1
				
			elif (i-(frame_count*3))%nAtoms == 0:
				#Account for periodic boundary conditions.
				# the box information - used in minimum image wrapping
				line_split = line.split()
				# assume the box to be cubic
				Lx = float(line_split[0])*10.0
				Ly = float(line_split[1])*10.0
				Lz = float(line_split[2])*10.0
		
				hmat[0][0] = Lx
				hmat[1][1] = Ly
				hmat[2][2] = Lz
		
				# load the rest of the box cell matrix (hmat) using gromacs form.
				if len(line_split) > 3:
					hmat[0][1] = float(line_split[3])*10.0
					hmat[0][2] = float(line_split[4])*10.0
					hmat[1][0] = float(line_split[5])*10.0
					hmat[1][2] = float(line_split[6])*10.0
					hmat[2][0] = float(line_split[7])*10.0
					hmat[2][1] = float(line_split[8])*10.0
				else:
					# for cubic system
					hmat[0][1] = 0
					hmat[0][2] = 0
					hmat[1][0] = 0
					hmat[1][2] = 0
					hmat[2][0] = 0
					hmat[2][1] = 0
		
				# let's determine the inverse of the hmat
				# first the determinant...
				determinant = 0
				for i in range(3):
				    determinant = determinant + (hmat[0][i] * (hmat[1][(i+1)%3] * hmat[2][(i+2)%3] - hmat[1][(i+2)%3] * hmat[2][(i+1)%3]))
				    invdeterminant = 1.0 / determinant
				    # now the double loop inverse...
				    for i in range(3):
				        for j in range(3):
				            invhmat[i][j] = ((hmat[(j+1)%3][(i+1)%3] * hmat[(j+2)%3][(i+2)%3]) - (hmat[(j+1)%3][(i+2)%3] * hmat[(j+2)%3][(i+1)%3])) * invdeterminant
	
				halfDistance_tol_array = np.array([distance_tol/2 for i in range(int(nMols))])
				dx = np.subtract.outer(oPosX,oPosX)
				dy = np.subtract.outer(oPosY,oPosY)
				dz = np.subtract.outer(oPosZ,oPosZ)
				distance_tol_array = np.add.outer(halfDistance_tol_array,halfDistance_tol_array)
				
				# do vector wrapping of periodic boundary conditions
				def pbc(oVec, k):
					scaledVec = oVec * invhmat[k][k]
					scaledVec -= np.round(scaledVec)
					oVec = scaledVec * hmat[k][k]
					return oVec
				dx = pbc(dx, 0)
				dy = pbc(dy, 1)
				dz = pbc(dz, 2)
				distance = np.array( np.sqrt( np.square(dx) + np.square(dy) + np.square(dz) ) )
				
				# to find the number of surrounding water molecules around central water molecule
				touching = np.where( distance < distance_tol_array, 1.0, 0.0) #elementwise conditional operator (condition ? 1 : 0 in this case)
				
				# building up the water cluster with proper identification
				hbondListIndex = [[] for i in range(int(nMols))]
				for i in range(int(nMols)):
				    for j in range(int(nMols)):
				        if i==j: continue
				        if (abs(touching[i,j]*distance_tol_array[i,j] - distance_tol_array[i,j])) < tolerance:
				            hbondListIndex[i].append(j)
				        else: continue
			
				# energy criteria of hbonding
				if hbond_energy:
					# energy value default to -2kcal/mol
					# etol = -2.0			# kcal/mol; energy criteria if require other that -2.0 kcal/mol (-4kcal/mol has higher hbond energy than -2kcal/mol)
					# this is to loop over water environment to test H-bonding
					for i, hbondWat in enumerate(hbondListIndex):
						# this is a reference water
						pos1 = water[i] 
						for j in hbondWat:
							pos2 = water[j]
							isHBond = energy_defn(pos1, pos2, i, j, isHBond=0)
							hbonds += isHBond
					
					# Speedy et al. criteria to take only 4 hbonds and prune additional hbonds if any
					hbond_pot = []
					hbond_trim = 0
					for i, neighbors in enumerate(hbond_ndx):
						if len(neighbors) > 4:
							pos1 = water[i] 
							for j in neighbors:
								pos2 = water[j]
								hbond_tetra(pos1, pos2, i, j, isHBond=0)
							for k in range(len(neighbors)-4):
								extra_hbond_ndx = hbond_pot.index(max(hbond_pot))
								hbond_pot.pop(extra_hbond_ndx)
								hbond_ndx[i].pop(extra_hbond_ndx)
								hbond_trim += 1
							hbond_pot = []
					hbonds -= hbond_trim

				# calling hbond geometry for ring connectivity
				else:
					# this is to loop over water environment to test H-bonding
					for i, hbondWat in enumerate(hbondListIndex):
						# this is a reference water
						pos1 = water[i] 
						for j in hbondWat:
							pos2 = water[j]
							isHBond = IsHBonding(pos1, pos2, i, j, isHBond=0)
							hbonds += isHBond

				# we test the rings if it crosses the pbc; rings should come back from the same wall
				# This check ensures that a chain like |-A-B-C|-A-B-C|-A- is not counted as a polygon 'A-B-C'
				def ring_pbc(member1, member2, outRing_x, outRing_y, outRing_z):
					xVecRing = water[member1][0] - water[member2][0]
					yVecRing = water[member1][1] - water[member2][1]
					zVecRing = water[member1][2] - water[member2][2]
					if abs(xVecRing) > 0.5 * Lx:
						outRing_x += 1
					if abs(yVecRing) > 0.5 * Ly:
						outRing_y += 1
					if abs(zVecRing) > 0.5 * Lz:
						outRing_z += 1
					# if 'outRing%2 == 0' is even then consider a ring else not a ring
					return outRing_x, outRing_y, outRing_z

				# now counting all possible non short-circuit rings
				# test for hbond donor in the hbond_ndx
				wat_loop = []
				count_x = 0
				count_y = 0
				count_z = 0
				for i in range(int(nMols)):
					temp_loop = []
					for j, element in enumerate(hbond_ndx):
						if i in element:
							for k, element2 in enumerate(hbond_ndx):
								if k == i:				# skip already counted rings; saves computational time
									continue
								elif j in element2:
									if k in hbond_ndx[i]:
										temp_loop.extend([i,j,k])
										# we do not want double counting of rings . . .
										wat_sort = sorted(list(set(temp_loop)))
										if (wat_sort not in wat_loop):
											# for forming rings, there should be equal outgoing and incoming hbonds
											for var_a in range(len(temp_loop)-1):
												tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[var_a], temp_loop[var_a+1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
												count_x += tempCount_x
												count_y += tempCount_y
												count_z += tempCount_z
											# test for initial and final members
											tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[0], temp_loop[-1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
											count_x += tempCount_x
											count_y += tempCount_y
											count_z += tempCount_z
											is_ring = count_x + count_y + count_z
											if is_ring%2 == 0:
												wat_loop.append(wat_sort)
											count_x = 0
											count_y = 0
											count_z = 0
											temp_loop = []
										temp_loop = []
									for l, element3 in enumerate(hbond_ndx):
										if l == i or l ==j:
											continue
										elif k in element3:
											if l in hbond_ndx[i]:
												temp_loop.extend([i,j,k,l])
												wat_sort = sorted(list(set(temp_loop)))
												if (wat_sort not in wat_loop):
													# we take 'temp_loop' for avoiding a test for hbonds
													# test for ring pbc, there should be equal outgoing and incoming hbonds
													for var_a in range(len(temp_loop)-1):
														tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[var_a], temp_loop[var_a+1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
														count_x += tempCount_x
														count_y += tempCount_y
														count_z += tempCount_z
													# test for initial and final members
													tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[0], temp_loop[-1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
													count_x += tempCount_x
													count_y += tempCount_y
													count_z += tempCount_z
													is_ring = count_x + count_y + count_z
													if is_ring%2 == 0:
														wat_loop.append(wat_sort)
													count_x = 0
													count_y = 0
													count_z = 0
													temp_loop = []
												temp_loop = []
											for m, element4 in enumerate(hbond_ndx):
												if m == i or m == j or m == k:
													continue
												elif l in element4:
													if m in hbond_ndx[i]:
														temp_loop.extend([i,j,k,l,m])
														wat_sort = sorted(list(set(temp_loop)))
														if (wat_sort not in wat_loop):
															# we take 'temp_loop' for avoiding a test for hbonds
															# test for ring pbc, there should be equal outgoing and incoming hbonds
															for var_a in range(len(temp_loop)-1):
																tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[var_a], temp_loop[var_a+1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																count_x += tempCount_x
																count_y += tempCount_y
																count_z += tempCount_z
															# test for initial and final members
															tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[0], temp_loop[-1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
															count_x += tempCount_x
															count_y += tempCount_y
															count_z += tempCount_z
															is_ring = count_x + count_y + count_z
															if is_ring%2 == 0:
																wat_loop.append(wat_sort)
															count_x = 0
															count_y = 0
															count_z = 0
															temp_loop = []
														temp_loop = []
													for n, element5 in enumerate(hbond_ndx):
														if n == i or n == j or n == k or n == l:
															continue
														elif m in element5:	
															if n in hbond_ndx[i]:
																temp_loop.extend([i,j,k,l,m,n])
																wat_sort = sorted(list(set(temp_loop)))
																if (wat_sort not in wat_loop):
																	# we take 'temp_loop' for avoiding a test for hbonds
																	# test for ring pbc, there should be equal outgoing and incoming hbonds
																	for var_a in range(len(temp_loop)-1):
																		tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[var_a], temp_loop[var_a+1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																		count_x += tempCount_x
																		count_y += tempCount_y
																		count_z += tempCount_z
																	# test for initial and final members
																	tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[0], temp_loop[-1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																	count_x += tempCount_x
																	count_y += tempCount_y
																	count_z += tempCount_z
																	is_ring = count_x + count_y + count_z
																	if is_ring%2 == 0:
																		wat_loop.append(wat_sort)
																	count_x = 0
																	count_y = 0
																	count_z = 0
																	temp_loop = []
																temp_loop = []
															for o, element6 in enumerate(hbond_ndx):
																if o == i or o == j or o == k or o == l or o == m:
																	continue
																elif n in element6:	
																	if o in hbond_ndx[i]:
																		temp_loop.extend([i,j,k,l,m,n,o])
																		wat_sort = sorted(list(set(temp_loop)))
																		if (wat_sort not in wat_loop):
																			# we take 'temp_loop' for avoiding a test for hbonds
																			# test for ring pbc, there should be equal outgoing and incoming hbonds
																			for var_a in range(len(temp_loop)-1):
																				tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[var_a], temp_loop[var_a+1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																				count_x += tempCount_x
																				count_y += tempCount_y
																				count_z += tempCount_z
																			# test for initial and final members
																			tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[0], temp_loop[-1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																			count_x += tempCount_x
																			count_y += tempCount_y
																			count_z += tempCount_z
																			is_ring = count_x + count_y + count_z
																			if is_ring%2 == 0:
																				wat_loop.append(wat_sort)
																			count_x = 0
																			count_y = 0
																			count_z = 0
																			temp_loop = []
																		temp_loop = []
																	for p, element7 in enumerate(hbond_ndx):
																		if p == i or p == j or p == k or p == l or p == m or p == n:
																			continue
																		elif o in element7:	
																			if p in hbond_ndx[i]:
																				temp_loop.extend([i,j,k,l,m,n,o,p])
																				wat_sort = sorted(list(set(temp_loop)))
																				if (wat_sort not in wat_loop):
																					# we take 'temp_loop' for avoiding a test for hbonds
																					# test for ring pbc, there should be equal outgoing and incoming hbonds
																					for var_a in range(len(temp_loop)-1):
																						tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[var_a], temp_loop[var_a+1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																						count_x += tempCount_x
																						count_y += tempCount_y
																						count_z += tempCount_z
																					# test for initial and final members
																					tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[0], temp_loop[-1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																					count_x += tempCount_x
																					count_y += tempCount_y
																					count_z += tempCount_z
																					is_ring = count_x + count_y + count_z
																					if is_ring%2 == 0:
																						wat_loop.append(wat_sort)
																					count_x = 0
																					count_y = 0
																					count_z = 0
																					temp_loop = []
																				temp_loop = []
																			if max_ring > 8:
																				for q, element8 in enumerate(hbond_ndx):
																					if q == i or q == j or q == k or q == l or q == m or q == n or q == o:
																						continue
																					elif p in element8:	
																						if q in hbond_ndx[i]:
																							temp_loop.extend([i,j,k,l,m,n,o,p,q])
																							wat_sort = sorted(list(set(temp_loop)))
																							if (wat_sort not in wat_loop):
																								# we take 'temp_loop' for avoiding a test for hbonds
																								# test for ring pbc, there should be equal outgoing and incoming hbonds
																								for var_a in range(len(temp_loop)-1):
																									tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[var_a], temp_loop[var_a+1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																									count_x += tempCount_x
																									count_y += tempCount_y
																									count_z += tempCount_z
																								# test for initial and final members
																								tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[0], temp_loop[-1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																								count_x += tempCount_x
																								count_y += tempCount_y
																								count_z += tempCount_z
																								is_ring = count_x + count_y + count_z
																								if is_ring%2 == 0:
																									wat_loop.append(wat_sort)
																								count_x = 0
																								count_y = 0
																								count_z = 0
																								temp_loop = []
																							temp_loop = []
																						if max_ring > 9:
																							for r, element9 in enumerate(hbond_ndx):
																								if r == i or r == j or r == k or r == l or r == m or r == n or r == o or r == p:
																									continue
																								elif q in element9:	
																									if r in hbond_ndx[i]:
																										temp_loop.extend([i,j,k,l,m,n,o,p,q,r])
																										wat_sort = sorted(list(set(temp_loop)))
																										if (wat_sort not in wat_loop):
																											# we take 'temp_loop' for avoiding a test for hbonds
																											# test for ring pbc, there should be equal outgoing and incoming hbonds
																											for var_a in range(len(temp_loop)-1):
																												tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[var_a], temp_loop[var_a+1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																												count_x += tempCount_x
																												count_y += tempCount_y
																												count_z += tempCount_z
																											# test for initial and final members
																											tempCount_x, tempCount_y, tempCount_z = ring_pbc(temp_loop[0], temp_loop[-1], outRing_x = 0, outRing_y = 0, outRing_z = 0)
																											count_x += tempCount_x
																											count_y += tempCount_y
																											count_z += tempCount_z
																											is_ring = count_x + count_y + count_z
																											if is_ring%2 == 0:
																												wat_loop.append(wat_sort)
																											count_x = 0
																											count_y = 0
																											count_z = 0
																											temp_loop = []
																										temp_loop = []
				# rings segregation
				ring3 = []
				ring4 = []
				ring5 = []
				ring6 = []
				ring7 = []
				ring8 = []
				ring9 = []
				ring10 = []

				for i in wat_loop:
					if len(i) == 3:
						ring3.append(i)
					elif len(i) ==4:
						ring4.append(i)
					elif len(i) ==5:
						ring5.append(i)
					elif len(i) ==6:
						ring6.append(i)
					elif len(i) ==7:
						ring7.append(i)
					elif len(i) ==8:
						ring8.append(i)
					elif len(i) ==9:
						ring9.append(i)
					elif len(i) ==10:
						ring10.append(i)

				# pruning down larger rings if smaller closed rings are available
				# now we do primitive/minimal ring counting
				if ring_closure == 0:
					primitiveRings = []
					minimal_rings = []
					# test if two molecules are hbonded; algorithm of 'common hbond' for pruning larger rings
					def common_edge(larger_ring, small_ring):
						common_elem = list(set(small_ring).intersection(set(larger_ring)))
						for i in common_elem:			# looping over the same list/water is fine because we do not see water index in its own vector
							for j in common_elem:
								# test for hbond formation
								if i in hbond_ndx[j]:
									minimal_rings.append(larger_ring)
									return minimal_rings

					# test if three molecules are bonded but not necessarily forming closed ring; algorithm of 'common angle' for pruning larger rings
					def common_angle(larger_ring, small_ring):
						common_elem = list(set(small_ring).intersection(set(larger_ring)))
						for i in common_elem:
							for j in common_elem:
								if i in hbond_ndx[j]:
									for k in common_elem:
										if k == i:			# 'i' is already included on the path
											continue
										# forms three molecules connected by hbond
										elif j in hbond_ndx[k]:
											minimal_rings.append(larger_ring)
											return minimal_rings

					# test if four molecules are bonded to form torsion but not necessarily forming closed ring
					# algorithm of 'common torsion' for pruning larger rings
					def common_torsion(larger_ring, small_ring):
						common_elem = list(set(small_ring).intersection(set(larger_ring)))
						for i in common_elem:
							for j in common_elem:
								if i in hbond_ndx[j]:
									for k in common_elem:
										if k == i:			# 'i' is already included on the path
											continue
										elif j in hbond_ndx[k]:
											for l in common_elem:
												if l == i or l == j:		# 'i' and 'j' are already included on the path
													continue
												# forms three molecules connected by hbond
												elif k in hbond_ndx[l]:
													minimal_rings.append(larger_ring)
													return minimal_rings

					if algorithm == "hbond":
						# 3 membered rings
						for ring_size in [ring4, ring5, ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring3:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 2:		# test if there are 2 or more common elements
										common_edge(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []
									
						# 4 membered rings
						for ring_size in [ring5, ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring4:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 2:
										common_edge(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 5 membered rings
						for ring_size in [ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring5:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 2:
										common_edge(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 6 membered rings
						for ring_size in [ring7, ring8, ring9, ring10]:
							for small_ring in ring6:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 2:
										common_edge(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 7 membered rings
						for ring_size in [ring8, ring9, ring10]:
							for small_ring in ring7:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 2:
										common_edge(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 8 membered rings
						for ring_size in [ring9, ring10]:
							for small_ring in ring8:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 2:
										common_edge(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 9 membered rings
						for small_ring in ring9:
							for larger_ring in ring10:
								if len(set(small_ring).intersection(set(larger_ring))) >= 2:
									common_edge(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

					elif algorithm == "torsion":
						# 3 membered rings lack torsion bonds

						# 4 membered rings
						for ring_size in [ring5, ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring4:
								for larger_ring in ring_size:
									if set(small_ring).issubset(larger_ring):		# test if 'small_ring' belongs to 'larger_ring'
										primitiveRings.append(larger_ring)
									
						# 5 membered rings
						for ring_size in [ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring5:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 4:		# test if there are 4 or more common elements
										common_torsion(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 6 membered rings
						for ring_size in [ring7, ring8, ring9, ring10]:
							for small_ring in ring6:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 4:
										common_torsion(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 7 membered rings
						for ring_size in [ring8, ring9, ring10]:
							for small_ring in ring7:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 4:
										common_torsion(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 8 membered rings
						for ring_size in [ring9, ring10]:
							for small_ring in ring8:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 4:
										common_torsion(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 9 membered rings
						for small_ring in ring9:
							for larger_ring in ring10:
								if len(set(small_ring).intersection(set(larger_ring))) >= 4:
									common_torsion(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

					# by default to bond angle for pruning larger rings; algorithm = "angle"
					elif algorithm == "hbondAngle":
						# 3 membered rings
						for ring_size in [ring4, ring5, ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring3:
								for larger_ring in ring_size:
									if set(small_ring).issubset(larger_ring):		# test if 'small_ring' belongs to 'larger_ring'
										primitiveRings.append(larger_ring)
									
						# 4 membered rings
						for ring_size in [ring5, ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring4:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 3:		# test if there are 3 or more common elements
										common_angle(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 5 membered rings
						for ring_size in [ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring5:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 3:
										common_angle(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 6 membered rings
						for ring_size in [ring7, ring8, ring9, ring10]:
							for small_ring in ring6:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 3:
										common_angle(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 7 membered rings
						for ring_size in [ring8, ring9, ring10]:
							for small_ring in ring7:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 3:
										common_angle(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 8 membered rings
						for ring_size in [ring9, ring10]:
							for small_ring in ring8:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 3:
										common_angle(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

						# 9 membered rings
						for small_ring in ring9:
							for larger_ring in ring10:
								if len(set(small_ring).intersection(set(larger_ring))) >= 3:
									common_angle(larger_ring, small_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []

					# test if single water molecule are found common vertex between two rings
					elif algorithm == "vertex":
						# 3 membered rings 
						for ring_size in [ring4, ring5, ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring3:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 1:		# test if there are 1 or more common elements
										minimal_rings.append(larger_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []
						primitiveRings = set(tuple(element) for element in primitiveRings)
						primitiveRings = list(map(list, primitiveRings))

						# 4 membered rings 
						for ring_size in [ring5, ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring4:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 1:	
										minimal_rings.append(larger_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []
						primitiveRings = set(tuple(element) for element in primitiveRings)
						primitiveRings = list(map(list, primitiveRings))

						# 5 membered rings 
						for ring_size in [ring6, ring7, ring8, ring9, ring10]:
							for small_ring in ring5:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 1:	
										minimal_rings.append(larger_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []
						primitiveRings = set(tuple(element) for element in primitiveRings)
						primitiveRings = list(map(list, primitiveRings))

						# 6 membered rings 
						for ring_size in [ring7, ring8, ring9, ring10]:
							for small_ring in ring6:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 1:	
										minimal_rings.append(larger_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []
						primitiveRings = set(tuple(element) for element in primitiveRings)
						primitiveRings = list(map(list, primitiveRings))

						# 7 membered rings 
						for ring_size in [ring8, ring9, ring10]:
							for small_ring in ring7:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 1:	
										minimal_rings.append(larger_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []
						primitiveRings = set(tuple(element) for element in primitiveRings)
						primitiveRings = list(map(list, primitiveRings))

						# 8 membered rings 
						for ring_size in [ring9, ring10]:
							for small_ring in ring8:
								for larger_ring in ring_size:
									if len(set(small_ring).intersection(set(larger_ring))) >= 1:	
										minimal_rings.append(larger_ring)
							if minimal_rings:
								primitiveRings.extend(minimal_rings)
							minimal_rings = []
						primitiveRings = set(tuple(element) for element in primitiveRings)
						primitiveRings = list(map(list, primitiveRings))

						# 9 membered rings 
						for small_ring in ring3:
							for larger_ring in ring10:
								if len(set(small_ring).intersection(set(larger_ring))) >= 1:	
									minimal_rings.append(larger_ring)
						if minimal_rings:
							primitiveRings.extend(minimal_rings)
						minimal_rings = []
						primitiveRings = set(tuple(element) for element in primitiveRings)
						primitiveRings = list(map(list, primitiveRings))

					else:
						usage()
						print("Non-existing algorithm")
						sys.exit()

					primitiveRings = set(tuple(element) for element in primitiveRings)		# removes repeatation
					primitiveRings = list(map(list, primitiveRings))						# converts 2d tuple to list

					# now we take only the minimal paths for polygons
					for pruneRings in primitiveRings:
						wat_loop.remove(pruneRings)

					# reinitialie closed rings to add primitive rings only
					ring3 = []
					ring4 = []
					ring5 = []
					ring6 = []
					ring7 = []
					ring8 = []
					ring9 = []
					ring10 = []

					# rings segregation
					for i in wat_loop:
						if len(i) == 3:
							ring3.append(i)
						elif len(i) ==4:
							ring4.append(i)
						elif len(i) ==5:
							ring5.append(i)
						elif len(i) ==6:
							ring6.append(i)
						elif len(i) ==7:
							ring7.append(i)
						elif len(i) ==8:
							ring8.append(i)
						elif len(i) ==9:
							ring9.append(i)
						elif len(i) ==10:
							ring10.append(i)

				# now writing the number of enumerated closed rings
				frame_count +=1
				outFile.write("{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}".format(frame_count, hbonds, len(ring3), len(ring4), len(ring5), len(ring6), len(ring7), len(ring8)))
				if max_ring > 8:
					outFile.write("{:>10d}".format(len(ring9)))
					if max_ring > 9:
						outFile.write("{:>10d}\n".format(len(ring10)))
					else:
						outFile.write("\n")
				else:
					outFile.write("\n")
	
				if ring_traj:
		 			# writing xyz file for the location of directional rings
					# wrt to oxygen (heavy atom), if take hydrogen atoms, then divide by total items taken
					outFile2 = open("rings_location.xyz", 'w')
					outFile2.write("%d\n" %(len(wat_loop)))
					outFile2.write("ring location in the system\n")
	
					# minimum image wrapping and finding center of the ring
					if Lx != Ly and Ly != Lz:
						print("Note: Not cubic! please check the system size for orthogonal lattice/edge vector")
					def min_image_wrap(ring):
						# we directly find the average position of the ring
						tempX_pbc = 0
						tempY_pbc = 0
						tempZ_pbc = 0
						# "water[ring[0]][0]" is the first member of the rings and taken as the reference molecule for "min_image_wrap"
						tempX_pbc += water[ring[0]][0]
						tempY_pbc += water[ring[0]][1]
						tempZ_pbc += water[ring[0]][2]
						# now we loop over all remaining waters and determine the pbc of other molecule with reference to the first molecule
						for i in ring[1:]:
							xVec = water[i][0] - water[ring[0]][0]
							yVec = water[i][1] - water[ring[0]][1]
							zVec = water[i][2] - water[ring[0]][2]
							# this is condition of pbc in the system
							# we need the CG of the ring, so directly outputing the position of the ring center
							if xVec > 0.5 * Lx:
								tempX_pbc += (water[i][0] - Lx)
							elif xVec <= 0.5 * (-Lx):
								tempX_pbc += (water[i][0] + Lx)
							else:
								tempX_pbc += water[i][0]
							if yVec > 0.5 * Ly:
								tempY_pbc += (water[i][1] - Ly)
							elif yVec <= 0.5 * (-Ly):
								tempY_pbc += (water[i][1] + Ly)
							else:
								tempY_pbc += water[i][1]
							if zVec > 0.5 * Lz:
								tempZ_pbc += (water[i][2] - Lz)
							elif zVec <= 0.5 * (-Lz):
								tempZ_pbc += (water[i][2] + Lz)
							else:
								tempZ_pbc += water[i][2]
						# now output the average position of the ring
						return tempX_pbc/(len(ring)), tempY_pbc/(len(ring)), tempZ_pbc/(len(ring))

					# we loop over all the segregated rings
					for ring in ring3:
						tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
						outFile2.write("Li  {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
					for ring in ring4:
						tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
						outFile2.write("Be  {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
					for ring in ring5:
						tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
						outFile2.write("B   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
					for ring in ring6:
						tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
						outFile2.write("C   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
					for ring in ring7:
						tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
						outFile2.write("N   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
					for ring in ring8:
						tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
						outFile2.write("S   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
					if max_ring > 8:
						for ring in ring9:
							tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
							outFile2.write("F   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
						if max_ring > 9:
							for ring in ring10:
								tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
								outFile2.write("Ne  {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
					
					outFile2.close()
					print("The location of the closed rings has been written to 'rings_location.xyz'\n")
					# for visualizer [chimera]; for testing code
					def visualizer(ring):
						for i in range(len(ring)):
							for j in range(len(ring[i])):
								ring[i][j] += 1
						print(ring)
						return ring
					# visualizer(ring4)			# P,Cl,Ar might be good for Li, B, and Be

				hbond_ndx = [[] for _ in range(int(nMols))]
				hbondIndex = []
				hbonds = 0
		
				# reinitialize variables
				oPosX = []
				oPosY = []
				oPosZ = []
				h1Posx = []
				h1Posy = []
				h1Posz = []
				h2Posx = []
				h2Posy = []
				h2Posz = []
				wat_loop = []
				water = [[0 for _ in range(9)] for _ in range(int(nMols))]
				pos = 0
				ref_ndx = 0
			else:
				continue
	outFile.close()

	if directionality:
		print("The closed directional rings have been output to 'ringsCount.dat\n")
	else:
		print("The closed rings have been output to 'ringsCount.dat\n")


def main(argv):
	global directionality, hbond_energy, ring_traj

	directionality = False
	hbond_energy = False
	ring_traj = False
	_haveinputFileName = 0
	max_ring = 0
	ring_closure = 0
	algorithm = "hbondAngle"
	try:
		opts, args = getopt.getopt(argv, "hexdf:r:c:m:", ["help", "energy_defn", "ring_traj", "directional_rings", "input-file=", "ring_size=", "minimal_ring=", "method_algorithm="])	# 'hed' do not take arguments, 'frcm' take argument
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-f", "--input-file"):
			inputFileName = arg
			_haveinputFileName = 1
		elif opt in ("-r", "--ring_size"):
			max_ring = int(arg)
		elif opt in ("-c", "--minimal_ring"):
			ring_closure = int(arg)
		elif opt in ("-m", "--method_algorithm"):
			algorithm = arg
		elif opt in ("-d", "--directional_rings"):
			directionality = True
		elif opt in ("-e", "--energy_defn"):
			hbond_energy = True
		elif opt in ("-x", "--ring_traj"):
			ring_traj = True

	if (_haveinputFileName != 1):
		usage()
		print("No input file was specified")
		sys.exit()

	beginCalc(inputFileName, max_ring, ring_closure, algorithm)

if __name__ == "__main__":
	if len(sys.argv) == 1:
		usage()
		sys.exit()
	main(sys.argv[1:])
