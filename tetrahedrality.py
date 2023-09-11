#!/usr/bin/env python3

"""
This program is used to calculate the tetrahedrality of water using Errington and Debenedetti equation of order paramter.\n

USAGE: ./tetrahedrality filename.gro\n
"""

# import multiprocessing
import sys
import getopt
import os
import os.path
from os import path
import shutil
import math
import numpy as np
import linecache

# input file
inputFileName = sys.argv[1]
filename, file_extension = os.path.splitext(inputFileName)

# list of order parameter and pdb file output files
outputFileName = open("%s_tetra.pdb"%(filename), 'w')
outputFileName2 = open("%s_tetra.dat"%(filename), 'w')
outputFileName2.write(" avg_<q>    stdev\n")

distance_tol = 3.5    # hbond distance tolerance (in Ångström)
distance_tol2 = distance_tol*distance_tol  # hbond square distance tolerance
distance_tol2_long = 2.0*distance_tol*distance_tol   # long square distance tolerance
angle_tol = 30        # hbond angle tolerance (in degrees)
angle_tol_rad = angle_tol * 3.1415926536 / 180.0  # hbond angle tolerance in radians
tolerance = 10e-3
neighbor_tol = 0.5

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
deviation = 0
tetrahedrality = []

# defining vectors with central water and a surrounding water
def defVec(pos1, pos2):
	# first, we do the O-O vector
	posVec[0] = pos2[0]-pos1[0]
	posVec[1] = pos2[1]-pos1[1]
	posVec[2] = pos2[2]-pos1[2]

	# do vector wrapping of periodic boundary conditions
	for k in range(3):
		scaledVec[k] = posVec[k] * invhmat[k][k]
		scaledVec[k] -= round(scaledVec[k])
		posVec[k] = scaledVec[k] * hmat[k][k]
	return posVec

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
			isHBond = 1
			hbondList[i].append(j)
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
				isHBond = 1
				hbondList[i].append(j)
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
					isHBond = 1
					hbondList[i].append(j)
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
						isHBond = 1
						hbondList[i].append(j)
	return isHBond

pos = 0
frame_count = 0
with open (inputFileName, 'r') as inFile:
	nAtoms = int(linecache.getline(inputFileName, 2))
	dummy_atom = linecache.getline(inputFileName, 6)
	five_site = linecache.getline(inputFileName, 8)
	if str(dummy_atom[10:15]).strip() == "MW":
		fourAtomWater = True
		nMols = int(nAtoms/4)
	elif str(five_site[10:15]).strip() in ["OW", "O"]:
		fiveAtomWater = True
		nMols = int(nAtoms/5)
	else:
		threeAtomWater = True
		nMols = nAtoms/3

	# the 9 positio:ns are sorted by O(x,y,z), H1(x,y,z), H2(x,y,z)
	water = [[0 for _ in range(9)] for _ in range(int(nMols))]	# '-1' because '0' is also part of index
	# assume a water can for hbonds with 10 neighbors in worst case senario in liquid phase
	# '-1' because '0' is also part of index
	hbondList= [[] for _ in range(int(nMols))]

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
			# print("test")
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

			halfDistance_tol_array = np.array([distance_tol2_long/2 for i in range(int(nMols))])
			distance_tol_array = np.add.outer(halfDistance_tol_array,halfDistance_tol_array)
			dx = np.subtract.outer(oPosX,oPosX)
			dy = np.subtract.outer(oPosY,oPosY)
			dz = np.subtract.outer(oPosZ,oPosZ)
			
			# do vector wrapping of periodic boundary conditions
			def pbc(oVec, k):
				scaledVec = oVec * invhmat[k][k]
				scaledVec -= np.round(scaledVec)
				oVec = scaledVec * hmat[k][k]
				return oVec
			dx = pbc(dx, 0)
			dy = pbc(dy, 1)
			dz = pbc(dz, 2)
			distance = np.array( np.square(dx) + np.square(dy) + np.square(dz) )
			
			# to find the number of surrounding water molecules around central water molecule
			touching = np.where( distance < distance_tol_array, 1.0, 0.0) #elementwise conditional operator (condition ? 1 : 0 in this case)
			
			# building up the water cluster with proper identification
			neighborList = [[] for i in range(int(nMols))]
			neighborDistance = [[] for i in range(int(nMols))]
			nearestNeighbor_dist = [[0 for _ in range(4)] for _ in range(int(nMols))]
			nearestNeighbor_index = [[0 for _ in range(4)] for _ in range(int(nMols))]
			for i in range(int(nMols)):
				for j in range(int(nMols)):
					if i==j: continue
					if (abs(touching[i,j]*distance_tol_array[i,j] - distance_tol_array[i,j])) < tolerance:
						neighborList[i].append(j)
						neighborDistance[i].append(math.sqrt(distance[i,j]))
					else: continue

			# this is to loop over water environment to test H-bonding
			for i, hbondWat in enumerate(neighborList):
				# this is a reference water
				pos1 = water[i] 
				for j in hbondWat:
					pos2 = water[j]
					isHBond = IsHBonding(pos1, pos2, i, j, isHBond=0)
					hbonds += isHBond

			# build the nearest neighbor list and calculate tetrahedrality
			for i in range(len(neighborList)):
				# let's sort our nearest neighbor list!
				# sort based on distances
				singleSort_dist = neighborDistance[i]
				singleSort_index = neighborList[i]
				# correspondingly rearragne index list based on distance in ascending order
				singleSort_dist, singleSort_index = (list(t) for t in zip(*sorted(zip(singleSort_dist, singleSort_index))))

				# check neighbors to see if we need to sort with H-bonding as a
				# priority . . .
				closest_neighbor = singleSort_dist[0]
				count = 0
				for j in range(1, len(singleSort_dist)):
					temp_dist = singleSort_dist[j] - closest_neighbor
					if (temp_dist < neighbor_tol):
						count += 1

				# print(singleSort_dist)
				# print(len(singleSort_dist))
				# print(count)
				if (count > 4):
					# load the 4 neighbors with a H-bonding+distance preference,
					# then distance only up to the 4 total neighbors
					temp_count = 0
					for j in range(len(hbondList[i])):
						for k in range(count):		# applicable because of sorted distance with central water
							if singleSort_index[k] == hbondList[i][j]:
								nearestNeighbor_index[i][temp_count] = singleSort_index[k]
								nearestNeighbor_dist[i][temp_count] = singleSort_dist[k]
								temp_count += 1
								break			# loop of 'count' is done

						if temp_count > 3:
							break
			
					# print(nearestNeighbor_index[i])
					if temp_count < 4:
						# load in closest neighbors that haven't already been
						# selected as H bonding
						for j in range(count):
							loaded = 0
							for k in range(temp_count):
								if nearestNeighbor_index[i][k] == singleSort_index[j]:
									loaded = 1
							if loaded == 0:
								nearestNeighbor_index[i][temp_count] = singleSort_index[j]
								nearestNeighbor_dist[i][temp_count] = singleSort_dist[j]
								temp_count += 1
							if temp_count > 3:
								break

					# creating a vector matrix of central water with four surrounding waters
					vec_OO = [[] for i in range(int(4))]
					# this is a reference water
					pos1 = water[i] 
					# this is to loop over water environment to test H-bonding
					ndx = 0
					for j in nearestNeighbor_index[i]:
						pos2 = water[j]
						tempVec = defVec(pos1, pos2)
						vec_OO[ndx].append(tempVec[0])
						vec_OO[ndx].append(tempVec[1])
						vec_OO[ndx].append(tempVec[2])
						ndx += 1

			
				else:
					# load the first 4 neighbors in the nearest neighbor lists
					for j in range(4):
						nearestNeighbor_dist[i][j] = singleSort_dist[j]

					# creating a vector matrix of central water with four surrounding waters
					vec_OO = [[] for i in range(int(4))]
					# this is a reference water
					pos1 = water[i] 
					# this is to loop over water environment to test H-bonding
					ndx = 0
					for j in singleSort_index[0:4]:
						pos2 = water[j]
						tempVec = defVec(pos1, pos2)
						vec_OO[ndx].append(tempVec[0])
						vec_OO[ndx].append(tempVec[1])
						vec_OO[ndx].append(tempVec[2])
						ndx += 1

				# tetrahedrality double loop over the nearest neighbors of each molecule
				cosine_sum = 0
				for j in range(3):
					inv_nearestNeighborMag1 = 1.0/nearestNeighbor_dist[i][j]
					unitVec1x = vec_OO[j][0] * inv_nearestNeighborMag1
					unitVec1y = vec_OO[j][1] * inv_nearestNeighborMag1
					unitVec1z = vec_OO[j][2] * inv_nearestNeighborMag1

					for k in range(j+1, 4):
						inv_nearestNeighborMag2 = 1.0/nearestNeighbor_dist[i][k]
						unitVec2x = vec_OO[k][0] * inv_nearestNeighborMag2
						unitVec2y = vec_OO[k][1] * inv_nearestNeighborMag2
						unitVec2z = vec_OO[k][2] * inv_nearestNeighborMag2

						dotProduct = (unitVec1x * unitVec2x) + (unitVec1y * unitVec2y) + (unitVec1z * unitVec2z)
						cosine_sum += pow((dotProduct + 1.0/3.0), 2)

				tetrahedralParam = 1 - (0.375 * cosine_sum)
				tetrahedrality.append(tetrahedralParam)
				
			# calculate the tetrahedrality for the frame
			avg_tetrahedrality = 0
			for i in range(len(tetrahedrality)):
				avg_tetrahedrality += tetrahedrality[i]

			avg_tetrahedrality /= len(tetrahedrality)

			for i in range(len(tetrahedrality)):
				deviation += (tetrahedrality[i] - avg_tetrahedrality)*(tetrahedrality[i] - avg_tetrahedrality)
			stdev_tetrahedrality = math.sqrt(deviation / (len(tetrahedrality)-1))

			# output a file with the list of order parameters
			outputFileName2.write(" %6.4f     %6.4f\n"%(avg_tetrahedrality, stdev_tetrahedrality))

			frame_count += 1
			# output a pdb file for molecular visualization
			outputFileName.write("COMPD    %s_%d\n" %(inputFileName, frame_count))
			outputFileName.write("{:6}{:9.3f}{:9.3f}{:9.3f}\n".format("CRYST1", Lx, Ly, Lz))
			outputFileName.write("MODEL {:>8d}\n".format(frame_count))
			m = 1
			o = 1
			occupancy = 1
			
			for i in range(len(oPosX)):
	 			outputFileName.write("ATOM {:>6d}   OW WAT A{:>4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12}\n".format(m, o, water[i][0], water[i][1], water[i][2], occupancy, tetrahedrality[i], "O"))
	 			m += 1
	 			outputFileName.write("ATOM {:>6d}  HW1 WAT A{:>4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12}\n".format(m, o, water[i][3], water[i][4], water[i][5], occupancy, tetrahedrality[i], "H"))
	 			m += 1
	 			outputFileName.write("ATOM {:>6d}  HW2 WAT A{:>4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12}\n".format(m, o, water[i][6], water[i][7], water[i][8], occupancy, tetrahedrality[i], "H"))
	 			m += 1
	 			o += 1
			outputFileName.write("MASTER     0     0     0     0     0     0     0     0     %d     0     0     0\n" %(3*len(oPosX)))
			outputFileName.write("END\n")
	
			tetrahedrality = []
			deviation = 0
			water = [[0 for _ in range(9)] for _ in range(int(nMols))]
			pos = 0
			hbonds = 0
			hbondList= [[] for _ in range(int(nMols))]
			oPosX = []
			oPosY = []
			oPosZ = []
	outputFileName.close()
print("The files have been output to '%s_tetra.dat' with enlisted order parameters \n and '%s_tetra.pdb' for molecular visualization"%(filename, filename))

