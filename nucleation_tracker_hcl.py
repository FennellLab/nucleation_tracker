#!/usr/bin/env python3

"""
The program is used to calculate the number of hydrogen bonded rings in the pure water system.\n

USAGE: ./nucleation_tracker.py -f filename.gro [-r 10 -c 0 -r 10 -m vertex -d]\n
       -f = input a gro file with correct box dimension
       -r = size of the rings; 9 or 10 membered closed rings; by default to 8 membered polygons
       -c = '0' for minimal ring counting and '1' for total non-self-intersecting rings; by default to minimal rings
       -m = algorithm for ring pruning [vertex, hbond, hbondAngle, torsion]; by default to hbondAngle
       -d = generates only directional rings tracking proton acceptor waters
       -x = Output an .xyz trajectory file containing the ring center locations and type (using atomic number of elements for size). (default=off)
       -b = binning rings for ring distribution; by default at 3.0 nm bin width
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
	fiveAtomWater = False

	if os.path.isfile("ringsCount.dat"):
		# shutil.copy("ringsCount.dat", "ringsCount2.dat")
		os.remove("ringsCount.dat")
	if os.path.isfile("rings_location.xyz"):
		# shutil.copy("rings_location.xyz", "rings_location2.xyz")	
		os.remove("rings_location.xyz")
	
	if ring_traj:
		outFile2 = open("rings_location.xyz", 'w')
		# for binning rings
		setRings_Li = []
		setRings_Be = []
		setRings_B = []
		setRings_C = []
		setRings_N = []
		setRings_O = []
		setRings_F = []
		setRings_Ne = []

	# for binning ring distribution; the 3 positions are sorted by polygons (x,y,z)
	if binning_rings:
		outFile3 = open("rings_dist.dat", 'w')
		outFile3.write("{:^10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}".format("bins", "rings_3", "rings_4", "rings_5", "rings_6", "rings_7", "rings_8"))
		if max_ring > 8:
			outFile3.write("{:>10s}".format("rings_9"))
			if max_ring > 9:
				outFile3.write("{:>10s}\n".format("rings_10"))
			else:
				outFile3.write("\n")
		else:
			outFile3.write("\n")

	print("\nEnumerating the closed rings in the system . . . \n")
	outFile = open("ringsCount.dat", 'w')
	outFile.write("{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}".format("Frames", "hbonds", "RSF_val", "rings_3", "rings_4", "rings_5", "rings_6", "rings_7", "rings_8"))
	if max_ring > 8:
		outFile.write("{:>10s}".format("rings_9"))
		if max_ring > 9:
			outFile.write("{:>10s}".format("rings_10"))
	outFile.write("\n")
	
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
			
			if i < cl_atoms:
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
					else:
						if cl_atoms <= j < 2*cl_atoms:
							# last level... I promise...
							# the O-H3
							posVec[0] = pos2[9]-pos2[0]
							posVec[1] = pos2[10]-pos2[1]
							posVec[2] = pos2[11]-pos2[2]
							
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
						if cl_atoms <= i < 2*cl_atoms:
							# continue down the rabbit hole
							# the O-H3
							posVec[0] = pos1[9]-pos1[0]
							posVec[1] = pos1[10]-pos1[1]
							posVec[2] = pos1[11]-pos1[2]
							
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
								else:
									if cl_atoms <= j < 2*cl_atoms:
										# last level... I promise...
										# the O-H3
										posVec[0] = pos2[9]-pos2[0]
										posVec[1] = pos2[10]-pos2[1]
										posVec[2] = pos2[11]-pos2[2]
										
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
	pos2 = 0
	pos3 = 0
	frame_count = 0
	cl_atoms = 0

	with open (inputFileName, 'r') as inFile:
		line_lists = inFile.readlines()
		last_atom = line_lists[-2]			# 'line_lists[-2].strip()' removes all blank spaces
		nMols = int(last_atom[0:5])

		for i in range(nMols):
			if str(line_lists[i][10:15]).strip() in ["CL"]:
				cl_atoms += 1

	# assume a water can for hbonds with 10 neighbors in worst case senario in liquid phase
	# '-1' because '0' is also part of index
	hbond_ndx = [[] for _ in range(nMols)]

	# the 3 positions are sorted by Cl(x,y,z)
	chloride = [[0 for _ in range(3)] for _ in range(cl_atoms)]	# '-1' because '0' is also part of index

	# the 12 positions are sorted by HY1(x,y,z), HY2(x,y,z), HY3(x,y,z), OY(x,y,z)
	hydronium = [[0 for _ in range(12)] for _ in range(cl_atoms)]

	# number of water molecules
	nWat = (nMols - (2*cl_atoms))		# consider TIP4P water model
	# the 9 positio:ns are sorted by O(x,y,z), H1(x,y,z), H2(x,y,z), H3(x,y,z)
	water = [[0 for _ in range(9)] for _ in range(int(nWat))]

	with open (inputFileName, 'r') as inFile:
		line = inFile.readline()					# we are incrementing the lines everytime we call this built-in function
		line = inFile.readline()
		while line != "":
			nAtoms = int(line)

			for i in range(nAtoms):			# this will finish reading lines equal to atom number
				line = inFile.readline()

				if str(line[10:15]).strip() in ["MW"]:
					continue

				elif str(line[10:15]).strip() in ["CL"]:
					chloride[pos][0] = float(line[20:28])*10 
					chloride[pos][1] = float(line[28:36])*10
					chloride[pos][2] = float(line[36:44])*10
					oPosX.append(float(line[20:28])*10)
					oPosY.append(float(line[28:36])*10)
					oPosZ.append(float(line[36:44])*10)
					pos += 1

				elif str(line[10:15]).strip() in ["OY"]:
					hydronium[pos2][0] = float(line[20:28])*10 
					hydronium[pos2][1] = float(line[28:36])*10
					hydronium[pos2][2] = float(line[36:44])*10
					oPosX.append(float(line[20:28])*10)
					oPosY.append(float(line[28:36])*10)
					oPosZ.append(float(line[36:44])*10)
					pos2 += 1

				elif str(line[10:15]).strip() in ["HY1"]:
					hydronium[pos2][3] = float(line[20:28])*10 
					hydronium[pos2][4] = float(line[28:36])*10
					hydronium[pos2][5] = float(line[36:44])*10

				elif str(line[10:15]).strip() in ["HY2"]:
					hydronium[pos2][6] = float(line[20:28])*10 
					hydronium[pos2][7] = float(line[28:36])*10
					hydronium[pos2][8] = float(line[36:44])*10

				elif str(line[10:15]).strip() in ["HY3"]:
					hydronium[pos2][9] = float(line[20:28])*10 
					hydronium[pos2][10] = float(line[28:36])*10
					hydronium[pos2][11] = float(line[36:44])*10

				elif str(line[10:15]).strip() in ["OW", "O"]:
					water[pos3][0] = float(line[20:28])*10 
					water[pos3][1] = float(line[28:36])*10
					water[pos3][2] = float(line[36:44])*10
					oPosX.append(float(line[20:28])*10)
					oPosY.append(float(line[28:36])*10)
					oPosZ.append(float(line[36:44])*10)
		
				elif str(line[10:15]).strip() in ["HW1"]:
					water[pos3][3] = float(line[20:28])*10 
					water[pos3][4] = float(line[28:36])*10
					water[pos3][5] = float(line[36:44])*10

				elif str(line[10:15]).strip() in ["HW2"]:
					water[pos3][6] = float(line[20:28])*10 
					water[pos3][7] = float(line[28:36])*10
					water[pos3][8] = float(line[36:44])*10
					pos3 += 1

			line = inFile.readline()
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
			
			# calling hbond geometry for ring connectivity
			# this is to loop over water environment to test H-bonding
			for i, hbondWat in enumerate(hbondListIndex):
				# this is a reference water
				if i < cl_atoms:
					pos1 = chloride[i]
				elif cl_atoms <= i < (2 * cl_atoms):
					pos1 = hydronium[i - cl_atoms]
				else:
					pos1 = water[i - (2*cl_atoms)]
				for j in hbondWat:
					if j < cl_atoms:
						pos2 = chloride[j]
					elif cl_atoms <= j < (2 * cl_atoms):
						pos2 = hydronium[j - cl_atoms]
					else:
						pos2 = water[j - (2*cl_atoms)]
					isHBond = IsHBonding(pos1, pos2, i, j, isHBond=0)
					hbonds += isHBond

			# we test the rings if it crosses the pbc; rings should come back from the same wall
			# This check ensures that a chain like |-A-B-C|-A-B-C|-A- is not counted as a polygon 'A-B-C'
			def ring_pbc(member1, member2, outRing_x, outRing_y, outRing_z):
				if member1 < cl_atoms:
					tempVec1 = chloride[member1]
				elif cl_atoms <= member1 < (2 * cl_atoms):
					tempVec1 = hydronium[member1 - cl_atoms]
				else:
					tempVec1 = water[member1 - (2*cl_atoms)]

				if member2 < cl_atoms:
					tempVec2 = chloride[member2]
				elif cl_atoms <= member2 < (2 * cl_atoms):
					tempVec2 = hydronium[member2 - cl_atoms]
				else:
					tempVec2 = water[member2 - (2*cl_atoms)]

				xVecRing = tempVec1[0] - tempVec2[0]
				yVecRing = tempVec1[1] - tempVec2[1]
				zVecRing = tempVec1[2] - tempVec2[2]
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
			for init_wat in range(int(nMols)):
				temp_loop = []
				for i in hbond_ndx[init_wat]:
					for j in hbond_ndx[i]:
						if j == init_wat:		# already a hbond; skip already counted rings; saves computational time
							continue
						else:
							for k in hbond_ndx[j]:
								if k == i or k == j:		# skip already counted rings; saves computational time
									continue
								elif k == init_wat:			# rings back to the first water
									temp_loop.extend([init_wat, i, j])
									# we do not want double counting of rings . . .
									wat_sort = sorted(list(set(temp_loop)))
									# we take 'temp_loop' for avoiding a test for hbonds
									if (wat_sort not in wat_loop):
										# test for ring pbc; for forming rings, there should be equal outgoing and incoming hbonds
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
										if count_x%2 == 0 and count_y%2 == 0 and count_z%2 == 0:
											wat_loop.append(wat_sort)
										count_x = 0
										count_y = 0
										count_z = 0
										temp_loop = []
									temp_loop = []
								else:
									for l in hbond_ndx[k]:
										if l == i or l == j or l == k:			# skip already counted rings; saves computational time
											continue
										elif l == init_wat:
											temp_loop.extend([init_wat, i, j, k])
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
												if count_x%2 == 0 and count_y%2 == 0 and count_z%2 == 0:
													wat_loop.append(wat_sort)
												count_x = 0
												count_y = 0
												count_z = 0
												temp_loop = []
											temp_loop = []
										else:
											for m in hbond_ndx[l]:
												if m == i or m == j or m == k or m == l:
													continue
												elif m == init_wat:
													temp_loop.extend([init_wat, i, j, k, l])
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
														if count_x%2 == 0 and count_y%2 == 0 and count_z%2 == 0:
															wat_loop.append(wat_sort)
														count_x = 0
														count_y = 0
														count_z = 0
														temp_loop = []
													temp_loop = []
												else:
													for n in hbond_ndx[m]:
														if n == i or n == j or n == k or n == l or n == m:
															continue
														elif n == init_wat:
															temp_loop.extend([init_wat, i, j, k, l, m])
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
																if count_x%2 == 0 and count_y%2 == 0 and count_z%2 == 0:
																	wat_loop.append(wat_sort)
																count_x = 0
																count_y = 0
																count_z = 0
																temp_loop = []
															temp_loop = []
														else:
															for o in hbond_ndx[n]:
																if o == i or o == j or o == k or o == l or o == m or o == n:
																	continue
																elif o == init_wat:
																	temp_loop.extend([init_wat, i, j, k, l, m, n])
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
																		if count_x%2 == 0 and count_y%2 == 0 and count_z%2 == 0:
																			wat_loop.append(wat_sort)
																		count_x = 0
																		count_y = 0
																		count_z = 0
																		temp_loop = []
																	temp_loop = []
																else:
																	for p in hbond_ndx[o]:
																		if p == i or p == j or p == k or p == l or p == m or p == n or p == o:
																			continue
																		elif p == init_wat:
																			temp_loop.extend([init_wat, i, j, k, l, m, n, o])
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
																				if count_x%2 == 0 and count_y%2 == 0 and count_z%2 == 0:
																					wat_loop.append(wat_sort)
																				count_x = 0
																				count_y = 0
																				count_z = 0
																				temp_loop = []
																			temp_loop = []
																		else:
																			if max_ring > 8:
																				for q in hbond_ndx[p]:
																					if q in [i, j, k, l, m, n, o, p]:
																						continue
																					elif q == init_wat:
																						temp_loop.extend([init_wat, i, j, k, l, m, n, o, p])
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
																							if count_x%2 == 0 and count_y%2 == 0 and count_z%2 == 0:
																								wat_loop.append(wat_sort)
																							count_x = 0
																							count_y = 0
																							count_z = 0
																							temp_loop = []
																						temp_loop = []
																					else:
																						if max_ring > 9:
																							for r in hbond_ndx[q]:
																								if r in [i, j, k, l, m, n, o, p, q]:
																									continue
																								elif r == init_wat:
																									temp_loop.extend([init_wat, i, j, k, l, m, n, o, p, q])
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
																										if count_x%2 == 0 and count_y%2 == 0 and count_z%2 == 0:
																											wat_loop.append(wat_sort)
																										count_x = 0
																										count_y = 0
																										count_z = 0
																										temp_loop = []
																									temp_loop = []
																			continue

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
			rsf_val = len(ring3) + len(ring4) + len(ring5) + len(ring6) + len(ring7) + len(ring8)
			if max_ring > 8:
				rsf_val += len(ring9)
				if max_ring > 9:
					rsf_val += len(ring10)
			rsf_val /= (2*nMols)
			outFile.write("{:>10d}{:>10d}{:>10.4f}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}".format(frame_count, hbonds, rsf_val, len(ring3), len(ring4), len(ring5), len(ring6), len(ring7), len(ring8)))
			if max_ring > 8:
				outFile.write("{:>10d}".format(len(ring9)))
				if max_ring > 9:
					outFile.write("{:>10d}\n".format(len(ring10)))
				else:
					outFile.write("\n")
			else:
				outFile.write("\n")
			outFile.flush()
	
			if ring_traj:
		 		# writing xyz file for the location of directional rings
				# wrt to oxygen (heavy atom), if take hydrogen atoms, then divide by total items taken
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

					if ring[0] < cl_atoms:
						tempVec1 = chloride[ring[0]]
					elif cl_atoms <= ring[0] < (2 * cl_atoms):
						tempVec1 = hydronium[ring[0] - cl_atoms]
					else:
						tempVec1 = water[ring[0] - (2*cl_atoms)]

					# "water[ring[0]][0]" is the first member of the rings and taken as the reference molecule for "min_image_wrap"
					tempX_pbc += tempVec1[0]
					tempY_pbc += tempVec1[1]
					tempZ_pbc += tempVec1[2]
					# now we loop over all remaining waters and determine the pbc of other molecule with reference to the first molecule
					for i in ring[1:]:
						if i < cl_atoms:
							tempVec2 = chloride[i]
						elif cl_atoms <= i < (2 * cl_atoms):
							tempVec2 = hydronium[i - cl_atoms]
						else:
							tempVec2 = water[i - (2*cl_atoms)]

						xVec = tempVec2[0] - tempVec1[0]
						yVec = tempVec2[1] - tempVec1[1]
						zVec = tempVec2[2] - tempVec1[2]
						# this is condition of pbc in the system
						# we need the CG of the ring, so directly outputing the position of the ring center
						if xVec > 0.5 * Lx:
							tempX_pbc += (tempVec2[0] - Lx)
						elif xVec <= 0.5 * (-Lx):
							tempX_pbc += (tempVec2[0] + Lx)
						else:
							tempX_pbc += tempVec2[0]
						if yVec > 0.5 * Ly:
							tempY_pbc += (tempVec2[1] - Ly)
						elif yVec <= 0.5 * (-Ly):
							tempY_pbc += (tempVec2[1] + Ly)
						else:
							tempY_pbc += tempVec2[1]
						if zVec > 0.5 * Lz:
							tempZ_pbc += (tempVec2[2] - Lz)
						elif zVec <= 0.5 * (-Lz):
							tempZ_pbc += (tempVec2[2] + Lz)
						else:
							tempZ_pbc += tempVec2[2]
					# now output the average position of the ring
					return tempX_pbc/(len(ring)), tempY_pbc/(len(ring)), tempZ_pbc/(len(ring))

				# we loop over all the segregated rings
				for ring in ring3:
					tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
					setRings_Li.append([tempX_pbc, tempY_pbc, tempZ_pbc])
					outFile2.write("Li  {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
				for ring in ring4:
					tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
					setRings_Be.append([tempX_pbc, tempY_pbc, tempZ_pbc])
					outFile2.write("Be  {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
				for ring in ring5:
					tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
					setRings_B.append([tempX_pbc, tempY_pbc, tempZ_pbc])
					outFile2.write("B   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
				for ring in ring6:
					tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
					setRings_C.append([tempX_pbc, tempY_pbc, tempZ_pbc])
					outFile2.write("C   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
				for ring in ring7:
					tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
					setRings_N.append([tempX_pbc, tempY_pbc, tempZ_pbc])
					outFile2.write("N   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
				for ring in ring8:
					tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
					setRings_O.append([tempX_pbc, tempY_pbc, tempZ_pbc])
					outFile2.write("O   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
				if max_ring > 8:
					for ring in ring9:
						tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
						setRings_F.append([tempX_pbc, tempY_pbc, tempZ_pbc])
						outFile2.write("F   {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
					if max_ring > 9:
						for ring in ring10:
							tempX_pbc, tempY_pbc, tempZ_pbc = min_image_wrap(ring)
							setRings_Ne.append([tempX_pbc, tempY_pbc, tempZ_pbc])
							outFile2.write("Ne  {:>10.5f} {:>10.5f} {:>10.5f}\n".format(tempX_pbc, tempY_pbc, tempZ_pbc))
				outFile2.flush()
				
				# for visualizer [chimera]; for testing code
				def visualizer(ring):
					for i in range(len(ring)):
						for j in range(len(ring[i])):
							ring[i][j] += 1
					print(ring)
					return ring
				# visualizer(ring4)			# P,Cl,Ar might be good for Li, B, and Be

			# Now, let us bin the ring distribution
			if binning_rings:			# please check for the box dimension to manage the bin width; we consider z-axis in this case
				bin_width = 30				# unit in Angstrom
				binNum = 0
				outFile3.write(";%f 0 0; 0 %f 0; 0 0 %f; hbonds = %d\n"%(Lx, Ly, Lz, hbonds))
				
				# final binning of the system
				for i in range(0, int(Lz), bin_width):
					binNum += 1
					# last binning of the system; less than double of the bin_width (60); compensation bin
					if (Lz - i) < 50:	
						def count_rings(ringList,countRings):
							for j in range(len(ringList)):
								if i < ringList[j][2] <= Lz:
									countRings += 1
								elif ringList[j][2] > Lz or ringList[j][2] <= 0:	# rings outside the pbc
									countRings += 1
							return countRings
						count_Li = count_rings(setRings_Li, 0)
						count_Be = count_rings(setRings_Be, 0)
						count_B = count_rings(setRings_B, 0)
						count_C = count_rings(setRings_C, 0)
						count_N = count_rings(setRings_N, 0)
						count_O = count_rings(setRings_O, 0)
						if max_ring > 8:
							count_F = count_rings(setRings_F, 0)
						if max_ring > 9:
							count_Ne = count_rings(setRings_Ne, 0)

						# writing a file after binning
						outFile3.write("  bin{:<2d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}".format(binNum, count_Li, count_Be, count_B, count_C, count_N, count_O))
						if max_ring > 8:
							outFile3.write("{:>10d}".format(count_F))
							if max_ring > 9:
								outFile3.write("{:>10d}\n".format(count_Ne))
							else:
								outFile3.write("\n")
						else:
							outFile3.write("\n")
						outFile3.write("\n")
						break
					else:
						# sequential binning
						def count_rings(ringList,countRings):
							for j in range(len(ringList)):
								if i < ringList[j][2] <= (i+bin_width):
									countRings += 1
							return countRings
						count_Li = count_rings(setRings_Li, 0)
						count_Be = count_rings(setRings_Be, 0)
						count_B = count_rings(setRings_B, 0)
						count_C = count_rings(setRings_C, 0)
						count_N = count_rings(setRings_N, 0)
						count_O = count_rings(setRings_O, 0)
						if max_ring > 8:
							count_F = count_rings(setRings_F, 0)
						if max_ring > 9:
							count_Ne = count_rings(setRings_Ne, 0)

						# writing a file after binning
						outFile3.write("  bin{:<2d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}".format(binNum, count_Li, count_Be, count_B, count_C, count_N, count_O))
						if max_ring > 8:
							outFile3.write("{:>10d}".format(count_F))
							if max_ring > 9:
								outFile3.write("{:>10d}\n".format(count_Ne))
							else:
								outFile3.write("\n")
						else:
							outFile3.write("\n")
				outFile3.flush()

				# for binning ring distribution; the 3 positions are sorted by polygons (x,y,z)
				setRings_Li = []
				setRings_Be = []
				setRings_B = []
				setRings_C = []
				setRings_N = []
				setRings_O = []
				setRings_F = []
				setRings_Ne = []

			# reinitialize variables
			hbond_ndx = [[] for _ in range(nMols)]
			hbondIndex = []
			hbonds = 0
			oPosX = []
			oPosY = []
			oPosZ = []
			wat_loop = []
			primitiveRings = []
			minimal_rings = []
			water = [[0 for _ in range(9)] for _ in range(nWat)]
			chloride = [[0 for _ in range(3)] for _ in range(cl_atoms)]	# '-1' because '0' is also part of index
			hydronium = [[0 for _ in range(12)] for _ in range(cl_atoms)]
			ref_ndx = 0
			pos = 0
			pos2 = 0
			pos3 = 0
			# rings segregation
			ring3 = []
			ring4 = []
			ring5 = []
			ring6 = []
			ring7 = []
			ring8 = []
			ring9 = []
			ring10 = []
			# for "readline" built-in functions
			line = inFile.readline()
			line = inFile.readline()

	outFile.close()
	if ring_traj:
		outFile2.close()
	if binning_rings:
		outFile3.close()

	if directionality:
		print("The closed directional rings have been output to 'ringsCount.dat\n")
	else:
		print("The closed rings have been output to 'ringsCount.dat\n")
	if ring_traj:
		print("The location of the closed rings has been written to 'rings_location.xyz'\n")
	if binning_rings:
		print("The binning of the trajectory has been written to 'rings_dist.dat'")

def main(argv):
	global directionality, ring_traj, binning_rings

	directionality = False
	ring_traj = False
	binning_rings = False
	_haveinputFileName = 0
	max_ring = 0
	ring_closure = 0
	algorithm = "hbondAngle"
	try:
		opts, args = getopt.getopt(argv, "hxbdf:r:c:m:", ["help", "ring_traj", "binning_rings", "directional_rings", "input-file=", "ring_size=", "minimal_ring=", "method_algorithm="])	# 'hed' do not take arguments, 'frcm' take argument
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-x", "--ring_traj"):
			ring_traj = True
		elif opt in ("-b", "--binning_rings"):
			binning_rings = True
			ring_traj = True		# requires ring trajectory file
		elif opt in ("-d", "--directional_rings"):
			directionality = True
		elif opt in ("-f", "--input-file"):
			inputFileName = arg
			_haveinputFileName = 1
		elif opt in ("-r", "--ring_size"):
			max_ring = int(arg)
		elif opt in ("-c", "--minimal_ring"):
			ring_closure = int(arg)
		elif opt in ("-m", "--method_algorithm"):
			algorithm = arg

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

