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
       -p = Build a pov_files directory containing .pov files for rendering of ring closures and locations. (default=off)
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
import itertools
import threading
import time

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

	# creating povray directory for saving pov files
	if povrayOutOpt:
		# create a pov_files directroy
		directory = "pov_files"
		current_directory = os.getcwd()
		path = os.path.join(current_directory, directory)
		# if not os.path.exists(final_directory):
   		# 	os.makedirs(final_directory)
		try: 
			os.makedirs(path, exist_ok = False) 
			# print("Directory '%s' created successfully" % directory) 
		except OSError as error: 
			print("Directory '%s' already existed" % directory) 

		# writing a single file for povray program execution
		fileName, file_extension = os.path.splitext(inputFileName)
		pov_info = open("%s/%s_pov.txt"%(path, fileName), 'w')

	outFile = open("ringsCount.dat", 'w')
	outFile.write("{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}".format("Frames", "hbonds", "RSF_val", "rings_3", "rings_4", "rings_5", "rings_6", "rings_7", "rings_8"))
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
		watAtoms = nAtoms
		for lc in range(3,watAtoms):
			if str(linecache.getline(inputFileName, lc))[10:15].strip() in ["OW", "O"]:
				break
		dummy_atom = linecache.getline(inputFileName, lc+3)
		five_site = linecache.getline(inputFileName, lc+5)
		if str(dummy_atom[10:15]).strip() == "MW":
			fourAtomWater = True
			nMols = int((watAtoms-(lc-3))/4)
		elif str(five_site[10:15]).strip() in ["OW", "O"]:
			fiveAtomWater = True
			nMols = int((watAtoms-(lc-3))/5)
		else:
			threeAtomWater = True
			nMols = (watAtoms-(lc-3))/3
	
		print("\nEnumerating the closed rings in the system . . .")
		#here is the animation
		done = False
		def animate():
			for c in itertools.cycle(['|', '/', '-', '\\']):
				if done or frame_count > 5:      # this is to limit the computational cost for large frames loading:
					break
				sys.stdout.write('\rloading ' + c)
				sys.stdout.flush()
				time.sleep(0.1)
			sys.stdout.write('\rDone!      \n ')
		
		t = threading.Thread(target=animate)
		t.start()
		
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
				
			elif (i-(frame_count*3))%nAtoms == 0 and i != 0:
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
						print("Note: Need to be orthogonal lattice/edge vector for minimum image wrapping")
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
					
					# for visualizer [chimera]; for testing code; increment by "+1" in chimera
					def visualizer(ring):
						for i in range(len(ring)):
							for j in range(len(ring[i])):
								ring[i][j] += 1
						print(ring)
						return ring
					# visualizer(ring5)			# P,Cl,Ar might be good for Li, B, and Be

				# Now, let us render the POV-Ray images
				class PovObjects:
					# count = 0
					def __init__(self, filename, boxLength, ring3_dot_size=0.5, ring4_dot_size=0.6, ring5_dot_size=0.7, ring6_dot_size=0.8, ring7_dot_size=0.9, ring8_dot_size=1.0, ring9_dot_size=1.1, ring10_dot_size=1.2):		# this is class constructor; __init__ is executed automatically on we call class. It is mostly used to initialize the class attributes. 
						self.filename = filename
						self.boxLength = boxLength
						self.ring3_dot_size = ring3_dot_size
						self.ring4_dot_size = ring4_dot_size
						self.ring5_dot_size = ring5_dot_size
						self.ring6_dot_size = ring6_dot_size
						self.ring7_dot_size = ring7_dot_size
						self.ring8_dot_size = ring8_dot_size
						self.ring9_dot_size = ring9_dot_size
						self.ring10_dot_size = ring10_dot_size

						# writing ring coordinates to pov files
						self.buffer_size = 2		# the imaging buffer size
						self.slab_thickness = 100	# 0.5*slab thickness
						self.scale_factor = 20		# scale the pixels!
						self.transparency = 0.5		# the transparency of ring dots

					def outFile(self):
						self.filename.write("""  //**************************************
  // generated by nucleation_tracker
  //**************************************
  
  //**************************************
  // Lights, camera, resolution!
  //**************************************
  global_settings{ max_trace_level 100 }\n\n""")
						self.filename.write("#declare Ratio = {};\n".format(self.boxLength[0]/self.boxLength[1]))
						self.filename.write("#declare zoom = {:>6.2f};\n".format(5*self.boxLength[2]))
						self.filename.write("""#declare RAD = off;
global_settings {
#if(RAD)
radiosity {
pretrace_start 0.08
pretrace_end   0.01
count 500
nearest_count 10
error_bound 0.02
recursion_limit 1
low_error_factor 0.2
gray_threshold 0.0
minimum_reuse 0.015
brightness 1.4
adc_bailout 0.01/2
}
#end
}

camera{
orthographic
location < 0, 0, zoom >
direction < 0, 0, 2 >\n""")
						self.filename.write("up < 0, {}, 0 >\n".format(self.boxLength[1]))
						self.filename.write("// Ratio is negative to switch povray to a right hand coordinate system.\n")
						self.filename.write("right < -{}, 0, 0 >\n".format(self.boxLength[0]))
						self.filename.write("""look_at < 0, 0, 0 >
}

background { color rgb 1 }

global_settings { ambient_light rgb 1 }\n""")
#         "    location < 0,0, 0.5333333*" << boxLength << "*zoom >
#         //"    location < 0,0, 0.8*" << boxLength << "*zoom >
#         "        direction < 0,0, zoom >
#         "            // Ratio is negative to switch povray to a right hand coordinate system.
#         "                right < -Ratio ,0 , 0 >
#         "                    look_at < 0, -3.5, 0 >
						self.filename.write("light_source { < 0, 0, 10*zoom > rgb < 1.0, 1.0, 1.0 > }\n")
						self.filename.write("light_source { < 0, 0, 10*zoom > rgb < 0.1, 0.1, 0.1 > }\n")
						return

					def printRing3(self):
						#### For sphere
						# self.filename.write("""\n#macro ring3 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)
  # intersection{
    # sphere{
      # < center_x, center_y, center_z >,\n""")
						# self.filename.write("        {}\n".format(self.ring3_dot_size))
						#### For triangle; Number of vertices = 3, Make sure the last vertex is the same as the first one
						self.filename.write("""\n#macro ring3 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)
  intersection{
    polygon{
      3,
      <center_x, center_y+1.125, center_z>,
      <center_x-0.487125, center_y-0.5625, center_z>,
      <center_x+0.487125, center_y-0.5625, center_z>\n""")
						self.filename.write("""        texture{
          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
          finish{
            ambient .2
	          diffuse .6
	          specular .1
              roughness .1
          }
        }
    }
   box{\n""")
						self.filename.write("      < -{}, -{}, center_z-0.01 >,< {}, {}, center_z+0.01 >".format(self.boxLength[0], self.boxLength[1], self.boxLength[0], self.boxLength[0])) 
						self.filename.write("""
        texture{
	      pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
	      finish{
	        ambient .2
	          diffuse .6
	          specular .1
	          roughness .1
	      }
	    }
     }
   }
#end\n""")
						return

					#def printRing3_loc(self, setRings_Li, red_val = "0.5", green_val = "0.5", blue_val = "1.0", zVal = -0.1):
					def printRing3_loc(self, setRings_Li, red_val = "1.0", green_val = "0.0", blue_val = "0.5", zVal = -0.1):
						for j in range(len(setRings_Li)):
							if (setRings_Li[j][2] <= self.slab_thickness and setRings_Li[j][2] >= -self.slab_thickness):
								self.red_val = red_val
								self.green_val = green_val
								self.blue_val = blue_val
								self.zVal = zVal
								self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Li[j][0], setRings_Li[j][1], self.zVal, self.red_val, green_val, blue_val, self.transparency))

								if (setRings_Li[j][0] < (-0.5*self.boxLength[0]+self.buffer_size)):
									if (setRings_Li[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]+self.boxLength[0]), (setRings_Li[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]+self.boxLength[0]), setRings_Li[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Li[j][0], (setRings_Li[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_Li[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]+self.boxLength[0]), (setRings_Li[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]+self.boxLength[0]), setRings_Li[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Li[j][0], (setRings_Li[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]+self.boxLength[0]), setRings_Li[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_Li[j][0] > (0.5*self.boxLength[0]-self.buffer_size)):
									if (setRings_Li[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]-self.boxLength[0]), (setRings_Li[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]-self.boxLength[0]), setRings_Li[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Li[j][0], (setRings_Li[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_Li[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]-self.boxLength[0]), (setRings_Li[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]-self.boxLength[0]), setRings_Li[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Li[j][0], (setRings_Li[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Li[j][0]-self.boxLength[0]), setRings_Li[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_Li[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
									self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Li[j][0], (setRings_Li[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_Li[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
									self.filename.write("ring3(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Li[j][0], (setRings_Li[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))

					def printRing4(self):
						self.filename.write("""\n#macro ring4 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)
  intersection{
    polygon{
      4,
      <center_x - 0.5261, center_y + 0.5261, center_z>,
      <center_x - 0.5261, center_y - 0.5261, center_z>,
      <center_x + 0.5261, center_y - 0.5261, center_z>,
      <center_x + 0.5261, center_y + 0.5261, center_z>\n""")
						self.filename.write("""        texture{
          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
          finish{
            ambient .2
	          diffuse .6
	          specular .1
              roughness .1
          }
        }
    }
   box{\n""")
						self.filename.write("      < -{}, -{}, center_z-0.01 >,< {}, {}, center_z+0.01 >".format(self.boxLength[0], self.boxLength[1], self.boxLength[0], self.boxLength[0])) 
						self.filename.write("""
        texture{
	      pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
	      finish{
	        ambient .2
	          diffuse .6
	          specular .1
	          roughness .1
	      }
	    }
     }
   }
#end\n""")
						return

					def printRing4_loc(self, setRings_Be, red_val = "0.5", green_val = "0.0", blue_val = "1.0", zVal = -0.1):
						for j in range(len(setRings_Be)):
							if (setRings_Be[j][2] <= self.slab_thickness and setRings_Be[j][2] >= -self.slab_thickness):
								self.red_val = red_val
								self.green_val = green_val
								self.blue_val = blue_val
								self.zVal = zVal
								self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Be[j][0], setRings_Be[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))

								if (setRings_Be[j][0] < (-0.5*self.boxLength[0]+self.buffer_size)):
									if (setRings_Be[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]+self.boxLength[0]), (setRings_Be[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]+self.boxLength[0]), setRings_Be[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Be[j][0], (setRings_Be[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_Be[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]+self.boxLength[0]), (setRings_Be[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]+self.boxLength[0]), setRings_Be[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Be[j][0], (setRings_Be[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]+self.boxLength[0]), setRings_Be[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_Be[j][0] > (0.5*self.boxLength[0]-self.buffer_size)):
									if (setRings_Be[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]-self.boxLength[0]), (setRings_Be[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]-self.boxLength[0]), setRings_Be[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Be[j][0], (setRings_Be[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_Be[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]-self.boxLength[0]), (setRings_Be[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]-self.boxLength[0]), setRings_Be[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Be[j][0], (setRings_Be[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Be[j][0]-self.boxLength[0]), setRings_Be[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_Be[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
									self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Be[j][0], (setRings_Be[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_Be[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
									self.filename.write("ring4(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Be[j][0], (setRings_Be[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))

					def printRing5(self):
						self.filename.write("""\n#macro ring5 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)
  intersection{
    polygon{
      5,
      <center_x, center_y + 0.825, center_z>,
      <center_x - 0.64658, center_y + 0.275, center_z>,
      <center_x - 0.39842, center_y - 0.6149, center_z>,
      <center_x + 0.39842, center_y - 0.6149, center_z>,
      <center_x + 0.64658, center_y + 0.275, center_z>""")
						self.filename.write("""        texture{
          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
          finish{
            ambient .2
	          diffuse .6
	          specular .1
              roughness .1
          }
        }
    }
   box{\n""")
						self.filename.write("      < -{}, -{}, center_z-0.01 >,< {}, {}, center_z+0.01 >".format(self.boxLength[0], self.boxLength[1], self.boxLength[0], self.boxLength[0])) 
						self.filename.write("""
        texture{
	      pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
	      finish{
	        ambient .2
	          diffuse .6
	          specular .1
	          roughness .1
	      }
	    }
     }
   }
#end\n""")
						return

					def printRing5_loc(self, setRings_B, red_val = "0.0", green_val = "0.5", blue_val = "1.0", zVal = -0.1):
						for j in range(len(setRings_B)):
							if (setRings_B[j][2] <= self.slab_thickness and setRings_B[j][2] >= -self.slab_thickness):
								self.red_val = red_val
								self.green_val = green_val
								self.blue_val = blue_val
								self.zVal = zVal
								self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_B[j][0], setRings_B[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))

								if (setRings_B[j][0] < (-0.5*self.boxLength[0]+self.buffer_size)):
									if (setRings_B[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]+self.boxLength[0]), (setRings_B[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]+self.boxLength[0]), setRings_B[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_B[j][0], (setRings_B[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_B[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]+self.boxLength[0]), (setRings_B[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]+self.boxLength[0]), setRings_B[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_B[j][0], (setRings_B[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]+self.boxLength[0]), setRings_B[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_B[j][0] > (0.5*self.boxLength[0]-self.buffer_size)):
									if (setRings_B[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]-self.boxLength[0]), (setRings_B[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]-self.boxLength[0]), setRings_B[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_B[j][0], (setRings_B[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_B[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]-self.boxLength[0]), (setRings_B[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]-self.boxLength[0]), setRings_B[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_B[j][0], (setRings_B[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_B[j][0]-self.boxLength[0]), setRings_B[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_B[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
									self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_B[j][0], (setRings_B[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_B[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
									self.filename.write("ring5(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_B[j][0], (setRings_B[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))


					def printRing6(self):
						self.filename.write("""\n#macro ring6 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)
  intersection{
    polygon{
      6,
      <center_x, center_y+0.75, center_z>,
      <center_x-0.6495, center_y+0.375, center_z>,
      <center_x-0.6495, center_y-0.375, center_z>,
      <center_x, center_y-0.75, center_z>,
      <center_x+0.6495, center_y-0.375, center_z>,
      <center_x+0.6495, center_y+0.375, center_z>\n""")
						self.filename.write("""        texture{
          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
          finish{
            ambient .2
	          diffuse .6
	          specular .1
              roughness .1
          }
        }
    }
   box{\n""")
						self.filename.write("      < -{}, -{}, center_z-0.01 >,< {}, {}, center_z+0.01 >".format(self.boxLength[0], self.boxLength[1], self.boxLength[0], self.boxLength[0])) 
						self.filename.write("""
        texture{
	      pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
	      finish{
	        ambient .2
	          diffuse .6
	          specular .1
	          roughness .1
	      }
	    }
     }
   }
#end\n""")
						return

					def printRing6_loc(self, setRings_C, red_val = "1.0", green_val = "0.0", blue_val = "0.0", zVal = 0.0):
						for j in range(len(setRings_C)):
							if (setRings_C[j][2] <= self.slab_thickness and setRings_C[j][2] >= -self.slab_thickness):
								self.red_val = red_val
								self.green_val = green_val
								self.blue_val = blue_val
								self.zVal = zVal
								self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_C[j][0], setRings_C[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))

								if (setRings_C[j][0] < (-0.5*self.boxLength[0]+self.buffer_size)):
									if (setRings_C[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]+self.boxLength[0]), (setRings_C[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]+self.boxLength[0]), setRings_C[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_C[j][0], (setRings_C[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_C[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]+self.boxLength[0]), (setRings_C[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]+self.boxLength[0]), setRings_C[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_C[j][0], (setRings_C[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]+self.boxLength[0]), setRings_C[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_C[j][0] > (0.5*self.boxLength[0]-self.buffer_size)):
									if (setRings_C[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]-self.boxLength[0]), (setRings_C[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]-self.boxLength[0]), setRings_C[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_C[j][0], (setRings_C[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_C[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]-self.boxLength[0]), (setRings_C[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]-self.boxLength[0]), setRings_C[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_C[j][0], (setRings_C[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_C[j][0]-self.boxLength[0]), setRings_C[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_C[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
									self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_C[j][0], (setRings_C[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_C[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
									self.filename.write("ring6(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_C[j][0], (setRings_C[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))

					def printRing7(self):
						self.filename.write("""\n#macro ring7 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)
  intersection{
    polygon{
      7,
      <center_x, center_y + 0.75, center_z>,
      <center_x - 0.5878, center_y + 0.4755, center_z>,
      <center_x - 0.9511, center_y - 0.1545, center_z>,
      <center_x - 0.5878, center_y - 0.7845, center_z>,
      <center_x + 0.5878, center_y - 0.7845, center_z>,
      <center_x + 0.9511, center_y - 0.1545, center_z>,
      <center_x + 0.5878, center_y + 0.4755, center_z>\n""")
						self.filename.write("""        texture{
          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
          finish{
            ambient .2
	          diffuse .6
	          specular .1
              roughness .1
          }
        }
    }
   box{\n""")
						self.filename.write("      < -{}, -{}, center_z-0.01 >,< {}, {}, center_z+0.01 >".format(self.boxLength[0], self.boxLength[1], self.boxLength[0], self.boxLength[0])) 
						self.filename.write("""
        texture{
	      pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
	      finish{
	        ambient .2
	          diffuse .6
	          specular .1
	          roughness .1
	      }
	    }
     }
   }
#end\n""")
						return

					def printRing7_loc(self, setRings_N, red_val = "1.0", green_val = "0.5", blue_val = "0.0", zVal = -0.2):
						for j in range(len(setRings_N)):
							if (setRings_N[j][2] <= self.slab_thickness and setRings_N[j][2] >= -self.slab_thickness):
								self.red_val = red_val
								self.green_val = green_val
								self.blue_val = blue_val
								self.zVal = zVal
								self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_N[j][0], setRings_N[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))

								if (setRings_N[j][0] < (-0.5*self.boxLength[0]+self.buffer_size)):
									if (setRings_N[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]+self.boxLength[0]), (setRings_N[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]+self.boxLength[0]), setRings_N[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_N[j][0], (setRings_N[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_N[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]+self.boxLength[0]), (setRings_N[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]+self.boxLength[0]), setRings_N[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_N[j][0], (setRings_N[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]+self.boxLength[0]), setRings_N[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_N[j][0] > (0.5*self.boxLength[0]-self.buffer_size)):
									if (setRings_N[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]-self.boxLength[0]), (setRings_N[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]-self.boxLength[0]), setRings_N[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_N[j][0], (setRings_N[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_N[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]-self.boxLength[0]), (setRings_N[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]-self.boxLength[0]), setRings_N[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_N[j][0], (setRings_N[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_N[j][0]-self.boxLength[0]), setRings_N[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_N[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
									self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_N[j][0], (setRings_N[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_N[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
									self.filename.write("ring7(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_N[j][0], (setRings_N[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))

					def printRing8(self):
						self.filename.write("""\n#macro ring8 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)
  intersection{
    polygon{
      8,
      <center_x, center_y+0.75, center_z>,
      <center_x-0.5303, center_y+0.5303, center_z>,
      <center_x-0.75, center_y, center_z>,
      <center_x-0.5303, center_y-0.5303, center_z>,
      <center_x, center_y-0.75, center_z>,
      <center_x+0.5303, center_y-0.5303, center_z>,
      <center_x+0.75, center_y, center_z>,
      <center_x+0.5303, center_y+0.5303, center_z>\n""")
						self.filename.write("""        texture{
          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
          finish{
            ambient .2
	          diffuse .6
	          specular .1
              roughness .1
          }
        }
    }
   box{\n""")
						self.filename.write("      < -{}, -{}, center_z-0.01 >,< {}, {}, center_z+0.01 >".format(self.boxLength[0], self.boxLength[1], self.boxLength[0], self.boxLength[0])) 
						self.filename.write("""
        texture{
	      pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
	      finish{
	        ambient .2
	          diffuse .6
	          specular .1
	          roughness .1
	      }
	    }
     }
   }
#end\n""")
						return

					def printRing8_loc(self, setRings_O, red_val = "0.0", green_val = "1.0", blue_val = "0.5", zVal = -0.2):
						for j in range(len(setRings_O)):
							if (setRings_O[j][2] <= self.slab_thickness and setRings_O[j][2] >= -self.slab_thickness):
								self.red_val = red_val
								self.green_val = green_val
								self.blue_val = blue_val
								self.zVal = zVal
								self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_O[j][0], setRings_O[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))

								if (setRings_O[j][0] < (-0.5*self.boxLength[0]+self.buffer_size)):
									if (setRings_O[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]+self.boxLength[0]), (setRings_O[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]+self.boxLength[0]), setRings_O[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_O[j][0], (setRings_O[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_O[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]+self.boxLength[0]), (setRings_O[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]+self.boxLength[0]), setRings_O[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_O[j][0], (setRings_O[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]+self.boxLength[0]), setRings_O[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_O[j][0] > (0.5*self.boxLength[0]-self.buffer_size)):
									if (setRings_O[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]-self.boxLength[0]), (setRings_O[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]-self.boxLength[0]), setRings_O[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_O[j][0], (setRings_O[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_O[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]-self.boxLength[0]), (setRings_O[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]-self.boxLength[0]), setRings_O[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_O[j][0], (setRings_O[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_O[j][0]-self.boxLength[0]), setRings_O[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_O[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
									self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_O[j][0], (setRings_O[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_O[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
									self.filename.write("ring8(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_O[j][0], (setRings_O[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))


					def printRing9(self):
						self.filename.write("""\n#macro ring9 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)
  intersection{
    polygon{
      9,
      <center_x, center_y+0.75, center_z>,
      <center_x-0.6495, center_y+0.375, center_z>,
      <center_x-0.9495, center_y-0.2165, center_z>,
      <center_x-0.6495, center_y-0.807, center_z>, 
      <center_x, center_y-1.125, center_z>,        
      <center_x+0.6495, center_y-0.807, center_z>, 
      <center_x+0.9495, center_y-0.2165, center_z>,
      <center_x+0.6495, center_y+0.375, center_z>,
      <center_x, center_y+0.75, center_z>\n""")
						self.filename.write("""        texture{
          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
          finish{
            ambient .2
	          diffuse .6
	          specular .1
              roughness .1
          }
        }
    }
   box{\n""")
						self.filename.write("      < -{}, -{}, center_z-0.01 >,< {}, {}, center_z+0.01 >".format(self.boxLength[0], self.boxLength[1], self.boxLength[0], self.boxLength[0])) 
						self.filename.write("""
        texture{
	      pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
	      finish{
	        ambient .2
	          diffuse .6
	          specular .1
	          roughness .1
	      }
	    }
     }
   }
#end\n""")
						return

					def printRing9_loc(self, setRings_F, red_val = "0.5", green_val = "1.0", blue_val = "0.0", zVal = -0.2):
						for j in range(len(setRings_F)):
							if (setRings_F[j][2] <= self.slab_thickness and setRings_F[j][2] >= -self.slab_thickness):
								self.red_val = red_val
								self.green_val = green_val
								self.blue_val = blue_val
								self.zVal = zVal
								self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_F[j][0], setRings_F[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))

								if (setRings_F[j][0] < (-0.5*self.boxLength[0]+self.buffer_size)):
									if (setRings_F[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]+self.boxLength[0]), (setRings_F[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]+self.boxLength[0]), setRings_F[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_F[j][0], (setRings_F[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_F[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]+self.boxLength[0]), (setRings_F[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]+self.boxLength[0]), setRings_F[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_F[j][0], (setRings_F[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]+self.boxLength[0]), setRings_F[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_F[j][0] > (0.5*self.boxLength[0]-self.buffer_size)):
									if (setRings_F[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]-self.boxLength[0]), (setRings_F[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]-self.boxLength[0]), setRings_F[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_F[j][0], (setRings_F[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_F[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]-self.boxLength[0]), (setRings_F[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]-self.boxLength[0]), setRings_F[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_F[j][0], (setRings_F[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									else:
										self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_F[j][0]-self.boxLength[0]), setRings_F[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_F[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
									self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_F[j][0], (setRings_F[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
								elif (setRings_F[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
									self.filename.write("ring9(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_F[j][0], (setRings_F[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))

					def printRing10(self):
						self.filename.write("""\n#macro ring10 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)
  intersection{
    polygon{
      10,
      <center_x, center_y+0.75, center_z>,
      <center_x-0.6125, center_y+0.5, center_z>,
      <center_x-0.75, center_y+0.1768, center_z>,
      <center_x-0.75, center_y-0.1768, center_z>,
      <center_x-0.6125, center_y-0.5, center_z>,
      <center_x, center_y-0.75, center_z>,
      <center_x+0.6125, center_y-0.5, center_z>,
      <center_x+0.75, center_y-0.1768, center_z>,
      <center_x+0.75, center_y+0.1768, center_z>,
      <center_x+0.6125, center_y+0.5, center_z>\n""")
						self.filename.write("""        texture{
          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
          finish{
            ambient .2
	          diffuse .6
	          specular .1
              roughness .1
          }
        }
    }
   box{\n""")
						self.filename.write("      < -{}, -{}, center_z-0.01 >,< {}, {}, center_z+0.01 >".format(self.boxLength[0], self.boxLength[1], self.boxLength[0], self.boxLength[0])) 
						self.filename.write("""
        texture{
	      pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }
	      finish{
	        ambient .2
	          diffuse .6
	          specular .1
	          roughness .1
	      }
	    }
     }
   }
#end\n""")
						return

					def printRing10_loc(self, setRings_Ne, red_val = "1.0", green_val = "0.0", blue_val = "0.5", zVal = -0.2):
							for j in range(len(setRings_Ne)):
								if (setRings_Ne[j][2] <= self.slab_thickness and setRings_Ne[j][2] >= -self.slab_thickness):
									self.red_val = red_val
									self.green_val = green_val
									self.blue_val = blue_val
									self.zVal = zVal
									self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Ne[j][0], setRings_Ne[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))

									if (setRings_Ne[j][0] < (-0.5*self.boxLength[0]+self.buffer_size)):
										if (setRings_Ne[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]+self.boxLength[0]), (setRings_Ne[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]+self.boxLength[0]), setRings_Ne[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Ne[j][0], (setRings_Ne[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										elif (setRings_Ne[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]+self.boxLength[0]), (setRings_Ne[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]+self.boxLength[0]), setRings_Ne[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Ne[j][0], (setRings_Ne[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										else:
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]+self.boxLength[0]), setRings_Ne[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_Ne[j][0] > (0.5*self.boxLength[0]-self.buffer_size)):
										if (setRings_Ne[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]-self.boxLength[0]), (setRings_Ne[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]-self.boxLength[0]), setRings_Ne[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Ne[j][0], (setRings_Ne[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										elif (setRings_Ne[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]-self.boxLength[0]), (setRings_Ne[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]-self.boxLength[0]), setRings_Ne[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Ne[j][0], (setRings_Ne[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
										else:
											self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%((setRings_Ne[j][0]-self.boxLength[0]), setRings_Ne[j][1], self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_Ne[j][1] < (-0.5*self.boxLength[1]+self.buffer_size)):
										self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Ne[j][0], (setRings_Ne[j][1]+self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))
									elif (setRings_Ne[j][1] > (0.5*self.boxLength[1]-self.buffer_size)):
										self.filename.write("ring10(%f, %f, %f, %s, %s, %s, %f)\n"%(setRings_Ne[j][0], (setRings_Ne[j][1]-self.boxLength[1]), self.zVal, red_val, green_val, blue_val, self.transparency))


				if povrayOutOpt:
					pov_out = open("%s/%s_%d.pov"%(path, fileName, frame_count), 'w')
					boxLength = [Lx, Ly, Lz]

					# calling a PovObjects class
					povray = PovObjects(pov_out, boxLength)
					povray.outFile()
					povray.printRing3()
					povray.printRing4()
					povray.printRing5()
					povray.printRing6()
					povray.printRing7()
					povray.printRing8()
					if max_ring > 8:
						povray.printRing9()
						if max_ring > 9:
							povray.printRing10()

					povray.printRing3_loc(setRings_Li)
					povray.printRing4_loc(setRings_Be)
					povray.printRing5_loc(setRings_B)
					povray.printRing6_loc(setRings_C)
					povray.printRing7_loc(setRings_N)
					povray.printRing8_loc(setRings_O)
					if max_ring > 8:
						povray.printRing9_loc(setRings_F)
						if max_ring > 9:
							povray.printRing10_loc(setRings_Ne)

				if (povrayOutOpt):
					pov_imageX = str(int(povray.scale_factor * povray.boxLength[0]))
					pov_imageY = str(int(povray.scale_factor * povray.boxLength[1]))
					pov_info.write("povray -w" + pov_imageX + " -h" + pov_imageY + " +a0.1 -D " + fileName + "_" + str(frame_count) + ".pov\n")


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

				hbond_ndx = [[] for _ in range(int(nMols))]
				hbondIndex = []
				hbonds = 0
		
				# reinitialize variables
				oPosX = []
				oPosY = []
				oPosZ = []
				wat_loop = []
				primitiveRings = []
				minimal_rings = []
				water = [[0 for _ in range(9)] for _ in range(int(nMols))]
				pos = 0
				ref_ndx = 0
				# rings segregation
				ring3 = []
				ring4 = []
				ring5 = []
				ring6 = []
				ring7 = []
				ring8 = []
				ring9 = []
				ring10 = []
			else:
				continue
	outFile.close()
	if ring_traj:
		outFile2.close()
	if binning_rings:
		outFile3.close()
	if povrayOutOpt:
		pov_info.close()
		pov_out.close()

	#long process here and end of animation
	time.sleep(5)
	done = True			# use this if used 'done' in the function
	time.sleep(1)

	if directionality:
		print("\nThe closed directional rings have been output to 'ringsCount.dat\n")
	else:
		print("\nThe closed rings have been output to 'ringsCount.dat\n")
	if ring_traj:
		print("The location of the closed rings has been written to 'rings_location.xyz'\n")
	if povrayOutOpt:
		print("""POVRay rendering command written to pov_files/conf_pov.txt
       ...potentially useful for rendering files in the NEWLY made pov_files directory\n""")
	if binning_rings:
		print("The binning of the trajectory has been written to 'rings_dist.dat'")

def main(argv):
	global directionality, hbond_energy, ring_traj, binning_rings, povrayOutOpt

	directionality = False
	hbond_energy = False
	ring_traj = False
	binning_rings = False
	povrayOutOpt = False
	_haveinputFileName = 0
	max_ring = 0
	ring_closure = 0
	algorithm = "hbondAngle"
	try:
		opts, args = getopt.getopt(argv, "hexbpdf:r:c:m:", ["help", "energy_defn", "ring_traj", "binning_rings", "povrayOutOpt", "directional_rings", "input-file=", "ring_size=", "minimal_ring=", "method_algorithm="])	# 'hed' do not take arguments, 'frcm' take argument
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-e", "--energy_defn"):
			hbond_energy = True
		elif opt in ("-x", "--ring_traj"):
			ring_traj = True
		elif opt in ("-b", "--binning_rings"):
			binning_rings = True
			ring_traj = True		# requires ring trajectory file
		elif opt in ("-p", "--povrayOutOpt"):
			povrayOutOpt = True
			ring_traj = True		# rendering povray also requires ring trajectory file
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

