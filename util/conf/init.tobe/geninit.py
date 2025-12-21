#!/usr/bin/env python3

import math, sys, argparse
from tqdm import tqdm
import numpy as np
import ase
from ase.io import read,write
from decimal import Decimal

def vector_to_sequential_id(indices, vprocs):
	return indices[0] + indices[1]*vprocs[0] + indices[2]*vprocs[0]*vprocs[1]

def element_mapper(elem_map_input_file):

	with open(elem_map_input_file, 'r') as fin:
		lines = fin.read().split('\n')
		data = [line.split() for line in lines]

	elem_map = {}
	for d in data:
		if len(d) == 0: continue
		if d[0] == 'allegro':
			for elem_id, elem_name in enumerate(d[3::2]):
				elem_map[elem_name] = elem_id + 1 # 1-indexed in Fortran
		
	print('element map: ', elem_map)

	return elem_map


class geninit:

	def __init__(self, inputxyz, repeat, vprocs, elem_map, outputbin="rxff.bin", outputxyz="output.xyz"):

		self.atoms = ase.io.read(inputxyz).repeat(repeat)
		self.repeat = repeat
		self.cell_l = self.atoms.cell.lengths()
		self.cell_a = self.atoms.cell.angles()
		self.elem_map = elem_map

		self.vprocs = vprocs
		self.nprocs = math.prod(vprocs)

		minxyz = [1e9]*3
		maxxyz = [-1e9]*3
		
		for d in self.atoms.get_scaled_positions():
			for i, p in enumerate(d): 
				if minxyz[i] > p: minxyz[i] = p
				if maxxyz[i] < p: maxxyz[i] = p

		shifted_positions = []
		for d in self.atoms.get_scaled_positions():
			fractional_pos = [(Decimal(p)-Decimal(m))%Decimal(1.0) for m, p in zip(minxyz, d)]
			shifted_positions.append(np.array(fractional_pos, dtype=np.float64))

		per_node_atoms = {}

		fractional_domain_dims = [1.0/v for v in self.vprocs]

		for gid, (sd, elem) in tqdm(enumerate(zip(shifted_positions, self.atoms.get_chemical_symbols()))):
			vprocs_per_atom = [int(d/fdd) for d,fdd in zip(sd, fractional_domain_dims)]
			sid = vector_to_sequential_id(vprocs_per_atom, self.vprocs)

			##NOTE## the global ID format used in RXMD
			elem_id = elem_map[elem] + gid*1e-13

			#print(elem, sd, elem_map[elem], elem_id)

			local_coords = [sd[i] - vprocs_per_atom[i]*fractional_domain_dims[i] for i in range(len(sd))]

			if sid not in per_node_atoms: per_node_atoms[sid] = []
			per_node_atoms[sid].append([local_coords, elem_id])

		#dump current system for inspection
		write(outputxyz, self.atoms)

		#print('per_node_atoms len ', len(per_node_atoms))

		with open(outputbin, 'wb') as fbin:
			np.array([self.nprocs], dtype=np.int32).tofile(fbin)
			np.array(self.vprocs, dtype=np.int32).tofile(fbin)

			pna_array = np.array([len(per_node_atoms[i]) for i in range(len(per_node_atoms))])
			np.array(pna_array, dtype=np.int32).tofile(fbin)

			np.array([0], dtype=np.int32).tofile(fbin) # current_step. initialized with 0

			np.array(self.cell_l, dtype=np.float64).tofile(fbin)
			np.array(self.cell_a, dtype=np.float64).tofile(fbin)

			for sid in tqdm(range(len(per_node_atoms))):
				for atom in per_node_atoms[sid]:
					p = atom[0]
					e = atom[1]
					# rr(3), vv(3), qq, dtype, qfsp0, qfsv0
					atom = np.array([p[0],p[1],p[2],0,0,0,0,e,0,0],dtype=np.float64)
					atom.tofile(fbin)
					
					#print(sid, atom)

		return 


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputxyz", default="input.xyz", 
                	help="Input file name. must be in extended XYZ format.")
parser.add_argument("-r", "--repeat", nargs=3, type=int, default=[1,1,1], 
	                help="Number of repetition of unit cell.")
parser.add_argument("-v", "--vprocs", nargs=3, type=int, default=[1,1,1], 
	                help="Number of MD domains.")
parser.add_argument("-em", "--element_map", type=str, default=None,
	                help="e.g. --element_mapping rxmdnn.in.")
parser.add_argument("-oxyz", "--outputxyz", default="output.xyz", 
	                help="Output XYZ file name. must be extended XYZ format.")
parser.add_argument("-o", "--outputbin", default="rxff.bin", 
	                help="Output binary file name. used to start RXMD simulation.")

args = parser.parse_args()

print(args)

elem_map = element_mapper(args.element_map)

g = geninit(
	inputxyz = args.inputxyz, 
	repeat = args.repeat,
	vprocs = args.vprocs, 
	elem_map = elem_map, 
	outputxyz = args.outputxyz,
	outputbin = args.outputbin,
	)

