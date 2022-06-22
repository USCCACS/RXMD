import os,sys,glob
import numpy as np
import pickle
import torch
import torch.nn as nn
import torch.nn.functional as F

ALL_ELEMENTS = ["Pb", "Ti", "O"]
type2int = {}
for i, e in enumerate(ALL_ELEMENTS):
	type2int[e] = i

RS=[1.0, 2.0, 3.0, 4.0, 5.0]
ETA=[1.0, 1.5, 2.0]
#RS=[1.0, 2.0]
#ETA=[2.0]
feature_size = len(RS)*len(ETA)

print(f'feature_size: {RS} {ETA}')

def get_mdframe(filename):

	mdframe = {}
	with open(filename, 'r') as fin:

		# read total energy
		data = fin.readline().split()
		total_energy = float(data[4])

		# skip two lines
		data = fin.readline()
		data = fin.readline()

		# get H-matrix
		lattice = []
		for i in range(3): 
			lattice.append([float(f) for f in fin.readline().split()])

		# skip one line
		data = fin.readline()

		# get number of atoms
		num_atoms = int(fin.readline().split()[0])

		#print(f'{filename} {num_atoms} {total_energy} {lattice}')

		mdframe['meta'] = {'filename':filename, 'num_atoms':num_atoms, 'total_energy':total_energy, 'lattice':lattice}
		mdframe['element'] = {}
		mdframe['position'] = {}
		mdframe['force'] = {}

		for n in range(num_atoms):
			data = fin.readline().split()
			element = data[0]
			position = [float(f) for f in data[1:4]]
			force = [float(f) for f in data[4:8]]
			mdframe['element'][n] = element
			mdframe['position'][n] = position
			mdframe['force'][n] = force

		for n in range(num_atoms):
			mdframe['position'][n] = torch.tensor(mdframe['position'][n], dtype=torch.float, requires_grad=True)
			mdframe['force'][n] = torch.tensor(mdframe['force'][n], dtype=torch.float, requires_grad=True)

	return mdframe

def pbc(dr, lattice):

	for ia in range(3):
		if dr[ia] >= 0.5*lattice[ia][ia]:
			dr[ia] -= lattice[ia][ia]
		elif dr[ia] <- 0.5*lattice[ia][ia]: 
			dr[ia] += lattice[ia][ia]

def train(mdframe, nets, optimizers):

	for opt in optimizers.values():
		opt.zero_grad()

	num_atoms = mdframe['meta']['num_atoms']
	elems = mdframe['element']
	pos = mdframe['position']
	ref_force = mdframe['force']
	ref_energy = torch.Tensor([mdframe['meta']['total_energy']])

	features = {}
	drtensors = {}
	gdrtensors = {}

	total_energy = torch.zeros(1)
	for i in range(num_atoms):
		dr_vec = []
		jtype_vec = []
		for j in range(num_atoms):
			if i == j: continue
			#dr = [pos[i][a] - pos[j][a] for a in range(3)]
			dr = pos[i] - pos[j]
			pbc(dr,mdframe['meta']['lattice'])
			drr = torch.linalg.norm(dr[0:3])

			dr_vec.append(drr)
			jtype_vec.append(elems[j])
			#print(i,j,dr,type2int[elems[j]])

		#dgdrvec = [None]*feature_size*len(type2int)

		#drtensor = torch.Tensor(dr_vec).requires_grad_()
		drtensor = dr_vec
		drtensors[i] = drtensor

		g2vec = [None]*feature_size*len(type2int)
		idx = 0
		for rs in RS:
			for eta in ETA:

				for jtype in ALL_ELEMENTS:

					#drtensor = torch.Tensor(dr_vec).requires_grad_()
					g2 = torch.zeros(1)

					for ib, dr in enumerate(drtensor):
						if jtype == jtype_vec[ib]:
							diff = dr - rs
							g2 = g2 + torch.exp(-eta*diff*diff)
					#g2.backward()

					index = feature_size*type2int[jtype] + idx


					#dgdrvec[index] = drtensor.grad
					g2vec[index]  = g2

				idx += 1

		features[i] = torch.cat(g2vec)

		itype = mdframe['element'][i]
		energy = nets[itype](features[i])
		total_energy += energy

	total_energy.backward(retain_graph=True)

	force_loss = torch.zeros([1])
	for i in range(num_atoms):
		diff = -pos[i].grad - ref_force[i]
		force_loss += torch.dot(diff, diff)
	force_loss /= num_atoms/3.0

	energy_loss = mse_loss(total_energy, ref_energy) 

	print('totalE,refE,Eloss,Floss: ', total_energy, ref_energy, energy_loss, force_loss)
	loss = energy_loss + force_loss
	loss.backward()	

	total_loss = loss.item()

	for opt in optimizers.values():
		opt.step()

	return total_loss

class NeuralNetwork(nn.Module):
	def __init__(self, input_size, output_size):
		super().__init__()
		self.fc1 = nn.Linear(input_size, 20, bias=True)
		self.fc2 = nn.Linear(20, 20, bias=True)
		self.fc3 = nn.Linear(20, output_size, bias=True)

	def forward(self, x):
		x = torch.tanh(self.fc1(x))
		x = torch.tanh(self.fc2(x))
		return self.fc3(x)

learning_rate = 1e-6

nets = {}
optimizers = {}
for e in ALL_ELEMENTS:
	nets[e] = NeuralNetwork(feature_size*len(type2int),1)
	optimizers[e] = torch.optim.SGD(nets[e].parameters(), lr=learning_rate)

mse_loss = F.mse_loss

dirname = sys.argv[1]

filenames = glob.glob(dirname+"/*.xsf")
filenames = filenames[:1]
mdframe = get_mdframe(filenames[0])

features_all = {}
drtensors_all = {}
gdrtensors_all = {}
num_epochs = 500
threshold = 1e-3

for epoch in range(num_epochs):

	total_energies = []
	ref_energies = []

	loss = 0.0
	for filename in filenames:
		mdframe = get_mdframe(filename)
		loss += train(mdframe,  nets, optimizers)

	print(f'epoch, loss : {epoch} {loss:8.3e}')

	feature = torch.rand(feature_size*len(type2int))
	for itype, net in nets.items():
		outfile =f'net_{itype}.pt'
		print(f'saving models {outfile}')
		traced_script_module = torch.jit.trace(net, feature)
		traced_script_module.save(outfile)

	if loss < threshold: break


