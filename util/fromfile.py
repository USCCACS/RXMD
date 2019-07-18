import struct
import numpy as np
import sys

# class to store model profile
class ModelProf:

	def __init__(self,profs):
	
		self.rad_types = profs[0]
		self.ang_types = profs[1]
		self.num_features_rad = profs[2]
		self.num_features_ang = profs[3]
		self.num_features_total = profs[4]
		self.feature_ptr_rad = profs[5]
		self.feature_ptr_ang = profs[6]


# binary file from rxmd code
filename = sys.argv[1]
fin = open(filename, 'rb')

# first two integers, number of atoms & number of models
num_atoms = np.fromfile(fin, dtype='int32', count=1)[0]
num_models = np.fromfile(fin, dtype='int32', count=1)[0]
print('num_atoms: ', num_atoms, ', num_models: ', num_models)

num_features = []
model_profs = []
for i in range(num_models):
	rad_types = np.fromfile(fin, dtype='int32', count=1)[0]
	ang_types = np.fromfile(fin, dtype='int32', count=1)[0]
	num_features_rad = np.fromfile(fin, dtype='int32', count=1)[0]
	num_features_ang = np.fromfile(fin, dtype='int32', count=1)[0]
	num_features_total = np.fromfile(fin, dtype='int32', count=1)[0]
	
	feature_ptr_rad = []
	for i in range(rad_types):
		feature_ptr_rad.append(np.fromfile(fin, dtype='int32', count=1)[0]-1)
	feature_ptr_ang = []
	for i in range(ang_types):
		feature_ptr_ang.append(np.fromfile(fin, dtype='int32', count=1)[0]-1)

	mp = ModelProf([rad_types,ang_types,num_features_rad,num_features_ang,num_features_total,feature_ptr_rad,feature_ptr_ang])
	model_profs.append(mp)

	num_features.append(num_features_total)

for i in range(num_models):
	print('---------------------------------')
	print('model Id: ', i)
	print('rad&ang type counts: ', model_profs[i].rad_types, model_profs[i].ang_types)
	print('rad&ang feature sizes: ', model_profs[i].num_features_rad, model_profs[i].num_features_ang)
	print('total feature sizes: ', model_profs[i].num_features_rad + model_profs[i].num_features_ang)
	print('rad feature ptr: ', model_profs[i].feature_ptr_rad)
	print('ang feature ptr: ', model_profs[i].feature_ptr_ang)
	print('---------------------------------')

features = {}
for i in range(num_atoms):
	gid = np.fromfile(fin, dtype='int32', count=1)[0] - 1 # change to zero-indexed
	ity = np.fromfile(fin, dtype='int32', count=1)[0] - 1 # change to zero-indexed

	#print('gid,ity,num_features[ity]: ',gid,ity,num_features[ity])

	feature_x = np.fromfile(fin, dtype='float32', count=num_features[ity])
	feature_y = np.fromfile(fin, dtype='float32', count=num_features[ity])
	feature_z = np.fromfile(fin, dtype='float32', count=num_features[ity])
	features[gid] = [feature_x,feature_y,feature_z]

sys.exit(0)
