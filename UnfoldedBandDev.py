#!/usr/bin/python2.6

import linecache
import os
import numpy as np
import matplotlib.pyplot as plt

# Setting some parameters for data-preprocessing
Height_BoVB = 13		## The energy height (bottom of valence band as the reference
				# points) below which the deviation calculation includes
intensity_threshold = 0.090	## Threshold for sampling (keeping) intensity peaks:
				# the spectral points with intensity less than maximum intensity
				# times the intensity_threshold will be discarded (treated as
				# numerical noise)
file_SC_BandUP = "Data.dat"	## Path to the BandUP output file.

# Do not change these
intensity_maximum = 0		# A random small value
energy_minimum = 100		# A random large value

# Removing the energy grid with zero intensity
with open(file_SC_BandUP) as f:
    content = f.readlines()
indices = 0, 1, 2, 3, 4, 5
content = [i for j, i in enumerate(content) if j not in indices]
content = [x.split() for x in content]
content = [[float(x) for x in y] for y in content]
content = np.asarray(content)
nonzero = []
for row in content:
    if row[2]!=0:
        nonzero.append(row)
nonzero = np.asarray(nonzero)
nonzero = np.delete(nonzero, [3,4,5,6,7],1)
l = nonzero.shape[0]
index = 0
tmp = -1
for i in range(l):
    if nonzero[i][0] != tmp:
        tmp = nonzero[i][0]
        index += 1
    nonzero[i][0] = index

# Labelling k points (index starts from 1)
l = nonzero.shape[0]
index = 0
tmp = 100
for i in range(l):
    if nonzero[i][1] < tmp:	# If the eigenvalue is smaller than the previous eigenvalue, label it a new k point
        index += 1
    nonzero[i][0] = index
    tmp = nonzero[i][1]

# Draw the supercell spectrum
tmp = np.delete(nonzero, [2], 1)
tmp = np.transpose(tmp)
plt.xticks(np.arange(0, 15, 1.0))
plt.plot(tmp[0],tmp[1],'+')
plt.xlabel('k point')
plt.ylabel('Energy (eV)')
plt.title('Supercell Spectrum')
plt.savefig('SC_raw.png')

# Find the max of intensity and the min of energy
for line in nonzero:
    if line[2] > intensity_maximum:
        intensity_maximum = line[2]
    if line[1] < energy_minimum:
        energy_minimum = line[1]
keep1 = []
# Use the BoVB as the reference point (E=0)
for line in nonzero:
    line[1] = line[1] - energy_minimum
    if line[1] <= Height_BoVB:
        keep1.append(line)
keep1 = np.asarray(keep1)
# Eliminate the points with low intensity
keep2 = []
for line in keep1:
    if line[2] >= intensity_threshold*intensity_maximum:		# Setting the desired magnitude of sampling
        keep2.append(line)
keep2 = np.asarray(keep2)

tmp = np.delete(keep2, [2], 1)
tmp = np.transpose(tmp)
plt.xticks(np.arange(0, 15, 1.0))
plt.plot(tmp[0],tmp[1],'+')
plt.xlabel('k point')
plt.ylabel('Energy (eV)')
plt.title('Supercell Spectrum (sampled and aligned)')
plt.savefig('SC_refined.png')

# Getting the PC k points from the IBZKPT of PC calculation
Nk = int(linecache.getline("IBZKPT",2))
print("Number of KPOINTS for PC: "+str(Nk))
kpts = []
for i in range(Nk):
    kpts.append(linecache.getline("IBZKPT", 4+i))
kpts = [x.split() for x in kpts]
kpts = [[float(x) for x in y] for y in kpts]
kpts = np.asarray(kpts)

# Getting the PC spectrum from the EIGENVAL of PC calculation
[_, _, Nvalues] = linecache.getline("EIGENVAL",6).split()
Nvalues = int(Nvalues)
print("Number of Eigenvalues per KPOINT: "+str(Nvalues))
start = 9
eigenvalues = np.empty([Nk, Nvalues])
for i in range(Nk):
    this_k = []
    for j in range(Nvalues):
        this_k.append(linecache.getline("EIGENVAL",start+i*(Nvalues+2)+j))
    this_k = [x.split() for x in this_k]
    this_k = [[float(x) for x in y] for y in this_k]
    this_k = np.asarray(this_k)
    this_k = np.delete(this_k,0,1)
    this_k = np.transpose(this_k)
    eigenvalues[i,:] = this_k

w = []
for i in range(Nk):
    for v in eigenvalues[i,:]:
        w.append([i+1,v])
w = np.asarray(w)
minimum = 100
for i in w:
    if i[1] < minimum:
        minimum = i[1]  
print("min="+str(minimum))
out = []
for i in w:
    i[1] = i[1] - minimum
    if i[1] <= Height_BoVB:
        out.append(i)
out = np.asarray(out)

# Plotting the PC spectrum
tmp = np.transpose(out)
plt.xticks(np.arange(0, Nk+1, 1.0))
plt.plot(tmp[0], tmp[1], 'o')
plt.xlim([0,Nk+1])
plt.title('Primitive cell spectrum')
plt.xlabel('k point')
plt.ylabel('Energy (eV)')
plt.savefig('PC.png')

# Getting the k points and their corresponding weights
Nk = int(linecache.getline("IBZKPT",2))
kpts = []
number = []
for i in range(Nk):
    number.append(i+1)
    kpts.append(linecache.getline("IBZKPT", 4+i))
kpts = [x.split() for x in kpts]
kpts = [[float(x) for x in y] for y in kpts]
kpts = np.asarray(kpts)
kpts = np.delete(kpts, [0,1,2],1)
kpts = np.reshape(kpts,(1,-1))
number = np.reshape(np.asarray(number),(1,-1))
kpts = np.transpose(np.concatenate((number,kpts),axis=0))

# Getting rid of the degeneracies
P_dat = out
store = []
tmp = -100
for row in P_dat:
    if np.absolute(row[1]-tmp) > 0.02:	# Adjacent eigenvalues with distance less than 0.02
					# eV will be treated as degeneracy.
        tmp = row[1]
        store.append(row)
P_dat = np.asarray(store)
del store

# P_dat pre-pocessing
tmp = 0
nlist = []
index = 1
for row in P_dat:
    if row[0] != tmp:
        tmp = row[0]
        index = 1
    else:
        index += 1
    nlist.append(index)
z = np.zeros(len(nlist))
z = np.reshape(z,(1,-1))
nlist = np.reshape(np.asarray(nlist),(1,-1))
P_dat = np.transpose(P_dat)
P_dat = np.transpose(np.concatenate((np.reshape(P_dat[0],(1,-1)),nlist,np.reshape(P_dat[1],(1,-1)),z,z,z),axis=0))
del nlist

target_dat = keep2
target_dat = np.delete(target_dat, [2],1)	# Getting rid of the intensity data

"""
# Searching for "peaks"
tmp = -100
peaks = []
for row in target_dat:
    if row[1]-tmp == 0.05:
        peaks.append(row[1])
    tmp = row[1]
peaks = np.asarray(peaks)
"""

# Assigning the unfolded peaks to PC k points
for row in target_dat:
    k = row[0]			# k point index
    e = row[1]			# energy
    for i in range(P_dat.shape[0]):
        if P_dat[i][0] == k:	# Searching for the corresponging k point, store the index in i
            break
    dist = np.absolute(e-P_dat[i][2])
    while (i+1<P_dat.shape[0]) and (P_dat[i+1][0] == k) and (np.absolute(e-P_dat[i+1][2]) < dist):
        dist = np.absolute(e-P_dat[i+1][2])
        i += 1
    P_dat[i][3] += 1
    P_dat[i][4] += dist ** 2
for row in P_dat:
    row[5] = row[4]/row[3]	# Averaged over every eigen energy

# Computing the total average
total_delta_E = 0
# total_number = 1	# Includes a degeneracy
for row in P_dat:
    number = kpts[row[0]-1][1]
    total_number += number
    total_delta_E += number * row[5]
# total_delta_E += P_dat[3][5]	# Degeneracy
average = (total_delta_E/total_number) ** 0.5
print("File: "+file_SC_BandUP)
print("Average variation per eigenvalue: "+str(average))
