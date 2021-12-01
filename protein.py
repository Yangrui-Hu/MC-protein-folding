import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numba import jit
import random

matplotlib.rcParams.update({'font.size': 15})

@jit(nopython=True)
def step_2D_jit():
	d = np.random.randint(1,5)
	if d == 1:
		x_s, y_s = 1, 1
	if d == 2:
		x_s, y_s = 1, -1
	if d == 3:
		x_s, y_s = -1, 1
	if d == 4:
		x_s, y_s = -1, -1
	return x_s, y_s

@jit(nopython=True)
def neighbor(x, y, x_test, y_test):
	if (x-x_test)**2+(y-y_test)**2 < 1.1:
		return 1
	if (x-x_test)**2+(y-y_test)**2 > 1.1:
		return 0

@jit(nopython=True)
def Calculate_delta_Energy_jit(N_aminoacids, R, J, A, index, xnew, ynew, temp):

	k = index
	x_new = xnew
	y_new = ynew
	R_new = np.zeros((N_aminoacids, 2))
	R_before = np.zeros((N_aminoacids, 2))
	T = temp
	beta = 1.0/T

	for i in range(N_aminoacids):
		R_new[i,0] = R[i,0]
		R_new[i,1] = R[i,1]
		R_before[i,0] = R[i,0]
		R_before[i,1] = R[i,1]
		if i == k :
			R_new[i,0] = x_new
			R_new[i,1] = y_new


	energy_before = 0
	energy_new = 0
	distance_before = 0
	distance_new = 0

	for i in range(N_aminoacids):
		for j in range(i):
			if i-j > 1.1:
				distance_before = (R_before[i,0]-R_before[j,0])**2 + (R_before[i,1]-R_before[j,1])**2
				distance_new = (R_new[i,0]-R_new[j,0])**2 + (R_new[i,1]-R_new[j,1])**2
				if distance_before < 1.1:
					energy_before = energy_before - J[A[i],A[j]]
				if distance_new < 1.1:
					energy_new = energy_new - J[A[i],A[j]]
		
	delta_energy = energy_new - energy_before
	if delta_energy <= 0:
		R[k,0] = x_new
		R[k,1] = y_new
	elif delta_energy > 0:
		p = random.random()
		if p < np.exp(-beta*delta_energy):
			R[k,0] = x_new
			R[k,1] = y_new
	return R

@jit(nopython=True)
def ProteinFolding_jit(N_aminoacids, J, A, Temp = 10, Nsteps = 1000000):

	R = np.zeros((N_aminoacids, 2))
	# Initialize the protein chain, R is the position array
	for i in range(N_aminoacids):
		R[i,0] = i # x_i
		R[i,1] = 0 # y_i

	FreeEnergy = []
	length = []
	step = []

	for i in range(Nsteps):
		# choose k-th amino acid randomly and make the movement
		k = np.random.randint(0, N_aminoacids)
		x = R[k,0]
		y = R[k,1]
		x_s, y_s = step_2D_jit()
		x_new = x + x_s
		y_new = y + y_s

		#judge whether the walking is self-avoiding
		SAW = 1
		for j in range(N_aminoacids):
			if (R[j,0]-x_new)**2+(R[j,1]-y_new)**2 < 0.1:
				SAW = 0
		if SAW == 0 : 
			continue

		#judge whether the length is not changed
		if k == 0:
			x_next = R[k+1,0]
			y_next = R[k+1,1]
			nbor = neighbor(x_new, y_new, x_next, y_next)
			if nbor == 1:
				R = Calculate_delta_Energy_jit(N_aminoacids, R, J, A, k, x_new, y_new, Temp)
			if nbor == 0:
				continue

		if k == N_aminoacids-1:
			x_pre = R[k-1,0]
			y_pre = R[k-1,1]
			nbor = neighbor(x_new, y_new, x_pre, y_pre)
			if nbor == 1:
				R = Calculate_delta_Energy_jit(N_aminoacids, R, J, A, k, x_new, y_new, Temp)
			if nbor == 0:
					continue

		if 0 < k < N_aminoacids-1:
			x_next = R[k+1,0]
			y_next = R[k+1,1]
			x_pre = R[k-1,0]
			y_pre = R[k-1,1]
			nbor_next = neighbor(x_new, y_new, x_next, y_next)
			nbor_pre = neighbor(x_new, y_new, x_pre, y_pre)
			if nbor_next == 1 and nbor_pre == 1:
				R = Calculate_delta_Energy_jit(N_aminoacids, R, J, A, k, x_new, y_new, Temp)
			if nbor_pre == 0:
				continue
			if nbor_next == 0:
				continue

		# calculate the free energy, r square
		F = 0
		for j in range(N_aminoacids):
			for l in range(j):
				if j-l > 1.1:
					distance = (R[j,0]-R[l,0])**2 + (R[j,1]-R[l,1])**2
					if distance < 1.1:
						F = F - J[A[j],A[l]]

		r2 = (R[N_aminoacids-1,0]-R[0,0])**2 + (R[N_aminoacids-1,1]-R[0,1])**2		
		step.append(i+1)
		FreeEnergy.append(F)
		length.append(np.sqrt(r2))

	return step, FreeEnergy, length, R

@jit(nopython=True)
def ProteinFolding_low_temp_jit(N_aminoacids, J, A, Nsteps = 1000000):

	R = np.zeros((N_aminoacids, 2))
	# Initialize the protein chain, R is the position array
	for i in range(N_aminoacids):
		R[i,0] = i # x_i
		R[i,1] = 0 # y_i

	FreeEnergy = []
	length = []
	step = []

	for i in range(Nsteps):

		if 0 <= i < 2000000:
			temp = 10.0
		if 2000000 <= i < 2500000:
			temp = 9.5
		if 2500000 <= i < 3000000:
			temp = 9.0
		if 3000000 <= i < 3500000:
			temp = 8.5
		if 3500000 <= i < 4000000:
			temp = 8.0
		if 4000000 <= i < 4500000:
			temp = 7.5
		if 4500000 <= i < 5000000:
			temp = 7.0
		if 5000000 <= i < 5500000:
			temp = 6.5
		if 5500000 <= i < 6000000:
			temp = 6.0
		if 6000000 <= i < 6500000:
			temp = 5.5
		if 6500000 <= i < 7000000:
			temp = 5.0
		if 7000000 <= i < 7500000:
			temp = 4.5
		if 7500000 <= i < 8000000:
			temp = 4.0
		if 8000000 <= i < 8500000:
			temp = 3.5
		if 8500000 <= i < 9000000:
			temp = 3.0
		if 9000000 <= i < 9500000:
			temp = 2.5
		if 9500000 <= i < 10000000:
			temp = 2.0
		if 10000000 <= i < 10500000:
			temp = 1.5
		if 10500000 <= i < 11000000:
			temp = 1.0
		if 11000000 <= i < 11500000:
			temp = 0.5
			
		# choose k-th amino acid randomly and make the movement
		k = np.random.randint(0, N_aminoacids)
		x = R[k,0]
		y = R[k,1]
		x_s, y_s = step_2D_jit()
		x_new = x + x_s
		y_new = y + y_s

		#judge whether the walking is self-avoiding
		SAW = 1
		for j in range(N_aminoacids):
			if (R[j,0]-x_new)**2+(R[j,1]-y_new)**2 < 0.1:
				SAW = 0
		if SAW == 0 : 
			continue

		#judge whether the length is not changed
		if k == 0:
			x_next = R[k+1,0]
			y_next = R[k+1,1]
			nbor = neighbor(x_new, y_new, x_next, y_next)
			if nbor == 1:
				R = Calculate_delta_Energy_jit(N_aminoacids, R, J, A, k, x_new, y_new, temp)
			if nbor == 0:
				continue

		if k == N_aminoacids-1:
			x_pre = R[k-1,0]
			y_pre = R[k-1,1]
			nbor = neighbor(x_new, y_new, x_pre, y_pre)
			if nbor == 1:
				R = Calculate_delta_Energy_jit(N_aminoacids, R, J, A, k, x_new, y_new, temp)
			if nbor == 0:
				continue

		if 0 < k < N_aminoacids-1:
			x_next = R[k+1,0]
			y_next = R[k+1,1]
			x_pre = R[k-1,0]
			y_pre = R[k-1,1]
			nbor_next = neighbor(x_new, y_new, x_next, y_next)
			nbor_pre = neighbor(x_new, y_new, x_pre, y_pre)
			if nbor_next == 1 and nbor_pre == 1:
				R = Calculate_delta_Energy_jit(N_aminoacids, R, J, A, k, x_new, y_new, temp)
			if nbor_pre == 0:
				continue
			if nbor_next == 0:
				continue

		# calculate the free energy, r square
		F = 0
		for j in range(N_aminoacids):
			for l in range(j):
				if j-l > 1.1:
					distance = (R[j,0]-R[l,0])**2 + (R[j,1]-R[l,1])**2
					if distance < 1.1:
						F = F - J[A[j],A[l]]

		r2 = (R[N_aminoacids-1,0]-R[0,0])**2 + (R[N_aminoacids-1,1]-R[0,1])**2		
		step.append(i+1)
		FreeEnergy.append(F)
		length.append(np.sqrt(r2))

	return temp, step, FreeEnergy, length, R


def initialize(N_aminoacids):
	A = np.zeros(N_aminoacids, int)
	J = np.zeros((20,20))
	# Set the sequence of the protein chain A[N_aminoacids]
	for i in range(N_aminoacids):
		A[i] = np.random.randint(0,20)

	# Set the interaction energy matrix J
	for i in range(20):
		for j in range(i):
			J[i,j] = 2 * random.random() + 2
			J[j,i] = J[i,j]
		J[i,i] = random.random()

	return A, J

def initialize_Uniform_J(N_aminoacids):
	A = np.zeros(N_aminoacids, int)
	J = np.zeros((20,20))+ 5
	for i in range(N_aminoacids):
		A[i] = np.random.randint(0,20)
	return A, J

def plot_folding_structure(N_aminoacids, R):

	for i in range(N_aminoacids-1):
		x1, y1 = R[i,0], R[i,1]
		x2, y2 = R[i+1,0], R[i+1,1]
		x = [x1, x2]
		y = [y1, y2]
		plt.plot(x,y,'r')
	for i in range(N_aminoacids):
		plt.scatter(R[i,0], R[i,1], c ='k')
	plt.xlabel('x', fontsize=20)
	plt.ylabel('y', fontsize=20)
	plt.title('Final Structure (T = 10, 1000000 steps)', fontsize=20)
	plt.show()
		

def equilibrium(N_aminoacids):

	Temp = 1
	A, J = initialize(N_aminoacids)
	
	for n in range(10):
		step, FreeEnergy, length, R = ProteinFolding_jit(N_aminoacids, J, A, Temp, 6000000)
		#t, step, FreeEnergy, length, R = ProteinFolding_low_temp_jit(N_aminoacids, J, A, 6000000)
		file = open('full_equilibrium_%d_T%.1f-beforeSA-%d.txt'%(N_aminoacids, Temp, n), 'w')
		for i in range(len(step)):
			file.write(str(step[i]) + " " + str(FreeEnergy[i]) + " " + str(length[i]) + "\n")
	
	#plot_folding_structure(N_aminoacids, R)


def phase_transition(N_aminoacids):

	A, J = initialize_Uniform_J(N_aminoacids)

	#high temperature simulation
	N_steps = 2000000
	ave_int = 1800000
	temp_int = 30
	temp_final = 20
	dtemp = 0.5
	temp_step = int((temp_int - temp_final)/dtemp) + 1
	temp_array = np.linspace(temp_int, temp_final, temp_step, endpoint=True)

	for t in temp_array:

		file = open('UniJ5_equilibrium_%d_T%.1f.txt'%(N_aminoacids, t), 'w')

		for n in range(10):

			step, FreeEnergy, L, R = ProteinFolding_jit(N_aminoacids, J, A, t, N_steps)

			for i in range(len(step)):
				if step[i] < ave_int:
					continue
				if step[i] > ave_int:
					file.write(str(step[i]) + " " + str(FreeEnergy[i]) + " " + str(L[i]) + "\n")

	#low temperature simulation
	#low_temp_array = np.linspace(9.5, 0.5, 19, endpoint=True)
	
	#for t in low_temp_array:

	#	if t > 1.5:
	#		continue

	#	Nsteps = 2500000 + int(500000 * ( 9.5 - t) * 2)
	#	ave_int_low = Nsteps - 200000
		
	#	file = open('UniJ5_equilibrium_%d_T%.1f(1).txt'%(N_aminoacids, t), 'w')

	#	for n in range(20):

	#		temp, step, FreeEnergy, L, R = ProteinFolding_low_temp_jit(N_aminoacids, J, A, Nsteps)
			#print(temp)
	#		for i in range(len(step)):
	#			if step[i] < ave_int_low:
	#				continue
	#			if step[i] > ave_int_low:
	#				file.write(str(step[i]) + " " + str(FreeEnergy[i]) + " " + str(L[i]) + "\n")


@jit(nopython=True)
def SAW_ProteinFolding_jit(N_aminoacids, J, A, Nsteps = 1000000):

	R = np.zeros((N_aminoacids, 2))
	# Initialize the protein chain, R is the position array
	for i in range(N_aminoacids):
		R[i,0] = i # x_i
		R[i,1] = 0 # y_i

	FreeEnergy = []
	length = []
	step = []

	for i in range(Nsteps):
		# choose k-th amino acid randomly and make the movement
		k = np.random.randint(0, N_aminoacids)
		x = R[k,0]
		y = R[k,1]
		x_s, y_s = step_2D_jit()
		x_new = x + x_s
		y_new = y + y_s

		#judge whether the walking is self-avoiding
		SAW = 1
		for j in range(N_aminoacids):
			if (R[j,0]-x_new)**2+(R[j,1]-y_new)**2 < 0.1:
				SAW = 0
		if SAW == 0 : 
			continue

		#judge whether the length is not changed
		if k == 0:
			x_next = R[k+1,0]
			y_next = R[k+1,1]
			nbor = neighbor(x_new, y_new, x_next, y_next)
			if nbor == 1:
				R[k,0], R[k,1] = x_new, y_new
			if nbor == 0:
				continue

		if k == N_aminoacids-1:
			x_pre = R[k-1,0]
			y_pre = R[k-1,1]
			nbor = neighbor(x_new, y_new, x_pre, y_pre)
			if nbor == 1:
				R[k,0], R[k,1] = x_new, y_new
			if nbor == 0:
					continue

		if 0 < k < N_aminoacids-1:
			x_next = R[k+1,0]
			y_next = R[k+1,1]
			x_pre = R[k-1,0]
			y_pre = R[k-1,1]
			nbor_next = neighbor(x_new, y_new, x_next, y_next)
			nbor_pre = neighbor(x_new, y_new, x_pre, y_pre)
			if nbor_next == 1 and nbor_pre == 1:
				R[k,0], R[k,1] = x_new, y_new
			if nbor_pre == 0:
				continue
			if nbor_next == 0:
				continue

		# calculate the free energy, r square
		F = 0
		for j in range(N_aminoacids):
			for l in range(j):
				if j-l > 1.1:
					distance = (R[j,0]-R[l,0])**2 + (R[j,1]-R[l,1])**2
					if distance < 1.1:
						F = F - J[A[j],A[l]]

		r2 = (R[N_aminoacids-1,0]-R[0,0])**2 + (R[N_aminoacids-1,1]-R[0,1])**2		
		step.append(i+1)
		FreeEnergy.append(F)
		length.append(np.sqrt(r2))

	return step, FreeEnergy, length


def SAW(N_aminoacids):

	#to confirm the assumption is correct
	N_steps = 2000000
	ave_int = 1800000
	A, J = initialize_Uniform_J(N_aminoacids)
	
	file = open('SAW_UniJ1_%d.txt'%N_aminoacids, 'w')
	
	for n in range(100):

		step, FreeEnergy, L = SAW_ProteinFolding_jit(N_aminoacids, J, A, N_steps)
		for i in range(len(step)):
				if step[i] < ave_int:
					continue
				if step[i] > ave_int:
					file.write(str(step[i]) + " " + str(FreeEnergy[i]) + " " + str(L[i]) + "\n")

def main():

	phase_transition(15)
	#equilibrium(100)
	#SAW(100)


if __name__ == '__main__':
	main()


