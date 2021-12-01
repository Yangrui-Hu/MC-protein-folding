import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import random

class ProteinFolding2D_low_temp(object):

	def __init__(self, N_aminoacids = 100, N_steps = 1000000):
		# the number of the amino acids in the chain or the length of the chain
		self.N_aminoacids = N_aminoacids
		# the total time = the number of steps
		self.N_steps = N_steps
		
		self.A = np.zeros(N_aminoacids, int)
		self.J = np.zeros((20,20))
		self.R = np.zeros((N_aminoacids, 2))


	#@jit(nopython=True)
	def step_2D(self):
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


	# if nbor = 0, [x,y] and [x_test, y_test] are not neighbor
	# if nbor = 1, [x,y] and [x_test, y_test] are neighbors
	def neighbor(self, x, y, x_test, y_test):
		x, y, x_test, y_test = x, y, x_test, y_test
		if (x-x_test)**2+(y-y_test)**2 < 1.1:
			return 1
		if (x-x_test)**2+(y-y_test)**2 > 1.1:
			return 0


	# calculate the change of the energy after the movement and judge whether receive this movement
	#@jit(nopython=True)
	def Calculate_delta_Energy(self, index, xnew, ynew, temp):
		#print('1')
		k = index
		x_new = xnew
		y_new = ynew
		T = temp
		beta = 1.0/T
		R_new = np.zeros((self.N_aminoacids, 2))
		R_before = np.zeros((self.N_aminoacids, 2))

		for i in range(self.N_aminoacids):
			R_new[i,0] = self.R[i,0]
			R_new[i,1] = self.R[i,1]
			R_before[i,0] = self.R[i,0]
			R_before[i,1] = self.R[i,1]
			if i == k :
				R_new[i,0] = x_new
				R_new[i,1] = y_new
		#print(R_new)

		energy_before = 0
		energy_new = 0
		distance_before = 0
		distance_new = 0

		for i in range(self.N_aminoacids):
			for j in range(i):
				if i-j > 1.1:
					distance_before = (R_before[i,0]-R_before[j,0])**2 + (R_before[i,1]-R_before[j,1])**2
					#print(distance_before)
					distance_new = (R_new[i,0]-R_new[j,0])**2 + (R_new[i,1]-R_new[j,1])**2
					if distance_before < 1.1:
						energy_before = energy_before - self.J[self.A[i],self.A[j]]
					if distance_new < 1.1:
						energy_new = energy_new - self.J[self.A[i],self.A[j]]
		
		delta_energy = energy_new - energy_before
		if delta_energy <= 0:
			#print('2')
			self.R[k,0] = x_new
			self.R[k,1] = y_new
		elif delta_energy > 0:
			p = random.random()
			if p < np.exp(-beta*delta_energy):
				#print('3')
				self.R[k,0] = x_new
				self.R[k,1] = y_new


	#@jit(nopython=True)
	def ProteinFolding_low_temp(self):

		# Set the sequence of the protein chain A[N_aminoacids]
		for i in range(self.N_aminoacids):
			self.A[i] = np.random.randint(0,20)

		# Initialize the protein chain, R is the position array
		for i in range(self.N_aminoacids):
			self.R[i,0] = i # x_i
			self.R[i,1] = 0 # y_i
		#print(self.R)

		# Set the interaction energy matrix J
		for i in range(20):
			for j in range(i):
				self.J[i,j] = 2 * random.random() + 2
				self.J[j,i] = self.J[i,j]
			self.J[i,i] = random.random()
		#print(self.J)	
		FreeEnergy = []
		length = []
		step = []

		for i in range(self.N_steps):

			if 0 <= i < 500000:
				temp = 5.0
			if 500000 <= i < 1000000:
				temp = 4.0
			if 1000000 <= i < 1500000:
				temp = 3.0
			if 1500000 <= i < 2000000:
				temp = 2.5
			if 2000000 <= i < 2500000:
				temp = 2.0
			if 2500000 <= i < 3000000:
				temp = 1.5
			if 3000000 <= i < 3500000:
				temp = 1.0
			if 3500000 <= i < 4000000:
				temp = 0.5
			
			# choose k-th amino acid randomly and make the movement
			k = np.random.randint(0, self.N_aminoacids)
			x = self.R[k,0]
			y = self.R[k,1]
			x_s, y_s = self.step_2D()
			x_new = x + x_s
			y_new = y + y_s

			#judge whether the walking is self-avoiding
			SAW = 1
			for j in range(self.N_aminoacids):
				if (self.R[j,0]-x_new)**2+(self.R[j,1]-y_new)**2 < 0.1:
					SAW = 0
			if SAW == 0 : 
				continue

			#judge whether the length is not changed
			if k == 0:
				x_next = self.R[k+1,0]
				y_next = self.R[k+1,1]
				nbor = self.neighbor(x_new, y_new, x_next, y_next)
				if nbor == 1:
					self.Calculate_delta_Energy(k, x_new, y_new, temp)
				if nbor == 0:
					continue

			if k == self.N_aminoacids-1:
				x_pre = self.R[k-1,0]
				y_pre = self.R[k-1,1]
				nbor = self.neighbor(x_new, y_new, x_pre, y_pre)
				if nbor == 1:
					self.Calculate_delta_Energy(k, x_new, y_new, temp)
				if nbor == 0:
					continue

			if 0 < k < self.N_aminoacids-1:
				x_next = self.R[k+1,0]
				y_next = self.R[k+1,1]
				x_pre = self.R[k-1,0]
				y_pre = self.R[k-1,1]
				nbor_next = self.neighbor(x_new, y_new, x_next, y_next)
				nbor_pre = self.neighbor(x_new, y_new, x_pre, y_pre)
				if nbor_next == 1 and nbor_pre == 1:
					self.Calculate_delta_Energy(k, x_new, y_new, temp)
				if nbor_pre == 0:
					continue
				if nbor_next == 0:
					continue
			# calculate the free energy, r square
			F = 0
			for j in range(self.N_aminoacids):
				for l in range(j):
					if j-l > 1.1:
						distance = (self.R[j,0]-self.R[l,0])**2 + (self.R[j,1]-self.R[l,1])**2
						if distance < 1.1:
							F = F - self.J[self.A[j],self.A[l]]

			r2 = (self.R[self.N_aminoacids-1,0]-self.R[0,0])**2 + (self.R[self.N_aminoacids-1,1]-self.R[0,1])**2		
			step.append(i+1)
			FreeEnergy.append(F)
			length.append(np.sqrt(r2))

		return temp, step, FreeEnergy, length

	def plot_folding_structure(self):

		for i in range(self.N_aminoacids-1):
			x1, y1 = self.R[i,0], self.R[i,1]
			x2, y2 = self.R[i+1,0], self.R[i+1,1]
			x = [x1, x2]
			y = [y1, y2]
			plt.plot(x,y,'r')
		for i in range(self.N_aminoacids):
			plt.scatter(self.R[i,0], self.R[i,1], c ='k')
		plt.xlabel('x')
		plt.ylabel('y')
		plt.title('Final Structure (T = 1, 3500000 steps)')
		plt.show()


def low_temp():

	P = ProteinFolding2D_low_temp(N_aminoacids = 30, N_steps = 3000000)
	temp, step, FreeEnergy, r = P.ProteinFolding_low_temp()
	energy_ave = 0
	r_ave = 0
	ave_step = 0
	ave_int = 2800000
	for i in range(len(step)):
		if step[i] < ave_int:
			continue
		if step[i] > ave_int:
			energy_ave = energy_ave + FreeEnergy[i]
			r_ave = r_ave + r[i]
			ave_step = ave_step + 1

	energy_ave = energy_ave / (1.0 * ave_step)
	r_ave = r_ave / (1.0 * ave_step)
	print(temp, energy_ave, r_ave)	
		
low_temp()

def plot_equilibrium():

	P = ProteinFolding2D_low_temp(N_aminoacids = 100, T = 2.0, N_steps = 3500000)
	temp, step, FreeEnergy, length = P.ProteinFolding_low_temp()

	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)
	ax1.plot(step, FreeEnergy, label='100 amino acids, T = 1')
	ax1.set_xlabel("time (Monte Carlo steps)")
	ax1.set_ylabel("Energy")
	ax1.legend()
	ax2.plot(step, length, label='100 amino acids, T = 1')
	ax2.set_xlabel("time (Monte Carlo steps)")
	ax2.set_ylabel("Length")
	ax2.legend()
	plt.show()

	P.plot_folding_structure()

#plot_equilibrium()
