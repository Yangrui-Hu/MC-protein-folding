import numpy as np
import matplotlib.pyplot as plt

def equilibrium_mean_SD(fname):
	a = np.loadtxt(fname)
	E = a[:,1]
	L = a[:,2]
	std_E = np.std(E)
	std_L = np.std(L)
	se_E = std_E/np.sqrt(len(E))
	se_L = std_L/np.sqrt(len(L))
	return np.mean(E), np.mean(L), std_E, se_E, std_L, se_L

def phase_data(N_aminoacids = 15):
	temp_int = 30
	temp_final = 0.5
	dtemp = 0.5
	temp_step = int((temp_int - temp_final)/dtemp) + 1
	temp_array = np.linspace(temp_int, temp_final, temp_step, endpoint=True)
	file = open('UniJ5_phase_%d.txt'%N_aminoacids, 'w')
	for t in temp_array:
		fname = "UniJ5_equilibrium_%d_T%.1f.txt"%(N_aminoacids, t)
		mean_E, mean_L, std_E, se_E, std_L, se_L = equilibrium_mean_SD(fname) 
		file.write(str(t) + " " + str(mean_E) + " " + str(se_E) + " " + 
			str(std_E) + " " + str(mean_L) +  " " + str(se_L)+ " " + str(std_L) 
			+ "\n")

def SAW_data(N_aminoacids = 15):
	temp_int = 10
	temp_final = 0.5
	dtemp = 0.5
	temp_step = int((temp_int - temp_final)/dtemp) + 1
	temp_array = np.linspace(temp_int, temp_final, temp_step, endpoint=True)
	
	fname = 'SAW_UniJ1_%d.txt'%N_aminoacids
	a = np.loadtxt(fname)
	E = a[:,1]
	L = a[:,2]
	mean_E = np.mean(E)
	mean_L = np.mean(L)
	std_E = np.std(E)
	std_L = np.std(L)
	se_E = std_E/np.sqrt(len(E))
	se_L = std_L/np.sqrt(len(L))

	file = open('SAW_UniJ1_phase_%d.txt'%N_aminoacids, 'w')
	for t in temp_array:
		file.write(str(t) + " " + str(mean_E) + " " + str(se_E) + " " + str(std_E) + " " + str(mean_L) +  " " + str(se_L)+ " " + str(std_L) + "\n")

def main():
	phase_data(N_aminoacids = 15)
	#SAW_data(100)

if __name__ == '__main__':
	main()
