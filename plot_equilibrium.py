import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 17})

def plot_equilibrium(fname, N_aminoacids, Temp):
	a = np.loadtxt(fname)
	step0 = a[:,0]
	E = a[:,1]
	L = a[:,2]

	step = step0[:]*0.00001

	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)

	ax1.plot(step, E, label='%d amino acids, T = %.1f, w/o SA'%(N_aminoacids, Temp))
	ax2.plot(step, L, label='%d amino acids, T = %.1f, w/o SA'%(N_aminoacids, Temp))
	ax1.set_xlabel("Time ($\\times 10^5$ Monte Carlo Step)", fontsize=20)
	ax1.set_ylabel("Energy", fontsize=20)
	ax1.legend(fontsize=15)
	ax2.set_xlabel("Time ($\\times 10^5$ Monte Carlo Step)", fontsize=20)
	ax2.set_ylabel("Length", fontsize=20)
	ax2.legend(fontsize=15)
	plt.show()

def plot_SA_compare(N_aminoacids, Temp):
	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)

	for n in range(10):
		fname = 'full_equilibrium_%d_T%.1f-beforeSA-%d.txt'%(N_aminoacids, Temp, n)
		a = np.loadtxt(fname)
		step0 = a[:,0]
		step = step0[:]*0.00001
		E = a[:,1]
		L = a[:,2]
		ax1.plot(step, E, label = '%d'%n)
		ax2.plot(step, L, label = '%d'%n)

	ax1.set_xlabel("Time ($\\times 10^5$ Monte Carlo Step)", fontsize=20)
	ax1.set_ylabel("Energy", fontsize=20)
	ax1.legend()
	ax2.set_xlabel("Time ($\\times 10^5$ Monte Carlo Step)", fontsize=20)
	ax2.set_ylabel("Length", fontsize=20)
	ax2.legend()
	plt.suptitle('%d amino acids, T = %.1f, w/o SA'%(N_aminoacids, Temp), fontsize=20)
	plt.show()


def main():
	N_aminoacids = 100
	Temp = 1.0
	#fname = 'full_equilibrium_%d_T%.1f-beforeSA-%d.txt'%(N_aminoacids, Temp, n)
	#plot_equilibrium(fname, N_aminoacids, Temp)
	plot_SA_compare(N_aminoacids, Temp)

if __name__ == '__main__':
	main()
