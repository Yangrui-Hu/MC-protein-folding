import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

matplotlib.rcParams.update({'font.size': 12})
#plt.subplots_adjust(wspace = 0.5)

def func_1(x, a, b):  
	return a*x+b

def analysis_phase_transition():
	a = np.loadtxt('UniJ5_phase_15.txt')

	temp_array = a[:,0]
	E = a[:,1]
	se_E = a[:,2]
	L = a[:,4]
	se_L = a[:,5]
	file = open('fitJ5_15.txt', 'w')
	for t in range(11):
		Tc = 1.5 + t * 0.1
		yE = []
		yL = []
		x = []
		for i in range(len(temp_array)):
			if Tc < temp_array[i] < 5:
				yE.append(np.log(abs(E[i])))
				yL.append(np.log(abs(L[i])))
				x.append(np.log(abs(temp_array[i]-Tc)))

		poptE, pcovE = curve_fit(func_1, x, yE)
		poptL, pcovL = curve_fit(func_1, x, yL)
		chi2_E = 0
		chi2_L = 0
		for i in range(len(x)):
			chi2_E += (yE[i]-poptE[0]*x[i]-poptE[1])**2
			chi2_L += (yL[i]-poptL[0]*x[i]-poptL[1])**2

		file.write(str(Tc) + " " + str(poptE[0]) + " " + str(chi2_E) + " " + str(poptL[0]) +  " " + str(chi2_L) + "\n")

def plot_fit_phasetransition():
	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)

	b = np.loadtxt('fitJ5_15.txt')
	Tc = b[:,0]
	kE = b[:,1]
	chi2E = b[:,2]
	kL = b[:,3]
	chi2L = b[:,4]
	ax1.plot(Tc, kE, label='critical exponent')
	ax1.set_xlabel("$T_c$", fontsize=12)
	ax1.set_ylabel("k", fontsize=12)
	ax1.legend(fontsize=12)
	ax3 = ax1.twinx()
	ax3.plot(Tc,chi2E, 'red', label = '$\chi^2$')
	ax3.set_xlabel("$T_c$", fontsize=12)
	ax3.set_ylabel('$\chi^2$',fontsize=12)
	ax3.legend(fontsize=12)
	ax1.set_title('Energy')
	
	ax2.plot(Tc, kL, label='critical exponent')
	ax2.set_xlabel("$T_c$", fontsize=12)
	ax2.set_ylabel("k", fontsize=12)
	ax2.legend(fontsize=12)
	ax4 = ax2.twinx()
	ax4.plot(Tc, chi2L, 'red', label='$\chi^2$')
	ax4.set_xlabel("$T_c$", fontsize=12)
	ax4.set_ylabel('$\chi^2$',fontsize=12)
	ax4.legend(fontsize=12)
	ax2.set_title('Length')
	plt.suptitle('15 amino acids, $J = -5$')
	plt.show()

def funcE(x,a,r):
	Tc = 2.0
	if x > Tc:
		return a*(x/Tc-1)**r-36
	if x <= Tc:
		return -36

def fit_funcE(x,a,r):
	f = np.zeros(len(x))
	for i in range(len(f)):
		f[i] = funcE(x[i],a,r)
	return f
 
def funcL(x,a,r):
	Tc = 2.0
	if x > Tc:
		return a*(x/Tc-1)**r-2.44
	if x <= Tc:
		return 2.44

def fit_funcL(x,a,r):
	f = np.zeros(len(x))
	for i in range(len(f)):
		f[i] = funcL(x[i],a,r)
	return f

def plot_15_phase():

	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)
	
	#a = np.loadtxt('phase_15.txt')
	#a1 = np.loadtxt('UniJ1_phase_15.txt')
	#a1_SAW = np.loadtxt('SAW_UniJ1_phase_15.txt')
	#a3 = np.loadtxt('UniJ3_phase_15.txt')
	#a3_SAW = np.loadtxt('SAW_UniJ3_phase_15.txt')
	a5 = np.loadtxt('UniJ5_phase_15_full.txt')
	#a5_SAW = np.loadtxt('SAW_UniJ5_phase_15.txt')

	#temp_array = a[:,0]
	#E = a[:,1]
	#se_E = a[:,2]
	#std_E = a[:,3]
	#L = a[:,4]
	#se_L = a[:,5]
	#std_L = a[:,6]

	#temp_array1 = a1[:,0]
	#E1 = a1[:,1]
	#se_E1 = a1[:,2]
	#std_E1 = a1[:,3]
	#L1 = a1[:,4]
	#se_L1 = a1[:,5]
	#std_L1 = a1[:,6]

	#temp_array1_SAW = a1_SAW[:,0]
	#E1_SAW = a1_SAW[:,1]
	#se_E1_SAW = a1_SAW[:,2]
	#std_E1_SAW = a1_SAW[:,3]
	#L1_SAW = a1_SAW[:,4]
	#se_L1_SAW = a1_SAW[:,5]
	#std_L1_SAW = a1_SAW[:,6]

	#temp_array3 = a3[:,0]
	#E3 = a3[:,1]
	#se_E3 = a3[:,2]
	#std_E3 = a3[:,3]
	#L3 = a3[:,4]
	#se_L3 = a3[:,5]
	#std_L3 = a3[:,6]

	#temp_array3_SAW = a3_SAW[:,0]
	#E3_SAW = a3_SAW[:,1]
	#se_E3_SAW = a3_SAW[:,2]
	#std_E3_SAW = a3_SAW[:,3]
	#L3_SAW = a3_SAW[:,4]
	#se_L3_SAW = a3_SAW[:,5]
	#std_L3_SAW = a3_SAW[:,6]

	temp_array5 = a5[:,0]
	E5 = a5[:,1]
	se_E5 = a5[:,2]
	std_E5 = a5[:,3]
	L5 = a5[:,4]
	se_L5 = a5[:,5]
	std_L5 = a5[:,6]

	temp_array5_fitE = []
	temp_array5_fitL = []
	E5_fit = []
	L5_fit = []
	for i in range(len(temp_array5)):
		if temp_array5[i] < 10:
			temp_array5_fitE.append(temp_array5[i])
			E5_fit.append(E5[i])
		if temp_array5[i] < 10:
			temp_array5_fitL.append(temp_array5[i])
			L5_fit.append(L5[i])
	poptE, pcovE = curve_fit(fit_funcE, temp_array5_fitE, E5_fit)
	yvalsE = fit_funcE(temp_array5_fitE, poptE[0],poptE[1])
	poptL, pcovL = curve_fit(fit_funcL, temp_array5_fitL, L5_fit)
	yvalsL = fit_funcL(temp_array5_fitL, poptL[0],poptL[1])

	#temp_array5_SAW = a5_SAW[:,0]
	#E5_SAW = a5_SAW[:,1]
	#se_E5_SAW = a5_SAW[:,2]
	#std_E5_SAW = a5_SAW[:,3]
	#L5_SAW = a5_SAW[:,4]
	#se_L5_SAW = a5_SAW[:,5]
	#std_L5_SAW = a5_SAW[:,6]


	
	#ax1.errorbar(temp_array, E, se_E, fmt='o' , mfc = 'red', label='$J\in[-2,-4]$')
	#ax1.errorbar(temp_array1, E1, se_E1, fmt='o' , mfc = 'black', label='$J = -1$')
	#ax1.plot(temp_array1_SAW, E1_SAW, 'orange', label='SAW, $J = -1$')
	#ax1.errorbar(temp_array3, E3, se_E3, fmt='o' , mfc = 'blue', label='$J = -3$')
	ax1.errorbar(temp_array5, E5, se_E5, fmt='o' , mfc = 'black', label='$J = -5$')
	ax1.plot(temp_array5_fitE,yvalsE,label='fitting curve, $\\beta=%.2f$'%poptE[1])
	#ax1.errorbar(temp_array, E, se_E, fmt='o' , mfc = 'black', label='$J = -7$')
	
	#ax2.errorbar(temp_array, L, se_L, fmt='o' , mfc = 'red', label='$J\in[-2,-4]$')
	#ax2.errorbar(temp_array1, L1, se_L1, fmt='o' , mfc = 'black', label='$J = -1$')
	#ax2.plot(temp_array1_SAW, L1_SAW, 'orange', label='SAW, $J = -1$')
	#ax2.errorbar(temp_array3, L3, se_L3, fmt='o' , mfc = 'blue', label='$J = -3$')
	ax2.errorbar(temp_array5, L5, se_L5, fmt='o' , mfc = 'black', label='$J = -5$')
	ax2.plot(temp_array5_fitL,yvalsL,label='fitting curve, $\\beta=%.2f$'%poptL[1])
	#ax2.errorbar(temp_array, L, se_L, fmt='o' , mfc = 'black', label='$J = -7$')
	
	#ax1.set_ylim(-30,0)
	#ax2.set_ylim(1,9)
	ax1.set_xlabel("Temperature", fontsize=20)
	ax1.set_ylabel("Energy", fontsize=20)
	ax1.legend(fontsize=12)
	ax2.set_xlabel("Temperature", fontsize=20)
	ax2.set_ylabel("Length", fontsize=20)
	ax2.legend(fontsize=12)
	plt.suptitle('15 amino acids')
	plt.show()

def plot_30_phase():

	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)
	
	a = np.loadtxt('phase_30.txt')
	a1 = np.loadtxt('UniJ1_phase_30.txt')
	a1_SAW = np.loadtxt('SAW_UniJ1_phase_30.txt')
	a3 = np.loadtxt('UniJ3_phase_30.txt')
	a5 = np.loadtxt('UniJ5_phase_30.txt')

	temp_array = a[:,0]
	E = a[:,1]
	se_E = a[:,2]
	std_E = a[:,3]
	L = a[:,4]
	se_L = a[:,5]
	std_L = a[:,6]

	temp_array1 = a1[:,0]
	E1 = a1[:,1]
	se_E1 = a1[:,2]
	std_E1 = a1[:,3]
	L1 = a1[:,4]
	se_L1 = a1[:,5]
	std_L1 = a1[:,6]

	temp_array1_SAW = a1_SAW[:,0]
	E1_SAW = a1_SAW[:,1]
	se_E1_SAW = a1_SAW[:,2]
	std_E1_SAW = a1_SAW[:,3]
	L1_SAW = a1_SAW[:,4]
	se_L1_SAW = a1_SAW[:,5]
	std_L1_SAW = a1_SAW[:,6]

	temp_array3 = a3[:,0]
	E3 = a3[:,1]
	se_E3 = a3[:,2]
	std_E3 = a3[:,3]
	L3 = a3[:,4]
	se_L3 = a3[:,5]
	std_L3 = a3[:,6]

	temp_array5 = a5[:,0]
	E5 = a5[:,1]
	se_E5 = a5[:,2]
	std_E5 = a5[:,3]
	L5 = a5[:,4]
	se_L5 = a5[:,5]
	std_L5 = a5[:,6]

	
	#ax1.errorbar(temp_array, E, se_E, fmt='o' , mfc = 'red', label='30 amino acids, $J\in[-2,-4]$')
	ax1.errorbar(temp_array1, E1, se_E1, fmt='o' , mfc = 'black', label='$J = -1$')
	ax1.plot(temp_array1_SAW, E1_SAW, 'orange', label='SAW, $J = -1$')
	#ax1.errorbar(temp_array3, E3, se_E3, fmt='o' , mfc = 'blue', label='30 amino acids, $J = -3$')
	#ax1.errorbar(temp_array5, E5, se_E5, fmt='o' , mfc = 'green', label='30 amino acids, $J = -5$')

	#ax2.errorbar(temp_array, L, se_L, fmt='o' , mfc = 'red', label='30 amino acids, $J\in[-2,-4]$')
	ax2.errorbar(temp_array1, L1, se_L1, fmt='o' , mfc = 'black', label='$J = -1$')
	ax2.plot(temp_array1_SAW, L1_SAW, 'orange', label='SAW, $J = -1$')
	#ax2.errorbar(temp_array3, L3, se_L3, fmt='o' , mfc = 'blue', label='30 amino acids, $J = -3$')
	#ax2.errorbar(temp_array5, L5, se_L5, fmt='o' , mfc = 'green', label='30 amino acids, $J = -5$')

	#ax1.set_ylim(-30,0)
	#ax2.set_ylim(1,9)
	ax1.set_xlabel("Temperature", fontsize=20)
	ax1.set_ylabel("Energy", fontsize=20)
	ax1.legend(fontsize=12)
	ax2.set_xlabel("Temperature", fontsize=20)
	ax2.set_ylabel("Length", fontsize=20)
	ax2.legend(fontsize=12)
	plt.suptitle('30 amino acids')
	plt.show()

def plot_100_phase():

	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)
	
	a = np.loadtxt('phase_100.txt')
	a1 = np.loadtxt('UniJ1_phase_100.txt')
	a1_SAW = np.loadtxt('SAW_UniJ1_phase_100.txt')
	a3 = np.loadtxt('UniJ3_phase_100.txt')
	a5 = np.loadtxt('UniJ5_phase_100.txt')

	temp_array = a[:,0]
	E = a[:,1]
	se_E = a[:,2]
	std_E = a[:,3]
	L = a[:,4]
	se_L = a[:,5]
	std_L = a[:,6]

	temp_array1 = a1[:,0]
	E1 = a1[:,1]
	se_E1 = a1[:,2]
	std_E1 = a1[:,3]
	L1 = a1[:,4]
	se_L1 = a1[:,5]
	std_L1 = a1[:,6]

	temp_array1_SAW = a1_SAW[:,0]
	E1_SAW = a1_SAW[:,1]
	se_E1_SAW = a1_SAW[:,2]
	std_E1_SAW = a1_SAW[:,3]
	L1_SAW = a1_SAW[:,4]
	se_L1_SAW = a1_SAW[:,5]
	std_L1_SAW = a1_SAW[:,6]

	temp_array3 = a3[:,0]
	E3 = a3[:,1]
	se_E3 = a3[:,2]
	std_E3 = a3[:,3]
	L3 = a3[:,4]
	se_L3 = a3[:,5]
	std_L3 = a3[:,6]

	temp_array5 = a5[:,0]
	E5 = a5[:,1]
	se_E5 = a5[:,2]
	std_E5 = a5[:,3]
	L5 = a5[:,4]
	se_L5 = a5[:,5]
	std_L5 = a5[:,6]

	
	#ax1.errorbar(temp_array, E, se_E, fmt='o' , mfc = 'red', label='100 amino acids, $J\in[-2,-4]$')
	ax1.errorbar(temp_array1, E1, se_E1, fmt='o' , mfc = 'black', label='$J = -1$')
	ax1.plot(temp_array1_SAW, E1_SAW, 'orange', label='SAW, $J = -1$')
	#ax1.errorbar(temp_array3, E3, se_E3, fmt='o' , mfc = 'blue', label='100 amino acids, $J = -3$')
	#ax1.errorbar(temp_array5, E5, se_E5, fmt='o' , mfc = 'green', label='100 amino acids, $J = -5$')

	#ax2.errorbar(temp_array, L, se_L, fmt='o' , mfc = 'red', label='100 amino acids, $J\in[-2,-4]$')
	ax2.errorbar(temp_array1, L1, se_L1, fmt='o' , mfc = 'black', label='$J = -1$')
	ax2.plot(temp_array1_SAW, L1_SAW, 'orange', label='SAW, $J = -1$')
	#ax2.errorbar(temp_array3, L3, se_L3, fmt='o' , mfc = 'blue', label='100 amino acids, $J = -3$')
	#ax2.errorbar(temp_array5, L5, se_L5, fmt='o' , mfc = 'green', label='100 amino acids, $J = -5$')

	#ax1.set_ylim(-30,0)
	#ax2.set_ylim(1,9)
	ax1.set_xlabel("Temperature", fontsize=20)
	ax1.set_ylabel("Energy", fontsize=20)
	ax1.legend(fontsize=12)
	ax2.set_xlabel("Temperature", fontsize=20)
	ax2.set_ylabel("Length", fontsize=20)
	ax2.legend(fontsize=12)
	plt.suptitle('100 amino acids')
	plt.show()

def plot_E():

	fig = plt.figure()
	ax1 = fig.add_subplot(131)
	ax2 = fig.add_subplot(132)
	ax3 = fig.add_subplot(133)
	
	a15 = np.loadtxt('phase_15.txt')
	a15_1 = np.loadtxt('UniJ1_phase_15.txt')
	a15_3 = np.loadtxt('UniJ3_phase_15.txt')
	a15_5 = np.loadtxt('UniJ5_phase_15.txt')
	a30 = np.loadtxt('phase_30.txt')
	a30_1 = np.loadtxt('UniJ1_phase_30.txt')
	a30_3 = np.loadtxt('UniJ3_phase_30.txt')
	a30_5 = np.loadtxt('UniJ5_phase_30.txt')
	a100 = np.loadtxt('phase_100.txt')
	a100_1 = np.loadtxt('UniJ1_phase_100.txt')
	a100_3 = np.loadtxt('UniJ3_phase_100.txt')
	a100_5 = np.loadtxt('UniJ5_phase_100.txt')

	temp_array15 = a15[:,0]
	E15 = a15[:,1]
	se_E15 = a15[:,2]

	temp_array15_1 = a15_1[:,0]
	E15_1 = a15_1[:,1]
	se_E15_1 = a15_1[:,2]

	temp_array15_3 = a15_3[:,0]
	E15_3 = a15_3[:,1]
	se_E15_3 = a15_3[:,2]

	temp_array15_5 = a15_5[:,0]
	E15_5 = a15_5[:,1]
	se_E15_5 = a15_5[:,2]

	ax1.errorbar(temp_array15, E15, se_E15, fmt='o' , mfc = 'red', label='$J\in[-2,-4]$')
	ax1.errorbar(temp_array15_1, E15_1, se_E15_1, fmt='o' , mfc = 'black', label='$J = -1$')
	ax1.errorbar(temp_array15_3, E15_3, se_E15_3, fmt='o' , mfc = 'blue', label='$J = -3$')
	ax1.errorbar(temp_array15_5, E15_5, se_E15_5, fmt='o' , mfc = 'green', label='$J = -5$')


	temp_array30 = a30[:,0]
	E30 = a30[:,1]
	se_E30 = a30[:,2]

	temp_array30_1 = a30_1[:,0]
	E30_1 = a30_1[:,1]
	se_E30_1 = a30_1[:,2]

	temp_array30_3 = a30_3[:,0]
	E30_3 = a30_3[:,1]
	se_E30_3 = a30_3[:,2]

	temp_array30_5 = a30_5[:,0]
	E30_5 = a30_5[:,1]
	se_E30_5 = a30_5[:,2]

	ax2.errorbar(temp_array30, E30, se_E30, fmt='o' , mfc = 'red', label='$J\in[-2,-4]$')
	ax2.errorbar(temp_array30_1, E30_1, se_E30_1, fmt='o' , mfc = 'black', label='$J = -1$')
	ax2.errorbar(temp_array30_3, E30_3, se_E30_3, fmt='o' , mfc = 'blue', label='$J = -3$')
	ax2.errorbar(temp_array30_5, E30_5, se_E30_5, fmt='o' , mfc = 'green', label='$J = -5$')

	temp_array100 = a100[:,0]
	E100 = a100[:,1]
	se_E100 = a100[:,2]

	temp_array100_1 = a100_1[:,0]
	E100_1 = a100_1[:,1]
	se_E100_1 = a100_1[:,2]

	temp_array100_3 = a100_3[:,0]
	E100_3 = a100_3[:,1]
	se_E100_3 = a100_3[:,2]

	temp_array100_5 = a100_5[:,0]
	E100_5 = a100_5[:,1]
	se_E100_5 = a100_5[:,2]

	ax3.errorbar(temp_array100, E100, se_E100, fmt='o' , mfc = 'red', label='$J\in[-2,-4]$')
	ax3.errorbar(temp_array100_1, E100_1, se_E100_1, fmt='o' , mfc = 'black', label='$J = -1$')
	ax3.errorbar(temp_array100_3, E100_3, se_E100_3, fmt='o' , mfc = 'blue', label='$J = -3$')
	ax3.errorbar(temp_array100_5, E100_5, se_E100_5, fmt='o' , mfc = 'green', label='$J = -5$')

	ax1.set_xlabel("Temperature", fontsize=20)
	ax1.set_ylabel("Energy", fontsize=20)
	ax1.set_title("15 amino acids")
	ax1.legend(fontsize=12)
	ax2.set_xlabel("Temperature", fontsize=20)
	#ax2.set_ylabel("Energy", fontsize=20)
	ax2.legend(fontsize=12)
	ax2.set_title("30 amino acids")
	ax3.set_xlabel("Temperature", fontsize=20)
	#ax3.set_ylabel("Energy", fontsize=20)
	ax3.legend(fontsize=12)
	ax3.set_title("100 amino acids")
	plt.show()

def plot_SAW_compare():
	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)

	a15_1 = np.loadtxt('UniJ1_phase_15.txt')
	a15_1_SAW = np.loadtxt('SAW_UniJ1_phase_15.txt')
	a30_1 = np.loadtxt('UniJ1_phase_30.txt')
	a30_1_SAW = np.loadtxt('SAW_UniJ1_phase_30.txt')
	a100_1 = np.loadtxt('UniJ1_phase_100.txt')
	a100_1_SAW = np.loadtxt('SAW_UniJ1_phase_100.txt')

	temp_array15_1 = a15_1[:,0]
	E15_1 = a15_1[:,1]
	se_E15_1 = a15_1[:,2]
	L15_1 = a15_1[:,4]
	se_L15_1 = a15_1[:,5]
	temp_array15_1_SAW = a15_1_SAW[:,0]
	E15_1_SAW = a15_1_SAW[:,1]
	se_E15_1_SAW = a15_1_SAW[:,2]
	L15_1_SAW = a15_1_SAW[:,4]
	se_L15_1_SAW = a15_1_SAW[:,5]

	temp_array30_1 = a30_1[:,0]
	E30_1 = a30_1[:,1]
	se_E30_1 = a30_1[:,2]
	L30_1 = a30_1[:,4]
	se_L30_1 = a30_1[:,5]
	temp_array30_1_SAW = a30_1_SAW[:,0]
	E30_1_SAW = a30_1_SAW[:,1]
	se_E30_1_SAW = a30_1_SAW[:,2]
	L30_1_SAW = a30_1_SAW[:,4]
	se_L30_1_SAW = a30_1_SAW[:,5]

	temp_array100_1 = a100_1[:,0]
	E100_1 = a100_1[:,1]
	se_E100_1 = a100_1[:,2]
	L100_1 = a100_1[:,4]
	se_L100_1 = a100_1[:,5]
	temp_array100_1_SAW = a100_1_SAW[:,0]
	E100_1_SAW = a100_1_SAW[:,1]
	se_E100_1_SAW = a100_1_SAW[:,2]
	L100_1_SAW = a100_1_SAW[:,4]
	se_L100_1_SAW = a100_1_SAW[:,5]

	ax1.errorbar(temp_array15_1, E15_1, se_E15_1, fmt='o' , mfc = 'black', label='15 amino acids')
	ax1.plot(temp_array15_1_SAW, E15_1_SAW, 'orange', label='15 amino acids, SAW')
	ax1.errorbar(temp_array30_1, E30_1, se_E30_1, fmt='o' , mfc = 'blue', label='30 amino acids')
	ax1.plot(temp_array30_1_SAW, E30_1_SAW, 'darkblue', label='30 amino acids, SAW')
	ax1.errorbar(temp_array100_1, E100_1, se_E100_1, fmt='o' , mfc = 'green', label='100 amino acids')
	ax1.plot(temp_array100_1_SAW, E100_1_SAW, 'darkgreen', label='100 amino acids, SAW')

	ax2.errorbar(temp_array15_1, L15_1, se_L15_1, fmt='o' , mfc = 'black', label='15 amino acids')
	ax2.plot(temp_array15_1_SAW, L15_1_SAW, 'orange', label='15 amino acids, SAW')
	ax2.errorbar(temp_array30_1, L30_1, se_L30_1, fmt='o' , mfc = 'blue', label='30 amino acids')
	ax2.plot(temp_array30_1_SAW, L30_1_SAW, 'darkblue', label='30 amino acids, SAW')
	ax2.errorbar(temp_array100_1, L100_1, se_L100_1, fmt='o' , mfc = 'green', label='100 amino acids')
	ax2.plot(temp_array100_1_SAW, L100_1_SAW, 'darkgreen', label='100 amino acids, SAW')

	ax1.set_xlabel("Temperature", fontsize=20)
	ax1.set_ylabel("Energy", fontsize=20)
	ax1.legend(fontsize=12)
	ax2.set_xlabel("Temperature", fontsize=20)
	ax2.set_ylabel("Length", fontsize=20)
	ax2.legend(fontsize=12)
	plt.suptitle('$J = -1$')
	plt.show()

def main():

	#plot_15_phase()
	#plot_30_phase()
	#plot_100_phase()
	#plot_E()
	#plot_SAW_compare()
	analysis_phase_transition()
	plot_fit_phasetransition()

if __name__ == '__main__':
	main()

