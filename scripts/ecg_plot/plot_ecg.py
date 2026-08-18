import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def read_csv_file (filename):
   with open(filename, 'r') as file:
        reader = csv.reader(file)
        lines = list(reader)
        return lines

def read_clinical_ecg (filename):
    ecg_data = read_csv_file(filename)
    nleads = len(ecg_data)
    x_arr, y_arr = [], []
    for k in range(nleads):
        x, y = [], []
        for i in range(len(ecg_data[k])):
            x.append(i)
            y.append(float(ecg_data[k][i]))
        x_arr.append(x)
        y_arr.append(y)
        #print("[Clinical ECG] {Lead %d} Maximum value = %g" % (k,max(y)))
    return np.array(x_arr[0]), np.array(y_arr), nleads

def read_ecg_readings (input_file):
	data = np.genfromtxt(input_file)
	nlin, ncol = np.shape(data)
	timesteps = data[:,0]
	currents = data[:,1:]
	num_leads = ncol-1
	return timesteps, currents, num_leads
	
def normalise_readings (ecg, leads):
    nb_leads = len(leads)
    ECG = np.zeros( (nb_leads, len(ecg[0])), dtype='float64')
    # Normalize by the max value
    max_val = np.amax(np.abs(ecg))	# Overwrite input value
    ECG = (ecg / max_val)
    # Normalize by the max value on each lead
    #for i in range(nb_leads):
    #    max_val = np.amax(np.abs(ecg[i, :]))
    #    ECG[i, :] = ecg[i, :] / max_val
    # Non-normalized
    #ECG = ecg
    return ECG

def post_process_clinical_ecg(t_0, ECG_0, ecg_qrs_shift):
    t_0 = t_0[ecg_qrs_shift:]-ecg_qrs_shift
    ECG_0 = ECG_0[:,ecg_qrs_shift:]
    return t_0, ECG_0


def read_monoalg3d_ecg (input_file):
	data = np.genfromtxt(input_file, skip_footer=1)
	t = data[:, 0]
	LA = data[:, 1]
	RA = data[:, 2]
	LL = data[:, 3]
	RL = data[:, 4]
	V1 = data[:, 5]
	V2 = data[:, 6]
	V3 = data[:, 7]
	V4 = data[:, 8]
	V5 = data[:, 9]
	V6 = data[:, 10]

	# Ealuate Wilson's central terminal
	VW = 1.0 / 3.0 * (RA + LA + LL)

	# Evaluate simulated ECG lead traces
	V1 = V1 - VW
	V2 = V2 - VW
	V3 = V3 - VW
	V4 = V4 - VW
	V5 = V5 - VW
	V6 = V6 - VW
	I = LA - RA
	II = LL - RA
	III = LL - LA
	aVL = LA - (RA + LL) / 2.0
	aVF = LL - (LA + RA) / 2.0
	aVR = RA - (LA + LL) / 2.0
	ecgs = np.vstack((I, II, V1, V2, V3, V4, V5, V6))
	return t, ecgs

def read_personalisation_ecg (filename):
	ecg_data = read_csv_file(filename)
	ecg_data = ecg_data[1:]	# Do not consider first line of the CSV
	nlin = len(ecg_data)
	ncol = len(ecg_data[0])
	timesteps = np.arange(0,nlin)
	x_arr, y_arr = [], []
	for k in range(ncol):
		y = []	
		for i in range(len(timesteps)):
			y.append(float(ecg_data[i][k]))
		y_arr.append(y)	
	for i in range(len(timesteps)):
		x_arr.append(timesteps[i])
	return np.array(x_arr), np.array(y_arr), ncol

def plot_ecg_readings (t, currents, nleads):
	for i in range(nleads):
		#plt.grid()
		plt.plot(t, currents[:,i], label="lead_%d" % (i), c="black", linewidth=3.0)
		plt.xlabel("t (ms)",fontsize=15)
		plt.ylabel("Current (mA)",fontsize=15)
		plt.title("ECG reading - Lead %d" % (i),fontsize=14)
		plt.legend(loc=0,fontsize=14)
		plt.show()
		#plt.savefig("ecg_lead_%d.pdf" % (i))

def plot_ecg (ecg_data, lead_names, output_filename):
	n_leads = ecg_data.shape[0]
	fig, axes = plt.subplots(2, 4, figsize=(20, 10))
	axes = np.reshape(axes, n_leads)
	for lead_i in range(n_leads):
		t = np.arange(ecg_data.shape[1])
		axes[lead_i].plot(t, ecg_data[lead_i], color='black', linewidth=1.)
		axes[lead_i].set_title(lead_names[lead_i])
		axes[lead_i].set_ylim([-1.5, 1.5])
		#for tick in axes[lead_i].xaxis.get_major_ticks():
        #    tick.label1.set_fontsize(14)
        #for tick in axes[lead_i].yaxis.get_major_ticks():
        #    tick.label1.set_fontsize(14)
    #axes[lead_i].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20)
	
	fig.savefig(output_filename, dpi=300, bbox_inches='tight')
	plt.close(fig)

def plot_all_ecgs (t_clinical, t_monoalg, t_personalisation, clinical_ecg, monoalg_ecg, personalisation_ecg, lead_names, output_filename):
	n_leads = clinical_ecg.shape[0]
	fig, axes = plt.subplots(1, n_leads, figsize=(20, 10))
	axes = np.reshape(axes, n_leads)
	for lead_i in range(n_leads):
		axes[lead_i].plot(t_clinical, clinical_ecg[lead_i], color='red', label='Clinical', linewidth=1.)
		axes[lead_i].plot(t_monoalg, monoalg_ecg[lead_i], color='blue', label='MonoAlg3D', linewidth=1.)
		axes[lead_i].plot(t_personalisation, personalisation_ecg[lead_i], color='green', label='Personalisation', linewidth=1.)
		axes[lead_i].set_title(lead_names[lead_i])
		axes[lead_i].set_ylim([-1.5, 1.5])
		#for tick in axes[lead_i].xaxis.get_major_ticks():
        #    tick.label1.set_fontsize(14)
        #for tick in axes[lead_i].yaxis.get_major_ticks():
        #    tick.label1.set_fontsize(14)
	axes[n_leads-1].legend(loc='lower right', fontsize=10)
	
	fig.savefig(output_filename, dpi=300, bbox_inches='tight')
	plt.close(fig)

def plot_monoalg3d_personalisation_ecgs (t_monoalg, t_personalisation, monoalg_ecg, personalisation_ecg, lead_names, output_filename):
	n_leads = monoalg_ecg.shape[0]
	fig, axes = plt.subplots(1, n_leads, figsize=(20, 10))
	axes = np.reshape(axes, n_leads)
	for lead_i in range(n_leads):
		axes[lead_i].plot(t_monoalg, monoalg_ecg[lead_i], color='blue', label='MonoAlg3D', linewidth=1.)
		axes[lead_i].plot(t_personalisation, personalisation_ecg[lead_i], color='green', label='Personalisation', linewidth=1.)
		axes[lead_i].set_title(lead_names[lead_i])
		#axes[lead_i].set_ylim([-2.5, 2.5])
		#for tick in axes[lead_i].xaxis.get_major_ticks():
        #    tick.label1.set_fontsize(14)
        #for tick in axes[lead_i].yaxis.get_major_ticks():
        #    tick.label1.set_fontsize(14)
	axes[n_leads-1].legend(loc='lower right', fontsize=10)
	
	fig.savefig(output_filename, dpi=300, bbox_inches='tight')
	plt.close(fig)


def main():
	
	if len(sys.argv) != 2:
		print("-"*100)
		print("Usage:> python %s <monoalg3d_ecg>" % sys.argv[0])
		print("-"*100)
		return 1

	monoalg3d_ecg_file = sys.argv[1]

	lead_names = ['I', 'II', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6']
	ecg_qrs_shift = 0
	# Clinical
	#t_0, c_0, nleads_0 = read_clinical_ecg(clinical_ecg_file)
	#ECG_0 = normalise_readings(c_0, lead_names)
	#t_0, ECG_0 = post_process_clinical_ecg(t_0, ECG_0, ecg_qrs_shift)
	#plot_ecg(ECG_0, lead_names, 'outputs/clinical_ecg.png')

	# Personalisation
	#t_2, c_2, nleads_2 = read_personalisation_ecg(personalisation_ecg_file)
	#ECG_2 = normalise_readings(c_2, lead_names)
	#t_2, ECG_2 = post_process_clinical_ecg(t_2, ECG_2, ecg_qrs_shift)
	#plot_ecg(ECG_2, lead_names, 'outputs/personalisation_ecg.png')

	# MonoAlg3D
	t_1, c_1 = read_monoalg3d_ecg(monoalg3d_ecg_file)
	ECG_1 = normalise_readings(c_1, lead_names)
	t_1, ECG_1 = post_process_clinical_ecg(t_1, ECG_1, ecg_qrs_shift)
	plot_ecg(ECG_1, lead_names, 'outputs/monoalg3d_ecg.png')

	#plot_monoalg3d_personalisation_ecgs(t_1, t_2, ECG_1, ECG_2, lead_names, 'outputs/ecg_monoalg3d_personalisation_comparison.png')
	#plot_all_ecgs(t_0, t_1, t_2, ECG_0, ECG_1, ECG_2, lead_names, 'outputs/ecg_comparison.png')

if __name__ == "__main__":
	main()
