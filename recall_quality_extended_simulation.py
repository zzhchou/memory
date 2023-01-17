from __future__ import division
from neuron import h
from scipy.stats import pearsonr
import time as cookie
import pickle
import random
import numpy as np
import cell_network
import matplotlib.pyplot as plt
h.load_file("stdrun.hoc")
np.random.seed(149)

def generate_spike_times(frequency, stim_length):
    intervals = []
    mu = 1000./frequency
    elapsed_time = 0
    flag = 1
    while flag:
        roll = np.random.uniform(0, 1)
        interval = -mu*np.log(roll)
        elapsed_time += interval
        if elapsed_time > stim_length:
            flag = 0
        else:
            intervals.append(interval)
    return list(intervals)

def generate_full_spike_times(spike_intervals, spike_times, stim_length, delay_interval, num_repetitions):
    true_spikes_full = []
    if len(spike_times) != 0:
        hanging_interval = stim_length - spike_times[-1]
        spike_intervals[0] = spike_intervals[0] + hanging_interval + delay_interval
         
        for i in range(num_repetitions):
            for spike in spike_intervals:
                true_spikes_full.append(spike)
        
        true_spikes_full[0] = true_spikes_full[0] - (hanging_interval  + delay_interval)

    return true_spikes_full

#################
# Creating cell #
#################
ID = 0
location = (0,0)

synvars = {}
#synvars['type'] = "E2"
#synvars['type'] = "Log_E2Syn"
#synvars['type'] = "Add_E2Syn"
synvars['type'] = "Mult_E2Syn"
#
fname_morph = 'output0_updated.swc'
modeltype = 'Single'
celltype = 'ca3pyramidalcell_Bakerb_network'

num_input = 50
num_cells = 50
sparsity = 0.25

radiatum_inputs = 16
oriens_inputs = 41000

num_patterns = 1
num_degradations = 10
num_bins = 20
degradation = 0.1
bins = np.linspace(0, 1, num_bins + 1)

# Network Topography
# initialize with no connections
connection_matrix = [[0 for ii in range(num_cells)] for ii in range(num_cells)]

# fully connected with self connections
connection_matrix = [[1 for ii in range(num_cells)] for ii in range(num_cells)]

## sparsely connected network
#connection_matrix = [[random.binomial(1, sparsity) for ii in range(num_cells)] for ii in range(num_cells)]

# remove self connections
for k in range(num_cells):
    connection_matrix[k][k] = 0
    
cell_list = []
for ii in range(num_cells):
    pc_ca3 = cell_network.Cell(ID,location,synvars,celltype,fname_morph,num_input,num_cells,modeltype)
    cell_list.append(pc_ca3)

# variables
a = 0.01
num_patterns = 1
num_exposures = 10
num_degradations = 10
num_neurons = num_cells
num_bins = 20     # bins for error bar groupings of correlations
degradation = 0.1
stim_length = 10
stim_delay = 40
max_frequency = 100

bins = np.linspace(0,1, num_bins + 1)

vecstims = [h.VecStim() for i in range(num_cells)]
evecs = [h.Vector() for i in range(num_cells)]


# Initialize data lists
initial_cue_correlations = []
final_state_correlations = []
initial_cue_correlations_binned = [[] for i in range(num_bins)]
final_state_correlations_binned = [[] for i in range(num_bins)]
y_equals_x = [0,1]

# Create unpatterned input
initializing_time = 150
tstop = initializing_time + (stim_length+stim_delay)*(num_exposures + num_degradations)
average_frequency = max_frequency/2
unpatterned_intervals = [0]*num_neurons
unpatterned_spikes = [0]*num_neurons
for i in range(num_neurons):
    unpatterned_intervals[i] = generate_spike_times(average_frequency, initializing_time)
    unpatterned_spikes[i] = np.cumsum(unpatterned_intervals[i])[:]

# Calculate the quality of retrieval of pattern
true_pattern = [0]*num_neurons
retrieved_pattern = [0]*num_neurons
input_pattern = [0]*num_neurons

true_spikes = [0]*num_neurons
true_spikes_full = [0]*num_neurons
degraded_spikes = [0]*num_neurons
complete_input_vector = [0]*num_neurons

# Define pattern of firing rates
for j in range(num_patterns):
    for i in range(num_neurons):
        true_pattern[i] = random.random()*max_frequency
        retrieved_pattern[i] = -99
        while retrieved_pattern[i] < 0:
            retrieved_pattern[i] = true_pattern[i] + random.gauss(0, 0.1)

        spike_intervals = generate_spike_times(true_pattern[i], stim_length)
        true_spikes[i] = np.cumsum(spike_intervals)[:]

        true_spikes_full[i] = generate_full_spike_times(spike_intervals, true_spikes[i], stim_length, stim_delay, num_exposures)
        if len(true_spikes[i]) != 0:
            true_hanging_interval = stim_length - true_spikes[i][-1]
        else:
            true_hanging_interval = stim_length*num_exposures
        degraded_spikes_full = []
        hanging_interval = 0
        for k in range(num_degradations):
            if k == 0:
                old_pattern = true_pattern[i]
            else:
                old_pattern = input_pattern[i]
            input_pattern_tmp = -99
            while input_pattern_tmp < 0:
                input_pattern_tmp = old_pattern + random.gauss(0, degradation)
                input_pattern[i] = input_pattern_tmp
            degraded_intervals = generate_spike_times(input_pattern[i], stim_length)
            degraded_spikes[i] = np.cumsum(degraded_intervals)[:]
            if len(degraded_spikes[i]) != 0:
                degraded_intervals[0] = degraded_intervals[0] + hanging_interval + stim_delay
                hanging_interval = stim_length - degraded_spikes[i][-1]
                for spike in degraded_intervals:
                    degraded_spikes_full.append(spike)
        if len(degraded_spikes_full) != 0:
            degraded_spikes_full[0] = degraded_spikes_full[0] - stim_delay + true_hanging_interval

        if len(true_spikes_full[i]) != 0:
            hanging_interval2 = initializing_time - unpatterned_spikes[i][-1]
            true_spikes_full[i][0] = true_spikes_full[i][0] + hanging_interval2
            complete_input_vector[i] = unpatterned_intervals[i] + true_spikes_full[i] + degraded_spikes_full
        else:
            complete_input_vector[i] = unpatterned_intervals[i]

        evecs[i] = h.Vector(complete_input_vector[i])
        vecstims[i].play(evecs[i])
        # Calculate Cue Correlations
        final_state_correlation, _ = pearsonr(true_pattern, retrieved_pattern)
        initial_cue_correlation, _ = pearsonr(true_pattern, input_pattern)
        initial_cue_correlations.append(initial_cue_correlation)
        final_state_correlations.append(final_state_correlation)

        if initial_cue_correlation >= 0 and initial_cue_correlation <= 1:
            for ii in range(num_bins):
                if initial_cue_correlation > bins[ii] and initial_cue_correlation <= bins[ii+1]:
                    initial_cue_correlations_binned[ii].append(initial_cue_correlation)
                    final_state_correlations_binned[ii].append(final_state_correlation)
                    break

# Calculate stats for binned correlations
initial_cue_correlations_means = []
initial_cue_correlations_stddev = []
final_state_correlations_means = []
final_state_correlations_stddev = []
for ii in range(num_bins):
    if initial_cue_correlations_binned[ii] != []:
        initial_cue_correlations_means.append(np.mean(initial_cue_correlations_binned[ii]))
        initial_cue_correlations_stddev.append(np.std(initial_cue_correlations_binned[ii]))
        final_state_correlations_means.append(np.mean(final_state_correlations_binned[ii]))
        final_state_correlations_stddev.append(np.std(final_state_correlations_binned[ii]))

# Plot correlations
fig = plt.figure()
plt.errorbar(initial_cue_correlations_means, final_state_correlations_means, yerr=final_state_correlations_stddev, linestyle='none', marker='o', capsize=3)
plt.plot(y_equals_x, y_equals_x, color='black')
plt.xlabel('Initial cue correlation', fontsize=16)
plt.ylabel('Final cue correlation', fontsize=16)
plt.show()

#####################
# Connecting inputs #
#####################
# Synaptic weight
#weight = 1.049365e-03
stim_weight = (1.049365e-03)*90
weight = (3.603577e-04)*0.01

stim_weight = (1.049365e-03)*1000
weight = (3.603577e-04)*1000
# NMDA parameters
ratio = 0.0775389348453*0.872
tau1 = 20.3877436154
tau2 = 26.6830234133
tau3 = 158.729359569
wtau2 = 0.963468100127
wtau3 = 1-wtau2

t0 = np.linspace(0,100,10000)
left = -np.log(tau1)-t0/tau1
right = np.log(wtau2/tau2*np.exp(-t0/tau2) + wtau3/tau3*np.exp(-t0/tau3))
tp = t0[np.argmin(np.abs(left-right))]
factor = -np.exp(-tp/tau1) + wtau2*np.exp(-tp/tau2) + wtau3*np.exp(-tp/tau3)
factor = 1/factor

stim_cons = []
net_cons = []

for jj in range(num_cells):
    pc_ca3 = cell_list[jj]
#    nc = h.NetCon(vecstims[jj], pc_ca3.synGroups['AMPA']['lucidum'][0])
#    pc_ca3.synGroups['AMPA']['lucidum'][0].tau1 = 0.1
#    pc_ca3.synGroups['AMPA']['lucidum'][0].tau2 = 144.031250
#    pc_ca3.synGroups['AMPA']['lucidum'][0].e = 0
#    nc.weight[0] = stim_weight
#    nc.delay = 0
#    pc_ca3.synGroups['AMPA']['lucidum'][0].k = weight*4/5
#    stim_cons.append(nc)
    for ii in range(num_input):
       if ii == jj:
           nc = h.NetCon(vecstims[ii], pc_ca3.synGroups['NMDA']['lucidum'][ii])
           pc_ca3.synGroups['NMDA']['lucidum'][ii].tau1 = tau1
           pc_ca3.synGroups['NMDA']['lucidum'][ii].tau2 = tau2
           pc_ca3.synGroups['NMDA']['lucidum'][ii].tau3 = tau3
           pc_ca3.synGroups['NMDA']['lucidum'][ii].wtau2 = wtau2
           pc_ca3.synGroups['NMDA']['lucidum'][ii].factor = factor
           pc_ca3.synGroups['NMDA']['lucidum'][ii].e = 0
           nc.weight[0] = stim_weight*ratio
           nc.delay = 0
           stim_cons.append(nc)
    
for kk in range(num_cells):
    post_cell = cell_list[kk]
    for jj in range(num_cells):
        if connection_matrix[jj][kk] == 1:
           pre_cell = cell_list[jj]
   
           nc = h.NetCon(pre_cell.soma(0.5)._ref_v, post_cell.synGroups['AMPA']['radiatum'][jj], sec=pre_cell.soma)
           post_cell.synGroups['AMPA']['radiatum'][jj].tau1 = 0.1
           post_cell.synGroups['AMPA']['radiatum'][jj].tau2 = 8.384766
           post_cell.synGroups['AMPA']['radiatum'][jj].e = 0
           nc.weight[0] = weight
           nc.delay = 0
           post_cell.synGroups['AMPA']['radiatum'][jj].k = weight*8
           net_cons.append(nc)
       
           nc = h.NetCon(pre_cell.soma(0.5)._ref_v, post_cell.synGroups['NMDA']['radiatum'][jj], sec=pre_cell.soma)
           post_cell.synGroups['NMDA']['radiatum'][jj].tau1 = tau1
           post_cell.synGroups['NMDA']['radiatum'][jj].tau2 = tau2
           post_cell.synGroups['NMDA']['radiatum'][jj].tau3 = tau3
           post_cell.synGroups['NMDA']['radiatum'][jj].wtau2 = wtau2
           post_cell.synGroups['NMDA']['radiatum'][jj].factor = factor
           post_cell.synGroups['NMDA']['radiatum'][jj].e = 0
           nc.weight[0] = weight*ratio
           nc.delay = 0
           net_cons.append(nc)        

################################
# Setting up vectors to record #
################################
v = []
t = []

in_nc_spikes = []
out_nc_spikes = []

for ii in range(num_cells):
    v.append(h.Vector())
    v[-1].record(cell_list[ii].soma(0.5)._ref_v)
    
    t.append(h.Vector())
    t[-1].record(h._ref_t)
    
for ii in range(num_input):
    in_nc_spikes.append(h.Vector())
    stim_cons[ii].record(in_nc_spikes[-1])
for ii in range(num_cells):
    out_nc_spikes.append(h.Vector())
    net_cons[ii].record(out_nc_spikes[-1])

#########################
# Setting up simulation #
#########################
h.v_init = -66
h.t = 0
h.dt = 0.025
h.celsius = 35.0
h("tstep = 0")
h("period = 2")
h.tstop = tstop
h("steps_per_ms = 10")
#h.cvode_active(1)
h.load_file("negative_init.hoc")

weight_changes = [[] for x in range(num_cells)]
for jj in range(num_cells):
    pc_ca3 = cell_list[jj]
    for ii in range(radiatum_inputs):
    	weight_changes[jj].append(h.Vector())
    #	print (net_cons[ii].syn())
    	weight_changes[jj][-1].record(pc_ca3.synGroups['AMPA']['radiatum'][ii]._ref_deltaw,200)

##################
# Run simulation #
##################
print "Starting...!"
ST = cookie.time()
h.run()
ET = cookie.time()-ST
print "Finished in %f seconds" % ET

########
# Plot #
########
for ii in range(num_cells):
    if ii <= 0:
       _=plt.plot(t[ii],v[ii])
       _=plt.xlabel('Time (ms)')
       _=plt.ylabel('Somatic Voltage (mV)')
   #    plt.xlim(0,2000)
       plt.show()
