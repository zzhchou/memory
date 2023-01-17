    ###################################################################
# 2019-05-08                                                      #
# Code to test implementation of log-normal STDP synapse in a CA3 #
# pyramidal cell                                                  #
###################################################################
from __future__ import division
from neuron import h
from scipy.stats import pearsonr
import time as cookie
import pickle
import numpy as np
import cell_network
import matplotlib.pyplot as plt
import numpy.random as random
h.load_file("stdrun.hoc")
np.random.seed(149)

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

random_input_period = 500

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
    

################################
# Create spike times for input #
################################


max_frequency = 5.5 # units: Hz
min_frequency = 0.5

vecstims = [h.VecStim() for ii in range(num_input*num_cells)]
evecs = [h.Vector() for ii in range(num_input*num_cells)]

true_firing_rates = [0 for ii in range(num_input*num_cells)]

# Rate Dependent Inputs
stim_length = 50 # units: ms
num_reps = 10
rep_interval_delay = 200 # units: ms

#Constants
num_patterns = 1
num_degradations = 10
num_bins = 20
degradation = 0.1   
bins = np.linspace(0, 1, num_bins + 1)

tstop = random_input_period + (num_reps + num_degradations)*(stim_length + rep_interval_delay)# unts: ms

stim_rates = []
spike_counts = []
for ii in range(num_input*num_cells):
    intervals= []
    frequency = np.random.uniform(min_frequency, max_frequency)
    true_firing_rates[ii] = frequency
    stim_rates.append(frequency)
    mu = 1000./frequency
    elapsed_time = 0
    flag = 1
    while flag:
        roll = np.random.uniform(0,1)
        interval = -mu*np.log(roll)
        
        elapsed_time += interval
        
        if elapsed_time > stim_length:
            flag = 0
        else:
            intervals.append(interval)
            
#            spikes = np.cumsum(intervals)[:-1]
#            for spike in spikes:
#                evecs[ii].append(spike)
            
#        intervals = [1]
    
    rep_intervals = list(intervals)
    last_interval = 0
    for j in range(num_reps - 1):
        if len(intervals) == 0:
            last_interval = last_interval - stim_length
        else:
            last_interval = intervals[-1]
            tmp_interval = rep_interval_delay - last_interval
            last_interval = 0
            elapsed_time += tmp_interval
    
            for interval_copy in intervals:
                elapsed_time += interval_copy + tmp_interval
                rep_intervals.append(interval_copy + tmp_interval)
                tmp_interval = 0

    if len(intervals) == 0:
        last_interval = 0
    else:
        last_interval = intervals[-1]
    for k in range(num_degradations):
        if k == 0:
           old_frequency = frequency
        else:
           old_frequency = degraded_frequency
        degraded_frequency = old_frequency + random.normal(0, degradation)
        while degraded_frequency < 0:
            degraded_frequency = old_frequency + random.normal(0, degradation)

        mu = 1000./frequency
        elapsed_time = 0
        flag = 1        
        intervals = []

        while flag:
            roll = np.random.uniform(0,1)
            interval = -mu*np.log(roll)
            
            elapsed_time += interval
            
            if elapsed_time > stim_length:
                flag = 0
            else:
                intervals.append(interval)
                
        if len(intervals) == 0:
            last_interval = last_interval - stim_length
        else:
            
            tmp_interval = rep_interval_delay - last_interval
            last_interval = 0
            for interval_copy in intervals:
                rep_intervals.append(interval_copy + tmp_interval)
                tmp_interval = 0
            last_interval = last_interval + intervals[-1]
                
            
            

    spikes = np.cumsum(rep_intervals)[:-1]
    for spike in spikes:
        evecs[ii].append(spike)
            
    vecstims[ii].play(evecs[ii])
            
    spike_counts.append(len(spikes))

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

true_input_patterns = []
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
           true_input_pattern = true_firing_rates[jj*num_input + ii]
           true_input_patterns.append(true_input_pattern)
           nc = h.NetCon(vecstims[jj*num_input + ii], pc_ca3.synGroups['NMDA']['lucidum'][ii])
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



training_reps = num_reps
training_delay = rep_interval_delay
combined_stim_time = (training_reps+num_degradations)*num_patterns*(stim_length+training_delay)
weight_array = np.array(weight_changes)

rep_response = []
out_spike_rate = []
for gg in range(training_reps):
    last_rep = gg*num_patterns*(stim_length+training_delay)
    last_rep_end = (gg + 1)*num_patterns*(stim_length + training_delay)
    
    trained_post_spikes = []
    t_array = np.linspace(0, combined_stim_time, np.ma.size(weight_array, 2))        
    last_rep_t = len(v[0])*last_rep/combined_stim_time
    last_rep_t_end = len(v[0])*last_rep_end/combined_stim_time
    
    for ii in range(len(weight_array)):
        spike_intervals = []
            
        post_spikes = []
        v_array = np.array(v[ii])
        t2_array = np.array(t[ii])
        for kk in range(int(last_rep_t), int(last_rep_t_end)-1):
            if v_array[kk] > 0 and v_array[kk] > v_array[kk+1] and v_array[kk] > v_array[kk-1]:
                post_spikes.append(t2_array[kk])

        trained_post_spikes.append(post_spikes)
    rep_response.append(trained_post_spikes)
out_spike_rate = ([[len(x)/stim_length for x in rep] for rep in rep_response])

#weight_changes = np.array(weight_changes)
#diff_weights = [0]*len(net_cons)
#for k in range(len(net_cons)):
#    diff_weights[k] = weight_changes[k][-1] - weight_changes[k][0]
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

#post_spikes = []
#v_array = np.array(v)
#t_array = np.array(t)
#for ii in range(len(v_array)):
#    if v_array[ii] > 0 and v_array[ii] > v_array[ii+1] and v_array[ii] > v_array[ii-1]:
#        post_spikes.append(t_array[ii])
#
#spike_delays = []
#times = np.linspace(0, tstop, num_reps + 1)
#for jj in range(len(times) - 1):
#    for ii in range(len(post_spikes) - 1):
#        if jj == 0:
#            spike_delays.append(post_spikes[ii])
#            break
#        elif post_spikes[ii] < times[jj] and post_spikes[ii+1] > times[jj]:
#            spike_delays.append(post_spikes[ii+1]- times[jj])
#            break

weight_array = np.array(weight_changes)
weight_array = weight_array + weight
w1 = weight_array[:,:,-1].ravel()
w2 = w1[w1 != 0]

#plt.figure()
#plt.hist(w2, 25)
#plt.show()


#t_array = np.linspace(0, tstop, np.ma.size(weight_array, 2))
#for ii in range(len(weight_array)):
#    plt.figure()
#    for jj in range(np.ma.size(weight_array, 1)):
#        plt.plot(t_array, weight_array[ii, jj, :])
#    plt.show()



# Initialize data lists
initial_cue_correlations = []
final_state_correlations = []
initial_cue_correlations_binned = [[] for i in range(num_bins)]
final_state_correlations_binned = [[] for i in range(num_bins)]

retrieved_pattern = np.array(out_nc_spikes)
degraded_input_pattern = [0]*num_cells

# Sort correlations into bins
j_offset = num_reps
for j in range(num_patterns):
   true_pattern = true_input_patterns
   for k in range(num_degradations):
      initial_cue_correlation, _ = pearsonr(true_pattern, degraded_input_pattern)
      final_state_correlation, _ = pearsonr(true_pattern, out_spike_rate[j_offset + j])
         
         
      if initial_cue_correlation >= 0 and initial_cue_correlation <= 1:
         for ii in range(num_bins):
            if initial_cue_correlation > bins[ii] and initial_cue_correlation <= bins[ii+1]:
               initial_cue_correlations_binned[ii].append(initial_cue_correlation)
               final_state_correlations_binned[ii].append(final_state_correlation)
               break
# Plot Retrieval Quality
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

y_equals_x = [0,1]
fig = plt.figure()
plt.errorbar(initial_cue_correlations_means, final_state_correlations_means, yerr=final_state_correlations_stddev, linestyle='none', marker='o', capsize=3)
plt.plot(y_equals_x, y_equals_x, color='black')
plt.xlabel('Initial cue correlation')
plt.ylabel('Final cue correlation')
#plt.title('alpha = {}'.format(a))
plt.show()

###############
# End of file #
###############
