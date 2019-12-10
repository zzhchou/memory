import random
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

# variables
a = 0.01
num_patterns = 10
num_degradations = 100
num_neurons = 50
num_bins = 20     # bins for error bar groupings of correlations
degredation = 0.1

bins = np.linspace(0,1, num_bins + 1)

# Initialize data lists
initial_cue_correlations = []
final_state_correlations = []
initial_cue_correlations_binned = [[] for i in range(num_bins)]
final_state_correlations_binned = [[] for i in range(num_bins)]
y_equals_x = [0,1]

# Calculate the quality of retrieval of pattern
true_pattern = [0]*num_neurons
retrieved_pattern = [0]*num_neurons
input_pattern = [0]*num_neurons

# Define pattern of firing rates
for j in range(num_patterns):
   for i in range(num_neurons):
      true_pattern[i] = random.random()
      retrieved_pattern[i] = true_pattern[i] + random.gauss(0, 0.1)
   for k in range(num_degradations):
      if k == 0:
         old_pattern = true_pattern
      else:
         old_pattern = input_pattern
      for i in range(num_neurons):
         input_pattern[i] = old_pattern[i] + random.gauss(0, degredation)
         b = random.random()
         retrieved_pattern[i] = (b*input_pattern[i] + (1-b)*true_pattern[i])/1 + random.gauss(0, 0.2)

      
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
#plt.title('alpha = {}'.format(a))
plt.show()
