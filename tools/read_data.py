import numpy as np
import csv
"""
read in a series of simulation output files and save the results as numpy arrays
with shapes:
{prefix}_data_lengths [EPS, RPT]
{prefix}_time_series [EPS, PRT, num_props, max_length]
{prefix}_data_avgs [EPS, RPT, num_props]

notes:
-average is of each time series not over repeated runs
-average does not include first skip_avg lines of each file
-time series still includes the first skip_avg lines
"""

output_dir = '..'
experiment_prefix = '2019aug01gpp'

# sequence identifier (number of monomers)
N=100

# # of repeated runs
RPT =240

# cohesive strength
EPS = [round(0.5*i, 1) for i in range(1,11)] + [round(5+0.1*i, 1) for i in range(1,11)]

max_length = 20000 # max # of lines in simulation output files
skip_avg = 2000 # lines to skip in averaging
num_props = 13 # columns in simulation output files

data_lengths = np.zeros((len(EPS), RPT), dtype = int)
time_series = np.empty((len(EPS), RPT, num_props, max_length))
time_series.fill(np.nan)
data_avgs = np.zeros((len(EPS), RPT, num_props))

for ieps in range(len(EPS)):
    eps = EPS[ieps]
    for rpt in range(RPT):
        # see settings_suffix in input_writer.py
        file_name = output_dir +'/'
        +experiment_prefix
        +'_{0}_0_{1:1.1f}_1_{2}.txt'.format(N,eps,rpt)
        
        with open(file_name, 'r') as f:
            print(file_name)
            csv_iter = csv.reader(f,
                           delimiter = ',')
            header = next(csv_iter)

            np_data = np.asarray([csv_line[:num_props] for csv_line in csv_iter],
                                 dtype = float)
            
            data_lengths[ieps, rpt] = np.shape(np_data)[0] # <= max_length
          
            time_series[ieps, rpt, :,:data_lengths[ieps,rpt]] = np.transpose(np_data)
            data_avgs[ieps, rpt, ...] = np.transpose(np.mean(np_data[skip_avg:,:], axis = 0))

np.save(experiment_prefix+'EPS', EPS)
np.save(experiment_prefix + '_data_lengths', data_lengths)
np.save(experiment_prefix + '_data_avgs', data_avgs)
np.save(experiment_prefix + '_time_series', time_series)
