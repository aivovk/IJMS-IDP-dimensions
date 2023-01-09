import os, sys

"""
To create settings and scripts:
$ python nup_model_input_writer.py {experiment_name}
To submit scripts:
$ for f in scripts/q_{experiment_name}* ; do sbatch $f ; done

- run this script from the same directory as SLURM_SUBMIT_DIR
- .conf settings files will be placed in ./settings
- job submission scripts will be placed in ./scripts
- output_dir should be set to your $SCRATCH directory

"""

script_name = sys.argv[0]
experiment_name = sys.argv[1]

output_dir = 'results/' + experiment_name
os.system ("mkdir %s" % (output_dir))
os.system ("cp %s %s" % (script_name, output_dir))

# rename to specific executable
exe = 'nup_model_'

steps= 2 * 10**9 
eqSteps = 1000

nupfiles = [100]
coulombstrength = [0]
debyelength = [2.225]
cohesivestrength = [round(0.5*i, 1) for i in range(1,11)] + [round(5+0.1*i, 1) for i in range(1,11)]
hydrodynamicradius = [1]
runrepeat = range(240)

# make a .conf settings file for each combination of parameters
settings = [nupfiles, coulombstrength, cohesivestrength, hydrodynamicradius, runrepeat]

settings_suffix = []
nupfile_settings = []
coulombstrength_settings = []
cohesivestrength_settings = []
hydrodynamicradius_settings = []
runrepeat_settings = []

for s0 in settings[0]:
    for s1 in settings[1]:
        for s2 in settings[2]:
            for s3 in settings[3]:
                for s4 in settings[4]:
                        settings_suffix.append(str(s0)+'_'+str(s1)+'_'+str(s2)+'_'+str(s3) +'_' +str(s4))
                        nupfile_settings.append(s0)
                        coulombstrength_settings.append(s1)
                        cohesivestrength_settings.append(s2)
                        hydrodynamicradius_settings.append(s3)
                        runrepeat_settings.append(s4)

              
for s in range(len(settings_suffix)):
    f = open('settings/'+experiment_name+'_'+settings_suffix[s]+'.conf', 'w')
    f.write('outSetting=FILE\n')
    f.write('outFile='+output_dir+'/'+experiment_name+'_'+settings_suffix[s]+'\n')
    f.write('checkpointFile='+output_dir+'/checkpoint_'+experiment_name+'_'+settings_suffix[s]+'\n')
    f.write('snapshotFile='+output_dir+'/'+experiment_name+'_'+settings_suffix[s]+'.bin\n')
    f.write('stepsBetweenSnapshot=10000\n')
    f.write('separationFile='+output_dir+'/'+experiment_name+'_'+settings_suffix[s]+'.sep\n')
    f.write('dt=0.001\n')
    f.write('avgBondLength=1.35\n')
    f.write('feneLength=2\n')
    f.write('MAX_FORCE=0.2\n')
    f.write('steps='+str(steps)+'\n')
    f.write('eqSteps='+str(eqSteps)+'\n')
    f.write('stepsBetweenFileOutput=1000\n')
    f.write('AAFile=Nups/GP.aa\n')
    f.write('eLJ=1\n')
    f.write('skip=0\n')
    f.write('hydrophobicStrength='+str(cohesivestrength_settings[s])+'\n')
    f.write('coulombStrength='+str(coulombstrength_settings[s])+'\n')
    f.write('debyeLength='+str(debyelength[0])+'\n')
    f.write('nupFile=GPP'+str(nupfile_settings[s])+'.txt\n')
    f.write('typeForceBond=TYPE_FORCE_BOND_FENE\n')
    f.write('typeTimeStep=TYPE_TIME_STEP_PERSISTENT\n')
    f.write('neighbourListBuffer=0.8\n')
    f.write('hydrophobicCutoff=4\n')
    f.write('skipNeighbourUpdate=4\n')
    f.write('includeHydrodynamics=TYPE_HYDRO_RPY\n')
    f.close()


# make a script file for each node which starts a process for each core
cores_per_node = 24

script_prefix_text = """#!/bin/bash
#SBATCH --time=18:00:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={0}
#SBATCH --output={1}/%j.out

cd $SLURM_SUBMIT_DIR

module load intel/2018.3 gsl/2.5
"""

for node_index in range(int(len(settings_suffix)/cores_per_node)):
    f = open('scripts/q_'+experiment_name+'_'+str(nupfiles[0])+'_'+str(node_index),'w')
    
    f.write(script_prefix_text.format(int(cores_per_node), output_dir))
    
    for core_index in range(cores_per_node):
    	f.write("./"+exe+" settings/"+experiment_name+"_"+settings_suffix[node_index*cores_per_node
                                                                          +core_index]+".conf &\n")
    f.write("wait\n")
    f.close()
