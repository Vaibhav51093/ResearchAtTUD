import tarfile
import subprocess
import os
import shutil
import csv
import os
import re


class read_OUTCAR_OSICAR():

    def time_vasp(job,filename='OUTCAR', skipSCFCycles = 5, get_energy=True, filename_2='OSZICAR'):

        # Works with Pyiron 
        # 1. Checking for specific lines in OUTCAR file:- 1. LOOP: CPU time, 2. Elapsed time, 3. Potenetail energy     
        _OUTCAR_TSCFRegex = re.compile("LOOP\:\s+cpu time\s+\d+\.\d+\:\s+real time\s+(?P<t_scf>\d+\.\d+)");
        _OUTCAR_TElapsedRegex = re.compile("Elapsed time \(sec\)\:\s+(?P<t_elapsed>\d+\.\d+)");
        _OSZICAR_TotalEnergyRegex = re.compile("E0= (?P<total_energy>[+-]?\d*\.\d+E[+-]?\d+)");

        numSCFSteps, tSCFAve, tElapsed = None, None, None;

        scfTimes = []

        job_file_compressed = os.path.join(job.working_directory, [f for f in os.listdir(job.working_directory) if ".tar.bz2" in f][0])
        os.system('tar -xvf %s %s'%(job_file_compressed,filename))
        os.system('mv OUTCAR' + ' ' + job.working_directory)
        file = job.working_directory + '/%s'%(filename)
    
        with open(file, 'r') as inputReader:
            for line in inputReader:
                match = _OUTCAR_TSCFRegex.search(line);
            
                if match:
                    scfTimes.append(
                        float(match.group('t_scf'))
                        );
                else:
                    match = _OUTCAR_TElapsedRegex.search(line);
                
                    if match:
                        tElapsed = float(match.group('t_elapsed'));

        if skipSCFCycles > 0:
            if len(scfTimes) > skipSCFCycles:
                scfTimes = scfTimes[skipSCFCycles:];
            else:
                print("WARNING: _ParseOUTCAR(): Number of SCF steps {0} <= skipSCFCycles {1}".format(len(scfTimes), skipSCFCycles));
                scfTimes = [];
    
        if len(scfTimes) > 0:                        
            numSCFSteps = len(scfTimes);
            tSCFAve = sum(scfTimes) / numSCFSteps;

        if get_energy == True:

            finalTotalEnergy = None;

            job_file_compressed_1 = os.path.join(job.working_directory, [f for f in os.listdir(job.working_directory) if ".tar.bz2" in f][0])
            os.system('tar -xvf %s %s'%(job_file_compressed_1,filename_2))
            os.system('mv OSZICAR' + ' ' + job.working_directory)
            file_1 = job.working_directory + '/%s'%(filename_2)
    
        with open(file_1, 'r') as inputReader:
            for line in inputReader:
                match = _OSZICAR_TotalEnergyRegex.search(line);
                if match:
                    finalTotalEnergy = float(match.group('total_energy'));
        print('Result contains:- 1. No of SCF steps, 2. Avg SCF time (sec), 3. Elapsed time (sec), 4. Final energy (eV)')
        return numSCFSteps, tSCFAve, tElapsed, finalTotalEnergy

class read_IBZKBT():

    # Get number of Irriducible k_points 
    def IBZ_vasp(job):
        job_file_compressed = os.path.join(job.working_directory, [f for f in os.listdir(job.working_directory) if ".tar.bz2" in f][0])
        tar = tarfile.open(job_file_compressed, "r:bz2")
        for member in tar.getmembers():
            if 'IBZKPT' == member.name:
                f = tar.extractfile(member)
                if f is not None:
                    lines = f.readlines()
                    return (int(float(lines[1])))

