temp_list = ['200','623','800','1400']
for temp in temp_list: 
  struct_nvt = job_mini.get_structure(iteration_step=-1)    # Last frame from minimization
  job_name = 'boilot_443K_2_05_NVT' + "_" + temp
  job_nvt = pr.create_job("Lammps", job_name, delete_existing_job=True)
  job_nvt.structure = struct_nvt
  job_nvt.potential = padone_potential
  job_nvt.calc_md(temperature=temp, pressure=0, n_ionic_steps=100, n_print=10, time_step=1.0)
  #print(job_mini.input.control)
  job_nvt.server.queue = "normal_l2" ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
  job_nvt.server.cores = 24 ### The number of cores you want for your job
  job_nvt.run()    
  # Check input file for NVT 
  #print(job_nvt.input.control)
  
# Wait and transfer all the jobs back to local 
count = 1
df = pr.job_table()
while job_name in df[df.status=="submitted"].job.to_list():
  count += 1
  job_mini.status.collect=True
  time.sleep(5)
#  pr.update_from_remote(recursive=True)
#  pr.refresh_job_status()
  if not job_name in df[df.status=="submitted"].job.to_list():
    job_nvt.transfer_from_remote()
    break
#  job_nvt.compress()   # Compress evrything to tar file

# 3. NPT equillibration 
# 4. Production NPT (5 Replicas each calculations)
# Note:- In analysis take the average of these properties as result. 

#for x in count:
job_status_lst = [
    "initialized",
    "appended",
    "created",
    "submitted",
    "running",
    "aborted",
    "collect",
    "suspended",
    "refresh",
    "busy",
    "finished",
    "not_converged",
    "warning"
]

#while True:
for x in count:
  print('start')
  job_mini.status.collect=True
  pr.update_from_remote(recursive=True)
  df = pr.job_table()
  print(df)
  #count += 1
  if all(df.status.isin(["finished", "aborted", "not_converged"])):
    print('here_2')
    finished = True
    job_mini.status.collect=True
    print('Vaibha')
    break  
  count.append(x+1)
  time.sleep(5)
    #pr.update_from_remote(recursive=True)
    #job_mini.refresh_job_status()
    

  
job_mini.status.collect=True
job_mini.compress() 
print('break')