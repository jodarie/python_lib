import subprocess

#idpart = 72839
#name = "job.pbs"

#command = ["qsub", "-W", "depend=afterok:", str(idpart), " ", name]

command = ["sort", "-nrk 1", "crap"]
    
p = subprocess.Popen(command, stdout=subprocess.PIPE)

print 'Result: '
print p.communicate()[0]
