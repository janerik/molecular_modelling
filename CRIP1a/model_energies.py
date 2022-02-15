from modeller import *
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

file = open("energies.txt","w")

for i in range(1, 11):
    for j in range(1, 6):
	# read model file
        if i < 10:
            if j < 10:
                code = "crip1a.BL000"+str(i)+"000" + str(j)+".pdb"
            if j >= 10:
                code = "crip1a.BL000"+str(i)+"00" + str(j)+".pdb"
        if i >= 10:
            if j < 10:
                code = "crip1a.BL00"+str(i)+"000" + str(j)+".pdb"
            if j >= 10:
                code = "crip1a.BL00"+str(i)+"00" + str(j)+".pdb"
        mdl = complete_pdb(env, code)
        s = selection(mdl)
        score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='TvLDH.profile',
                  normalize_profile=True, smoothing_window=15)
        file.write(code + " " + str(score) + "\n")
file.close()
