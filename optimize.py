# Example for: conjugate_gradients(), molecular_dynamics(), model.switch_trace()

# This will optimize stereochemistry of a given model, including
# non-bonded contacts.

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions

env = environ()
env.io.atom_files_directory = ['.']
env.edat.dynamic_sphere = True

env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

for i in range(1,9):
    for j in range(1,6):
        if i < 10:
            if j < 10:
                code = "crip.BL000"+str(i)+"000" + str(j)+".pdb"
            if j >= 10:
                code = "crip.BL000"+str(i)+"00" + str(j)+".pdb"
        if i >= 10:
            if j < 10:
                code = "crip.BL00"+str(i)+"000" + str(j)+".pdb"
            if j >= 10:
                code = "crip.BL00"+str(i)+"00" + str(j)+".pdb"
        mdl = complete_pdb(env, code)
        mdl.write(file=code+'.ini')

        # Select all atoms:
        atmsel = selection(mdl)

        # Generate the restraints:
        mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
        mdl.restraints.add(secondary_structure.strand(mdl.residue_range('7:', '16:')))
        mdl.restraints.add(secondary_structure.strand(mdl.residue_range('23:', '27:')))
        mdl.restraints.add(secondary_structure.strand(mdl.residue_range('36:', '41:')))

        mdl.restraints.add(secondary_structure.strand(mdl.residue_range('45:', '53:')))
        mdl.restraints.add(secondary_structure.strand(mdl.residue_range('57:', '64:')))

        mdl.restraints.add(secondary_structure.strand(mdl.residue_range('67:', '70:')))
        mdl.restraints.add(secondary_structure.strand(mdl.residue_range('82:', '90:')))

        mdl.restraints.add(secondary_structure.strand(mdl.residue_range('105:', '111:')))
        mdl.restraints.add(secondary_structure.strand(mdl.residue_range('119:', '125:')))

        mdl.restraints.write(file=code+'.rsr')

        mpdf = atmsel.energy()

        # Create optimizer objects and set defaults for all further optimizations
        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')

        # Open a file to get basic stats on each optimization
        trcfil = open(code+'.D00000001', 'w')

        # Run CG on the all-atom selection; write stats every 5 steps
        cg.optimize(atmsel, max_iterations=20, actions=actions.trace(5, trcfil))

        mpdf = atmsel.energy()

        mdl.write(file=code+'_opt.pdb')