# Example of changing the default optmization schedule
from modeller import *
from modeller.automodel import *
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions

log.verbose()
env = environ()

env.io.water = True
env.io.hydrogen = True

# Give less weight to all soft-sphere restraints:
#env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)
#env.io.atom_files_directory = ['.']

class MyModel(dope_loopmodel):
    def special_restraints(self, aln):
        rsr = self.restraints

        rsr.add(secondary_structure.strand(self.residue_range('7:', '16:')))

        rsr.add(secondary_structure.strand(self.residue_range('23:', '27:')))
        rsr.add(secondary_structure.strand(self.residue_range('36:', '41:')))

        rsr.add(secondary_structure.strand(self.residue_range('45:', '53:')))
        rsr.add(secondary_structure.strand(self.residue_range('57:', '64:')))

        rsr.add(secondary_structure.strand(self.residue_range('67:', '70:')))
        rsr.add(secondary_structure.strand(self.residue_range('82:', '90:')))

        rsr.add(secondary_structure.strand(self.residue_range('105:', '111:')))
        rsr.add(secondary_structure.strand(self.residue_range('119:', '125:')))


a = MyModel(env,
              alnfile  = 'alignment_4.ali',     # alignment filename
              knowns   = '1ds6',              # codes of the templates
	      #assess_methods=(assess.DOPE),
              sequence = 'crip')              # code of the target


a.starting_model= 1                 # index of the first model
a.ending_model  = 5                 # index of the last model

a.library_schedule = autosched.slow
a.max_var_iterations = 300

#a.md_level = refine.slow  

a.repeat_optimization = 4
a.max_molpdf = 1e6


a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 8           # Last loop model
#a.loop.md_level       = refine.slow # Loop model refinement level

a.make()                            # do comparative modeling
