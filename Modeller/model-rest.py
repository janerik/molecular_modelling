# Example of changing the default optmization schedule
from modeller import *
from modeller.automodel import *
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions

#log.verbose()
env = environ()

# Give less weight to all soft-sphere restraints:
env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)
env.io.atom_files_directory = ['.']

class MyModel(loopmodel):
    def special_restraints(self, aln):
        rsr = self.restraints

        rsr.add(secondary_structure.strand(self.residue_range('6:', '15:')))

        rsr.add(secondary_structure.strand(self.residue_range('22:', '26:')))
        rsr.add(secondary_structure.strand(self.residue_range('35:', '40:')))

        rsr.add(secondary_structure.strand(self.residue_range('44:', '52:')))
        rsr.add(secondary_structure.strand(self.residue_range('56:', '63:')))

        rsr.add(secondary_structure.strand(self.residue_range('66:', '69:')))
        rsr.add(secondary_structure.strand(self.residue_range('81:', '89:')))

        rsr.add(secondary_structure.strand(self.residue_range('104:', '110:')))
        rsr.add(secondary_structure.strand(self.residue_range('118:', '124:')))


a = MyModel(env,
              alnfile  = 'alignment_4.ali',     # alignment filename
              knowns   = '1ds6',              # codes of the templates
	      assess_methods=(assess.DOPE),
              sequence = 'crip')              # code of the target


a.starting_model= 1                 # index of the first model
a.ending_model  = 10                 # index of the last model

a.library_schedule = autosched.slow
a.max_var_iterations = 300

#a.md_level = refine.slow  

a.repeat_optimization = 2
a.max_molpdf = 1e6
                                    # (determines how many models to calculate)
#a.md_level = refine.slow                   # No refinement of model

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 5           # Last loop model
a.loop.md_level       = refine.slow # Loop model refinement level

a.make()                            # do comparative modeling
