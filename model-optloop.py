# Example of changing the default optmization schedule
from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()

# Give less weight to all soft-sphere restraints:
env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)
env.io.atom_files_directory = ['.']

a = loopmodel(env,
              alnfile  = 'alignment_3.ali',     # alignment filename
              knowns   = '1ds6',              # codes of the templates
              sequence = 'crip')              # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model
                                    # (determines how many models to calculate)
a.md_level = None                   # No refinement of model

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 10           # Last loop model
a.loop.md_level       = refine.slow # Loop model refinement level

# Very thorough VTFM optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 300

# Thorough MD optimization:
a.md_level = refine.slow

# Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
a.repeat_optimization = 10
a.max_molpdf = 1e6

a.make()
