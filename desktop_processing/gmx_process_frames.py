import os
import sys

ref_folder = "reference_index "

for folder in os.listdir("."):
    if folder[-9:] == "ancestors":
        os.cwd(folder)
        numstring = folder[0:14]
        os.system(f"gmx trjcat -f *-traj_comp.xtc -o {numstring}-trj.xtc -sort -cat")
        os.system(f"echo   1 0 | gmx trjconv -f {numstring}-trj.xtc                     -o {numstring}-trj-pbcmol-centered-tmd.xtc     -s topology/eq19.tpr -n index.ndx -center -pbc mol")
        os.system(f"echo 1 1 0 | gmx trjconv -f {numstring}-trj-pbcmol-centered-tmd.xtc -o {numstring}-trj-pbcmol-centered-tmd-rot.xtc -s topology/eq19.tpr -n index.ndx -center -fit rot+trans")
