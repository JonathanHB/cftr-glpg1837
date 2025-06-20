import os
import sys


for simulation in os.listdir():
	if simulation == "nonlip_glpg_1":
		continue

	indep_trjs = [s for s in os.listdir(simulation) if s[-18:] == "ancestors-2.5A-20A"]
	print(indep_trjs)
	for indtrj in indep_trjs:

		numstring = indtrj[0:13]
		if True: #not os.path.exists(f"./{simulation}/{indtrj}/{numstring}-trj.xtc"):

			os.chdir(f"./{simulation}/{indtrj}")
		        
			os.system(f"gmx trjcat -f *-traj_comp.xtc -o {numstring}-trj.xtc -sort -cat")
			os.system(f"echo   1 0 | gmx trjconv -f {numstring}-trj.xtc                     -o {numstring}-trj-pbcmol-centered-tmd.xtc     -s ../topology/eq19.tpr -n ../topology/index_withtmd.ndx -center -pbc mol")
			
			if simulation[0:3] == "lip":
				os.system(f"echo 20 0 | gmx trjconv -f {numstring}-trj-pbcmol-centered-tmd.xtc -o {numstring}-trj-pbcmol-centered-tmd-rot.xtc -s ../topology/eq19.tpr -n ../topology/index_withtmd.ndx -fit rot+trans")
                        
			elif simulation[0:3] == "non":
                                os.system(f"echo 21 0 | gmx trjconv -f {numstring}-trj-pbcmol-centered-tmd.xtc -o {numstring}-trj-pbcmol-centered-tmd-rot.xtc -s ../topology/eq19.tpr -n ../topology/index_withtmd.ndx -fit rot+trans")			
			
			os.chdir("../../")

