delete all
reset

load /home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_1/topology/input.gro, input_ref

load /home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_1/topology/input.gro, input1
load /home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_1/topology/input.gro, input2
load /home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_1/topology/input.gro, input3

load_traj /home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_1/well-001354-000029-traj_comp.xtc, object=input1, start=3
load_traj /home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_1/well-001417-000021-traj_comp.xtc, object=input2, start=3
load_traj /home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_2/well-000715-000025-traj_comp.xtc, object=input3, start=3

align input1, input_ref and poly
align input2, input_ref and poly
align input3, input_ref and poly

hide everything

show cart
show sticks, (resn LJP or poly) and not elem H

util.cbag
util.cbay input1 and resn LJP
util.cbac input2 and resn LJP
util.cbam input3 and resn LJP

center resn LJP
