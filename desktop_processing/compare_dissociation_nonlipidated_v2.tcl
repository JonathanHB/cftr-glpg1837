#Jonathan Borowsky
#Grabe lab
#5/6/25

#run this script with source ../cftr-glpg1837/desktop_processing/compare_dissociation_v1.tcl

#note these lists are not comma-delimited
set topologies {nonlip_glpg_1/topology/input.gro nonlip_glpg_1/topology/input.gro nonlip_glpg_1/topology/input.gro nonlip_glpg_2/topology/input.gro}
set trajectories {nonlip_glpg_1/001119-000101-trj-pbcmol-centered-tmd-rot.xtc nonlip_glpg_1/001655-000156-trj-pbcmol-centered-tmd-rot.xtc nonlip_glpg_1/001913-000187-trj-pbcmol-centered-tmd-rot.xtc nonlip_glpg_2/001042-000195-trj-pbcmol-centered-tmd-rot.xtc}

mol new [lindex $topologies 0] type gro 
mol new [lindex $topologies 1] type gro 


for {set i 2} {$i < 3} {incr i} {

    mol new [lindex $topologies $i] type gro 
    #lip_glpg_1/topology/input.gro
    mol addfile [lindex $trajectories $i] type xtc step 10
    #lip_glpg_1/002013-000132-trj-pbcmol-centered-tmd-rot.xtc step 10
    
    mol addrep $i
    mol modselect 0 $i protein
    mol modstyle 0 $i quicksurf
    mol modcolor 0 $i name
    mol smoothrep $i 0 5 
    #note that the molecule and representation id arguments go in reverse order for smoothrep

    #mol modselect 0 $i all and not name \"C.*\" and not name \"H.*\"
    #mol modstyle 0 $i Lines
    #mol modcolor 0 $i name

    mol addrep $i
    mol modselect 1 $i resname LJP and not name \"H.*\"
    mol modstyle 1 $i Licorice 0.300000 12.000000 12.000000
    mol modcolor 1 $i name #ColorId $i
    mol smoothrep $i 1 5

    mol addrep $i
    mol modselect 2 $i protein and not name \"H.*\" and resid 305 236 232 233 229 309 308 312 313 316 315 304 925 926 927 928 929 930 931 932 933 934 935 936 861 862 863 864 865 866 867 868 869 870 871 872 873 874 875 876 877
    mol modstyle 2 $i Licorice 0.300000 12.000000 12.000000
    mol modcolor 2 $i name
    mol smoothrep $i 2 5

    mol addrep $i
    mol modselect 3 $i protein
    mol modstyle 3 $i Ribbons
    mol modcolor 3 $i ColorID 7
    mol smoothrep $i 3 5

    #mol modselect 4 $i protein and name \"C.*\" and resid 305 236 232 233 229 309 308 312 313 316 315 304 925 926 927 928 929 930 931 932 933 934 935 936 861 862 863 864 865 866 867 868 869 870 871 872 873 874 875 876 877
    #mol modstyle 4 $i Licorice 0.300000 12.000000 12.000000
    #mol modcolor 4 $i ColorID 13

    mol addrep $i
    mol modselect 4 $i name P31
    mol modstyle 4 $i VDW
    mol modcolor 4 $i name
    mol smoothrep $i 4 5

}

#to make movie:
#ffmpeg -start_number 6 -framerate 20 -i frame_file_name.%05d.ppm movie_file_name.mp4