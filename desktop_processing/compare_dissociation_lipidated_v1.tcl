#Jonathan Borowsky
#Grabe lab
#5/6/25

#run this script with source ../cftr-glpg1837/desktop_processing/compare_dissociation_v1.tcl

#note these lists are not comma-delimited
set topologies {lip_glpg_1/topology/input.gro lip_glpg_1/topology/input.gro lip_glpg_2/topology/input.gro}
set trajectories {lip_glpg_1/001798-000087-trj-pbcmol-centered-tmd-rot.xtc lip_glpg_1/002013-000132-trj-pbcmol-centered-tmd-rot.xtc lip_glpg_2/001986-000211-trj-pbcmol-centered-tmd-rot.xtc}

for {set i 0} {$i < 3} {incr i} {

    mol new [lindex $topologies $i] type gro 
    #lip_glpg_1/topology/input.gro
    mol addfile [lindex $trajectories $i] type xtc step 10
    #lip_glpg_1/002013-000132-trj-pbcmol-centered-tmd-rot.xtc step 10


    #mol modselect 0 $i protein
    #mol modstyle 0 $i quicksurf
    #mol modcolor 0 $i name

    #mol modselect 0 $i all and not name \"C.*\" and not name \"H.*\"
    #mol modstyle 0 $i Lines
    #mol modcolor 0 $i name

    mol addrep $i
    mol modselect 1 $i resname LJP and not name \"H.*\"
    mol modstyle 1 $i Licorice 0.300000 12.000000 12.000000
    mol modcolor 1 $i ColorId $i

    mol addrep $i
    mol modselect 2 $i protein and not name \"H.*\" and resid 305 236 232 233 229 309 308 312 313 316 315 304 925 926 927 928 929 930 931 932 933 934 935 936 861 862 863 864 865 866 867 868 869 870 871 872 873 874 875 876 877
    mol modstyle 2 $i Licorice 0.300000 12.000000 12.000000
    mol modcolor 2 $i name

    mol addrep $i
    mol modselect 3 $i protein
    mol modstyle 3 $i Ribbons
    mol modcolor 3 $i ColorID 7

    #mol modselect 4 $i protein and name \"C.*\" and resid 305 236 232 233 229 309 308 312 313 316 315 304 925 926 927 928 929 930 931 932 933 934 935 936 861 862 863 864 865 866 867 868 869 870 871 872 873 874 875 876 877
    #mol modstyle 4 $i Licorice 0.300000 12.000000 12.000000
    #mol modcolor 4 $i ColorID 13

    mol addrep $i
    mol modselect 4 $i name P31
    mol modstyle 4 $i VDW
    mol modcolor 4 $i name

}


#mol addrep 0
#mol modselect 1 0 protein
#mol modstyle 1 0 quicksurf
#mol modcolor 1 0 name

#mol addrep 0
#mol modselect 2 0 name P31
#mol modstyle 2 0 VDW
#mol modcolor 2 0 name

#terminal command to make movie:
#ffmpeg -framerate 2 -i movie_frame_name.%05d.ppm movie_name.mp4
#see https://askubuntu.com/questions/610903/how-can-i-create-a-video-file-from-a-set-of-jpg-images
#and https://superuser.com/questions/820134/why-cant-quicktime-play-a-movie-file-encoded-by-ffmpeg 
#    It suggests  -pix_fmt yuv420p -vcodec libx264 to make mac quicktime player recognize mp4 files, but this does not work for my mac