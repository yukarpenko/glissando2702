#!/bin/bash

# you may replace the postscript viewer with your favorite one
VIEW=gv

# init
function pause(){
   read -p "$*"
}

# complile with standard precompiler options

make clean
make

# running

: <<'END'

ECHO this block is commented out

END

pause "Press ENTER to see the help"
./glissando2 -h

pause "Press ENTER to see the version"
./glissando2 -v

pause "Press ENTER to run the fastest fun example: a single event in the hot-spot model, Pb+Pb @ LHC"
./glissando2 input/input_snap.dat 2> err
cd output
root -b -l -q -x "../macro/density.C(\"glissando.root\")"
$VIEW 3D-density.eps 2> err &
$VIEW 2D-density.eps 2> err &
cd ..

pause "Press ENTER to run central deuteron-Au collisions (2 min.)"
./glissando2 input/input_d_Au.dat output/dAu.root 2> err
cd output
root -b -l -q -x "../macro/info.C(\"dAu.root\")" 
root -b -l -q -x "../macro/density.C(\"dAu.root\")" 
$VIEW 2D-density.eps 2> err &
cd ..

pause "Press ENTER to run minimum bias proton-Pb collisions (1 min. - times for i7-4700 @ 2.4GHz)"

make clean
make 'PREPROCESS = -D_nnwp_=2'

./glissando2 input/input_p_Pb.dat output/pPb.root 2> err
cd output
root -b -l -q -x "../macro/centrality2.C(\"pPb.root\")" 
$VIEW centrality_b.eps 2> err &
$VIEW centrality_nw.eps 2> err &
$VIEW centrality_RDS.eps 2> err &
cd ..

pause "Press ENTER to run a minimum-bias Pb-Pb run for LHC (14 min.)"
./glissando2 input/input_minbias_LHC.dat output/lhc.root 2> err
cd output
root -b -l -q -x "../macro/centrality2.C(\"lhc.root\")"
$VIEW centrality_b.eps 2> err &
$VIEW centrality_nw.eps 2> err &
$VIEW centrality_RDS.eps 2> err &
root -b -l -q -x "../macro/epsilon.C(\"lhc.root\")" 
$VIEW epsilon.eps 2> err &
$VIEW sigma_epsilon.eps 2> err &
root -b -l -q -x "../macro/fourier.C(\"lhc.root\")" 
$VIEW epsn.eps 2> err &
root -b -l -q -x "../macro/size.C(\"lhc.root\")" 
$VIEW sigma_r.eps 2> err &
root -b -l -q -x "../macro/mult.C(\"lhc.root\")" 
$VIEW om_RDS.eps 2> err &
root -b -l -q -x "../macro/hydro.C(\"lhc.root\",\"hydro_lhc.dat\")" 
cd ..

pause "Press ENTER to run the centrality 30-40% run for RHIC at 200 GeV (3 min.)"
make clean
make

./glissando2 input/input_30_40_RHIC.dat output/rhic_30_40.root 2> err
cd output        
root -b -l -q -x "../macro/density.C(\"rhic_30_40.root\")" 
$VIEW 3D-density.eps 2> err &
$VIEW 2D-density.eps 2> err &
root -b -l -q -x "../macro/core_mantle.C(\"rhic_30_40.root\")" 
$VIEW 3D-cm-density.eps 2> err &
$VIEW 2D-cm-density.eps 2> err &
root -b -l -q -x "../macro/profile2.C(\"rhic_30_40.root\")" 
$VIEW f0_fig.eps 2> err &
$VIEW f2_fig.eps 2> err &
$VIEW f4_fig.eps 2> err &
cd ..

pause "Press ENTER to run asymmetric collisions (1 min.)"
./glissando2 input/input_S_Pb_SPS.dat output/SPb.root 2> err
cd output
root -b -l -q -x "../macro/epsilon.C(\"SPb.root\")" 
$VIEW epsilon.eps 2> err &
root -b -l -q -x "../macro/mult.C(\"SPb.root\")" 
$VIEW om_RDS.eps 2> err &
cd ..


pause "Press ENTER to run deformed and spherical U+U, Au+Au and Cu+Cu collisions (5 min.)"
make clean
make 'PREPROCESS = -D_nnwp_=1 -D_profile_=1'

./glissando2 input/input_minbias_UU_mixed_deformed_norot.dat output/output_minbias_UU_mixed_deformed_norot.root 2> err
./glissando2 input/input_minbias_UU_mixed_spherical_norot.dat output/output_minbias_UU_mixed_spherical_norot.root 2> err
./glissando2 input/input_minbias_AuAu_mixed_deformed_norot.dat output/output_minbias_AuAu_mixed_deformed_norot.root 2> err
./glissando2 input/input_minbias_AuAu_mixed_spherical_norot.dat output/output_minbias_AuAu_mixed_spherical_norot.root 2> err
./glissando2 input/input_minbias_Cu63Cu63_mixed_deformed_norot.dat output/output_minbias_Cu63Cu63_mixed_deformed_norot.root 2> err
./glissando2 input/input_minbias_Cu63Cu63_mixed_spherical_norot.dat output/output_minbias_Cu63Cu63_mixed_spherical_norot.root 2> err

cd output 
root -b -l -q -x "../macro/profile2_deformation_U.C(\"output_minbias_UU_mixed_deformed_norot.root\",\"output_minbias_UU_mixed_spherical_norot.root\")"
$VIEW rcostheta_spherical_deformed_U.eps 2> err &
root -b -l -q -x "../macro/profile2_deformation_Au.C(\"output_minbias_AuAu_mixed_deformed_norot.root\",\"output_minbias_AuAu_mixed_spherical_norot.root\")"
$VIEW rcostheta_spherical_deformed_Au.eps 2> err &
root -b -l -q -x "../macro/profile2_deformation_63Cu.C(\"output_minbias_Cu63Cu63_mixed_deformed_norot.root\",\"output_minbias_Cu63Cu63_mixed_spherical_norot.root\")"
$VIEW rcostheta_spherical_deformed_63Cu.eps 2> err &
cd ..

pause "Press ENTER to see more functionality (4 min.)"
make clean
make 'PREPROCESS = -D_nnwp_=2 -D_profile_=1 -D_weight_=1 -D_rapidity_=1'

./glissando2 input/input_aux.dat output/glissando_aux.root 2> err    

cd output   
root -b -l -q -x "../macro/wounding_profile.C(\"glissando_aux.root\")" 
$VIEW wounding_pr.eps 2> err &
root -b -l -q -x "../macro/fitr.C(\"glissando_aux.root\")" 
$VIEW fitr.eps 2> err &
root -b -l -q -x "../macro/corr.C(\"glissando_aux.root\")" 
$VIEW corr.eps 2> err &
root -b -l -q -x "../macro/overlay.C(\"glissando_aux.root\")" 
$VIEW overlay_pr.eps 2> err &
root -b -l -q -x "../macro/tilted.C(\"glissando_aux.root\")" 
$VIEW tilted_pr.eps  2> err &
root -b -l -q -x "../macro/angles.C(\"glissando_aux.root\")" 
$VIEW angles.eps  2> err &
cd ..

pause "Press ENTER for a sample run with nuclear correlations read from external files (1 min.)"
pause "(First, you must get and save files from http://www.phys.psu.edu/∼malvioli/eventgenerator/ and create one big file, for instance running cat o16-1.dat o16-2.dat o16-3.dat [more files] > o16.dat, which must be placed in subdirectory 'nucl' ) Press ENTER"

if [ -f "nucl/o16.dat" ] ;
then
make clean
make 'PREPROCESS = -D_nnwp_=1 -D_files_=1 -D_profile_=1 -D_weight_=1'

./glissando2 input/input_O16.dat output/O16.root nucl/o16.dat nucl/o16.dat 2> err
cd output
root -b -l -q -x "../macro/corr.C(\"O16.root\")" 
$VIEW corr.eps 2> err &
root -b -l -q -x "../macro/wounding_profile.C(\"O16.root\")" 
$VIEW wounding_pr.eps 2> err &
cd ..

else
echo "file nucl/o16.dat does not exist!"
fi

pause "Press ENTER to run the interpolation code"
cd addons
make -f interpolation.mk
./interpolation ../output/SPb.root
make -f interpolation.mk clean
cd ..

pause "Press ENTER to output the event-by-event info to file events.dat (LARGE OUTPUT DATA!, 5 min.)"

# compile with standard precompiler options
make clean
make 'PREPROCESS = -D_evout_=1'
./glissando2


pause "Press ENTER to generate and retrieve the full event info (LARGE OUTPUT DATA!)"
pause "(FULL=1 flag must be set in the GLISSANDO input file) Press ENTER"
pause "This method is alternative to the method with _evout_=1 flag (preferred)"

# compile with standard precompiler options
make clean
make

./glissando2 input/input_full.dat output/full.root 2> err
cd addons
make -f retrieve.mk
./retrieve ../output/full.root
make -f retrieve.mk clean
cd ..

