#!/bin/bash
#! postprocess.sh
#! Script to handle the NEMO outputs

#INPUT PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#Number of domains (nodes) and number of CPUS (cores) 
export NUM_DOM=$XIOSCORES
export NUM_CPU=$OCEANCORES

#Resubmission parameters
#export SUBMDIR=/projects/nexcs-n02/astyles/subm_scripts/WSExp/testing/W10/R1/ #Path to submission script directory
#export SUBMSCRIPT=NEMO_subm_script.pbs #Submission script file name
#export Nresub_spup=0 #Total number of resubmissions for spinup
#export Nresub_exp=0 #Total number of resubmissions experiment
#export SPINUP_INT=4320 #Number of time steps in a single spinup submission
#export EXP_INT=4320 #Number of time steps in final resubmission (post spinup)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#Make sure your namelist is correct:
#nn_it000 = 1
#nn_itend = spinup interval
#nn_stock = spinup interval

#ln_rstart=.false.
#cn_ocersr_in="restart_tmp"


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1) Limited recombination of files using rebuild_nemo. Grid files should be recombined somewhere off of NEXCS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#Extract the time stamp for the restart files from MODEL_*_restart_0000.nc

export RES_TIMESTAMP=$(echo $(ls -d ${MODEL}_*_restart_0000.nc) | awk -F _ '{print $2 }')
echo "RES_TIMESTAMP = "${RES_TIMESTAMP}



# Write a series of namelists to define how to rebuild data files

echo "&nam_rebuild"             > nam_rebuild_grid_T
echo "filebase='${MODEL}_grid_T'" >> nam_rebuild_grid_T
echo "ndomain=${NUM_DOM}"      >> nam_rebuild_grid_T
echo "l_maskout=.false "       >> nam_rebuild_grid_T
echo "nslicesize=0"            >> nam_rebuild_grid_T
echo "deflate_level=1"         >> nam_rebuild_grid_T
echo "nc4_xchunk=50"           >> nam_rebuild_grid_T
echo "nc4_ychunk=50"           >> nam_rebuild_grid_T
echo "nc4_zchunk=1"            >> nam_rebuild_grid_T
echo "nc4_tchunk=1"            >> nam_rebuild_grid_T
echo "fchunksize=32000000"     >> nam_rebuild_grid_T
echo "/"                       >> nam_rebuild_grid_T

echo "&nam_rebuild"             > nam_rebuild_grid_U
echo "filebase='${MODEL}_grid_U'" >> nam_rebuild_grid_U
echo "ndomain=${NUM_DOM}"      >> nam_rebuild_grid_U
echo "l_maskout=.false "       >> nam_rebuild_grid_U
echo "nslicesize=0"            >> nam_rebuild_grid_U
echo "deflate_level=1"         >> nam_rebuild_grid_U
echo "nc4_xchunk=50"           >> nam_rebuild_grid_U
echo "nc4_ychunk=50"           >> nam_rebuild_grid_U
echo "nc4_zchunk=1"            >> nam_rebuild_grid_U
echo "nc4_tchunk=1"            >> nam_rebuild_grid_U
echo "fchunksize=32000000"     >> nam_rebuild_grid_U
echo "/"                       >> nam_rebuild_grid_U

echo "&nam_rebuild"             > nam_rebuild_grid_V
echo "filebase='${MODEL}_grid_V'" >> nam_rebuild_grid_V
echo "ndomain=${NUM_DOM}"      >> nam_rebuild_grid_V
echo "l_maskout=.false "       >> nam_rebuild_grid_V
echo "nslicesize=0"            >> nam_rebuild_grid_V
echo "deflate_level=1"         >> nam_rebuild_grid_V
echo "nc4_xchunk=50"           >> nam_rebuild_grid_V
echo "nc4_ychunk=50"           >> nam_rebuild_grid_V
echo "nc4_zchunk=1"            >> nam_rebuild_grid_V
echo "nc4_tchunk=1"            >> nam_rebuild_grid_V
echo "fchunksize=32000000"     >> nam_rebuild_grid_V
echo "/"                       >> nam_rebuild_grid_V

echo "&nam_rebuild"             > nam_rebuild_grid_W
echo "filebase='${MODEL}_grid_W'" >> nam_rebuild_grid_W
echo "ndomain=${NUM_DOM}"      >> nam_rebuild_grid_W
echo "l_maskout=.false "       >> nam_rebuild_grid_W
echo "nslicesize=0"            >> nam_rebuild_grid_W
echo "deflate_level=1"         >> nam_rebuild_grid_W
echo "nc4_xchunk=50"           >> nam_rebuild_grid_W
echo "nc4_ychunk=50"           >> nam_rebuild_grid_W
echo "nc4_zchunk=1"            >> nam_rebuild_grid_W
echo "nc4_tchunk=1"            >> nam_rebuild_grid_W
echo "fchunksize=32000000"     >> nam_rebuild_grid_W
echo "/"                       >> nam_rebuild_grid_W

echo "&nam_rebuild"             > nam_rebuild_grid_transport
echo "filebase='${MODEL}_grid_transport'" >> nam_rebuild_grid_transport
echo "ndomain=${NUM_DOM}"      >> nam_rebuild_grid_transport
echo "l_maskout=.false "       >> nam_rebuild_grid_transport
echo "nslicesize=0"            >> nam_rebuild_grid_transport
echo "deflate_level=1"         >> nam_rebuild_grid_transport
echo "nc4_xchunk=50"           >> nam_rebuild_grid_transport
echo "nc4_ychunk=50"           >> nam_rebuild_grid_transport
echo "nc4_zchunk=1"            >> nam_rebuild_grid_transport
echo "nc4_tchunk=1"            >> nam_rebuild_grid_transport
echo "fchunksize=32000000"     >> nam_rebuild_grid_transport
echo "/"                       >> nam_rebuild_grid_transport

echo "&nam_rebuild"             > nam_rebuild_restart
echo "filebase='${MODEL}_${RES_TIMESTAMP}_restart'" >> nam_rebuild_restart
echo "ndomain=${NUM_CPU}"      >> nam_rebuild_restart
echo "l_maskout=.false "       >> nam_rebuild_restart
echo "nslicesize=0"            >> nam_rebuild_restart
echo "deflate_level=1"         >> nam_rebuild_restart
echo "nc4_xchunk=50"           >> nam_rebuild_restart
echo "nc4_ychunk=50"           >> nam_rebuild_restart
echo "nc4_zchunk=1"            >> nam_rebuild_restart
echo "nc4_tchunk=1"            >> nam_rebuild_restart
echo "fchunksize=32000000"     >> nam_rebuild_restart
echo "/"                       >> nam_rebuild_restart

echo "&nam_rebuild"             > nam_rebuild_mesh_mask
echo "filebase=mesh_mask"      >> nam_rebuild_mesh_mask
echo "ndomain=${NUM_CPU}"      >> nam_rebuild_mesh_mask
echo "l_maskout=.false "       >> nam_rebuild_mesh_mask
echo "nslicesize=0"            >> nam_rebuild_mesh_mask
echo "deflate_level=1"         >> nam_rebuild_mesh_mask
echo "nc4_xchunk=50"           >> nam_rebuild_mesh_mask
echo "nc4_ychunk=50"           >> nam_rebuild_mesh_mask
echo "nc4_zchunk=1"            >> nam_rebuild_mesh_mask
echo "nc4_tchunk=1"            >> nam_rebuild_mesh_mask
echo "fchunksize=32000000"     >> nam_rebuild_mesh_mask
echo "/"                       >> nam_rebuild_mesh_mask


#Combine restart files and mesh_mask files
if $BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo nam_rebuild_restart ; then
    echo "Success restart files - removing separated files"
    rm -v ${MODEL}_${RES_TIMESTAMP}_restart_????.nc
else
    echo "Failure restart files - keeping separated files"
fi

#if $BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo nam_rebuild_mesh_mask ; then
#    echo "Success mesh_mask files - removing separated files"
#    rm -v mesh_mask_????.nc
#else
#    echo "Failure mesh_mask files - keeping separated files"
#fi

#if $BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo mesh_mask $NUM_CPU ; then
#    echo "Success mesh_mask files - removing separated files"
#    rm -v mesh_mask_????.nc
#else
#    echo "Failure mesh_mask files - keeping separated files"
#fi

#if $BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo nam_rebuild_grid_T ; then
#    echo "Success grid_T files - removing separated files"
#    rm -v ${MODEL}_grid_T_????.nc 
#else
#    echo "Failure grid_T files - keeping separated files"
#fi

#if $BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo nam_rebuild_grid_U ; then
#    echo "Success grid_U files - removing separated files"
#    rm -v ${MODEL}_grid_U_????.nc
#else
#    echo "Failure grid_U files - keeping separated files"
#fi

#if $BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo nam_rebuild_grid_V ; then
#    echo "Success grid_V files - removing separated files"
#    rm -v ${MODEL}_grid_V_????.nc
#else
#    echo "Failure grid_V files - keeping separated files"
#fi

#if $BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo nam_rebuild_grid_W ; then
#    echo "Success grid_W files - removing separated files"
#    rm -v ${MODEL}_grid_W_????.nc
#else
#    echo "Failure grid_W files - keeping separated files"
#fi

#if $BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo nam_rebuild_grid_transport ; then
#    echo "Success grid_W files - removing separated files"
#    rm -v ${MODEL}_grid_transport_????.nc
#else
#    echo "Failure grid_W files - keeping separated files"
#fi

# Export freq/start/end numbers fro grid output files (use T_0000.nc)
#export OUT_FREQ=$(echo $(ls -d ${MODEL}_*_grid_T_0000.nc) | awk -F _ '{print $2 }')
#export OUT_START=$(echo $(ls -d ${MODEL}_*_grid_T_0000.nc) | awk -F _ '{print $3 }')
#export OUT_END=$(echo $(ls -d ${MODEL}_*_grid_T_0000.nc) | awk -F _ '{print $4 }')

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2) Organise files
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



n=$(< resub_count)

export TIME=$(echo $(date +date_%d-%m-%y_time_%H-%M-%S))
export OUTDIR=OUTPUTS/${MODEL}_${TIME}_${label}_${n}
mkdir $OUTDIR

#Also make temporary restart file for resuming model if needed
cp -pv ${MODEL}_${RES_TIMESTAMP}_restart.nc ./RESTARTS/restart_tmp.nc

#Move restart files to output folder (combined and uncombined)
mv -v ${MODEL}_${RES_TIMESTAMP}_restart*.nc ./${OUTDIR}/

#Move mesh_mask files to output folder (combined and uncombined)
mv -v mesh_mask* ./${OUTDIR}/

#Move grid files to output folder
mv -v ${MODEL}*grid_*.nc ./${OUTDIR}/

#Copy run info into the output folder
cp -pv ./ocean.output ./${OUTDIR}/ocean.output.${RES_TIMESTAMP}
cp -pv ./solver.stat ./${OUTDIR}/solver.stat.${RES_TIMESTAMP}
cp -pv ./stdouterr ./${OUTDIR}/stdouterr.${RES_TIMESTAMP}
cp -pv ./namelist_cfg ./${OUTDIR}/namelist_cfg.${RES_TIMESTAMP}
cp -pv ./output.namelist.dyn ./${OUTDIR}/output.namelist.dyn.${RES_TIMESTAMP}
cp -pv ./nam_rebuild_grid_T  ./${OUTDIR}/nam_rebuild_grid_T.${RES_TIMESTAMP}
cp -pv ./nam_rebuild_grid_U  ./${OUTDIR}/nam_rebuild_grid_U.${RES_TIMESTAMP}
cp -pv ./nam_rebuild_grid_V  ./${OUTDIR}/nam_rebuild_grid_V.${RES_TIMESTAMP}
cp -pv ./nam_rebuild_grid_W  ./${OUTDIR}/nam_rebuild_grid_W.${RES_TIMESTAMP}
cp -pv ./nam_rebuild_restart ./${OUTDIR}/nam_rebuild_restart.${RES_TIMESTAMP}
cp -pv ./file_def_nemo-oce.xml ./${OUTDIR}/file_def_nemo-oce.xml.${RES_TIMESTAMP}
cp -pv ./field_def_nemo-oce.xml ./${OUTDIR}/field_def_nemo-oce.xml.${RES_TIMESTAMP}

#Move output foder to save location

if [ $n -le $Nresub_spup ]; then
   mv -v ./${OUTDIR} ${OUTSAVE}/SPINUP/
elif [ $n -gt $Nresub_spup ]; then
   mv -v ./${OUTDIR} ${OUTSAVE}/EXP/
fi

if [ $n -eq 1 ]; then
   mv -v ${OUTSAVE}/SPINUP/mesh_mask* ${OUTSAVE}
fi

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2) Updating the namelist_cfg file
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if [ $n -eq 1 ]; then
    #Save the original namelist
    cp -pv namelist_cfg namelist_cfg.original

    #Save the original file_def_nemo-oce.xml
    cp -pv file_def_nemo-oce.xml file_def_nemo-oce.xml.original
fi

if [ $n -lt $Nresub_spup ]; then
    export NEW_stock_NUM=$SPINUP_INT
elif [ $n -ge $Nresub_spup ]; then
    export NEW_stock_NUM=$EXP_INT
    #
    if [ $n -gt $(expr $Nresub_spup + $Nresub_exp) ]; then
       #Restore original namelist
       mv -v namelist_cfg.original namelist_cfg

       #Restore original file_def_nemo-oce.xml
       mv -v file_def_nemo-oce.xml.original file_def_nemo-oce.xml

      #Reset resub count to 1
      echo 1 > resub_count

      #Stop here
      exit
   fi
fi

echo "resubmission number " $n

# pull out strings from the namelist_cfg
#nn_itend>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
export OLD_itend_STR=$(grep -ri "nn_itend" namelist_cfg)
echo "Old itend String" $OLD_itend_STR

#remove 'nn_it000 =' part of string (variable name contains numbers)
export OLD_itend_STR_CLIP=$(echo ${OLD_itend_STR} | sed -e 's/.*=//g')
echo "Old itend String" $OLD_itend_STR

#extract just the number from the nn_it000
export OLD_itend_NUM=$(echo ${OLD_itend_STR_CLIP} | sed -e 's/[^0-9 ]//g' | awk '{print $NF}')
echo "Old itend Number" $OLD_itend_NUM
#Calculate new end time step
export NEW_itend_NUM=$((OLD_itend_NUM + $NEW_stock_NUM))
echo "New itend Number" $NEW_itend_NUM

#Create new replacement string
export NEW_itend_STR=$(echo ${OLD_itend_STR} | sed -r "s/${OLD_itend_NUM}/${NEW_itend_NUM}/g")
echo "NEW itend string " $NEW_itend_STR

#---------------------------------------

#nn_it000>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
export OLD_it000_STR=$(grep -ri "first time step" namelist_cfg)
echo "Old it000 String" $OLD_it000_STR

#remove 'nn_it000 =' part of string (variable name contains numbers)
export OLD_it000_STR_CLIP=$(echo ${OLD_it000_STR} | sed -e 's/.*=//g')
echo "Old it000 String" $OLD_it000_STR_CLIP

#extract just the number from the nn_it000
export OLD_it000_NUM=$(echo ${OLD_it000_STR_CLIP} | sed -e 's/[^0-9 ]//g' | awk '{print $NF}')
echo "Old it000 number " $OLD_it000_NUM

#Calculate new starting time step
export NEW_it000_NUM=$((OLD_itend_NUM + 1))
echo "New it000 number " $NEW_it000_NUM

#Create new replacement string
export NEW_it000_STR=$(echo ${OLD_it000_STR} | sed -r "s/${OLD_it000_NUM}/${NEW_it000_NUM}/g")
echo "NEW it000 string " $NEW_it000_STR
#---------------------------------------

##nn_stock>>>>>>>>>>>>>>>>>>>>>>>>>>>>
export OLD_stock_STR=$(grep -ri "frequency of creation of a restart file" namelist_cfg)
echo "Old stock String" $OLD_stock_STR

#remove comment part of string (contains numbers)
export OLD_stock_STR_CLIP=$(echo ${OLD_stock_STR} | sed -e 's/!.*//g')


#extract just the number from the nn_it000
export OLD_stock_NUM=$(echo ${OLD_stock_STR_CLIP} | sed -e 's/[^0-9 ]//g' | awk '{print $NF}')
echo "Old stock number: " $OLD_stock_NUM

echo "New stock number: " $NEW_stock_NUM

#Create new replacement string
export NEW_stock_STR=$(echo ${OLD_stock_STR} | sed -r "s/${OLD_stock_NUM}/${NEW_stock_NUM}/g")
echo "NEW stock string " $NEW_stock_STR

#-------------------------------------

#ln_start>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
export OLD_rstart_STR=$(grep -ri "start from rest (F) or from a restart file (T)" namelist_cfg)
echo $OLD_rstart_STR

#If false, set to true (essential for first restart)
export NEW_rstart_STR=$(echo ${OLD_rstart_STR} | sed -r 's/false/true/g')
echo $NEW_rstart_STR
#---------------------------------------

#ln_meshmask, set to false (essential for first restart)
export OLD_meshmask_STR=$(grep -ri "ln_meshmask" namelist_cfg)
export NEW_meshmask_STR=$(echo ${OLD_meshmask_STR} | sed -r 's/true/false/g')
echo $NEW_meshmask_STR

#Edit the namelist_cfg file
sed -i "s/${OLD_it000_STR}/   ${NEW_it000_STR}/g" namelist_cfg 
sed -i "s/${OLD_itend_STR}/   ${NEW_itend_STR}/g" namelist_cfg 
sed -i "s/${OLD_stock_STR}/   ${NEW_stock_STR}/g" namelist_cfg
sed -i "s/${OLD_rstart_STR}/   ${NEW_rstart_STR}/g" namelist_cfg
sed -i "s/${OLD_meshmask_STR}/   ${NEW_meshmask_STR}/g" namelist_cfg

#If going into the experimental window. Alter the output frequency in file_def_nemo-oce.xml
if [ $n -eq $Nresub_spup ]; then
    echo "Output frequency being changed from ${SPINUP_OUT_FREQ} to ${EXP_OUT_FREQ}"
    sed -i "s/output_freq=\"${SPINUP_OUT_FREQ}\"/output_freq=\"${EXP_OUT_FREQ}\"/g" file_def_nemo-oce.xml
fi

echo $(($n + 1))
echo $(($n + 1)) > resub_count

cd $SUBMDIR

qsub $SUBMSCRIPT

exit
