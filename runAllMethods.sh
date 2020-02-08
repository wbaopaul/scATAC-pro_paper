#!/bin/bash

methods=( seurat seurat_correct cisTopic scABC SCRAT chromVAR LSI ) 
#dtypes=( HSC_GSE96769 )
dtypes=( resample10k_hsc_clean resample10k_hsc_noisy_q2 resample10k_hsc_noisy_q4 )
#dtypes=( resample10k_hsc resample3k_hsc )
#dtypes=( bonemarrow_clean bonemarrow_noisy_p2 bonemarrow_noisy_p4 )
#dtypes=( erythropoiesis_clean erythropoiesis_noisy_p2 erythropoiesis_noisy_p4 )


clength=${#methods[@]}
rlength=${#dtypes[@]}

job0='job_run'
echo "#!/bin/bash" > $job0
echo "#$ -cwd" >> $job0
echo "#$ -j y" >> $job0
echo "#$ -pe smp 2"  >> $job0
echo "#$ -l h_vmem=16G"  >> $job0
echo "#$ -l mem_free=16G"  >> $job0

for (( i=0; i<${clength}; i++ ));
do
    method0=${methods[$i]}
    for (( j=0; j<${rlength}; j++ ));
        do
        dtype0=${dtypes[$j]}
        jobk=${job0}_${method0}_${dtype0} 
        cp $job0 $jobk
        if [ "$method0" == "cisTopic" ]; then
          echo "#$ -pe smp 8"  >> $jobk
        fi
        #echo "R --vanilla --args $dtype0 $method0 < runMethodDefaultSetting_diffComposition.R ">> $jobk 
        echo "R --vanilla --args $dtype0 $method0 < runMethodDefaultSetting.R ">> $jobk 
        qsub $jobk
        #nohup bash $jobk &
    done
done
