# Run HiC pro analysis
HiC-Pro -i 9311/data/ -o 9311/results/ -c config_9311.txt
HiC-Pro -i Nip/data/ -o Nip/results/ -c config_Nip.txt

# matrix format convert
binSize=10000
dir=`pwd`
mkdir $binSize
for i in HS_rep1  HS_rep2  Normal_rep1  Normal_rep2 ; do
    $python $bin/sparseToDense.py -b ${dir}/${i}/raw/$binSize/${i}_${binSize}_abs.bed ${dir}/${i}/iced/${binSize}/${i}_${binSize}_iced.matrix &
done


# TAD analysis with binSize=10000

binSize=10000
dir=`pwd`

for i in HS_rep1 HS_rep2 Normal_rep1 Normal_rep2 ; do
 cd ${dir}/${i}/iced/${binSize}
 hicConvertFormat -m ${i}_${binSize}_iced.matrix --bedFileHicpro ${dir}/${i}/raw/${binSize}/${i}_${binSize}_abs.bed  --inputFormat hicpro --outputFormat h5 -o ${i}_${binSize}_iced.h5 &
 hicFindTADs -m ${i}_${binSize}_iced.h5 --outPrefix ${i}_${binSize}_iced --minDepth 15000 --maxDepth 150000 --step 5000 --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr  -p 50 &
done

# compartment analysis 
dir=`pwd`
cd $dir
for i in Nip_HS_rep{1,2} Nip_Normal_rep{1,2} ; do
    runHiCpca.pl auto ${i}/ -res 50000 -window 100000 -genome Nip -cpu 10 -rpath /usr/bin/R  
done