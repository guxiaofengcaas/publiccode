DIR=`pwd`
REF="/9311/9311_star"
GTF="/9311_genome_gene_final.gff3"
cd $DIR/rawdata/
for NA in `ls *.gz |awk 'gsub("_good.+","")'|sort -u ` ; do
    echo $NA
    FQ=`ls $DIR/rawdata/${NA}*.gz`
    echo $FQ
    OUT=$NA\_STAR
    cd $DIR/matching/
    if [ ! -d $OUT ] ; then
        mkdir $OUT
    fi  
    cd $DIR/matching/$OUT/

    STAR --genomeDir $REF --genomeLoad LoadAndKeep --runThreadN 55 --readFilesIn $FQ --outFileNamePrefix $DIR/matching/$OUT/$NA --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --limitBAMsortRAM 20000000000
done
STAR --genomeDir $REF --genomeLoad Remove


featureCounts -T 30  -a $te -t CDS  -g Parent -o TE.featureCount.res matching/*STAR/*bam

