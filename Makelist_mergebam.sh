t="cgas"
for i in {1..4}
do
        grep $t"_tp"$i 8LGAxjADB9Lgukt2w0Fp.list > "$t"_tp"$i".list
done



#source /mnt/common/epfl/etc/bbcf_bashrc
#hts_minilims.py -m mapseq -e $jobkey -l > $jobkey".list"


##############mergebam

############ hay que generar un archivo .list por cada grupo:
#t="wt"
#for i in {1..4}
#do
#       grep $t"_tp"$i 4PPrqSjEWvQFlSde6U5u.list > "$t"_tp"$i".list
#done




#samtools merge out.bam in1.bam in2.bam in3.bam
for x in `ls *list`; do
echo $x;
awk '{if($3 ~ /bam$/){n=split($3,c,"_");for(j=1; j<=n;j++){if(c[j] ~ /rep[0-9]+/){rep=c[j];}};a[rep]=a[rep]" "$4}}END{pref=FILENAME; gsub(".list","",pref);for(i in a){outfile=pref"_"i"_merged.bam";print "*****"outfile;cmd="samtools merge "outfile" "a[i]" &"; print cmd; system(cmd)}; }' $x

done
