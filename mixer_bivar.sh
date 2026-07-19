sigfile="test/amygdala_phenotupes"

for phe in `cat $sigfile`;do
for i in $(seq 1 20);do
cmd="python3.8 mixer/mixer-master/precimed/mixer.py fit2 --trait1-file test/cFDR/traitfolder/"$phe".csv_noMHC.csv --trait2-file test/cFDR/traitfolder/MDD2025_noMHC.csv --trait1-params-file test/mixer/result/"$phe".fit.rep"$i".json --trait2-params-file test/mixer/result/MDD2025.fit.rep"$i".json --out test/mixer/bivar_res/MDD2025_vs_"$phe".fit.rep"$i" --extract mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep"$i".snps --bim-file mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --ld-file mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --lib mixer/mixer-master/src/build/lib/libbgmg.so"
echo $cmd > "test/mixer/script/MDD2025_vs_"$phe"_fit2_"$i".sh"
qsub -q cluster5 -e "test/mixer/script/MDD2025_vs_"$phe"_fit2_"$i".err" -o "test/mixer/script/MDD2025_vs_"$phe"_fit2_"$i".out" "test/mixer/script/MDD2025_vs_"$phe"_fit2_"$i".sh"
done
done


sigfile="test/amygdala_phenotupes"

for phe in `cat $sigfile`;do
for i in $(seq 1 20);do
cmd="python3.8 mixer/mixer-master/precimed/mixer.py test2 --trait1-file test/cFDR/traitfolder/"$phe".csv_noMHC.csv --trait2-file test/cFDR/traitfolder/MDD2025_noMHC.csv --load-params-file test/mixer/bivar_res/MDD2025_vs_"$phe".fit.rep"$i".json --out test/mixer/bivar_res/MDD2025_vs_"$phe".test2.rep"$i" --bim-file mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --ld-file mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --lib mixer/mixer-master/src/build/lib/libbgmg.so"
echo $cmd > "test/mixer/script/MDD2025_vs_"$phe"_test2_"$i".sh"
qsub -q cluster5 -e "test/mixer/script/MDD2025_vs_"$phe"_test2_"$i".err" -o "test/mixer/script/MDD2025_vs_"$phe"_test2_"$i".out" "test/mixer/script/MDD2025_vs_"$phe"_test2_"$i".sh"
done
done


for phe in `cat test/amygdala_phenotupes`;do
python3.8 mixer/mixer-master/precimed/mixer_figures.py combine --json test/mixer/bivar_res/MDD2025_vs_"$phe".fit.rep@.json --out test/mixer/bivar_res/MDD2025_vs_"$phe".fit

python3.8 mixer/mixer-master/precimed/mixer_figures.py combine --json test/mixer/bivar_res/MDD2025_vs_"$phe".test2.rep@.json --out test/mixer/bivar_res/MDD2025_vs_"$phe".test2

python3.8 mixer/mixer-master/precimed/mixer_figures.py two --json-fit test/mixer/bivar_res/MDD2025_vs_"$phe".fit.json --json-test test/mixer/bivar_res/MDD2025_vs_"$phe".test2.json --out test/mixer/bivar_res/MDD2025_vs_"$phe"_bivarres --statistic mean std
done