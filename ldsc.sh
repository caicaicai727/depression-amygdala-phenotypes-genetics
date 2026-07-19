sigfile="test/amygdala_phenotypes"


# Munge amygdala GWAS summary statistics
for i in `cat $sigfile`;do
python3.8 ldsc/munge_sumstats.py \
--sumstats test/ldsc/IDP/"$i"_noMHC.temp \
--N 33224 \
--snp rsid \
--a1 a1 \
--a2 a2 \
--p P \
--out test/ldsc/IDP/"$i" \
--merge-alleles ldsc/eur_w_ld_chr/w_hm3.snplist
done


# Munge MDD GWAS summary statistics
python3.8 ldsc/munge_sumstats.py \
--sumstats test/ldsc/mdd2025/pgcmdd2025_no23andMenoUKBB_eur_v3492411_noMHC.ma \
--N-col N \
--snp SNP \
--a1 A1 \
--a2 A2 \
--p P \
--out test/ldsc/mdd2025/MDD2025_noMHC \
--merge-alleles ldsc/eur_w_ld_chr/w_hm3.snplist


# Genetic correlation analysis
for i in `cat $sigfile`;do
python3.8 ldsc/ldsc.py \
--rg test/ldsc/mdd2025/MDD2025_noMHC.sumstats.gz,test/ldsc/IDP/"$i".sumstats.gz \
--ref-ld-chr ldsc/eur_w_ld_chr/ \
--w-ld-chr ldsc/eur_w_ld_chr/ \
--out test/ldsc/IDP/MDD2025_"$i"
done