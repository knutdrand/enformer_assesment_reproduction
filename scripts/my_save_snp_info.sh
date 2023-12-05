# From vcf, save variant information by gene for each SNP within 100k of the TSS

#base_path=/data/mostafavilab/bng/rosmapAD/projects/insilicoMutagenesis/extractSequence/results/
#vcf_path=/data/mostafavilab/bng/rosmapAD/data/wholeGenomeSeq/
#save_dir=/data/aspiro17/enformer_res/variant_info_100k/
#all_rel_genes_path=/data/aspiro17/enformer_res/gene_lists/ # save information for these genes

base_path=data/
vcf_path=data/
save_dir=data/debug/
all_rel_genes_path=data/ # save information for these genes
mkdir -p ${save_dir}
# readarray -t all_rel_genes < ${all_rel_genes_path}to_process.txt
gene_path='enformer_assesment_reproduction/tmp.csv'
ensg=(`cut -f1 ${gene_path}`)
chr=(`cut -f2 ${gene_path}`)
winStart=(`cut -f3 ${gene_path}`)
winEnd=(`cut -f4 ${gene_path}`)
for i in ${!ensg[@]}
do

start=`date +%s`

echo ${ensg[$i]}
echo ${all_rel_genes}
#if [[ " ${all_rel_genes[*]} " =~ " ${ensg[$i]} " ]]; then
    echo "yes"
    bcftools view -v snps -r ${chr[$i]}:${winStart[$i]}-${winEnd[$i]} ../bionumpy-example-data/big_phased.vcf.gz | bcftools query -f '%CHROM,%REF,%POS[,%TGT]\n' > ${save_dir}${ensg[$i]}.csv

    end=`date +%s`
    runtime=$((end-start))
    echo ${runtime}
#fi
done

