
sampleinfo="sample_info_example_ESrun.txt"
norm="TMM"
corr="RUVs"
enrich="stouffer"
dea="DESeq2shrink"
deaeq="0+W_1+group"


scriptdir="../scripts"
listdir="../gene_lists"
geneinfo="$listdir"/rareseq_targeted_coding_genes.txt
negctrlgenes="$listdir"/platelet_genes.txt
sig="$listdir"/SupplementaryTable10_gene_signatures.txt

indir="../../raw_RSEM_files"
outdir="output"


Rscript "$scriptdir"/importExpression.R $indir $sampleinfo
Rscript "$scriptdir"/normalizeExpression.R "$outdir"/tximport.rds $geneinfo $sampleinfo $norm $corr
Rscript "$scriptdir"/calcZscores.R "$outdir"/"$norm"_"$corr"/"$norm"_"$corr".rds
Rscript "$scriptdir"/calcEnrichment.R "$outdir"/"$norm"_"$corr"/"$norm"_"$corr"_withZscores.rds $enrich $sig test
Rscript "$scriptdir"/differentialExpression.R "$outdir"/"$norm"_"$corr"/"$norm"_"$corr".rds $dea $deaeq



