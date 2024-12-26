
sampleinfo="sample_info_example_ENrun.txt"
norm="TMM"
corr="RUVs"


scriptdir="../scripts"
listdir="../gene_lists"
geneinfo="$listdir"/rareseq_targeted_coding_genes.txt

indir="../../raw_RSEM_files"
outdir="output"


Rscript "$scriptdir"/importExpression.R $indir $sampleinfo
Rscript "$scriptdir"/normalizeExpression.R "$outdir"/tximport.rds $geneinfo $sampleinfo $norm $corr


calls="../../RARESeq-targeted-caller/example_output/variant_call_table.txt"
Rscript "$scriptdir"/setupEN.R "$outdir"/"$norm"_"$corr"/"$norm"_"$corr".rds $calls
Rscript "$scriptdir"/trainEN.R "$outdir"/"$norm"_"$corr"/caret/"$norm"_"$corr"_forCaret.rds
Rscript "$scriptdir"/validateEN.R "$outdir"/"$norm"_"$corr"/caret/"$norm"_"$corr"_forCaret.rds "$outdir"/"$norm"_"$corr"/caret/model.rds
