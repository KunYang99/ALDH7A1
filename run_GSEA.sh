pathway=$1
input=$2
output=$3

mkdir -p $output

java -cp /Users/kunyang/Applications/gsea-3.0.jar -Xmx5000m xtools.gsea.GseaPreranked -gmx $pathway -norm meandiv -nperm 1000 -rnk $input -scoring_scheme classic -rpt_label myAnal -create_svgs false -make_sets true -plot_top_x 10 -rnd_seed timestamp -set_max 500 -set_min 10 -zip_report false -out $output -gui false
