
module load StdEnv/2020 gcc/9.3.0 mugqic/R_Bioconductor/4.1.0_3.13 r/4.1.2

path="/lustre06/project/6004736/alvann/from_narval/221128_K36E/ChIP/mm/out/01A_genpipes_chipseq/HPC_10W"
path_dm="/lustre06/project/6004736/alvann/from_narval/221128_K36E/ChIP/mm/out/02A_genpipes_chipseq_dm6/HPC_10W"
out="/lustre06/project/6004736/alvann/from_narval/221128_K36E/ChIP/mm/out/02B_rx_values/HPC_10W"

batch=("batch1")

mkdir -p ${out}

for b in "${batch[@]}"
do
   Rscript -e "rmarkdown::render(
            '/lustre06/project/6004736/alvann/from_narval/221128_K36E/ChIP/mm/scripts/src/calculate_rx_norm.Rmd',
            params = list(path_mm = '${path}/${b}/metrics/SampleMetrics.tsv',
                          path_dm = '${path_dm}/${b}/metrics/SampleMetrics.tsv',
                          path_out = '${out}/',
                          prefix = '${b}'),
            output_file = '${out}/rx_norm.${b}.html')"
done
