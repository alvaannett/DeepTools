import pandas as pd

# checks final output files
rule all:
    input:
        expand(config["OUT_PATH"] + "/03_plot_qc_corr/{norm_meth}.scale_{scal_fact}.bin_{bin_size}.correlation.pdf",
               norm_meth = config["NORM_METHOD"],
               scal_fact = config["SCALE_FACTOR"],
               bin_size = config["BIN_SIZE"]),
        expand(config["OUT_PATH"] + "/03_plot_qc_pca/{norm_meth}.scale_{scal_fact}.bin_{bin_size}.pca_{pc}.pdf",
	       norm_meth = config["NORM_METHOD"],
	       scal_fact = config["SCALE_FACTOR"],
               pc = config["PC"],
               bin_size = config["BIN_SIZE"]),
        expand(config["OUT_PATH"] + "/03_plot_qc_corr/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.correlation.pdf",
               norm_meth = config["NORM_METHOD"],
               scal_fact = config["SCALE_FACTOR"],
               bed_region = config["SUMMARY_REGIONS"]),
        expand(config["OUT_PATH"] + "/03_plot_qc_pca/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.pca_{pc}.pdf",
               pc = config["PC"],
               norm_meth = config["NORM_METHOD"],
               scal_fact = config["SCALE_FACTOR"],
               bed_region = config["SUMMARY_REGIONS"]),
        expand(config["OUT_PATH"] + "/05_plot_heatmap/{norm_meth}.scale_{scal_fact}/{matrix}.pdf",
               norm_meth = config["NORM_METHOD"],
               scal_fact = config["SCALE_FACTOR"],
               matrix = config["PLOT_HEATMAP_PARAM"]),
        expand(config["OUT_PATH"] + "/06_plot_profile/{norm_meth}.scale_{scal_fact}/{matrix}.pdf",
               norm_meth = config["NORM_METHOD"],
               scal_fact = config["SCALE_FACTOR"],
               matrix = config["PLOT_PROFILE_PARAM"]),
        expand(config["OUT_PATH"] + "/07_z_score/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.{condition}.z_score.tab",
               norm_meth = config["NORM_METHOD"],
               scal_fact = config["SCALE_FACTOR"],
               bed_region = config["SUMMARY_REGIONS"],
               condition = config["Z_SCORE_PARAM"]["condition"]),
        expand(config["OUT_PATH"] + "/08_log2fc/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.{condition}.log2fc.tab",
               norm_meth = config["NORM_METHOD"],
               scal_fact = config["SCALE_FACTOR"],
	       bed_region = config["SUMMARY_REGIONS"],
	       condition = config["LOG2FC_PARAM"]["condition"])


# function to get combinations of samples and names that excist
def get_sample_mark_combo(paths, samples, marks):
  sample_mark = []
  for s in samples:
      for m in marks:
         for p in paths:
           bam_path = p + "/alignment/" + s + "/" + m + "/" + s + "." + m + ".sorted.dup.filtered.bam"
	   if os.path.exists(bam_path):
              sample_mark.append(s + "." + m)
  return(sample_mark)

# function to return bam path from sample and mark name
def get_bam_path_from_sample_mark(wildcards):
  sm = wildcards.sample_mark.split('.')
  for p in config["GENPIPES"]:
     path = p + "/alignment/" + sm[0] + "/" + sm[1] + "/" + sm[0] + "." + sm[1] + ".sorted.dup.filtered.bam"
     if os.path.exists(path):
        return(path)

# function to get sample bw for compute matrix and heatmaps
def get_bw_samples_for_compute_matrix(wildcards, norm_meth, scal_fact):
  bw = [config["OUT_PATH"] + "/01_normalize_tracks/" + norm_meth + ".scale_" + scal_fact + "/" + s + "." + norm_meth + ".scale_" + scal_fact + ".bw" for s in config["COMPUT_MATRIX_PARAM"][wildcards.matrix]["samples"]]
  return bw 

SAMPLE_MARK = get_sample_mark_combo(config["GENPIPES"], config["SAMPLES"], config["MARKS"])

if(os.path.exists(config['SCALE_FACTOR_FILE_RX'])):
    rx = pd.read_csv(config['SCALE_FACTOR_FILE_RX'])
    rx.columns = rx.columns.str.replace(' ', '_')
else:
    rx = pd.DataFrame(columns=['Sample', 'Mark_Name', 'rx'])


def get_rx_factor(sample_mark, rx):
    if not(rx.empty):
        sample_mark = sample_mark.split('.')
        return(rx.loc[(rx['Sample'] == sample_mark[0]) & (rx['Mark_Name'] == sample_mark[1]), 'rx'].values[0])
    else:
        return(1)


rule normalize_tracks:
    input:
        bam = get_bam_path_from_sample_mark
    output:
        bw = config["OUT_PATH"] + "/01_normalize_tracks/{norm_meth}.scale_{scal_fact}/{sample_mark}.{norm_meth}.scale_{scal_fact}.bw"
    params:
        rx_factor = lambda wildcards: get_rx_factor(wildcards.sample_mark, rx)
    threads: config["CPU"]
    resources:
        mem_mb=config["CPU"]*5000
    shell:
        '''
        module load mugqic/deepTools/3.3.1
        module load mugqic/python/2.7.14

        SCALE="{wildcards.scal_fact}"
        MATCH="rx"

        if [ "$SCALE" == "$MATCH" ]; then
        
        bamCoverage \
         --bam {input.bam} \
         --outFileName {output.bw} \
         --outFileFormat bigwig \
         --binSize 10 \
         --normalizeUsing {wildcards.norm_meth} \
         --scaleFactor {params.rx_factor} \
         --numberOfProcessors {config[CPU]} \
         --minMappingQuality 30 \
         --extendReads \
         --verbose

         else
         bamCoverage \
          --bam {input.bam} \
          --outFileName {output.bw} \
          --outFileFormat bigwig \
          --binSize 10 \
          --normalizeUsing {wildcards.norm_meth} \
          --scaleFactor {wildcards.scal_fact} \
          --numberOfProcessors {config[CPU]} \
          --minMappingQuality 30 \
          --extendReads \
          --verbose

          fi
        '''

# takes big wig tracks and computes average scores in genomic bins
rule bigwig_summary_bin:
   input:
      bw = expand(config["OUT_PATH"] + "/01_normalize_tracks/{{norm_meth}}.scale_{{scal_fact}}/{sample_mark}.{{norm_meth}}.scale_{{scal_fact}}.bw",
                  sample_mark = SAMPLE_MARK)
   output:
      out = config["OUT_PATH"] + "/02_bigwig_summary_bin/avrage_bin.{norm_meth}.scale_{scal_fact}.bin_{bin_size}.npz",
      raw_counts = config["OUT_PATH"] + "/02_bigwig_summary_bin/avrage_bin.{norm_meth}.scale_{scal_fact}.bin_{bin_size}.tab"
   threads: config["CPU"]
   resources:
      mem_mb=config["CPU"]*5000
   shell:
       '''
       module load mugqic/deepTools/3.3.1
       module load mugqic/python/2.7.14

       multiBigwigSummary bins \
        --bwfiles {input.bw} \
        --outFileName {output.out} \
        --outRawCounts {output.raw_counts}  \
        --binSize {wildcards.bin_size} \
        --numberOfProcessors {config[CPU]} \
        --smartLabels \
        --verbose
       '''

# takes big wig tracks and computes average scores in genomic regions
rule bigwig_summary_bed:
    input:
       bw = expand(config["OUT_PATH"] + "/01_normalize_tracks/{{norm_meth}}.scale_{{scal_fact}}/{sample_mark}.{{norm_meth}}.scale_{{scal_fact}}.bw",
                   sample_mark = SAMPLE_MARK)
    output:
       out = config["OUT_PATH"] + "/02_bigwig_summary_bed/avrage_bed.{norm_meth}.scale_{scal_fact}.bed_{bed_region}.npz",
       raw_counts = config["OUT_PATH"] + "/02_bigwig_summary_bed/avrage_bed.{norm_meth}.scale_{scal_fact}.bed_{bed_region}.tab"
    params:
       region_bed = lambda wildcards: config["SUMMARY_REGIONS"][wildcards.bed_region]
    threads: config["CPU"]
    resources:
        mem_mb=config["CPU"]*5000
    shell:
        '''
        module load mugqic/deepTools/3.3.1
        module load mugqic/python/2.7.14

        multiBigwigSummary BED-file \
         --bwfiles {input.bw} \
         --outFileName {output.out} \
         --outRawCounts {output.raw_counts} \
         --BED {params.region_bed} \
         --numberOfProcessors {config[CPU]} \
         --smartLabels \
         --verbose
        '''

# takes big wig tracks and computes average scores in genomic bins
rule plot_correlation_bin:
    input:
        bigwig_summary_bin = config["OUT_PATH"] + "/02_bigwig_summary_bin/avrage_bin.{norm_meth}.scale_{scal_fact}.bin_{bin_size}.npz"
    output:
        out_pdf = config["OUT_PATH"] + "/03_plot_qc_corr/{norm_meth}.scale_{scal_fact}.bin_{bin_size}.correlation.pdf",
        out_matrix = config["OUT_PATH"] + "/03_plot_qc_corr/{norm_meth}.scale_{scal_fact}.bin_{bin_size}.correlation.tab"
    threads: 2
    resources:
        mem_mb=2000
    shell:
        '''
        module load mugqic/deepTools/3.3.1
        module load mugqic/python/2.7.14

        plotCorrelation \
         --corData {input.bigwig_summary_bin} \
         --plotFile {output.out_pdf} \
         --outFileCorMatrix {output.out_matrix} \
         --plotFileFormat pdf \
         --corMethod spearman \
         --whatToPlot heatmap \
         --colorMap coolwarm \
         --plotNumbers \
         --skipZeros
        '''

rule plot_pca_bin:
    input:
        bigwig_summary_bin = config["OUT_PATH"] + "/02_bigwig_summary_bin/avrage_bin.{norm_meth}.scale_{scal_fact}.bin_{bin_size}.npz"
    output:
        out_pdf = config["OUT_PATH"] + "/03_plot_qc_pca/{norm_meth}.scale_{scal_fact}.bin_{bin_size}.pca_{pc}.pdf",
        out_matrix = config["OUT_PATH"] + "/03_plot_qc_pca/{norm_meth}.scale_{scal_fact}.bin_{bin_size}.pca_{pc}.tab"
    params:
        pc = lambda wildcards: config["PC"][wildcards.pc]
    threads: 2
    resources:
        mem_mb=2000
    shell:
        '''
        module load mugqic/deepTools/3.3.1
        module load mugqic/python/2.7.14

        plotPCA \
         --corData {input.bigwig_summary_bin} \
         --plotFile {output.out_pdf} \
         --outFileNameData {output.out_matrix} \
         --plotFileFormat pdf \
         --PCs {params.pc} \
         --rowCenter
        '''

rule plot_correlation_bed:
    input:
        bigwig_summary_bed = config["OUT_PATH"] + "/02_bigwig_summary_bed/avrage_bed.{norm_meth}.scale_{scal_fact}.bed_{bed_region}.npz"
    output:
        out_pdf = config["OUT_PATH"] + "/03_plot_qc_corr/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.correlation.pdf",
        out_matrix = config["OUT_PATH"] + "/03_plot_qc_corr/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.correlation.tab"
    threads: 2
    resources:
        mem_mb=2000
    shell:
        '''
        module load mugqic/deepTools/3.3.1
        module load mugqic/python/2.7.14

        plotCorrelation \
         --corData {input.bigwig_summary_bed} \
         --plotFile {output.out_pdf} \
         --outFileCorMatrix {output.out_matrix} \
         --plotFileFormat pdf \
         --corMethod spearman \
         --whatToPlot heatmap \
         --colorMap coolwarm \
         --plotNumbers \
         --skipZeros
        '''

rule plot_pca_bed:
    input:
        bigwig_summary_bed = config["OUT_PATH"] + "/02_bigwig_summary_bed/avrage_bed.{norm_meth}.scale_{scal_fact}.bed_{bed_region}.npz"
    output:
        out_pdf = config["OUT_PATH"] + "/03_plot_qc_pca/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.pca_{pc}.pdf",
        out_matrix = config["OUT_PATH"] + "/03_plot_qc_pca/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.pca_{pc}.tab"
    params:
        pc = lambda wildcards: config["PC"][wildcards.pc]
    threads: 2
    resources:
        mem_mb=2000
    shell:
        '''
        module load mugqic/deepTools/3.3.1
        module load mugqic/python/2.7.14

        plotPCA \
         --corData {input.bigwig_summary_bed} \
         --plotFile {output.out_pdf} \
         --outFileNameData {output.out_matrix} \
         --plotFileFormat pdf \
         --PCs {params.pc} \
         --rowCenter
        '''

rule compute_matrix:
    input:
        bw = lambda wildcards: get_bw_samples_for_compute_matrix(wildcards, norm_meth = wildcards.norm_meth, scal_fact = wildcards.scal_fact)
    output:
        matrix = config["OUT_PATH"] + "/04_compute_matrix/{norm_meth}.scale_{scal_fact}/{matrix}.gzip"
    params:
        type = lambda wildcards: config["COMPUT_MATRIX_PARAM"][wildcards.matrix]["type"],
        region_bed = lambda wildcards: config["COMPUT_MATRIX_PARAM"][wildcards.matrix]["region_bed"],
        length = lambda wildcards: config["COMPUT_MATRIX_PARAM"][wildcards.matrix]["length"],
        up = lambda wildcards: config["COMPUT_MATRIX_PARAM"][wildcards.matrix]["up"],
        down = lambda wildcards: config["COMPUT_MATRIX_PARAM"][wildcards.matrix]["down"],
        ref_point = lambda wildcards: config["COMPUT_MATRIX_PARAM"][wildcards.matrix]["ref_point"]
    resources:
        mem_mb = config["CPU_COMPUTE_MATRIX"]*5000
    threads: config["CPU_COMPUTE_MATRIX"]
    shell:
        '''
        module load mugqic/deepTools/3.3.1
        module load mugqic/python/2.7.14

        TYPE="{params.type}"
        MATCH1="scale"
        MATCH2="ref_point"

        if [ "$TYPE" == "$MATCH1" ]; then

        computeMatrix scale-regions \
         --scoreFileName {input.bw} \
         --regionsFileName {params.region_bed} \
         --outFileName {output.matrix} \
         --regionBodyLength {params.length} \
         --upstream {params.up}  \
         --downstream {params.down} \
         --numberOfProcessors {config[CPU_COMPUTE_MATRIX]} \
         --skipZeros

        elif [ "$TYPE" == "$MATCH2" ]; then

        computeMatrix reference-point \
         --scoreFileName {input.bw} \
         --regionsFileName {params.region_bed} \
         --outFileName {output.matrix} \
         --referencePoint "{params.ref_point}" \
         --beforeRegionStartLength {params.up}  \
         --afterRegionStartLength {params.down} \
         --numberOfProcessors {config[CPU_COMPUTE_MATRIX]} \
         --skipZeros \
         --missingDataAsZero 
        fi
        '''


rule plot_heatmap:
    input:
        matrix = config["OUT_PATH"] + "/04_compute_matrix/{norm_meth}.scale_{scal_fact}/{matrix}.gzip"
    output:
        heatmap = config["OUT_PATH"] + "/05_plot_heatmap/{norm_meth}.scale_{scal_fact}/{matrix}.pdf"
    params:
        color = lambda wildcards: config["PLOT_HEATMAP_PARAM"][wildcards.matrix]["color"],
        lables = lambda wildcards: config["PLOT_HEATMAP_PARAM"][wildcards.matrix]["lables"],
        startLabel = lambda wildcards: config["PLOT_HEATMAP_PARAM"][wildcards.matrix]["startLabel"],
        endLabel = lambda wildcards: config["PLOT_HEATMAP_PARAM"][wildcards.matrix]["startLabel"],
	title = lambda wildcards: config["PLOT_HEATMAP_PARAM"][wildcards.matrix]["title"],
        plotType = lambda wildcards: config["PLOT_HEATMAP_PARAM"][wildcards.matrix]["plotType"],
	k = lambda wildcards: config["PLOT_HEATMAP_PARAM"][wildcards.matrix]["k"],
        height = lambda wildcards: config["PLOT_HEATMAP_PARAM"][wildcards.matrix]["height"]
    threads: config["CPU"]/2
    resources:
        mem_mb=2000*5
    shell:
        '''
        module load mugqic/deepTools/3.3.1
        module load mugqic/python/2.7.14
        
        plotHeatmap \
         --matrixFile {input.matrix} \
         --outFileName {output.heatmap} \
         --plotTitle "{params.title}" \
         --samplesLabel {params.lables} \
         --colorList "{params.color}" \
         --colorNumber 100 \
         --startLabel "{params.startLabel}" \
         --endLabel "{params.endLabel}" \
         --plotType "{params.plotType}" \
         --heatmapHeight {params.height} \
         --heatmapWidth 3 \
         --kmeans {params.k}
         '''
        
rule plot_profile:
    input:
        matrix = config["OUT_PATH"] + "/04_compute_matrix/{norm_meth}.scale_{scal_fact}/{matrix}.gzip"
    output:
       profile = config["OUT_PATH"] + "/06_plot_profile/{norm_meth}.scale_{scal_fact}/{matrix}.pdf"
    params:
       lables = lambda wildcards: config["PLOT_PROFILE_PARAM"][wildcards.matrix]["lables"],
       startLabel = lambda wildcards: config["PLOT_PROFILE_PARAM"][wildcards.matrix]["startLabel"],
       endLabel = lambda wildcards: config["PLOT_PROFILE_PARAM"][wildcards.matrix]["startLabel"],
       title = lambda wildcards: config["PLOT_PROFILE_PARAM"][wildcards.matrix]["title"],
       plotType = lambda wildcards: config["PLOT_PROFILE_PARAM"][wildcards.matrix]["plotType"],
       colors = lambda wildcards: config["PLOT_PROFILE_PARAM"][wildcards.matrix]["colors"]
    threads: config["CPU"]
    resources:
       mem_mb=config["CPU"]*15000
    shell:
       '''
       module load mugqic/deepTools/3.3.1
       module load mugqic/python/2.7.14

       plotProfile \
        --matrixFile {input.matrix} \
        --outFileName {output.profile} \
        --samplesLabel {params.lables} \
        --startLabel "{params.startLabel}" \
        --endLabel "{params.endLabel}" \
	    --plotTitle "{params.title}" \
        --plotType "{params.plotType}" \
        --colors {params.colors} \
        --averageType mean \
        --perGroup \
	    --plotWidth 10 \
	    --plotHeight 10 
      '''

rule calculate_z_score_bed:
  input:
    counts = config["OUT_PATH"] + "/02_bigwig_summary_bed/avrage_bed.{norm_meth}.scale_{scal_fact}.bed_{bed_region}.tab"
  output:
    z_score = config["OUT_PATH"] + "/07_z_score/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.{condition}.z_score.tab"
  params:
    control = config["Z_SCORE_PARAM"]["control"],
    condition = lambda wildcards: config["Z_SCORE_PARAM"]["condition"][wildcards.condition]
  shell:
    '''
    module load StdEnv/2020
    module load gcc/9.3.0
    module load mugqic/R_Bioconductor/4.1.0_3.13
    module load r/4.1.2

    Rscript --vanilla /lustre06/project/6004736/alvann/from_narval/PIPELINES/DeepTools/scripts/calculate_z-score.R \
                      {input.counts} \
                      {output.z_score} \
                      "{params.control}"  \
                      "{params.condition}" \
                      "{config[MARKS]}"
    '''
 
rule calculate_log2fc_bed:
  input:
    counts = config["OUT_PATH"] + "/02_bigwig_summary_bed/avrage_bed.{norm_meth}.scale_{scal_fact}.bed_{bed_region}.tab"
  output:
    log2fc = config["OUT_PATH"] + "/08_log2fc/{norm_meth}.scale_{scal_fact}.bed_{bed_region}.{condition}.log2fc.tab"
  params:
    control = config["LOG2FC_PARAM"]["control"],
    condition = lambda wildcards: config["LOG2FC_PARAM"]["condition"][wildcards.condition]
  shell:
    '''
    module load StdEnv/2020
    module load gcc/9.3.0
    module load mugqic/R_Bioconductor/4.1.0_3.13
    module load r/4.1.2

    Rscript --vanilla /lustre06/project/6004736/alvann/from_narval/PIPELINES/DeepTools/scripts/calculate_log2fc.R \
                      {input.counts} \
                      {output.log2fc} \
                      "{params.control}"  \
                      "{params.condition}" \
                      "{config[MARKS]}"
    '''
