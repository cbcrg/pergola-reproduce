#!/usr/bin/env nextflow

/*
 *  Copyright (c) 2014-2017, Centre for Genomic Regulation (CRG).
 *  Copyright (c) 2014-2017, Jose Espinosa-Carrasco and the respective authors.
 *
 *  This file is part of Pergola.
 *
 *  Pergola is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Pergola is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Pergola.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Jose Espinosa-Carrasco. CB-CRG. March 2017
 *
 * Script to reproduce Pergola paper Jaaba annotated data analysis
 */

params.scores       = "$baseDir/small_data/scores/scores_chase.mat"
params.var_dir      = "$baseDir/small_data/perframe/"
params.var_dir_test = "$baseDir/"
params.variables    = "velmag"
params.mappings     = "$baseDir/small_data/jaaba2pergola.txt"
params.output       = "results/"
params.image_format = "tiff"
params.complete     = false

log.info "melanogaster_GAL4 - Pergola - Reproduce  -  version 0.1"
log.info "====================================="
log.info "annotated scores       : ${params.scores}"
log.info "variables directory    : ${params.var_dir}"
log.info "variables to run       : ${params.variables}"
log.info "mappings               : ${params.mappings}"
log.info "output                 : ${params.output}"
log.info "image format           : ${params.image_format}"
log.info "complete analysis      : ${params.complete}"
log.info "\n"

// Example command to run the script with the toy data provided in the repository
/*
nextflow run melanogaster_GAL4-Pergola-Reproduce.nf \
  --scores='small_data/scores/scores_chase_*.mat' \
  --var_dir='small_data/perframe_*' \
  --variables="velmag" \
  --mappings='small_data/mappings/jaaba2pergola.txt' \
  --output='results' \
  --image_format='png' \
  -with-docker
*/

/*
 * Input parameters validation
 */
mapping_file = file(params.mappings)

if( !mapping_file.exists() ) exit 1, "Missing mapping file: ${mapping_file}"

/*
 * Create a channel for scores
 * File naming is scores_behavior_strain.mat, thus tag corresponds to behavior and strain
 */
Channel
  .fromPath( params.scores )
  .ifEmpty { error "Cannot find any mat file with Jaaba annotated scores" }
  .set { score_files }

score_files_tag = score_files.map {
  def content = it
  def name = it.baseName.replaceAll('scores_','').split("\\.")[0].split("_")[-1]
  def behavior = it.baseName.replaceAll('scores_','').split("_")[0]
  [ content, name, behavior ]
}

score_files_tag.into { score_files_tag_bed; score_files_tag_comp }

/*
 * Create a channel for directory containing variables
 */
Channel.fromPath( params.var_dir, type: 'dir' )
                .ifEmpty { error "Cannot find variables directory" }
                .map {
                    def content = it
                    def name = it.name.split("\\_")[-1]
                    [ content, name ]
                }
                .set { variable_dir_tag }

variable_dir_tag.into { variable_dir_bg; variable_dir_scores }

/*
 * List of variable to extract from the folder. If set to "all", all variables are extracted.
 */
if ( params.variables == "all" ) {
  variables_list = "a absangle2wall absanglefrom1to2_anglesub absanglefrom1to2_nose2ell absdangle2wall absdtheta absdv_cor absphidiff_anglesub absphidiff_nose2ell abssmoothdtheta absthetadiff_anglesub absthetadiff_nose2ell absyaw accmag angle2wall anglefrom1to2_anglesub anglefrom1to2_nose2ell angleonclosestfly anglesub area areasmooth arena_angle arena_r b closestfly_anglesub closestfly_center closestfly_ell2nose closestfly_nose2ell closestfly_nose2ell_angle_30tomin30 closestfly_nose2ell_angle_min20to20 closestfly_nose2ell_angle_min30to30 closestfly_nose2tail corfrac_maj corfrac_min da dangle2wall danglesub darea db dcenter ddcenter ddell2nose ddist2wall ddnose2ell decc dell2nose dist2wall dnose2ell dnose2ell_angle_30tomin30 dnose2ell_angle_min20to20 dnose2ell_angle_min30to30 dnose2tail dphi dt dtheta du_cor du_ctr du_tail dv_cor dv_ctr dv_tail ecc flipdv_cor magveldiff_anglesub magveldiff_nose2ell phi phisideways signdtheta smoothdtheta theta timestamps velmag velmag_ctr velmag_nose velmag_tail veltoward_anglesub veltoward_nose2ell xnose_mm yaw ynose_mm".split(" ")
}
else {
  variables_list = params.variables.split(" ")
}

/*
 * Read image format
 */
image_format = "${params.image_format}"

/*
 * We use genome coverage to obtain the fraction of time performing the behavior
 */
process fract_time_behavior {
    publishDir = [path: "results/igv", mode: 'copy', pattern: "*_igv*"]

    input:
    set file (scores), val (tag_group), val (behavior) from score_files_tag_bed
    file mapping_file

    output:
    stdout into fraction_time
    set "results_score_${behavior}_${tag_group}", tag_group, behavior into results_bed_score
    set "results_score_igv_${behavior}_${tag_group}", tag_group into results_bed_score_igv
    file '*.fa' into out_fasta

    """
    cov_fraction_time_behavior.py -s ${scores} -m ${mapping_file} -t ${tag_group} -b ${behavior}
    mkdir results_score_${behavior}_${tag_group}
    mkdir results_score_igv_${behavior}_${tag_group}

    for f in *.bed
    do
        cat \${f} | sed 's/chase/\"\"/g' > \${f}.tmp
        cp \${f}.tmp results_score_${behavior}_${tag_group}/\${f}
        mv \${f}.tmp results_score_igv_${behavior}_${tag_group}/"`echo \${f} | sed s/L4//g | sed s/[a-z]*// | sed s/_//g | sed s/.bed//g | sed s/[a-zA-Z]//g`".bed
    done
    """
}

longest_fasta = out_fasta
                   .max { it.size() }

result_dir_IGV = file("results/igv")


longest_fasta.subscribe {
    fasta_file = it
    fasta_file.copyTo ( result_dir_IGV.resolve ( "drosophila.fa" ) )
}

fraction_time_comb_gr = fraction_time.collectFile()

/*
 * Compares the fraction of time spent of a given behavior on a boxplot
 */
process comparison_fract_time {
    publishDir = [path: "results/fraction_time", mode: 'copy']

    input:
    file (fraction_time_tbl) from fraction_time_comb_gr

    output:
    file "boxplot_fract_time*.${image_format}" into boxplot_fractime

    """
    boxplot_fract_time.R --path2tbl_fr_time=${fraction_time_tbl} \
        --image_format=${image_format}
    """
}

process variables_to_bedGraph {
    input:
    set file ('variable_d'), val (tag_group) from variable_dir_bg
    each var from variables_list
    file mapping_file

    output:
    set "results_var_${tag_group}_${var}", var, tag_group into results_bedg_var, bedGr_to_gviz

    when:
    params.complete

    script:
    """
    jaaba_to_pergola fp -i ${variable_d} -jf ${var} -m ${mapping_file} -f bedGraph -nt
    mkdir results_var_${tag_group}_${var}
    mv *.bedGraph results_var_${tag_group}_${var}/
    """
}

results_bed_score.into { results_bed_score_1; results_bed_score_2; results_bed_score_3 }

process sushi_plot_behavior_annot {
    publishDir = [path: "results/sushi", mode: 'copy']

    input:
    set scores_bed_dir, tag_group, behavior from results_bed_score_1

    output:
    file "*.${image_format}" into sushi_plot_annot

    """
    sushi_pergola_bed.R --path2scores=${scores_bed_dir} \
        --image_format=${image_format}

    mv "sushi_jaaba_annot.${image_format}" "sushi_${behavior}_${tag_group}.${image_format}"
    """
}

process gviz_plot_behavior_annot {
    publishDir = [path: "results/gviz", mode: 'copy']

    input:
    set scores_bed_dir, tag_group, behavior from results_bed_score_2

    output:
    file "*.${image_format}" into gviz_plot_var

    """
    melanogaster_gviz_annotations.R --path_bed_files=${scores_bed_dir} \
        --image_format=${image_format}

    mv "gviz_jaaba_annot.${image_format}" "gviz_${behavior}_${tag_group}.${image_format}"
    """
}

strain_output = "${params.output}"

score_files_tag_comp_var_dir = score_files_tag_comp
                                .spread ( variables_list )
                                .spread (variable_dir_scores)
                                .filter { it[1] == it[4] }

process jaaba_scores_vs_variables {
  	input:
  	set file (scores), val (behavior_strain), val(behavior), val (var), file ('variable_d'), val (strain) from score_files_tag_comp_var_dir

    file mapping_file

  	output:
    set "results_annot_${strain}_${var}", var, behavior_strain into annot_vs_non_annot_result
    set "results_bedGr_${strain}_${var}", var, behavior_strain into bedGr_to_sushi

    when:
    params.complete

    script:
  	"""
  	jaaba_scores_vs_variables.py -s ${scores} -t ${behavior_strain} -d ${variable_d} -v ${var} -m  ${mapping_file}
  	mkdir results_annot_${strain}_${var}
  	mkdir results_bedGr_${strain}_${var}

  	mv *.txt  results_annot_${strain}_${var}/
    mv *.bedGraph results_bedGr_${strain}_${var}/
  	"""
}

process sushi_plot_highlight_bg {
    input:
    set file (bedGr_dir), val (var), val (behavior_strain) from bedGr_to_sushi

    output:
    set "*.${image_format}", var, behavior_strain into sushi_plot_highlight

    when:
    params.complete

    script:
    """
    sushi_pergola_BedGraph_highlight.R --path2variables=${bedGr_dir} \
        --variable_name=${var} \
        --behavior_strain=${behavior_strain} \
        --image_format=${image_format}
    """
}

score_variables_vs_scores = results_bedg_var
                                .spread (results_bed_score_3)
                                .filter { it[2] == it[4] }

process sushi_plot_behavior_var {
    input:
    set var_bedg_dir, var, tag_group, scores_bed_dir, tag_group_score, behavior from score_variables_vs_scores

    output:
    file "*.${image_format}" into sushi_plot

    when:
    params.complete

    script:
    """
    sushi_pergola_bedAndBedGraph.R --path2variables=${var_bedg_dir} \
         --path2scores=${scores_bed_dir} \
         --variable_name=${var} \
         --image_format=${image_format}
    mv "sushi_jaaba_scores_annot_${var}.${image_format}" "sushi_${behavior}_${var}_${tag_group}.${image_format}"
    """
}

process gviz_plot_behavior_var {
    input:
    set var_bedg_dir, var_name, tag_group from bedGr_to_gviz

    output:
    file "*.${image_format}" into gviz_plot_annot

    when:
    params.complete

    script:
    """
    melanogaster_gviz_var.R  --path2variables=${var_bedg_dir} \
        --behavior_strain=${tag_group} \
        --variable_name=${var_name} \
        --image_format=${image_format}

    mv "gviz_jaaba_var.${image_format}" "gviz_jaaba_var_${var_name}_${tag_group}.${image_format}"
    """
}

process significance_variable_annotation {
    input:
    set file(dir_annot_vs_non_annot), var, strain from annot_vs_non_annot_result

    output:
    file "boxplot_${var}_${strain}.${image_format}"
    stdout into FC_pvalue

    when:
    params.complete

    script:
    """
    ttest_var_annotated_jaaba.R --path2files=${dir_annot_vs_non_annot} \
        --variable_name=${var} \
        --image_format=${image_format}

    mv ${var}.${image_format} boxplot_${var}_${strain}.${image_format}
    """
}

FC_pvalues_collected = FC_pvalue
                        .collectFile(name: 'FC_pvalue.csv', newLine: false)

process plot_volcano {
    input:
    file pvalues_FC from FC_pvalues_collected

    output:
    file "*.pdf"
    file 'tbl_fc_pvalues.txt'

    when:
    params.complete

    script:
    """
    volcano_plot_jaaba.R --path2file=${pvalues_FC}
    """
}
