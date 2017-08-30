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

params.scores      = "$baseDir/small_data/scores/scores_chase.mat"
params.var_dir     = "$baseDir/small_data/perframe/"
params.var_dir_test = "$baseDir/"
params.variables   = "velmag"
params.mappings    = "$baseDir/small_data/jaaba2pergola.txt"
params.output      = "results/"

log.info "drosophila_jaaba - Pergola - Reproduce  -  version 0.1"
log.info "====================================="
log.info "annotated scores       : ${params.scores}"
log.info "variables directory    : ${params.var_dir}"
log.info "variables to run       : ${params.variables}"
log.info "mappings               : ${params.mappings}"
log.info "output                 : ${params.output}"
log.info "\n"

// Example command to run the script with the toy data provided in the repository
/*
nextflow run melanogaster_GAL4-Pergola-Reproduce.nf \
  --scores='small_data/scores/scores_chase_*.mat' \
  --var_dir='small_data/perframe/' \
  --variables="dnose2ell dtheta velmag" \
  --mappings='small_data/jaaba2pergola.txt' \
  -with-docker

// by strain
nextflow run melanogaster_GAL4-Pergola-Reproduce.nf \
  --scores='small_data/scores/scores_chase_case_TrpA.mat' \
  --var_dir='small_data/perframe_TrpA/' \
  --variables="dnose2ell dtheta velmag" \
  --mappings='small_data/jaaba2pergola.txt' \
  --output='TrpA'
  -with-docker

nextflow run melanogaster_GAL4-Pergola-Reproduce.nf \
  --scores='small_data/scores/scores_chase_ctrl_pBDPGAL4.mat' \
  --var_dir='small_data/perframe_pBDPGAL4/' \
  --variables="dnose2ell dtheta velmag" \
  --mappings='small_data/jaaba2pergola.txt' \
  --output='pBDPGAL4'
  -with-docker
*/

/*
 * Input parameters validation
 */
mapping_file = file(params.mappings)

if( !mapping_file.exists() ) exit 1, "Missing mapping file: ${mapping_file}"
if( !file(params.var_dir).exists() ) exit 1, "Missing variable directory: ${params.var_dir}"

/*
 * Create a channel for scores
 * File naming is scores_behavior_strain.map, thus tag corresponds to behavior and strain
 */
Channel
  .fromPath( params.scores )
  .ifEmpty { error "Cannot find any mat file with Jaaba annotated scores" }
  .set { score_files }

score_files_tag = score_files.map {
  def content = it
  def name = it.baseName.replaceAll('scores_','').split("\\.")[0]
  [ content, name ]
}

score_files_tag.into { score_files_tag_bed; score_files_tag_comp; score_files_print }
score_files_print.println()

/*
 * Create a channel for directory containing variables
 */
variable_dir = Channel.fromPath( params.var_dir )
//variable_dir_test = Channel
//    	       	.fromPath( "${params.var_dir_test}/*.txt"  )
//    		.set { var_dir_test_ch }

//var_dir_test_ch.println()

variable_dir.into { variable_dir_bg; variable_dir_scores; variable_dir_print}

variable_dir_print.println()

/*
 * List of variable to extract from the folder. If set to "all", all variables are extracted.
 */
//if( params.variables == "" ) exit 1, "Variables to extract not provided: ${params.variables}"

if ( params.variables == "all" ) {
  variables_list = "a absangle2wall absanglefrom1to2_anglesub absanglefrom1to2_nose2ell absdangle2wall absdtheta absdv_cor absphidiff_anglesub absphidiff_nose2ell abssmoothdtheta absthetadiff_anglesub absthetadiff_nose2ell absyaw accmag angle2wall anglefrom1to2_anglesub anglefrom1to2_nose2ell angleonclosestfly anglesub area areasmooth arena_angle arena_r b closestfly_anglesub closestfly_center closestfly_ell2nose closestfly_nose2ell closestfly_nose2ell_angle_30tomin30 closestfly_nose2ell_angle_min20to20 closestfly_nose2ell_angle_min30to30 closestfly_nose2tail corfrac_maj corfrac_min da dangle2wall danglesub darea db dcenter ddcenter ddell2nose ddist2wall ddnose2ell decc dell2nose dist2wall dnose2ell dnose2ell_angle_30tomin30 dnose2ell_angle_min20to20 dnose2ell_angle_min30to30 dnose2tail dphi dt dtheta du_cor du_ctr du_tail dv_cor dv_ctr dv_tail ecc flipdv_cor magveldiff_anglesub magveldiff_nose2ell phi phisideways signdtheta smoothdtheta theta timestamps velmag velmag_ctr velmag_nose velmag_tail veltoward_anglesub veltoward_nose2ell xnose_mm yaw ynose_mm".split(" ")
}
else {
  variables_list = params.variables.split(" ")             
}

/*
 * We use genome coverage to obtain the fraction of time performing the behavior
 */
process frac_time_behavior {
    input:
    set file (scores), val (tag_group) from score_files_tag_bed
    file mapping_file

    output:
    stdout into fraction_time
    file 'tr*.bed' into bed_score_cov
    file 'chrom.sizes' into chrom_sizes
    set 'results_score', tag_group into results_bed_score, results_bed_score_2, results_bed_score_3

    """
    cov_fraction_time_behavior.py -s ${scores} -m  ${mapping_file} -t ${tag_group}
    mkdir results_score
    # mv *.bed results_score/
    cp *.bed results_score/
    """
}

fraction_time_comb_gr = fraction_time.collectFile()
//fraction_time_comb_gr.println()

/*
 * Compares the fraction of time spent of a given behavior on a boxplot
 */
process comparison_fract_time {
    input:
    file ('fraction_time_tbl') from fraction_time_comb_gr

    output:
    file "boxplot_fract_time*.pdf" into boxplot_fractime

    """
    boxplot_fract_time.R --path2tbl_fr_time=${fraction_time_tbl}
    """
}


process variables_to_bedGraph {
    input:
    file ('variable_d') from variable_dir_bg.first()
    each var from variables_list
    file mapping_file

    output:
    set 'results_var', var into results_bedg_var, bedGr_to_gviz

    """
    jaaba_to_pergola fp -i ${variable_d} -jf ${var} -m ${mapping_file} -f bedGraph -nt
    mkdir results_var
    mv *.bedGraph results_var/
    """
}

/*
The resulting plot is not used by the moment in the paper, eventually delete
*/

process sushi_plot {
    input:
    set var_bedg_dir, var from results_bedg_var
    set scores_bed_dir, tag_group from results_bed_score.first()    

    output:
    file "sushi_jaaba_scores_annot_${var}.png" into sushi_plot

    """
    sushi_pergola_bedAndBedGraph.R --path2variables=${var_bedg_dir} --path2scores=${scores_bed_dir} --variable_name=${var}
    """
}

//score_files_tag_comp.into { score_files_tag_comp; score_files_tag_comp2 } //del
score_files_tag_comp_var = score_files_tag_comp
    .spread ( variables_list )


process jaaba_scores_vs_variables {

  	input:
  	//set file (scores), val (behavior_strain) from score_files_tag_comp
  	set file (scores), val (behavior_strain), val (var) from score_files_tag_comp_var
  	file ('variable_d') from variable_dir_scores.first()
  	//each var from variables_list
        file mapping_file

  	output:
  	set file('results_annot'), var into annot_vs_non_annot_result
    set file('results_bedGr'), var, behavior_strain into bedGr_to_sushi

  	"""
  	jaaba_scores_vs_variables.py -s ${scores} -t ${behavior_strain} -d ${variable_d} -v ${var} -m  ${mapping_file}
  	mkdir results_annot
  	mkdir results_bedGr

  	mv *.txt  results_annot/
    mv *.bedGraph results_bedGr/
  	"""
}

bedGr_to_sushi.into { bedGr_to_sushi; bedGr_to_sushi2 }
//bedGr_to_sushi2.println()

process sushi_plot_highlight_bg {
    input:
    set file (bedGr_dir), val (var), val (behavior_strain) from bedGr_to_sushi

    output:
    set "*.pdf", var, behavior_strain into sushi_plot2

    //script:
    //behavior_strain.println()

    """
    sushi_pergola_BedGraph_highlight.R --path2variables=${bedGr_dir} --variable_name=${var} --behavior_strain=${behavior_strain}
    """
}

process sushi_plot_behavior_annot {
    input:
    //set scores_bed_dir, tag_group from results_bed_score_2.first()
    set scores_bed_dir, tag_group from results_bed_score_2

    output:
    file "*.pdf" into sushi_plot_annot

    """
    sushi_pergola_bed.R --path2scores=${scores_bed_dir}
    mv sushi_jaaba_annot.pdf "sushi_jaaba_annot_"${tag_group}".pdf"
    """
}


process gviz_plot_behavior_var {
    input:
    set var_bedg_dir, var_name from bedGr_to_gviz

    output:
    file "*.tiff" into gviz_plot_annot

    """
    melanogaster_gviz_var.R  --path2variables=${var_bedg_dir} --variable_name=${var_name}
    mv gviz_jaaba_var.tiff "gviz_jaaba_var_"${var_name}".tiff"
    """
}


process gviz_plot_behavior_annot {
    input:
    set scores_bed_dir, tag_group from results_bed_score_3

    output:
    file "*.tiff" into gviz_plot_var

    """
    melanogaster_gviz_annotations.R --path_bed_files=${scores_bed_dir}
    mv gviz_jaaba_annot.tiff "gviz_jaaba_annot_"${tag_group}".tiff"
    """
}

process sign_variable_annotation {
    input:
    set file(dir_annot_vs_non_annot), var from annot_vs_non_annot_result

    output:
    file "${var}.pdf"
    stdout into FC_pvalue

    """
    ttest_var_annotated_jaaba.R --path2files=${dir_annot_vs_non_annot} --variable_name=${var}
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

    """
    volcano_plot_jaaba.R --path2file=${pvalues_FC}
    """

}
