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
 * Script to reproduce Pergola paper figures
 */

params.scores      = "$baseDir/data/scores/.mat"
params.variables   = "$baseDir/data/perframe/*.mat"
params.mappings    = "$baseDir/data/b2p.txt"
params.output      = "results/"

log.info "drosophila_jaaba - Pergola - Reproduce  -  version 0.1"
log.info "====================================="
log.info "annotated scores       : ${params.scores}"
log.info "variables directory    : ${params.var_dir}"
log.info "variables folder       : ${params.variables}"
log.info "mappings               : ${params.mappings}"
log.info "output                 : ${params.output}"
log.info "\n"

// Example command to run the script
//nextflow run drosophila_jaaba-Pergola-Reproduce.nf \
//	--scores='data/scores/*.mat' \
//  --variables='data/perframe/*.mat' \
//  -with-docker

// Extract all names of variables from a folder to pass them as arguments
// for f in *.mat; do     printf '%s ' "${f%.mat}"; done

/*
nextflow run melanogaster_GAL4-Pergola-Reproduce.nf --scores='data/scores/*.mat' --var_dir='data/perframe/' --variables="velmag dtheta" --mappings='data/jaaba2pergola.txt'

*/

/*
 * Input parameters validation
 */

mapping_file = file(params.mappings)

/*
 * Create a channel for scores
 */
Channel
	.fromPath( params.scores )
    .ifEmpty { error "Cannot find any mat file with Jaaba annotated scores" }
	.set { score_files }

score_files_tag = score_files.map {
	def content = it
	def name = it.name.replaceAll('scores_',' ').split("\\.")[0]
	//println ">>>>>>>>>>>>>>>>>>" + name
	[ content, name ]
}

score_files_tag.into { score_files_tag_bed; score_files_tag_comp }

/*
 * Create a channel for directory containing variables
 */
variable_dir = file( params.var_dir )
//variable_dir = Channel.fromPath( params.var_dir )
                      //.println ()

/*
 * Variable list to extract from the folder
 */
variables_list = params.variables.split(" ")
                  //.println (  )

process scores_to_bed {
    input:
    set file (scores), val (annotated_behavior) from score_files_tag_bed
    file mapping_file

    output:
    file 'results_score' into results_bed_score

    """
    jaaba_to_pergola sp -i ${scores} -m ${mapping_file} -f bed -bl -nt
    mkdir results_score
    mv *.bed results_score/
    """
}

process variables_to_bedGraph {
    input:
    file (variable_d) from variable_dir
    each var from variables_list
    file mapping_file

    output:
    set 'results_var', var into results_bedg_var

    """
    jaaba_to_pergola fp -i ${variable_d} -jf ${var}  -m ${mapping_file} -f bedGraph -nt
    mkdir results_var
    mv *.bedGraph results_var/
    """
}

process sushi_plot {
    input:
    set var_bedg_dir, var from results_bedg_var
    file scores_bed_dir from results_bed_score.first()    

    output:
    file "sushi_jaaba_scores_annot_${var}.png" into sushi_plot

    """
    sushi_pergola_bedAndBedGraph.R --path2variables=${var_bedg_dir} --path2scores=${scores_bed_dir} --variable_name=${var}
    """
}

process jaaba_scores_vs_variables {

  	input:
  	set file (scores), val (annotated_behavior) from score_files_tag_comp
  	val variable_d from variable_dir
  	each var from variables_list
    file mapping_file

  	output:
  	set file('results_annot'), var into annot_vs_non_annot_result

  	"""
  	jaaba_scores_vs_variables.py -s ${scores} -t ${annotated_behavior} -d ${variable_d} -v ${var} -m  ${mapping_file}
  	mkdir results_annot
  	mv *.txt  results_annot/

  	"""
}

process sign_variable_annotation {
    input:
    set file(dir_annot_vs_non_annot), var from annot_vs_non_annot_result

    output:
    file "${var}.png"
    set stdout into FC_pvalue

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
    file 'volcano_plot.png'

    """
    volcano_plot_jaaba.R --path2file=${pvalues_FC}
    """

}
