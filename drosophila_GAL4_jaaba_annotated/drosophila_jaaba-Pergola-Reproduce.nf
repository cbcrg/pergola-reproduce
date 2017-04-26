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

/*
nextflow run drosophila_jaaba-Pergola-Reproduce.nf --scores='/Users/jespinosa/2017_sushi_pergola/data/drosophila_jaaba/scores/*.mat' \
                                                   --var_dir='/Users/jespinosa/2017_sushi_pergola/data/drosophila_jaaba/perframe/' \
                                                   --variables="velmag dtheta" \
                                                   --mappings='/Users/jespinosa/2017_sushi_pergola/lib/nxf/data/jaaba2pergola.txt'

nextflow run drosophila_jaaba-Pergola-Reproduce.nf --scores='/Users/jespinosa/2017_sushi_pergola/data/drosophila_jaaba/scores/*.mat' \
                                                   --var_dir='/Users/jespinosa/2017_sushi_pergola/data/drosophila_jaaba/perframe/' \
                                                   --variables="velmag" \
                                                   --mappings='/Users/jespinosa/2017_sushi_pergola/lib/nxf/data/jaaba2pergola.txt'

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


 /*
Channel
	.fromPath( params.variables )
    .ifEmpty { error "Cannot find any mat file with Jaaba derived variables from flies trajectory" }
	.set { variables_files }

variables_files_tag = variables_files.map {
	def content = it
	def name = it.name.split("\\.")[0]
	//println ">>>>>>>>>>>>>>>>>>" + name
	[ content, name ]
}
*/

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

// jaaba_to_pergola fp -i "/Users/jespinosa/JAABA_MAC_0.5.1/sampledata_v0.5/Chase1_TrpA_Rig1Plate15BowlA_20120404T141155/perframe/" -jf velmag dtheta  -m "/Users/jespinosa/git/pergola/test/jaaba2pergola.txt" -dd /Users/jespinosa/2017_sushi_pergola/data/ -f bedGraph -nt

process variables_to_bedGraph {
    input:
    val variable_d from variable_dir
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
    file scores_bed_dir from results_bed_score
    set var_bedg_dir, var from results_bedg_var

    output:
    file 'sushi_jaaba_scores_annot.png' into sushi_plot

    """
    sushi_pergola_bedAndBedGraph.R --path2variables=${var_bedg_dir} --path2scores=${scores_bed_dir} --var_name=${var}
    """
}

process jaaba_scores_vs_variables {

  	input:
  	set file (scores), val (annotated_behavior) from score_files_tag_comp
  	val variable_d from variable_dir
  	each var from variables_list
    file mapping_file

  	output:
  	set 'results_annot', var into annot_vs_non_annot_result

  	"""
  	jaaba_scores_vs_variables.py -s ${scores} -t ${annotated_behavior} -d ${variable_d} -v ${var} -m  ${mapping_file}
  	mkdir results_annot
  	mv *.txt  results_annot/

  	"""
}

process sign_variable_annotation {
    input:
    set dir_annot_vs_non_annot, var from annot_vs_non_annot_result

    output:
    file "${var}.png"

    """
    ttest_var_annotated_jaaba.R --path2files=${dir_annot_vs_non_annot} --variable_name=${var}
    """
}
