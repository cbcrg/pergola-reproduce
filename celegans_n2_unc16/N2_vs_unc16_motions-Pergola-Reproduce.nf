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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Pergola.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Jose Espinosa-Carrasco. CB/CSN-CRG. April 2016
 *
 * Script to reproduce Pergola paper figures using C.elegans trackings downloaded from the
 * Wormbehavior DB (http://wormbehavior.mrc-lmb.cam.ac.uk/). 
 * 
 */    

params.strain1_trackings = "$baseDir/small_data/unc_16/*.mat"
params.strain2_trackings = "$baseDir/small_data/N2/*.mat"
params.mappings_speed    = "$baseDir/small_data/mappings/worms_speed2p.txt"
params.mappings_bed      = "$baseDir/small_data/mappings/bed2pergola.txt"
params.mappings_motion   = "$baseDir/small_data/mappings/worms_motion2p.txt"
params.output            = "results/"
params.tag_results       = ""
params.image_format      = "tiff"

log.info "C. elegans locomotion phenotypes comparison - N F  ~  version 0.1"
log.info "========================================================="
log.info "C. elegans strain 1 data    : ${params.strain1_trackings}"
log.info "C. elegans strain 2 data    : ${params.strain2_trackings}"
log.info "mappings speed              : ${params.mappings_speed}"
log.info "mappings bed                : ${params.mappings_bed}"
log.info "mappings motion             : ${params.mappings_motion}"
log.info "output                      : ${params.output}"
log.info "tag results                 : ${params.tag_results}"
log.info "image format                : ${params.image_format}"
log.info "\n"

/*
nextflow run N2_vs_unc16_motions-Pergola-Reproduce.nf --strain1_trackings 'small_data/unc_16/*.mat' --strain2_trackings 'small_data/N2/*.mat' \
	--mappings_speed 'small_data/mappings/worms_speed2p.txt' \
	--mappings_bed 'small_data/mappings/bed2pergola.txt' \
	--mappings_motion small_data/mappings/worms_motion2p.txt \
	-with-docker
*/

/*
 * Input parameters validation
 */
map_speed = file(params.mappings_speed)
map_motion = file(params.mappings_motion)

/*
 * Input files validation
 */
if( !map_speed.exists() ) exit 1, "Missing speed mapping file: ${map_speed}"
if( !file(params.mappings_bed).exists() ) exit 1, "Missing bed mapping file: ${params.mappings_bed}" 
if( !map_motion.exists() ) exit 1, "Missing motion mapping file: ${map_motion}" 

/*
 * Create a channel for strain 1 worm trackings 
 */
Channel
    .fromPath( params.strain1_trackings )
    .ifEmpty { error "Cannot find any mat file with strain 1 data" }
	.set { strain1_files }
	
/*
 * Create a channel for strain 2 worm trackings 
 */
Channel
	 .fromPath( params.strain2_trackings )
    .ifEmpty { error "Cannot find any mat file with strain 2 data" }
	 .set { strain2_files }

/*
 * Create a channel for pergola mappings for bed files
 */
Channel
	.fromPath( params.mappings_bed )
    .ifEmpty { error "Missing mapping file: ${params.mappings_bed}" }
	.set { map_bed }

/*
 * Read tag for results if exists 
 */
tag_res = "${params.tag_results}"
tag_str1 = "strain1_worms"
tag_str2 = "strain2_worms"

/*
 * Read image format
 */
image_format = "${params.image_format}"

/*
 * Creates a channel with file content and name of input file without spaces
 * Substitutes spaces by "_" in file name - strain 1
 */
strain1_files_name = strain1_files.flatten().map { strain1_files_file ->      
  	def content = strain1_files_file
  	def name = strain1_files_file.name.replaceAll(/ /,'_')
  	def tag = tag_str1
    [ content, name, tag ]
}

/*
 * Creates a channel with file content and name of input file without spaces
 * Substitutes spaces by "_" in file name - strain 2
 */
strain2_files_name = strain2_files.flatten().map { strain2_files_file ->      
  	def content = strain2_files_file
  	def name = strain2_files_file.name.replaceAll(/ /,'_')
  	def tag = tag_str2
    [ content, name, tag ]
}

// Files tagged joined in a single channel
trackings_files_name = strain2_files_name.mix ( strain1_files_name )
trackings_files_name.into { trackings_loc; trackings_motion }

/*
 * Get locomotion phenotypic features from mat files (speed)
 */ 
process get_feature {
  	
  	input:
  	set file ('file_worm'), val (name_file_worm), val (exp_group) from trackings_loc
  
  	output:  
  	set '*_speed.csv', name_file_worm, exp_group into locomotions_files, locomotion_files_wr
  	
  	script:
  	println "Matlab file containing worm behavior processed: $name_file_worm"

  	"""
  	extract_worm_speed.py -i $file_worm
  	"""
}

/*
 * Transform locomotion files into bed format files
 */ 
body_parts =  ['midbody']

process feature_to_pergola {

  	input:
  	set file ('speed_file'), val (name_file), val (exp_group) from locomotions_files 
  	file worms_speed2p from map_speed
  	each body_part from body_parts
  
  	output: 
  	set '*.no_na.bed', body_part, name_file into bed_loc_no_nas
  	//set '*.no_na.bed', body_part, name_file, 'chrom.sizes' into bed_cov
  	set '*.zeros.bedGraph', body_part, name_file into bedGraph_loc_no_nas
  	set '*.no_tr.bedGraph', body_part, name_file into bedGraph_heatmap
  	
  	set name_file, body_part, '*.no_tr.bed', exp_group into bed_loc_no_track_line, bed_loc_no_track_line_cp
  	set name_file, body_part, '*.no_tr.bedGraph', exp_group into bedGraph_loc_no_track_line
  	
  	set '*.fa', body_part, name_file, val(exp_group) into out_fasta  	
  	
  	"""  	
  	cat $worms_speed2p | sed 's/behavioral_file:$body_part > pergola:dummy/behavioral_file:$body_part > pergola:data_value/g' > mod_map_file
  	pergola -i $speed_file -m mod_map_file
  	pergola -i $speed_file -m mod_map_file -f bedGraph -w 1 -min 0 -max 29000
  	
  	# This is done just because is easy to see in the display of the genome browsers
  	cat tr*.bed | sed 's/track name=\"1_a\"/track name=\"${body_part}\"/g' > bed_file.tmp
  	
  	cat tr*.bedGraph | sed 's/track name=\"1_a\"/track name=\"${body_part}\"/g' > bedGraph_file.tmp
  	
  	# delete values that were assigned as -10000 to skip na of the original file
  	# to avoid problems if a file got a feature with a feature always set to NA I add this code (short files for examples)
  	cat bed_file.tmp | grep -v "\\-10000" > ${name_file}".no_na.bed"  
  	cat ${name_file}".no_na.bed" | grep -v "track name" > ${name_file}.no_tr.bed || echo -e "chr1\t0\t100\t.\t-10000\t+\t0\t100\t135,206,250" > ${name_file}".no_tr.bed"
  	rm bed_file.tmp
  	
  	# delete values that were assigned as -10000 to skip na of the original file
  	# to avoid problems if a file got a feature with a feature always set to NA I add this code (short files for 
  	cat bedGraph_file.tmp | sed 's/-10000/0/g' > ${name_file}".zeros.bedGraph"
  	cat bedGraph_file.tmp > ${name_file}".no_na.bedGraph"  
  	cat ${name_file}".no_na.bedGraph" | grep -v "track name" > ${name_file}".no_tr.bedGraph" || echo -e echo -e "chr1\t0\t100\t1" > ${name_file}".no_tr.bedGraph"
  	rm bedGraph_file.tmp 
  	"""
}

/*
 * Transform motion intervals from mat files (forward, backward and paused)
 */ 
process get_motion {
	  
  	input:
  	set file ('file_worm'), val (name_file_worm) from trackings_motion
  
  	output: 
  	set name_file_worm, '*.csv' into motion_files, motion_files_wr
    
  	script:
  	println "Matlab file containing worm motion processed: $name_file_worm"

  	"""
  	extract_worm_motion.py -i \"$file_worm\"
  	"""
}

/*
process motion_to_pergola {

  	input:
  	set file (bed_file), val (body_part), val (name_file) from bed_cov
    falta chrom sizes
    file worms_motion_map from map_motion

  	output:
  	set name_file, 'tr*.bed', name_file_motion into bed_motion, bed_motion_wr

  	"""
  	pergola -i $motion_file -m $worms_motion_map
  	cat tr_1_dt_${motion}.bed | sed 's/track name=\"1_a\"/track name=\"${motion}\"/g' > tr_1_dt_${motion}.bed.tmp
  	cat tr_1_dt_${motion}.bed.tmp | grep -v 'track name' > tr_1_dt_${motion}.bed
  	rm tr_1_dt_${motion}.bed.tmp
  	"""
}
*/

/*
 * From one mat file 3 motion (forward, paused, backward) files are obtained
 * A channel is made with matfile1 -> forward
 *                        matfile1 -> backward
 *                        matfile1 -> paused
 *                        matfile2 -> forward ...
 */
motion_files_flat_p = motion_files.map { name_mat, motion_f ->
                          motion_f.collect { 
        	                def motion = it.name.split("\\.")[1]   
                          [ it, name_mat, it.name, motion ]
                          }
                      }
                      .flatMap()

motion_files_flat_p.into { motion_files_flat; motion_files_flat_wr }

process motion_to_pergola { 

  	input:
  	set file ('motion_file'), val (name_file), val (name_file_motion), val (motion) from motion_files_flat
    file worms_motion_map from map_motion
    
  	output:
  	set name_file, 'tr*.bed', name_file_motion into bed_motion, bed_motion_wr
  	  	
  	"""
  	pergola -i $motion_file -m $worms_motion_map
  	cat tr_1_dt_${motion}.bed | sed 's/track name=\"1_a\"/track name=\"${motion}\"/g' > tr_1_dt_${motion}.bed.tmp
  	cat tr_1_dt_${motion}.bed.tmp | grep -v 'track name' > tr_1_dt_${motion}.bed
  	rm tr_1_dt_${motion}.bed.tmp  	
  	"""
} 

map_bed.into { map_bed_loc; map_bed_bG; map_bed_turn} //del remove turn

/*
 * Filter is used to delete pairs that do not come from the same original mat file
 */
bed_loc_motion = bed_loc_no_track_line
	.spread (bed_motion)
	.filter { it[0] == it[4] }

/*
 * Using bedtools intersect motion with phenotypic feature bed files
 */
process intersect_loc_motion {

	  input:
	  set val (mat_file_loc), val (pheno_feature), file ('bed_loc_no_tr'), val (exp_group), val (mat_motion_file), file (motion_file), val (name_file_motion), val (direction) from bed_loc_motion
	  file bed2pergola from map_bed_loc.first()
	
	  output:
	  set '*.mean.bed', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into bed_mean_speed_motion	
	  //set '*.mean.bedGraph', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into bedGr_mean_loc_motion
	  set '*.intersect.bed', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into bed_intersect_loc_motion, bed_intersect_loc_motion2p, bed_intersect_l_m
	  //set '*.mean_file.bed', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into mean_intersect_loc_motion

	  """
	  intersectBed -a $bed_loc_no_tr -b $motion_file  > mean_speed_i_motionDir.intersect.bed

      awk '{ sum += \$5; n++ } END { if (n > 0) print sum / n; }' mean_speed_i_motionDir.intersect.bed > mean_speed_i_motionDir.mean.bed

	  # celegans_feature_i_motion.py -p $bed_loc_no_tr -m $motion_file -b $bed2pergola
	  """
}

/*
 * Intersected bed files transformed to bedgraph format for heatmaps visualizations
 */
process inters_to_bedGr {
	
	  input:
	  set file (file_bed_inters), val (pheno_feature), val (mat_file_loc), val (mat_motion_file), val (name_file_motion), val (exp_group) from bed_intersect_l_m
	  file bed2pergola from map_bed_bG.first()
	
	  output:
	  set '*.bedGraph', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into  bedGraph_intersect_loc_motion, bedGraph_intersect_loc_motion_str1, bedGraph_intersect_loc_motion_str2, bedGraph_intersect_loc_motion_str3, bedGraph_intersect_loc_motion_str4, bedGraph_intersect_loc_motion_str5, bedGraph_intersect_loc_motion_str6
	
	  """		
	  if [ -s $file_bed_inters ]
	  then		
		  pergola -i $file_bed_inters -m $bed2pergola -nh -s chrm start end nature value strain start_rep end_rep color -f bedGraph -w 1
	  else
		  touch tr_chr1_d.bedGraph
      fi
	  """
}

/*
 * Grouping (collect) bed files in order to plot the distribution by strain, motion direction and body part 
 */
bed_intersect_loc_motion_plot = bed_intersect_loc_motion2p.collectFile(newLine: false, sort:'none') { 
	  def name = it[3].split("_on_")[0] + "." + it[1] + "." +  it[4].tokenize(".")[1] + "." +  it[5]
	  [ name, it[0].text]	
}.map {  
	  def strain =  it.name.split("\\.")[0]	
	  def pheno_feature =  it.name.split("\\.")[1]	
	  def direction =  it.name.split("\\.")[2]
	  def exp_group =  it.name.split("\\.")[3]
      def name_file = it.name.replaceAll("." + tag_str2, '').replaceAll("." + tag_str1, '')

 	  [ it, strain, pheno_feature, direction, name_file, exp_group ]
}

/*
 * Tagging files for plotting
 */
process tag_bed_mean_files {
  	input: 
  	set file ('bed_file'), val (strain), val (pheno_feature), val (direction), val (strain_beh_dir), val (exp_group) from bed_intersect_loc_motion_plot
  	
  	output:
  	set '*.bed', strain, pheno_feature, direction, exp_group into bed_tagged
  	
  	"""
  	# Adds to the bed file a tag for being used inside the R dataframe
  	awk '{ print \$0, \"\\t$pheno_feature\\t$direction\\t$strain\" }' ${bed_file} > ${strain_beh_dir}.bed  
  	"""
}

bed_tagged.into { bed_tagged_for_str1; bed_tagged_for_str1_f; bed_tagged_for_str1_b; bed_tagged_for_str1_p; bed_tagged_for_str2; bed_tagged_for_str2_f; bed_tagged_for_str2_b; bed_tagged_for_str2_p}

str1_bed_tagged = bed_tagged_for_str1.filter { it[4] == tag_str1 }
str2_bed_tagged = bed_tagged_for_str2.filter { it[4] == tag_str2 }

str1_bed_for = bed_tagged_for_str1_f
    .filter { it[4] == tag_str1 }
    .filter { it[3] == "forward" }
    .map { it[0] }

str1_bed_back = bed_tagged_for_str1_b
    .filter { it[4] == tag_str1 }
    .filter { it[3] == "backward" }
    .map { it[0] }
				
str1_bed_paused = bed_tagged_for_str1_p
	  .filter { it[4] == tag_str1 }
	  .filter { it[3] == "paused" }
	  .map { it[0] }

str2_bed_for = bed_tagged_for_str2_f
    .filter { it[4] == tag_str2 }
    .filter { it[3] == "forward" }
    .map { it[0] }

str2_bed_back = bed_tagged_for_str2_b
    .filter { it[4] == tag_str2 }
    .filter { it[3] == "backward" }
    .map { it[0] }
				
str2_bed_paused = bed_tagged_for_str2_p
    .filter { it[4] == tag_str2 }
    .filter { it[3] == "paused" }
    .map { it[0] }

bedGraph_heatmap.into { bedGraph_heatmap_intro; bedGraph_heatmap_str1; bedGraph_heatmap_str2}

str1_bedGraph_heatmap = bedGraph_heatmap_str1
                        		.filter { it[2] =~ /^unc-16.*/  }
                        		.map { it[0] }
str2_bedGraph_heatmap = bedGraph_heatmap_str2
              							.filter { it[2] =~ /^N2.*/  }
              							.map { it[0] }						

/*
 * Plotting heat map and density plots
 */
process heat_and_density_plot {

  	input:  	
  	file (str1_f) from str1_bed_for
    file (str1_b) from str1_bed_back
    file (str1_p) from str1_bed_paused
    file (str2_f) from str2_bed_for
    file (str2_b) from str2_bed_back
    file (str2_p) from str2_bed_paused

    // This might be simplied as in the melanogaster example
    file (str1_bedGraph_heatmap_list) from str1_bedGraph_heatmap.toSortedList()
    file (str2_bedGraph_heatmap_list) from str2_bedGraph_heatmap.toSortedList()
    
  	output:
    file "heatmap_str1_str2.${image_format}" into heatmap
    file ('results_bedgr1') into results_bedgr1
    file ('results_bedgr2') into results_bedgr2
    file ('results_bedgr1_igv') into results_bedgr1_igv
    file ('results_bedgr2_igv') into results_bedgr2_igv

  	"""  	
  	mkdir str1
  	mkdir str2
  	mkdir igv_str1
  	mkdir igv_str2
  	a=1
  	
  	for bedg in ${str1_bedGraph_heatmap_list}
      do
      	  cp \${bedg} str1/bedg.str1.\${a}.bedGraph
  		  { echo -e "track name=\${a} description=\${a} visibility=full"; cat \$bedg; } > igv_str1/\${a}.\$bedg.new
  		  cat igv_str1/\${a}.\$bedg.new | sed 's/-10000/0/g' > igv_str1/\${a}.\$bedg
  		  rm igv_str1/*.new
  		  # mv igv_str1/\${a}.\$bedg{.new,}
  		  let a=a+1
  	  done
  	
  	a=1
  	for bedg in ${str2_bedGraph_heatmap_list}
      	do
      	  cp \${bedg} str2/bedg.str2.\${a}.bedGraph
  	      { echo -e "track name=\${a} description=\${a} visibility=full"; cat \$bedg; } > igv_str2/\${a}.\$bedg.new
  	      cat igv_str2/\${a}.\$bedg.new | sed 's/-10000/0/g' > igv_str2/\${a}.\$bedg
  	      rm igv_str2/*.new
  	      # mv igv_str2/\${a}.\$bedg{.new,}

  	      let a=a+1
      	done
  	
  	path_str1=`pwd`"/str1"
  	path_str2=`pwd`"/str2"
  	echo \$path_str2 > file.txt
  	
  	heatmap_density_pheno_feature.R --path_str1=\$path_str1 --path_str2=\$path_str2 \
  		   --file_f_str1=${str1_f} --file_f_str2=${str2_f} \
  		   --file_b_str1=${str1_b} --file_b_str2=${str2_b} \
  		   --file_p_str1=${str1_p} --file_p_str2=${str2_p} \
  		   --image_format=${image_format}

  	mv \$path_str1 results_bedgr1/
  	mv \$path_str2 results_bedgr2/

  	mv `pwd`"/igv_str1" results_bedgr1_igv/
  	mv `pwd`"/igv_str2" results_bedgr2_igv/
  	"""
}

results_bedgr_1 = results_bedgr1.map { [ it, 'unc-16' ] }
results_bedgr_2 = results_bedgr2.map { [ it, 'N2' ] }

results_bedgr_sushi = results_bedgr_1.concat (results_bedgr_2)

process heatmap_sushi {
    input:
    set var_bedg_dir, tag_group from results_bedgr_sushi

    output:
    file "*.${image_format}" into sushi_heatmap

    """
    heatmap_sushi.R --path_bedgr=${var_bedg_dir} --image_format=${image_format}
    mv "sushi_var.${image_format}" "sushi_heatmap_${tag_group}.${image_format}"
    """
}

result_dir_heatmap = file("$baseDir/heatmap$tag_res")
 
result_dir_heatmap.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_heatmap"
}

heatmap.subscribe {	
	  it.copyTo( result_dir_heatmap.resolve ( it.name ) )   
}

/*
result_dir_bedg_str1 = file("$baseDir/results$tag_str1")
result_dir_bedg_str2 = file("$baseDir/results$tag_str2")

result_dir_bedg_str1.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_bedg_str1"
}

str1_bedg_tagged.subscribe {   
    bedGraph_file = it[0]
    bedGraph_file.copyTo (result_dir_bedg_str1.resolve ( it[1] + "." + it[2] + ".bedGraph" ) )
}

result_dir_bedg_str2.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_bedg_str2"
}

str2_bedg_tagged.subscribe {   
    bedGraph_file = it[0]
    bedGraph_file.copyTo (result_dir_bedg_str2.resolve ( it[1] + "." + it[2] + ".bedGraph" ) )
}

dir_bedg_str1 = file("$baseDir/results$tag_str1")
dir_bedg_str2 = file("$baseDir/results$tag_str2")
*/

/*
 * Matching the control group for each strain in the data set
 */
str1_str2_bed_tag = str1_bed_tagged
    .spread (str2_bed_tagged)
	  //same feature and motion direction
	  .filter { it[2] == it[7] && it[3] == it[8] }
	  .map { [ it[0], it[1], it[2], it[3], it[5] ] }	

/*
 * Plots the distribution of bed containing all intervals by strain, motion and body part compairing the distro of ctrl and case strain
 */
process plot_distro {

  	input:
  	set file (intersect_feature_motion_str1), strain, pheno_feature, direction, file (intersect_feature_motion_strain2) from str1_str2_bed_tag
  
  	output:
  	// R creates a Rplots.pdf that is way we have to specify the tag "out" 
  	//set '*.png', strain, pheno_feature, direction into plots_pheno_feature_str1_str2
    set "*out.${image_format}", strain, pheno_feature, direction into plots_pheno_feature_str1_str2

  	"""
  	plot_pheno_feature_distro_paper_version.R --bed_file_str1=${intersect_feature_motion_str1} \
  	    --bed_file_str2=${intersect_feature_motion_strain2} \
  	    --image_format=${image_format}
  	"""
}

longest_fasta = out_fasta     
                   .max { it.size() }

result_dir_distro = file("$baseDir/plots_distro$tag_res")
 
result_dir_distro.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_distro"
}

plots_pheno_feature_str1_str2.subscribe {
	  it[0].copyTo( result_dir_distro.resolve ( it[1] + "." + it[2] + "." + it[3] + ".pdf" ) )   
}

/*
 * Creating folder to keep all files for data visualization on IGV
 */
result_dir_IGV = file("results_IGV$tag_res")

result_dir_IGV.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_IGV"
} 

longest_fasta.subscribe { 
    fasta_file = it[0]    
    fasta_file.copyTo ( result_dir_IGV.resolve ( "celegans.fa" ) )
}


result_dir_bed = file("results_bed$tag_res")

/*
result_dir_bed.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_bed"
} 

bed_loc_no_nas.subscribe {  
    bed_file = it[0]
    bed_file.copyTo ( result_dir_bed.resolve ( it[1] + "." + it[2] + ".bed" ) )
}
*/

/*
result_dir_bedGraph = file("results_bedGraph$tag_res")

result_dir_bedGraph.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_bedGraph"
} 
*/
/*
bedGraph_loc_no_nas.subscribe {
    bedGraph_file = it[0]
    bedGraph_file.copyTo (result_dir_IGV.resolve ( it[1] + "." + it[2] + ".bedGraph" ) )
}


results_bedgr1_igv.subscribe {
    it.copyTo (result_dir_IGV.resolve ( it[1] + "." + it[2] + ".bedGraph" ) )
}
*/

results_bedgr1_igv.subscribe {
    it.copyTo (result_dir_IGV)
}

results_bedgr2_igv.subscribe {
    it.copyTo (result_dir_IGV)
}

/*
results_bedgr1_igv.subscribe {
    bedGraph_file = it[0]
    bedGraph_file.copyTo (result_dir_IGV)
}

results_bedgr2_igv.subscribe {
    bedGraph_file = it[0]
    bedGraph_file.copyTo (result_dir_IGV)
}
*/

result_dir_IGV_intersect = file("$result_dir_IGV/motion_intersected")

result_dir_IGV_intersect.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created:   $result_dir_IGV_intersect"
} 

bedGraph_intersect_loc_motion.subscribe {
    bedGraph_file = it[0]
    bedGraph_file.copyTo ( result_dir_IGV_intersect.resolve ( "intersect." + it[1] + "." + it[3] + "." + it[4] + ".bedGraph" ) )
}

bed_motion_wr.subscribe {
    bed_direction = it[1]
    bed_direction.copyTo ( result_dir_bed.resolve ( it[0] + "_" + it[2] + "_direction" + ".bed" ))
}

/*
result_dir_csv = file("results_csv$tag_res")

locomotion_files_wr.subscribe {
    file_speed = it[0]
    file_speed.copyTo ( result_dir_csv.resolve ( it[1] + ".csv" ) )
}
*/

/*

i=1
for f in *N2*.csv
do
     echo -e "id\tframe_start\tframe_end\thead\theadTip\tmidbody\ttail\ttailTip\tforaging_speed\ttail_motion\tcrawling" > "N2_$i.csv"
     cat $f | grep -v "#"  | sed 's/-10000/0/g' | awk -v i="$i" 'NR>1 {print "N2_"i"\t"$0}' >> "N2_$i.csv"
     echo $i
     i=$[$(echo $i) + 1]
done

i=1
for f in *unc*.csv
do
     echo -e "id\tframe_start\tframe_end\thead\theadTip\tmidbody\ttail\ttailTip\tforaging_speed\ttail_motion\tcrawling" > "unc16_$i.csv"
     cat $f | grep -v "#" | sed 's/-10000/0/g' | awk -v i="$i" 'NR>1 {print "unc16_"i"\t"$0}' >> "unc16_$i.csv"
     echo $i
     i=$[$(echo $i) + 1]
done

pergola -i ./*.csv -m ./worm_speed2pergola.txt -f bedGraph -w 1 -min 0 -max 29000

awk '{ sum += $5; n++ } END { if (n > 0) print sum / n; }' mean_speed_i_motionDir.intersect.bed

*/