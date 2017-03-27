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
 * Wormbehavior DB (http://wormbehavior.mrc-lmb.cam.ac.uk/) processed by pergola for paper
 * Process mat files downloaded from the DB to extract locomotion phenotypes characterized 
 * as significantly different between a given strain and its N2 control subset
 * TODO  explain what pergola does
 */    

params.strain1_trackings = "$baseDir/small_data/unc_16/*.mat"
params.strain2_trackings = "$baseDir/small_data/N2/*.mat"
params.mappings_speed    = "$baseDir/small_data/mappings/worms_speed2p.txt"
params.mappings_bed      = "$baseDir/small_data/mappings/bed2pergola.txt"
params.mappings_motion   = "$baseDir/small_data/mappings/worms_motion2p.txt"

params.output      = "results/"

log.info "C. elegans locomotion phenotypes comparison - N F  ~  version 0.1"
log.info "========================================================="
log.info "C. elegans strain 1 data    : ${params.strain1_trackings}"
log.info "C. elegans strain 2 data    : ${params.strain2_trackings}"
log.info "mappings speed              : ${params.mappings_speed}"
log.info "mappings bed                : ${params.mappings_bed}"
log.info "mappings motion             : ${params.mappings_motion}"
log.info "output                      : ${params.output}"
log.info "\n"

/*
./N2_vs_unc16_motions-Pergola-Reproduce.nf --strain1_trackings 'data/unc_16/*.mat' --strain2_trackings 'data/N2/*.mat' \
	--mappings_speed 'data/mappings/worms_speed2p.txt' \
	--mappings_bed 'data/mappings/bed2pergola.txt' \
	--mappings_motion data/mappings/worms_motion2p.txt \
	-with-docker -resume
*/

/*
 * Input parameters validation
 */
map_speed = file(params.mappings_speed)
map_motion = file(params.mappings_motion)

/*
 * Input files validation
 */
if( !map_speed.exists() ) exit 1, "Missing mapping file: ${map_speed}"
if( !map_motion.exists() ) exit 1, "Missing mapping file: ${map_motion}" 

/*
 * Create a channel for strain 1 worm trackings 
 */
Channel
	.fromPath( params.strain1_trackings )
//    .println ()
    .ifEmpty { error "Cannot find any mat file with strain 1 data" }
	.set { strain1_files }
	
/*
 * Create a channel for strain 2 worm trackings 
 */
Channel
	.fromPath( params.strain2_trackings )
//    .println ()
    .ifEmpty { error "Cannot find any mat file with strain 2 data" }
	.set { strain2_files }

/*
 * Create a channel for strain 2 worm trackings 
 */
Channel
	.fromPath( params.mappings_bed )
//    .println ()
    .ifEmpty { error "Missing mapping file: ${map_motion}" }
	.set { map_bed }

params.tag_results = "tag"
tag_res = "${params.tag_results}"
tag_str1 = "strain1_worms"
tag_str2 = "strain2_worms"

/*
 * Creates a channel with file content and name of input file without spaces
 * Substitutes spaces by "_" in file name
 */ 
strain1_files_name = strain1_files.flatten().map { strain1_files_file ->      
	def content = strain1_files_file
	def name = strain1_files_file.name.replaceAll(/ /,'_')
	def tag = tag_str1
    [ content, name, tag ]
}

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
 * Get locomotion phenotypic features from mat files
 */ 
process get_feature {
	//container 'ipython/scipyserver'
  	
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
//body_parts =  ['head', 'headTip', 'midbody', 'tail', 'tailTip']
//body_parts =  ['head', 'headTip', 'midbody', 'tail', 'tailTip', 'foraging_speed', 'tail_motion', 'crawling']
body_parts =  ['midbody']
process feature_to_pergola {
	//container 'joseespinosa/pergola:celegans'
  
  	input:
  	set file ('speed_file'), val (name_file), val (exp_group) from locomotions_files 
  	file worms_speed2p from map_speed
  	each body_part from body_parts
  
  	output: 
  	set '*.no_na.bed', body_part, name_file into bed_loc_no_nas
  	set '*.no_na.bedGraph', body_part, name_file into bedGraph_loc_no_nas
  	set '*.no_tr.bedGraph', body_part, name_file into bedGraph_heatmap
  	
  	set name_file, body_part, '*.no_tr.bed', exp_group into bed_loc_no_track_line, bed_loc_no_track_line_cp
  	set name_file, body_part, '*.no_tr.bedGraph', exp_group into bedGraph_loc_no_track_line
  	
  	set '*.fa', body_part, name_file, val(exp_group) into out_fasta  	
  	
  	"""  	
  	cat $worms_speed2p | sed 's/behavioural_file:$body_part > pergola:dummy/behavioural_file:$body_part > pergola:data_value/g' > mod_map_file
  	pergola_rules.py -i $speed_file -m mod_map_file
  	pergola_rules.py -i $speed_file -m mod_map_file -f bedGraph -w 1 -min 0 -max 29000
  	
  	# This is done just because is easy to see in the display of the genome browsers
  	cat tr*.bed | sed 's/track name=\"1_a\"/track name=\"${body_part}\"/g' > bed_file.tmp
  	
  	cat tr*.bedGraph | sed 's/track name=\"1_a\"/track name=\"${body_part}\"/g' > bedGraph_file.tmp
  	
  	# delete values that were assigned as -10000 to skip na of the original file
  	# to avoid problems if a file got a feature with a feature always set to NA I add this code (short files for examples)
  	cat bed_file.tmp | grep -v "\\-10000" > ${name_file}".no_na.bed"  
  	cat ${name_file}".no_na.bed" | grep -v "track name" > ${name_file}.no_tr.bed || echo -e "chr1\t0\t100\t.\t-10000\t+\t0\t100\t135,206,250" > ${name_file}".no_tr.bed"
  	rm bed_file.tmp
  	
  	# delete values that were assigned as -10000 to skip na of the original file
  	# to avoid problems if a file got a feature with a feature always set to NA I add this code (short files for examples)
  	# cat bedGraph_file.tmp | grep -v "\\-10000" > ${name_file}".no_na.bedGraph"
  	cat bedGraph_file.tmp > ${name_file}".no_na.bedGraph"  
  	cat ${name_file}".no_na.bedGraph" | grep -v "track name" > ${name_file}".no_tr.bedGraph" || echo -e echo -e "chr1\t0\t100\t1" > ${name_file}".no_tr.bedGraph"
  	rm bedGraph_file.tmp 
  	"""
}

/*
 * Transform motion intervals from mat files (forward, backward and paused)
 */ 
process get_motion {
	//container 'ipython/scipyserver'
	  
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

//motion_files_wr.subscribe { println (it) }

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
	//container 'cbcrg/pergola:latest'  
  	
  	input:
  	set file ('motion_file'), val (name_file), val (name_file_motion), val (motion) from motion_files_flat
  	//set worms_motion_map from map_motion
    file worms_motion_map from map_motion
    
  	output:
  	set name_file, 'tr*.bed', name_file_motion into bed_motion, bed_motion_wr
  	  	
  	"""
  	pergola_rules.py -i $motion_file -m $worms_motion_map
  	cat tr_1_dt_${motion}.bed | sed 's/track name=\"1_a\"/track name=\"${motion}\"/g' > tr_1_dt_${motion}.bed.tmp
  	cat tr_1_dt_${motion}.bed.tmp | grep -v 'track name' > tr_1_dt_${motion}.bed
  	rm tr_1_dt_${motion}.bed.tmp  	
  	"""
} 

map_bed.into { map_bed_loc; map_bed_bG; map_bed_turn} //del remove turn

/*
 * Filter is used to delete pairs that do not come from the same original mat file
 */
//bed_motion.subscribe { println ("=========" + it) }  
//bed_loc_no_track_line.subscribe { println ("********" + it) }

bed_loc_motion = bed_loc_no_track_line
	.spread (bed_motion)
	.filter { it[0] == it[4] }

/*
 * Using bedtools intersect motion with phenotypic feature bed files
 */
process intersect_loc_motion {
	//container 'cbcrg/pergola:latest'
	
	input:
	set val (mat_file_loc), val (pheno_feature), file ('bed_loc_no_tr'), val (exp_group), val (mat_motion_file), file (motion_file), val (name_file_motion), val (direction) from bed_loc_motion
	file bed2pergola from map_bed_loc.first()
	
	output:
	set '*.mean.bed', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into bed_mean_speed_motion	
	set '*.mean.bedGraph', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into bedGr_mean_loc_motion
	set '*.intersect.bed', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into bed_intersect_loc_motion, bed_intersect_loc_motion2p, bed_intersect_l_m
	set '*.mean_file.bed', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into mean_intersect_loc_motion

	"""
	celegans_feature_i_motion.py -p $bed_loc_no_tr -m $motion_file -b $bed2pergola
	"""
}

/*
 * Intersected bed files transformed to bedgraph format for heatmaps visualizations
 */
process inters_to_bedGr {
	//container 'cbcrg/pergola:latest'
	
	input:
	set file (file_bed_inters), val (pheno_feature), val (mat_file_loc), val (mat_motion_file), val (name_file_motion), val (exp_group) from bed_intersect_l_m
	file bed2pergola from map_bed_bG.first()
	
	output:
	set '*.bedGraph', pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group into  bedGraph_intersect_loc_motion, bedGraph_intersect_loc_motion_str1, bedGraph_intersect_loc_motion_str2, bedGraph_intersect_loc_motion_str3, bedGraph_intersect_loc_motion_str4, bedGraph_intersect_loc_motion_str5, bedGraph_intersect_loc_motion_str6
	
	"""		
	if [ -s $file_bed_inters ]
	then		
		pergola_rules.py -i $file_bed_inters -m $bed2pergola -nh -s chrm start end nature value strain start_rep end_rep color -f bedGraph -w 1
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
	def exp_group =  it.name.split("\\.")[3] //OJO
	//def name_file = it.name.replaceAll('.strain2_worms', '').replaceAll('.strain1_worms', '')
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

bed_tagged.into { bed_tagged_for_str1_intro; bed_tagged_for_str1; bed_tagged_for_str1_f; bed_tagged_for_str1_b; bed_tagged_for_str1_p; bed_tagged_for_str2; bed_tagged_for_str2_f; bed_tagged_for_str2_b; bed_tagged_for_str2_p}

str1_bed_tagged = bed_tagged_for_str1.filter { it[4] == tag_str1 }
str2_bed_tagged = bed_tagged_for_str2.filter { it[4] == tag_str2 }

bed_tagged_for_str1_intro.subscribe { println ("=========" + it) }

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
bedGraph_heatmap_intro.subscribe { println ("=========" + it) }

str1_bedGraph_heatmap = bedGraph_heatmap_str1
							.filter { it[2] =~ /^unc-16.*/  }
//							.subscribe { println ("====str1===" + it) }	
							.map { it[0] }
str2_bedGraph_heatmap = bedGraph_heatmap_str2
							.filter { it[2] =~ /^N2.*/  }
//							.subscribe { println ("====str2===" + it) }	
							.map { it[0] }						

//bedGraph_intersect_loc_motion_str1

//str1_bedg_back = bedGraph_intersect_loc_motion_str2
//	.filter { it[5] == tag_str1 }
//	.filter { it[4]  =~ /^*.backward.*/ }
//	.map { it[0] }

/*
str1_bedg
	.toSortedList()
	.view()


str1_bedg	
	.toList()
    .subscribe onNext: { println it }, onComplete: 'Done'    
*/
/*
str1_bedg.toSortedList()
  .view()
*/


process heat_and_density_plot {

  	input:  	
  	file (str1_f)  from str1_bed_for
    file (str1_b) from str1_bed_back
    file (str1_p) from str1_bed_paused
    file (str2_f)  from str2_bed_for
    file (str2_b) from str2_bed_back
    file (str2_p) from str2_bed_paused
    
	file (str1_bedGraph_heatmap_list) from str1_bedGraph_heatmap.toSortedList()
    file (str2_bedGraph_heatmap_list) from str2_bedGraph_heatmap.toSortedList()
    
  	output:
  	// R creates a Rplots.pdf that is way we have to specify the tag "out" 
  	//set '*.png', strain, pheno_feature, direction into plots_pheno_feature_str1_str2
//    file 'file.txt' into result_d
    file 'heatmap_str1_str2.tiff' into heatmap
  	"""  	
  	mkdir str1
  	mkdir str2
  	
  	a=1
  	
  	for bedg in ${str1_bedGraph_heatmap_list}
      do
      	  cp \${bedg} str1/bedg.str1.\$a.bedGraph
  		  
  		  let a=a+1
  	  done
  	
  	a=1
  	for bedg in ${str2_bedGraph_heatmap_list}
      do
      	  cp \${bedg} str2/bedg.str2.\$a.bedGraph
  		  
  		  let a=a+1
  	  done
  	
  	path_str1=`pwd`"/str1"
  	path_str2=`pwd`"/str2"
  	echo \$path_str2 > file.txt
  	
  	heatmap_density_pheno_feature.R --path_str1=\$path_str1 --path_str2=\$path_str2 \
  		   --file_f_str1=${str1_f} --file_f_str2=${str2_f} \
  		   --file_b_str1=${str1_b} --file_b_str2=${str2_b} \
  		   --file_p_str1=${str1_p} --file_p_str2=${str2_p} 
  	"""
}

//result_d.subscribe { println ("=========" + it) }	
result_dir_heatmap = file("$baseDir/heatmap_$tag_res")
 
result_dir_heatmap.with {
     if( !empty() ) { deleteDir() }
     mkdirs()
     println "Created: $result_dir_heatmap"
}

heatmap.subscribe {	
	it.copyTo( result_dir_heatmap.resolve ( it.name ) )   
}



/*
process heat_and_density_plot {
  input:
  file bedg1_list from str1_bedg.collect()
  //file bedg2_list from str2_bedg.collect()
  
  output:
  set'*.tmp' into culo
  
  """
  mkdir bedg_str1
  mkdir bedg_str2
  
  a=1
  b=1
  for every bedg in ${bedg1_list}
    do
    	cp \$bedg "bedg_str1/${name_file_motion}"\$a".bedGraph"
    	let a=a+1
    done
  
  for every bedg in ${bedg2_list}
    do
    	cp $bedg "bedg_str1/${name_file_motion}"\$b".bedGraph"
    	let b=b+1
    done
    
  echo "culo" > cc.tmp
  """
}
*/

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
process dir_str1  {
	input:  	
  	set file (bedg), pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, tag_gr_str1 from str1_bedg_tagged
  	
  	output:
	file "bedg_${tag_gr_str1}" into bedg_out_dir_str1
  	
  	"""
  	mkdir bedg_${tag_gr_str1}
  	
  	cp $bedg "bedg_${tag_gr_str1}/${name_file_motion}.bedGraph" 
  	"""
}

process dir_str2  {
	input:  	
  	set file (bedg), pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, tag_gr_str2 from str2_bedg_tagged
  	
  	output:
  	file "bedg_${tag_gr_str2}" into bedg_out_dir_str2
    
  	"""
  	mkdir bedg_${tag_gr_str2}
  	
  	cp $bedg "bedg_${tag_gr_str2}/${name_file_motion}.bedGraph"  	
  	"""
}
*/

//set file (bedg), pheno_feature, mat_file_loc, mat_motion_file, name_file_motion, exp_group from str2_bedg_tagged
  
//bedg_out_dir_str2.subscribe { println "======================: $it"}

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
	//container 'joseespinosa/docker-r-ggplot2:v0.1'

  	input:
  	//set file (intersect_feature_motion_str1), strain, pheno_feature, direction, file (intersect_feature_motion_strain2) from str1_str2_bed_tag
  	set file (intersect_feature_motion_str1), strain, pheno_feature, direction, file (intersect_feature_motion_strain2) from str1_str2_bed_tag
  
  	output:
  	// R creates a Rplots.pdf that is way we have to specify the tag "out" 
  	//set '*.png', strain, pheno_feature, direction into plots_pheno_feature_str1_str2
    set '*out.pdf', strain, pheno_feature, direction into plots_pheno_feature_str1_str2
    
  	"""
  	plot_pheno_feature_distro_paper_version.R --bed_file_str1=${intersect_feature_motion_str1} --bed_file_str2=${intersect_feature_motion_strain2} 
  	"""
}

result_dir_distro = file("$baseDir/plots_distro_$tag_res")
 
result_dir_distro.with {
     if( !empty() ) { deleteDir() }
     mkdirs()
     println "Created: $result_dir_distro"
}

plots_pheno_feature_str1_str2.subscribe {
	/*it[0].copyTo( result_dir_distro.resolve ( it[1] + "." + it[2] + "." + it[3] + ".png" ) )*/
	it[0].copyTo( result_dir_distro.resolve ( it[1] + "." + it[2] + "." + it[3] + ".pdf" ) )   
}

/*
 * Creating folder to keep bed files to visualize data
 */
/*
result_dir_fasta = file("results_fasta_$tag_res")

result_dir_fasta.with {
     if( !empty() ) { deleteDir() }
     mkdirs()
     println "Created: $result_dir_fasta"
} 

out_fasta.subscribe {  
  fasta_file = it[0]
  fasta_file.copyTo( result_dir_fasta.resolve ( it[2] + ".fa" ) )
}

result_dir_bed = file("results_bed_$tag_res")

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

result_dir_bedGraph = file("results_bedGraph_$tag_res")

result_dir_bedGraph.with {
     if( !empty() ) { deleteDir() }
     mkdirs()
     println "Created: $result_dir_bedGraph"
} 

bedGraph_loc_no_nas.subscribe {   
  bedGraph_file = it[0]
  bedGraph_file.copyTo (result_dir_bedGraph.resolve ( it[1] + "." + it[2] + ".bedGraph" ) )
}

/*
bed_motion_wr.subscribe {
  bed_file = it[1]
  bed_file.copyTo ( result_dir_bed.resolve ( it[0] + it[2] + ".bed" ) )
}

result_dir_bed_intersect = file("$result_dir_bed/motion_intersected")

result_dir_bed_intersect.with {
     if( !empty() ) { deleteDir() }
     mkdirs()
     println "Created: $result_dir_bed_intersect"
} 

bed_intersect_loc_motion.subscribe {   
  bed_file = it[0]
  bed_file.copyTo ( result_dir_bed_intersect.resolve ( "intersect." + it[1] + "." + it[3] + "." + it[4] + ".bed" ) )
}
*/

result_dir_bedGraph_intersect = file("$result_dir_bedGraph/motion_intersected")

result_dir_bedGraph_intersect.with {
     if( !empty() ) { deleteDir() }
     mkdirs()
     println "Created: $result_dir_bedGraph_intersect"
} 

bedGraph_intersect_loc_motion.subscribe {
  bedGraph_file = it[0]
  bedGraph_file.copyTo ( result_dir_bedGraph_intersect.resolve ( "intersect." + it[1] + "." + it[3] + "." + it[4] + ".bedGraph" ) )
}
