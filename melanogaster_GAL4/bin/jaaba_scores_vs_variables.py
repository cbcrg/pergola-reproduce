#!/usr/bin/env python

#  Copyright (c) 2014-2016, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2014-2016, Jose Espinosa-Carrasco and the respective authors.
#
#  This file is part of Pergola.
#
#  Pergola is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Pergola is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Pergola.  If not, see <http://www.gnu.org/licenses/>.
############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Dec 2016                        ###
############################################################################
### Script process jaaba derived scores (chasing) and var derived from   ###
### the trajectory using Pergola. Then uses pybedtools to intersect      ###
### annotated periods and scores and complement of the annotated periods.###
### Then uses R for compare means                                        ###
############################################################################

from argparse   import ArgumentParser
from sys        import stderr
from os.path    import dirname, abspath, basename, splitext
from tempfile   import NamedTemporaryFile
from pergola    import parsers, jaaba_parsers, mapping, intervals
from pybedtools import BedTool

parser = ArgumentParser(description='File input arguments')
parser.add_argument('-s','--score_file', help='Jaaba file containing scores annotated for a given behavior', required=True)
parser.add_argument('-t','--type_annotated', help='Annotated type in score Jaaba file', required=True)
parser.add_argument('-d','--variables_dir', help='Variables file containing variables derived from the trajectory using Jaaba', required=True)
parser.add_argument('-v','--variable_name', help='Annotated type in score Jaaba file', required=True)
parser.add_argument('-m','--mapping_file', help='Mapping file to read Jaaba files', required=True)

args = parser.parse_args()

print >> stderr, "Score file: %s" % args.score_file
print >> stderr, "Annotated type: %s" % args.type_annotated
print >> stderr, "Variable directory: %s" % args.variables_dir
print >> stderr, "Variable name: %s" % args.variable_name
print >> stderr, "Mapping file: %s" % args.mapping_file

chase_score_f = args.score_file
annotated_type = args.type_annotated
var_traj = args.variable_name
path_var_jaaba = args.variables_dir
# map_file_jaaba = "/Users/jespinosa/git/pergola/test/jaaba2pergola.txt"
mappings_jaaba = mapping.MappingInfo(args.mapping_file)

# mapping_bed = mapping.MappingInfo("/Users/jespinosa/git/pergola/test/c_elegans_data_test/bed2pergola.txt")

# tmp file out
tmp_track = NamedTemporaryFile(prefix='jaaba_csv', suffix='.csv', delete=True)
name_tmp = splitext(basename(tmp_track.name))[0]
path_out = dirname(abspath(tmp_track.name))

jaaba_parsers.jaaba_scores_to_csv(input_file=chase_score_f, path_w=path_out, name_file=name_tmp, norm=True, data_type=annotated_type)

scores_chase_int = intervals.IntData (tmp_track.name, map_dict=mappings_jaaba.correspondence).read()
# mapping.write_chr_sizes (scores_chase_int, path_w=out_dir, file_n=chr_file_n)
chr_file_n = "chrom"
mapping.write_chr_sizes (scores_chase_int, file_n=chr_file_n)
chr_file = chr_file_n + ".sizes"

#jaaba_to_pergola sp -i "/Users/jespinosa/git/pergola/sample_data/jaaba_example/scores_chase.mat" -m /Users/jespinosa/git/pergola/test/jaaba2pergola.txt -f bed -bl -nt

bed_annotated_int = scores_chase_int.convert(mode="bed")

## Get the variables that Jaaba processes from the trajectory
# path_var_jaaba = '/Users/jespinosa/JAABA_MAC_0.5.1/sampledata_v0.5/Chase1_TrpA_Rig1Plate15BowlA_20120404T141155/perframe/'
#jaaba_to_pergola fp -i "/Users/jespinosa/JAABA_MAC_0.5.1/sampledata_v0.5/Chase1_TrpA_Rig1Plate15BowlA_20120404T141155/perframe/" -jf velmag dtheta  -m "/Users/jespinosa/git/pergola/test/jaaba2pergola.txt" -dd /Users/jespinosa/2017_sushi_pergola/data/ -f bedGraph -nt 


# var_traj = ["velmag", "dtheta"]
# var_traj = "dtheta"
# var_traj = "velmag"
# var_traj = "velmag_tail"
# var_traj = "dtheta"
# var_traj = "dell2nose"
# var_traj = "angle2wall"

# var_traj = "phi"

dict_bedGraph_features = dict()
# statistic="mean"
statistic="distinct"

dict_bedGraph_var = jaaba_parsers.extract_jaaba_features(dir_perframe=path_var_jaaba, map_jaaba=args.mapping_file, delimiter="\t",
                                                             feature=var_traj, output="IntData").read().convert(mode="bedGraph", window=1)

list_no_intervals = [("chr1", 0, 1, 0)]

for tr, bedGraph_var in dict_bedGraph_var.iteritems():

    id_worm, var = tr

    bedGraph_var_bt = bedGraph_var.create_pybedtools()

    ## only if the track has annotations it is possible to map annotations into the variable
    ## generates a bedGraph with the variable values intersecting annotations if possible and those out of annotations
    if (id_worm, annotated_type) in bed_annotated_int:

        bed_annotated_int_bt = bed_annotated_int[(id_worm, annotated_type)].create_pybedtools()

        if bedGraph_var_bt.count() != 0 and bed_annotated_int_bt.count() != 0:
            # bedtools map - a 3chase.bed - b 3chase.bedGraph - c 4 - o mean
            bed_annotated_int_comp = bed_annotated_int_bt.complement(g=chr_file)

            bed_annotated_int_bt.map(bedGraph_var_bt, c=4, o=statistic, null=0).saveas('values_' + id_worm + '_' + var_traj + '.txt')
            bed_annotated_int_comp.map(bedGraph_var_bt, c=4, o=statistic, null=0).saveas('values_' + id_worm + '_' + var_traj + '.comp.txt')

            bedGraph_var_bt.intersect(bed_annotated_int_bt).saveas('values_' + id_worm + '_' + var_traj + '.bedGraph')
            bedGraph_var_bt.intersect(bed_annotated_int_comp).saveas('values_' + id_worm + '_' + var_traj + '.comp.bedGraph')

            # bedGraph_var_bt.(bed_annotated_int_bt).saveas('values_' + id_worm + '_' + var_traj + '.comp.bedGraph')
            #
            # .intersect(motion_bt).saveas(tag_file + ".intersect.bed")
            #
            # bed_annotated_int_bt.intersect()
            # (bedGraph_var_bt, c=4, o=statistic, null=0).saveas('values_' + id_worm + '_' + var_traj + '.txt')
            # bed_annotated_int_comp.map(bedGraph_var_bt, c=4, o=statistic, null=0).saveas('values_' + id_worm + '_' + var_traj + '.comp.txt')

    ## if no annotations present then returns bedGraph for its plotting and and empty bedgraph for annotated regions
    else:
        bedGraph_var_bt.saveas('values_' + id_worm + '_' + var_traj + '.comp.bedGraph')
        bed_no_intervals = BedTool(list_no_intervals).saveas('values_' + id_worm + '_' + var_traj + '.bedGraph')

