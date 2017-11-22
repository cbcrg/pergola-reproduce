#!/usr/bin/env python

#  Copyright (c) 2014-2017, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2014-2017, Jose Espinosa-Carrasco and the respective authors.
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

from argparse import ArgumentParser
from os.path  import dirname, abspath, basename, splitext
from tempfile import NamedTemporaryFile
from pergola  import jaaba_parsers, mapping, intervals

parser = ArgumentParser(description='File input arguments')
parser.add_argument('-s','--score_file', help='Jaaba file containing scores annotated for a given behavior', required=True)
parser.add_argument('-m','--mapping_file', help='Mapping file to read JAABA files', required=True)
parser.add_argument('-t','--tag_group', help='Tag to distinguish data sets', required=True)
parser.add_argument('-b','--behavior', help='Annotated behavior using JAABA', required=True)

args = parser.parse_args()

chase_score_f = args.score_file
mappings_jaaba = mapping.MappingInfo(args.mapping_file)
tag_group = args.tag_group
annotated_behavior = args.behavior

tmp_track = NamedTemporaryFile(prefix='jaaba_csv', suffix='.csv', delete=True)
name_tmp = splitext(basename(tmp_track.name))[0]
path_out = dirname(abspath(tmp_track.name))

jaaba_parsers.jaaba_scores_to_csv(input_file=chase_score_f, path_w=path_out, name_file=name_tmp, norm=True, data_type=annotated_behavior)

scores_chase_int = intervals.IntData(tmp_track.name, map_dict=mappings_jaaba.correspondence).read()

data_type_col = {annotated_behavior: 'blue'}

if tag_group == "pBDPGAL4":
#    data_type_col = {annotated_behavior: 'orange'}
    data_type_col = {annotated_behavior: 'red'}

dict_bed_annotated_int = scores_chase_int.convert(mode="bed", color_restrictions=data_type_col)

chr_file_n = "chrom"
mapping.write_chr_sizes (scores_chase_int, file_n=chr_file_n)
mapping.write_chr(scores_chase_int)

chr_file = chr_file_n + ".sizes"

for tr, bed_tr in dict_bed_annotated_int.iteritems():
    id_worm, var = tr
    pybed_tr = bed_tr.create_pybedtools()
    pybed_tr.saveas('tr_' + id_worm + '_dt_' + tag_group + '.bed')
    df_by_tr = pybed_tr.genome_coverage(g=chr_file).to_dataframe()# .saveas('tr_' + id_worm + '.bed')
    print "%f\t%s" % (df_by_tr.score[(df_by_tr.chrom == 'chr1') & (df_by_tr.start == 1)].values[0], annotated_behavior + '_' + tag_group )
    # pybedtools.BedTool.to_dataframe()