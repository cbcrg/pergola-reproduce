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
#############################################################################
### Jose Espinosa-Carrasco NPMMD/CB-CRG Group. May 2016                   ###
#############################################################################
### Reads mat lab files downloaded from the Wormbehavior DB               ###
### (http://wormbehavior.mrc-lmb.cam.ac.uk/). Each file contains several  ###
### features extracted from C.elegans video trackings.                    ### 
### Read features are foraging_speed, crawling, tail_motion and speeds    ###
### measured on different body parts                                      ###
#############################################################################

# Loading libraries
from argparse import ArgumentParser
from scipy.io  import loadmat
import h5py
from time import strptime
from calendar import timegm
from sys import stderr
import numpy as np
from csv import writer
from os.path import basename

parser = ArgumentParser(description='File input arguments')
parser.add_argument('-i','--input', help='Worms data hdf5 format matlab file', required=True)

args = parser.parse_args()

print >> stderr, "Input file: %s" % args.input

# Input files
input_file =  args.input

file_name = basename(input_file).split('.')[0]
file_name = file_name.replace (" ", "_")

f = h5py.File(input_file)

### INFO
# sex_r = f['info']['experiment']['worm']['sex']
# sex = str(''.join(unichr(c) for c in sex_r))

# habituation_r = f['info']['experiment']['worm']['habituation']
# habituation = str(''.join(unichr(c) for c in habituation_r))

# annotations (empty)
# annotations_r = f['info']['experiment']['environment']['annotations']
# annotations = str(''.join(c.astype(str) for c in annotations_r))

# info/experiment/worm/genotype
# genotype_r = f['info']['experiment']['worm']['genotype'] #type u2
# genotype = str(''.join(unichr(c) for c in genotype_r))

# /info/experiment/worm/strain
strain_r = f['info']['experiment']['worm']['strain']
strain = str(''.join(unichr(c) for c in strain_r))

# age worm
# /info/experiment/worm/age
# age_r = f['info']['experiment']['worm']['age'] #type u2
# age = str(''.join(unichr(c) for c in age_r))

# /info/experiment/environment/food
# food_r = f['info']['experiment']['environment']['food'] #type u2
# food = str(''.join(unichr(c) for c in food_r))

# /info/experiment/environment/timestamp
# timestamp_r = f['info']['experiment']['environment']['timestamp'] #type u2
# timestamp = str(''.join(unichr(c) for c in timestamp_r))

# HH:MM:SS.mmmmmm
# my_date_object = strptime(timestamp, '%Y-%m-%d %H:%M:%S.%f')
# unix_time = timegm(my_date_object) # utc based # correct!!!

# /info/video/length/time
# time_recorded_r = f['info']['video']['length']['time']
# time_recorded = time_recorded_r[0][0]

# /info/video/length/frames
frames_r = f['info']['video']['length']['frames']
frames = frames_r[0][0] 

fps_r = f['info']['video']['resolution']['fps']
fps = fps_r[0][0]

# extracted phenotypic features (speeds)
velocity_keys = ['head', 'headTip', 'midbody', 'tail', 'tailTip']

fh = open(file_name + "_speed.csv",'wb')

# fh.write("#genotype;%s\n" % genotype)
fh.write("#strain;%s\n" % strain)
# fh.write("#age;%s\n" % age)
# fh.write("#habituation;%s\n" % habituation)
# fh.write("#food;%s\n" % food)
# fh.write("#unix_time;%s\n" % unix_time)
# fh.write("#time_recorded;%s\n" % time_recorded)
fh.write("#frames;%s\n" % frames)
fh.write("#fps;%s\n" % fps)
# fh.write("#annotations;%s\n" % annotations)

writer_out = writer(fh, dialect = 'excel-tab')

writer_out.writerow(['frame_start', 'frame_end']  + sorted(velocity_keys) + ['foraging_speed', 'tail_motion', 'crawling'])

# foraging angle speed
try:
    foraging_speed = f['worm']['locomotion']['bends']['foraging']['angleSpeed']
except KeyError:
    raise KeyError ("Foraging angle speed is corrupted and can not be retrieved from hdf5 file")
                            
# tail motion
try:
    tail_motion = f['worm']['locomotion']['velocity']['tail']['direction']
except KeyError:
    raise KeyError ("Tail motion is corrupted and can not be retrieved from hdf5 file")

# Crawling
try:
    crawling = f['worm']['locomotion']['bends']['midbody']['amplitude']
except KeyError:
    raise KeyError ("Crawling is corrupted and can not be retrieved from hdf5 file")
 
# range already substract one to frames 
for frame in range(0, int(frames)):
    list_v = list()
    list_v.extend ([frame, frame+1])
    
    for velocity_k in sorted(velocity_keys):
        try:
            v = f['worm']['locomotion']['velocity'][velocity_k]['speed'][frame][0] 
        except KeyError:
            raise KeyError ("Velocity field %s is corrupted and can not be retrieved from hdf5 file"
                            % (velocity_k, frame))
                            
        if np.isnan(v) : v = -10000
        
        list_v.append (v)    
    
    v_fs = abs(foraging_speed[frame][0]) 
    v_tm = abs(tail_motion[frame][0])
    v_c = abs(crawling[frame][0])
    
    if np.isnan(v_fs) : v_fs = -10000
    if np.isnan(v_tm) : v_tm = -10000
    if np.isnan(v_c) : v_c = -10000
    
    list_v.append (v_fs)
    list_v.append (v_tm)
    list_v.append (v_c)
    
    writer_out.writerows([list_v])

fh.close()
