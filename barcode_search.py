#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 10:10:36 2020

@author: u5501
"""


import genetic_library as gl
import os

parameters = gl.read_parameters('parameters_BCS.json')

if parameters['merge_fastqs']:
 	new_fastq_name            = os.path.join(parameters['outfolder'], 'merged_file.fastq')
 	gl.merge_fastq(parameters['fastq_names'], new_fastq_name)
 	parameters['fastq_names'] = [new_fastq_name] #from now on we work with the merged fastq

if parameters['analyze_fastqs']:
 	for f in parameters['fastq_names']:
		 gl.analyze_fastq(parameters, f)

if parameters['search_lib']:
 	gl.BCS(parameters)
 	gl.transversal_BCS_analysis(parameters, parameters["N_barplot"])

if parameters['find_coincidences']:
	gl.find_coincidences(parameters)

