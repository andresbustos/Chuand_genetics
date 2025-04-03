#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 10:23:43 2020

@author: u5501
"""


import genetic_library as gl
import os
import re
import sys

	
##################################################
parameters = gl.read_parameters('parameters_PMuS.json')

if parameters['merge_fastqs']:
	for pp in parameters['dict_conf'].keys():
		new_fastq_name            = os.path.join(parameters['outfolder'], pp+'_merged_file.fastq')
		gl.merge_fastq(parameters['dict_conf'][pp]['files'], new_fastq_name)

if parameters['analyze_fastqs']:
	for pp in parameters['dict_conf'].keys():
		for f in parameters['dict_conf'][pp]['files']:
			gl.analyze_fastq(parameters, f)	
		
gl.PMuS(parameters)
gl.transversal_PMuS_analysis(parameters)
