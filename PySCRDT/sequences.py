"""
Data and sequence holder class for PySCRDT 
"""
import numpy as np
from pathlib import Path
import json 
import xtrack as xt

sequence_folder = Path(__file__).resolve().parent.joinpath('../test_data').absolute()

class Sequences():
	"""Data containter for PS and SPS Pb ion sequences""" 
	def __init__(self, get_ps=True, get_sps=True):
	
		if get_ps:
			ps_path = '{}/PS_2022_Pb_ions_matched_with_RF.json'.format(sequence_folder)
			with open(ps_path, 'r') as fid:
				input_data_ps = json.load(fid)
			self.PS_Pb_line = xt.Line.from_dict(input_data_ps)

		if get_sps: 
			sps_path ='{}/SPS_2021_Pb_ions_matched_with_RF.json'.format(sequence_folder)
			with open(sps_path, 'r') as fid:
				input_data_sps = json.load(fid)   
			self.SPS_Pb_line = xt.Line.from_dict(input_data_sps)
