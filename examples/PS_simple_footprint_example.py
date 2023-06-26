"""
#  PySCRDT - simple example to generate tune footprint with the PS
#  
#  example by Elias Waagaard 
"""
import matplotlib.pyplot as plt
import xtrack as xt
import json 

from PySCRDT import PySCRDT, tune_footprint_maker

#### Load the sequence ####
fname_line = 'example_data/PS_2022_Pb_ions_matched_with_RF.json'
with open(fname_line, 'r') as fid:
     input_data = json.load(fid)    
line = xt.Line.from_dict(input_data)
particle_ref = line.particle_ref
line.build_tracker()
twiss_xtrack = line.twiss() 

#### PS default Pb beam settings #### 
bunch_intensity = 3.5e8 
sigma_z = 4.74
nemitt_x= 0.8e-6
nemitt_y= 0.5e-6 

#### Analytical tune footprint settings ###
Qh = 6.210000511153454
Qv = 6.24499998913536
plot_range  =   [[5.95,6.3],[5.95,6.3]]   # range in Qh & Qv for the plot
plot_order  =   5   # order of resonances to plot
periodicity =   16  # periodicity of ring for the colorcode of the plot


#### Create instance of PySCRDT, taking normalized emittances as input ####
s = PySCRDT()
s.setParameters(
    intensity = bunch_intensity,
    bunchLength = sigma_z,
    emittance_x = nemitt_x,
    emittance_y = nemitt_y, 
    dpp_rms = 1e-3,  # very small contribution anyway to beam size
    bF=None,
    ro = particle_ref.get_classical_particle_radius0() 
)
s.loadTwissFromXsuite(twissTableXsuite=twiss_xtrack)

#### Initiate tune footprint maker ####
PS_tune_footprint = tune_footprint_maker(Qh, Qv, plot_range = plot_range)

#### Plot the tune footprint ####
fig = plt.figure(figsize=(8, 8))
PS_tune_footprint.generate_tune_footprint(fig, s)