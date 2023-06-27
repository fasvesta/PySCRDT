"""
#  Simple tune footprint maker example with PySCRDT to generate tune footprint
#  with Pb ions in the PS
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
Qh = twiss_xtrack['qx']
Qv = twiss_xtrack['qy']
#plot_range  =   [[5.95,6.3],[5.95,6.3]]   # plot range in Qh & Qv, can be custom-provided to tune_footprint_maker

#### Create instance of PySCRDT, taking beam parameters as input ####
s = PySCRDT()
s.setParameters(
    intensity = bunch_intensity,
    bunchLength = sigma_z,
    emittance_x = nemitt_x,
    emittance_y = nemitt_y, 
    dpp_rms = 1e-3, 
    bF=None,
    ro = particle_ref.get_classical_particle_radius0() 
)
s.loadTwissFromXsuite(twissTableXsuite=twiss_xtrack)

#### Initiate tune footprint maker ####
PS_tune_footprint = tune_footprint_maker(Qh, Qv)

#### Plot the tune footprint ####
fig = plt.figure(figsize=(8, 8))
PS_tune_footprint.generate_tune_footprint(fig, s)