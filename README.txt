Geochemical modelling of chemical evolution of groundwater in the Pra Basin of Ghana

Authors: Evans Manu, Marco De Lucia and Michael Kühn
Institution: GFZ German Research Centre for Geoscience (Section 3.4 Fluid Systems Modelling)
For further inquiries, please contact Marco at delucia@gfz-potsdam.de 

---Combinatorial Inverse modelling with PHREEQC coupled with R

Description

1. This directory contains three PHREEQC script files representing groundwater input solutions from predefined zones: 'northern = Sol1', 'central = Sol2', and 'southern = Sol3', respectively. It is assumed that groundwater hypothetically flows from north to south, mirroring the regional topography.

2. The source code utilized for combinatorial inverse modelling is divided into two parts: one script containing all underlying functions (rfun) and separate codes for the northern (Northrtn_zone_good_data), central (Northern_Central_Good_data), and southern (Central_Southern_good_data) zones.

3. For clarity, each zone's script includes a starting solution: evaporated rainwater for the northern zone, groundwater from the northern zone for the central zone, and groundwater from the central zone for the southern zone.

4. The thermodynamic database used for simulation is a modified version of phreeqc.dat, incorporating minerals such as Phlogopite, Biotite, Plagioclase, and Organic matter, presented as phreeqc_invPra.dat.

5. Major ion concentrations for all groundwater samples in each zone are stored in an Excel file named Pra_data_M_Good.

---Reaction path modelling

Description

Based on the outcomes of combinatorial inverse modelling, batch reactions were incorporated to consistently replicate the observed groundwater composition. This process aimed to fine-tune the identified mineral assemblages and develop additional hypotheses.

The PHREEQC script for running the batching simulation is in 'Reaction_path_model_code'.

For a comprehensive understanding of combinatorial inverse and reaction path modelling, readers are encouraged to review the following papers:

a. https://doi.org/10.3390/w15213760 

b. https://doi.org/10.3390/w15071325

c. https://doi.org/10.3390/min13070899

The raw data has also been published in the GFZ Data services repository as open access

d. https://doi.org/10.5880/GFZ.3.4.2023.002
  
