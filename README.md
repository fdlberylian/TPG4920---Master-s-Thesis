# TPG4920-Masters-Thesis
The title of the thesis is "Calculation and Visualization of Energy Dissipation and Energy Balance in Reservoir Models". You can see the full thesis at the file named "Master's Thesis - Fadhil Berylian.pdf". Section 3.1 of the thesis details how the scripts work, and the workflow is presented in the thesis as Figure 3.2.

The scripts are written in Python. The objective of the scripts is to extract necessary data from the reservoir visualization software ResInsight, then process the data to get energy changes in the reservoir system during simulation of field production. Obtaining visualization of how energy is used in reservoir is important, especially during simulation phase prior to starting production. Results from these scripts could provide considerations to optimize production scenario by modifying well placement/well control/etc.

List of code scripts:
1) energy_reservoir_norne.py : script to fetch required reservoir data from ResInsight, calculate energy changes in reservoir based on the data, and transfer results back to ResInsight.
2) energy_well_norne.py : script to fetch required well data from ResInsight to calculate energy changes in well.
3) energy_balance_norne.py : script to calculate total energy changes in reservoir and well, and produce plots of the results.
