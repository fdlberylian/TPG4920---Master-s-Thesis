import rips
import time
import grpc
import math
import os
from operator import itemgetter
from ecl.eclfile import EclFile
from ecl.grid import EclGrid
from ecl.summary import EclSum
import matplotlib.pyplot as plt

# Define the function to calculate energy changes in well bottomholes
def energywell():
    # Connect to ResInsight
    resinsight = rips.Instance.find()
    case = resinsight.project.cases()[0]
    num_tsteps = len(case.time_steps())
    name = case.name
    grids = case.grids()
    for grid in grids:
        dimension = grid.dimensions()
    Nx = dimension.i
    Ny = dimension.j
    Nz = dimension.k

    class well:
        def __init__(self, name, idx, welltype):
            self.name = name
            self.idx = idx
            self.type = welltype      


    # Read EGRID, RST and INIT files
    summary_file = EclSum("%s.UNSMRY" % name)
    egrid_file = EclFile("%s.EGRID" % name)
    rst_file = EclFile("%s.UNRST" % name)
    timestep_width = []
    days = []
    for tstep in range(num_tsteps):
        if tstep==0:
            width = rst_file.iget_restart_sim_days(tstep)
        else:
            width = rst_file.iget_restart_sim_days(tstep) - rst_file.iget_restart_sim_days(tstep-1)
        
        timestep_width.append(width)
        days.append(rst_file.iget_restart_sim_days(tstep))

    # Read ACTNUM numbers from EGRID file
    actnum = egrid_file["ACTNUM"][0]
    active_cells = []
    for i in range(len(actnum)):
        if actnum[i] == 1:
            active_cells.append(i)

    # Convert summary file timesteps to restart file timesteps
    summary_days = summary_file.days
    idx_trim = []
    for d in days[2:]:
        add = summary_days.index(d)
        idx_trim.append(add)

    energy_balance = []
    energy_external = []
    energy_internal = []
    energy_dissipated = []
    for tstep in range(2,num_tsteps):
        print("Timestep", tstep, "of", num_tsteps)

        # List active wells in the timestep
        zwel = rst_file.iget_named_kw("ZWEL",tstep)
        nwells = int(len(zwel)/3)
        well_list = []
        for wel in range(nwells):
            welname = rst_file["ZWEL"][tstep][wel*3]
            welname = welname.rstrip()
            niwelz = int(len(rst_file["IWEL"][tstep]) / nwells)

            if rst_file["IWEL"][tstep][wel*niwelz + 6] == 1:
                weltype = "PROD"
            else:
                weltype = "INJE"
            
            ncwmax = int(len(rst_file["ICON"][tstep]) / (25*nwells))
            welidx = []
            for ncw in range(ncwmax):
                numicon = 25*(wel*ncwmax + ncw)
                weli = rst_file["ICON"][tstep][numicon+1]
                welj = rst_file["ICON"][tstep][numicon+2]
                welk = rst_file["ICON"][tstep][numicon+3]
                if weli != 0:
                    welidx.append(int(welk-1)*(Nx*Ny) + int(welj-1)*(Nx) + int(weli-1))

            well_list.append(well(welname,welidx,weltype))

        # Read results into list
        porv = case.active_cell_property('STATIC_NATIVE', 'PORV', 0)
        pres = case.active_cell_property('DYNAMIC_NATIVE', 'PRESSURE', tstep)
        bo = case.active_cell_property('DYNAMIC_NATIVE', 'BO', tstep)
        bw = case.active_cell_property('DYNAMIC_NATIVE', 'BW', tstep)
        bg = case.active_cell_property('DYNAMIC_NATIVE', 'BG', tstep)

        # Fetch results from the reservoir energy calculation
        edis = case.active_cell_property('GENERATED', 'Energy Dissipation', tstep)
        eint = case.active_cell_property('GENERATED', 'Internal Energy Change', tstep)

        # Calculate energy changes in well bottomholes
        e_external = []
        e_internal_well = []
        ed_well = []
        for wel in well_list:
            idx = idx_trim[tstep-2]
            pcel = case.active_cell_property('DYNAMIC_NATIVE', 'PRESSURE', tstep)
            for wel_idx in wel.idx:
                wel_idx_act = active_cells.index(wel_idx)

                pwel = pcel[wel_idx_act]
                bo_wel = bo[wel_idx_act]
                bw_wel = bw[wel_idx_act]
                bg_wel = bg[wel_idx_act]
                bhp = summary_file["WBHP:%s" % wel.name][idx].value

                if wel.type == "PROD":
                    wpr = summary_file["CWPR:%s:%i" % (wel.name,wel_idx+1)][idx].value
                    opr = summary_file["COPR:%s:%i" % (wel.name,wel_idx+1)][idx].value
                    gpr = summary_file["CGPR:%s:%i" % (wel.name,wel_idx+1)][idx].value
                    e_external.append(-bhp * (wpr*bw_wel+opr*bo_wel+gpr*bg_wel) * 1e5/86400)
                    e_internal_well.append(-pwel * (wpr*bw_wel+opr*bo_wel+gpr*bg_wel) * 1e5/86400)
                    ed_well.append((wpr*bw_wel+opr*bo_wel+gpr*bg_wel) * (pwel-bhp) * 1e5/86400)
                else:
                    wir = summary_file["CWIR:%s:%i" % (wel.name,wel_idx+1)][idx].value
                    gir = summary_file["CGIR:%s:%i" % (wel.name,wel_idx+1)][idx].value
                    e_external.append(bhp * (wir*bw_wel+gir*bg_wel) * 1e5/86400)
                    e_internal_well.append(pwel * (wir*bw_wel+gir*bg_wel) * 1e5/86400)
                    ed_well.append((wir*bw_wel+gir*bg_wel) * (bhp-pwel) * 1e5/86400)

        edis_res = sum(edis)
        edis_well = sum(ed_well)
        e_dissipated = edis_res + edis_well
        e_internal = sum(eint) + sum(e_internal_well)
        
        # Calculate each component in energy balance
        e_balance = (sum(e_external) - (e_dissipated + e_internal))
        energy_balance.append(e_balance)
        energy_dissipated.append(e_dissipated)
        energy_external.append(sum(e_external))
        energy_internal.append(e_internal)
    days = days[2:]

    return(days,energy_balance,energy_external,energy_internal,energy_dissipated)