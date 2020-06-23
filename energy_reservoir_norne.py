import rips
import time
import grpc
import math
from operator import itemgetter
from ecl.eclfile import EclFile
from ecl.grid import EclGrid

# Connect to ResInsight
resinsight = rips.Instance.find()
start = time.time()
case = resinsight.project.cases()[0]
num_tsteps = len(case.time_steps())
name = case.name
grids = case.grids()
for grid in grids:
    dimension = grid.dimensions()
Nx = dimension.i
Ny = dimension.j
Nz = dimension.k

# Read output files
egrid_file = EclFile("%s.EGRID" % name)
init_file = EclFile("%s.INIT" % name)
egrid = EclGrid("%s.EGRID" % name)

# Read ACTNUM numbers from EGRID file
actnum = egrid_file["ACTNUM"][0]
active_cells = []
for i in range(len(actnum)):
    if actnum[i] == 1:
        active_cells.append(i+1)
num_active_cells = len(active_cells)

# Read NNC pairs from EGRID and INIT file
nnc_list = []
if egrid_file.has_kw("NNC1"):
    nnc1 = egrid_file["NNC1"][0]
    nnc2 = egrid_file["NNC2"][0]
    tran = init_file["TRANNNC"][0]
    for g1,g2,t in zip(nnc1,nnc2,tran):
        nnc_list.append((g1,g2,t))

# Calculation of energy dissipation for each cell
total_ed = 0.0
for tstep in range(219,num_tsteps):
    print("Timestep", tstep, "of", num_tsteps)
    
    # Read results into list
    tranx_results = case.active_cell_property('STATIC_NATIVE', 'TRANX', 0)
    trany_results = case.active_cell_property('STATIC_NATIVE', 'TRANY', 0)
    tranz_results = case.active_cell_property('STATIC_NATIVE', 'TRANZ', 0)
    multx_results = case.active_cell_property('STATIC_NATIVE', 'MULTX', 0)
    multy_results = case.active_cell_property('STATIC_NATIVE', 'MULTY', 0)
    multz_results = case.active_cell_property('STATIC_NATIVE', 'MULTZ', 0)
    pres_results = case.active_cell_property('DYNAMIC_NATIVE', 'PRESSURE', tstep)
    krw_results = case.active_cell_property('DYNAMIC_NATIVE', 'WATKR', tstep)
    kro_results = case.active_cell_property('DYNAMIC_NATIVE', 'OILKR', tstep)
    krg_results = case.active_cell_property('DYNAMIC_NATIVE', 'GASKR', tstep)
    muw_results = case.active_cell_property('DYNAMIC_NATIVE', 'WAT_VISC', tstep)
    muo_results = case.active_cell_property('DYNAMIC_NATIVE', 'OIL_VISC', tstep)
    mug_results = case.active_cell_property('DYNAMIC_NATIVE', 'GAS_VISC', tstep)

    # Determination of upstream/downstream cells in each direction
    prevx_pres_results = [0.0] * num_active_cells
    nextx_pres_results = [0.0] * num_active_cells
    prevy_pres_results = [0.0] * num_active_cells
    nexty_pres_results = [0.0] * num_active_cells
    prevz_pres_results = [0.0] * num_active_cells
    nextz_pres_results = [0.0] * num_active_cells

    for c in range(num_active_cells):
        print("pres calculation: checking cell", c+1, "of", num_active_cells, end="\r")
        idx = active_cells[c]
        plane_num = idx%(Nx*Ny)
        layer_num = math.floor(idx/(Nx*Ny))
        if idx-1 in active_cells and idx%Nx != 0:
            prevx_pres_results[c] = pres_results[active_cells.index(idx-1)]
        
        if idx+1 in active_cells and (idx+1)%Nx != 0:
            nextx_pres_results[c] = pres_results[active_cells.index(idx+1)]

        if idx-Nx in active_cells:
            prevy_pres_results[c] = pres_results[active_cells.index(idx-Nx)]
            if plane_num < Nx:
                prevy_pres_results[c] = 0.0
        
        if idx+Nx in active_cells:
            nexty_pres_results[c] = pres_results[active_cells.index(idx+Nx)]
            if plane_num+Nx >= (Nx*Ny):
                nexty_pres_results[c] = 0.0
        
        if idx-(Nx*Ny) in active_cells:
            prevz_pres_results[c] = pres_results[active_cells.index(idx-Nx*Ny)]
            if layer_num == 0:
                prevz_pres_results[c] = 0.0
        
        if idx+(Nx*Ny) in active_cells:
            nextz_pres_results[c] = pres_results[active_cells.index(idx+Nx*Ny)]
            if layer_num == Nz-1:
                nextz_pres_results[c] = 0.0
    print(end='\r')
    
    # Calculate interblock flowrates in + direction(s)
    flow_x_results = []
    flow_y_results = []
    flow_z_results = []
    zip_results_f = list(zip(tranx_results, trany_results, tranz_results, multx_results, multy_results, multz_results, pres_results, prevx_pres_results, nextx_pres_results, prevy_pres_results, nexty_pres_results, prevz_pres_results, nextz_pres_results, krw_results, kro_results, krg_results, muw_results, muo_results, mug_results))
    for (tranx, trany, tranz, multx, multy, multz, pres, prevx_pres, nextx_pres, prevy_pres, nexty_pres, prevz_pres, nextz_pres, krw, kro, krg, muw, muo, mug) in zip_results_f:
        if nextx_pres > 0.0 and tranx > 0.0:
            flow_x = tranx * multx * ((krw/muw)+(kro/muo)+(krg/mug)) * (pres-nextx_pres)
        else:
            flow_x = 0.0
        
        if nexty_pres > 0.0 and trany > 0.0:
            flow_y = trany * multy * ((krw/muw)+(kro/muo)+(krg/mug)) * (pres-nexty_pres)
        else:
            flow_y = 0.0
        
        if nextz_pres > 0.0 and tranz > 0.0:
            flow_z = tranz * multz * ((krw/muw)+(kro/muo)+(krg/mug)) * (pres-nextz_pres)
        else:
            flow_z = 0.0
        
        flow_x_results.append(flow_x)
        flow_y_results.append(flow_y)
        flow_z_results.append(flow_z)

    # Determination of downstream flowrates in each direction
    prevx_flow_results = [0.0] * num_active_cells
    prevy_flow_results = [0.0] * num_active_cells
    prevz_flow_results = [0.0] * num_active_cells
    for c in range(num_active_cells):
        print("flow calculation: checking cell", c+1, "of", num_active_cells, end="\r")
        idx = active_cells[c]
        plane_num = idx%(Nx*Ny)
        layer_num = math.floor(idx/(Nx*Ny))
        if idx-1 in active_cells and idx%Nx != 0:
            prevx_flow_results[c] = flow_x_results[active_cells.index(idx-1)]

        if idx-Nx in active_cells:
            prevy_flow_results[c] = flow_y_results[active_cells.index(idx-Nx)]
            if plane_num < Nx:
                prevy_flow_results[c] = 0.0
        
        if idx-(Nx*Ny) in active_cells:
            prevz_flow_results[c] = flow_z_results[active_cells.index(idx-Nx*Ny)]
            if layer_num == 0:
                prevz_flow_results[c] = 0.0
    print(end='\r')

    # Generate energy change lists as results
    e_dis = []
    e_int = []
    # 1: Calculate energy changes in reservoir
    zip_results_e = list(zip(pres_results, prevx_pres_results, nextx_pres_results, prevy_pres_results, nexty_pres_results, prevz_pres_results, nextz_pres_results, flow_x_results, flow_y_results, flow_z_results, prevx_flow_results, prevy_flow_results, prevz_flow_results))
    for (pres, prevx_pres, nextx_pres, prevy_pres, nexty_pres, prevz_pres, nextz_pres, flow_x, flow_y, flow_z, prevx_flow, prevy_flow, prevz_flow) in zip_results_e:
        ed_x = flow_x * (pres - nextx_pres)
        ed_y = flow_y * (pres - nexty_pres)
        ed_z = flow_z * (pres - nextz_pres)

        ei_x = pres * (prevx_flow - flow_x)
        ei_y = pres * (prevy_flow - flow_y)
        ei_z = pres * (prevz_flow - flow_z)

        e_dis.append((ed_x + ed_y + ed_z) * (1e5/86400))     # Convert unit to J/s
        e_int.append((ei_x + ei_y + ei_z) * (1e5/86400))

    #2: Calculate energy changes from flow between NNC pairs    
    for n in range(len(nnc_list)):
        print("Checking NNC", n+1, "of", len(nnc_list), end="\r")
        orig = active_cells.index(nnc_list[n][0])
        dest = active_cells.index(nnc_list[n][1])
        param1 = zip_results_f[orig]
        param2 = zip_results_f[dest]
        tran = nnc_list[n][2]
        flow_nnc = tran * ((param1[13]/param1[16])+(param1[14]/param1[17])+(param1[15]/param1[18])) * (param1[6]-param2[6])
        ed_nnc = flow_nnc * (param1[6]-param2[6]) * (1e5/86400)
        ei_nnc_1 = -param1[6] * flow_nnc * (1e5/86400)
        ei_nnc_2 = param2[6] * flow_nnc * (1e5/86400)
        
        e_dis[orig] += ed_nnc
        e_int[orig] += ei_nnc_1
        e_int[dest] += ei_nnc_2
    print("Total energy dissipation:",sum(e_dis),"J/s")
    print("Total internal energy change:",sum(e_int),"J/s")

    try:
      # Send back output result
        case.set_active_cell_property(e_dis, 'GENERATED', 'Energy Dissipation', tstep)
        case.set_active_cell_property(e_int, 'GENERATED', 'Internal Energy Change', tstep)
    except grpc.RpcError as e:
        print("Exception Received: ", e)

end = time.time()
print("Time elapsed: ", end - start)
print("Transferred all results back\n")

view = case.views()[0].apply_cell_result('GENERATED', 'Energy Dissipation')