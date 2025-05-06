# -*- coding: utf-8 -*-
"""
Main Script for Fault Analysis (batch mode)
Generates results_B.txt with all fault types.
"""
import numpy as np
import io
from contextlib import redirect_stdout

import FaultAnalysis_46705 as fa
import LoadNetworkData4FA as lnd4fa
import ReadNetworkData as rd   # <-- new

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------
DisplayResults   = 1      # 0 = no console output, 1 = print each result to screen
filename         = "Nordic32_4FA.txt"
FaultBus         = 4011
FaultImpedance   = 0      # per-unit fault impedance
PrefaultVoltage  = 1.0    # per-unit prefault voltage

fault_types = {
    0: "Three-phase balanced fault",
    1: "Single Line-to-Ground fault",
    2: "Line-to-Line fault",
    3: "Double Line-to-Ground fault"
}

# ------------------------------------------------------------------------------
# READ BUS DATA & SYSTEM BASE (for Display function)
# ------------------------------------------------------------------------------
bus_data, load_data, gen_data, line_data, tran_data, \
mva_base, bus_to_ind, ind_to_bus = rd.read_network_data_from_file(filename)

# ------------------------------------------------------------------------------
# LOAD NETWORK DATA INTO Ybus/Zbus
# ------------------------------------------------------------------------------
lnd4fa.LoadNetworkData4FA(filename)

# ------------------------------------------------------------------------------
# OPEN OUTPUT FILE (UTF-8)
# ------------------------------------------------------------------------------
with open("results_B.txt", "w", encoding="utf-8") as outfile:
    outfile.write("********** Fault Analysis Batch Results **********\n\n")

    for ft_code, ft_desc in fault_types.items():
        Iph, Vph_mat = fa.FaultAnalysis(
            lnd4fa.Zbus0,
            lnd4fa.Zbus1,
            lnd4fa.Zbus2,
            lnd4fa.bus_to_ind,
            FaultBus,
            ft_code,
            FaultImpedance,
            PrefaultVoltage
        )

        # Header
        outfile.write(f"--- {ft_desc} (Type {ft_code}) at Bus {FaultBus} ---\n\n")

        # Capture display output (pu + kA)
        buf = io.StringIO()
        with redirect_stdout(buf):
            fa.DisplayFaultAnalysisResults(
                Iph,
                Vph_mat,
                FaultBus,
                ft_code,
                FaultImpedance,
                PrefaultVoltage,
                lnd4fa.bus_to_ind,
                mva_base,            # <-- new
                bus_data             # <-- new
            )
        outfile.write(buf.getvalue())
        outfile.write("\n\n")

        if DisplayResults:
            fa.DisplayFaultAnalysisResults(
                Iph,
                Vph_mat,
                FaultBus,
                ft_code,
                FaultImpedance,
                PrefaultVoltage,
                lnd4fa.bus_to_ind,
                mva_base,            # <-- new
                bus_data             # <-- new
            )

print("All fault cases processed. Results exported to results_B.txt")
