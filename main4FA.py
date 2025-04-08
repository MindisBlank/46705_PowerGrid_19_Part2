# -*- coding: utf-8 -*-
"""
Main Script for Fault Analysis
"""
import numpy as np  # Ensure NumPy is imported
import FaultAnalysis_46705 as fa    # import Fault Analysis functions
import LoadNetworkData4FA as lnd4fa # load the network data to global variables
filename = "Nordic32_SA.txt" # TestSystem4FA.txt  Nordic32_SA.txt
lnd4fa.LoadNetworkData4FA(filename) # makes Zbus0 available as lndfa.Zbus0 etc.

print('**********Start of Fault Analysis**********')

p=0 # Set p=1 to display Ybus and Zbus matrices, set p=0 to carry out fault analysis.
if p==1:
    print(f'Zbus0 = \n {np.array2string(np.round(lnd4fa.Zbus0, 2), separator=", ")}')
    print(f'Zbus1 = \n {np.array2string(np.round(lnd4fa.Zbus1, 2), separator=", ")}')
    print(f'Zbus2 = \n {np.array2string(np.round(lnd4fa.Zbus2, 2), separator=", ")}')
    print(f'Ybus0 = \n {np.array2string(np.round(lnd4fa.Ybus0, 2), separator=", ")}')
    print(f'Ybus1 = \n {np.array2string(np.round(lnd4fa.Ybus, 2), separator=", ")}')
    print(f'Ybus2 = \n {np.array2string(np.round(lnd4fa.Ybus2, 2), separator=", ")}')
    print('Bus to Index = \n', lnd4fa.bus_to_ind)
else :
    # Carry out the fault analysis ... 
    FaultBus = 4011
    # FaultType: 0 = 3-phase balanced fault; 1 = Single Line-to-Ground fault;
    #            2 = Line-to-Line fault;     3 = Double Line-to-Ground fault.
    FaultType = 3
    FaultImpedance = 0 # (in pu) 
    PrefaultVoltage = 1.000 # (in pu)
    # Iph: phase current array (0: phase a; 1: phase b; 2: phase c). 
    # Vph_mat: phase line-to-ground voltages (rows: busses; columns: phases a, b, c).
    Iph,Vph_mat = fa.FaultAnalysis(lnd4fa.Zbus0,lnd4fa.Zbus1,lnd4fa.Zbus2,lnd4fa.bus_to_ind, 
                                FaultBus,FaultType,FaultImpedance,PrefaultVoltage)
    # Display results
    fa.DisplayFaultAnalysisResults(Iph,Vph_mat,FaultBus,FaultType,FaultImpedance,PrefaultVoltage,lnd4fa.bus_to_ind)
    print('**********End of Fault Analysis**********')