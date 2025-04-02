# -*- coding: utf-8 -*-
"""
Main Script for Fault Analysis
"""

import FaultAnalysis_46705 as fa    # import Fault Analysis functions
import LoadNetworkData4FA as lnd4fa # load the network data to global variables
filename = "TestSystem4FA.txt" # TestSystem4FA.txt  Nordic32_SA.txt
lnd4fa.LoadNetworkData4FA(filename) # makes Zbus0 available as lndfa.Zbus0 etc.
print('**********Start of Fault Analysis**********')

# print('Zbus0 = \n', lnd4fa.Zbus0)
# print('Zbus1 = \n', lnd4fa.Zbus1)
# print('Zbus2 = \n', lnd4fa.Zbus2)
# print('Ybus0 = \n', lnd4fa.Ybus0)
# print('Ybus1 = \n', lnd4fa.Ybus)
# print('Ybus2 = \n', lnd4fa.Ybus2)
#print('Bus to Index = \n', lnd4fa.bus_to_ind)


# Carry out the fault analysis ... 
FaultBus = 3
# FaultType: 0 = 3-phase balanced fault; 1 = Single Line-to-Ground fault;
#            2 = Line-to-Line fault;     3 = Double Line-to-Ground fault.
FaultType = 1
FaultImpedance = 0 # (in pu) 
PrefaultVoltage = 1.000 # (in pu)
# Iph: phase current array (0: phase a; 1: phase b; 2: phase c). 
# Vph_mat: phase line-to-ground voltages (rows: busses; columns: phases a, b, c).
Iph,Vph_mat = fa.FaultAnalysis(lnd4fa.Zbus0,lnd4fa.Zbus1,lnd4fa.Zbus2,lnd4fa.bus_to_ind, 
                               FaultBus,FaultType,FaultImpedance,PrefaultVoltage)
# Display results
fa.DisplayFaultAnalysisResults(Iph,Vph_mat,FaultBus,FaultType,FaultImpedance,PrefaultVoltage)
print('**********End of Fault Analysis**********')