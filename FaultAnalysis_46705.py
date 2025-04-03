"""
46705 - Power Grid Analysis
This file contains the definitions of the functions needed to
carry out Fault Analysis calculations in python.
"""

import numpy as np

# 1. the FaultAnalysis() function
def FaultAnalysis(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,fault_type,Zf,Vf):
    ''' String here with purpose '''
    # calculate sequence fault currents
    Iseq = Calculate_Sequence_Fault_Currents(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,fault_type,Zf,Vf)
    # calculate sequence fault voltages
    Vseq_mat = Calculate_Sequence_Fault_Voltages(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,Vf,Iseq)
    # convert sequence currents to phase (fault) currents
    Iph = Convert_Sequence2Phase_Currents(Iseq)
    # convert sequence voltages to phase line-to-ground (fault) voltages
    Vph_mat = Convert_Sequence2Phase_Voltages(Vseq_mat)    
    return Iph, Vph_mat

# 1.1. the Calculate_Sequence_Fault_Currents() function
def Calculate_Sequence_Fault_Currents(Zbus0, Zbus1, Zbus2, bus_to_ind, fault_bus, fault_type, Zf, Vf):
    """
    Calculate the sequence fault currents for a fault at a given bus.
    
    Parameters:
        Zbus0 : np.array
            Zero-sequence bus impedance matrix.
        Zbus1 : np.array
            Positive-sequence bus impedance matrix.
        Zbus2 : np.array
            Negative-sequence bus impedance matrix.
        bus_to_ind : dict
            Mapping from bus numbers to their indices in the matrices.
        fault_bus : int
            Bus number at which the fault occurs.
        fault_type : int
            Fault type coded as:
              0 = three-phase fault,
              1 = single line-to-ground fault,
              2 = line-to-line fault,
              3 = double line-to-ground fault.
        Zf : complex
            Fault impedance (applied in the appropriate sequence network(s)).
        Vf : complex
            Prefault bus voltage (typically 1.0 pu unless specified otherwise).
            
    Returns:
        Iseq : list of complex
            List of sequence fault currents [I0, I1, I2],
            where I0 is the zero–sequence, I1 is the positive–sequence,
            and I2 is the negative–sequence fault current.
    """
    # Get the matrix index corresponding to the faulted bus.
    idx = bus_to_ind[fault_bus]
    
    # Extract the Thevenin equivalent impedances at the faulted bus.
    Z0_f = Zbus0[idx, idx]
    Z1_f = Zbus1[idx, idx]
    Z2_f = Zbus2[idx, idx]
    
    # Initialize sequence currents.
    Iseq = [0, 0, 0]  # [I0, I1, I2]
    
    if fault_type == 0:
        # Three-phase fault:
        # Only the positive-sequence network is active.
        I1 = Vf / (Z1_f + Zf)
        Iseq[0] = 0
        Iseq[1] = I1
        Iseq[2] = 0
        
    elif fault_type == 1:
        # Single line-to-ground fault:
        # All three sequence networks are connected in series.
        I_fault = Vf / (Z0_f + Z1_f + Z2_f + 3 * Zf)
        Iseq[0] = I_fault
        Iseq[1] = I_fault
        Iseq[2] = I_fault
        
    elif fault_type == 2:
        # Line-to-line fault:
        # The zero-sequence network is not involved (I0 = 0).
        # Typically, I1 = Vf / (Z1_f + Z2_f + Zf) and I2 = -I1.
        I1 = Vf / (Z1_f + Z2_f + Zf)
        Iseq[0] = 0
        Iseq[1] = I1
        Iseq[2] = -I1
        
    elif fault_type == 3:
        # Double line-to-ground fault:
        # For a DLG fault, the sequence networks are connected as follows:
        #     Z1 is in series with a parallel combination of (Z0+3Zf) and (Z2+3Zf).
        # The equivalent impedance is:
        #     Z_eq = Z1_f + [ (Z0_f+3Zf)*(Z2_f+3Zf) / (Z0_f+Z2_f+3Zf) ]
        # Then,
        #     I1 = Vf / Z_eq
        # and
        #     I0 = ((Z2_f+3Zf)/(Z0_f+Z2_f+3Zf)) * I1,
        #     I2 = ((Z0_f+3Zf)/(Z0_f+Z2_f+3Zf)) * I1.
        Z_parallel = ((Z0_f + 3 * Zf) * (Z2_f + 3 * Zf)) / (Z0_f + Z2_f + 3 * Zf)
        Z_eq = Z1_f + Z_parallel
        I1 = Vf / Z_eq
        I0 = ((Z2_f + 3 * Zf) / (Z0_f + Z2_f + 3 * Zf)) * I1
        I2 = ((Z0_f + 3 * Zf) / (Z0_f + Z2_f + 3 * Zf)) * I1
        Iseq[0] = I0
        Iseq[1] = I1
        Iseq[2] = I2
        
    else:
        raise ValueError("Unsupported fault type: {}. Supported types are 0, 1, 2, 3.".format(fault_type))
    
    return Iseq





# 1.2 the Calculate_Sequence_Fault_Voltages() function
def Calculate_Sequence_Fault_Voltages(Zbus0, Zbus1, Zbus2, bus_to_ind, fault_bus, Vf, Iseq):
    """
    Calculate the sequence fault voltages at all buses.
    
    Parameters:
        Zbus0 : np.array
            Zero-sequence bus impedance matrix (NxN).
        Zbus1 : np.array
            Positive-sequence bus impedance matrix (NxN).
        Zbus2 : np.array
            Negative-sequence bus impedance matrix (NxN).
        bus_to_ind : dict
            Mapping from bus numbers to their indices.
        fault_bus : int
            Bus number where the fault occurs.
        Vf : complex
            Prefault voltage (applied only to the positive-sequence network).
        Iseq : list of complex
            List of sequence fault currents [I0, I1, I2] (calculated at the faulted bus).
    
    Returns:
        Vseq_mat : np.array
            An Nx3 matrix of sequence voltages. For each bus i:
                Vseq_mat[i][0] = zero-sequence voltage,
                Vseq_mat[i][1] = positive-sequence voltage,
                Vseq_mat[i][2] = negative-sequence voltage.
    """
    # Number of buses
    N = Zbus0.shape[0]
    
    # Initialize the matrix (rows: buses, columns: sequence networks)
    Vseq_mat = np.zeros((N, 3), dtype=complex)
    
    # Get the index of the faulted bus.
    fault_idx = bus_to_ind[fault_bus]
    
    # For each bus, calculate the sequence voltages.
    # Note: It is typical to assume that pre-fault zero and negative sequence voltages are zero.
    for i in range(N):
        # Zero-sequence voltage (prefault = 0)
        Vseq_mat[i, 0] = 0 - Zbus0[i, fault_idx] * Iseq[0]
        # Positive-sequence voltage (prefault = Vf)
        Vseq_mat[i, 1] = Vf - Zbus1[i, fault_idx] * Iseq[1]
        # Negative-sequence voltage (prefault = 0)
        Vseq_mat[i, 2] = 0 - Zbus2[i, fault_idx] * Iseq[2]
    
    return Vseq_mat

# 1.3. the Convert_Sequence2Phase_Currents() function
def Convert_Sequence2Phase_Currents(Iseq):
    """
    Convert the sequence current array Iseq to the phase current array Iph.
    
    Parameters:
        Iseq : list or array-like
            A list [I0, I1, I2] containing the zero-, positive-, and negative-
            sequence currents.
    
    Returns:
        Iph : list
            A list [Ia, Ib, Ic] containing the phase currents for phases a, b, and c.
    """
    # Define the operator 'a' (120° rotation)
    a = np.exp(1j * 2 * np.pi / 3)  # a = exp(j*120°)
    a2 = a ** 2                    # a^2 = exp(j*240°)
    
    # Unpack the sequence currents
    I0, I1, I2 = Iseq[0], Iseq[1], Iseq[2]
    
    # Apply the transformation to obtain phase currents.
    Ia = I0 + I1 + I2
    Ib = I0 + a2 * I1 + a * I2
    Ic = I0 + a * I1 + a2 * I2
    
    # Create the phase current array.
    Iph = [Ia, Ib, Ic]
    return Iph


# 1.4 the Convert_Sequence2Phase_Voltages() function
def Convert_Sequence2Phase_Voltages(Vseq_mat):
    """
    Convert the sequence voltage matrix Vseq_mat to the phase line-to-ground voltage matrix Vph_mat.
    
    Parameters:
        Vseq_mat : np.array
            An N x 3 matrix where each row contains the sequence voltages [V0, V1, V2] for a bus.
    
    Returns:
        Vph_mat : np.array
            An N x 3 matrix where each row contains the phase voltages [Va, Vb, Vc] for a bus.
            The ordering is such that Vph_mat[i][0] = Va at bus i+1,
            Vph_mat[i][1] = Vb, and Vph_mat[i][2] = Vc.
    """
    # Define the operator a = exp(j*120°) and a^2 = exp(j*240°)
    a = np.exp(1j * 2 * np.pi / 3)
    a2 = a ** 2

    # Number of buses (rows in Vseq_mat)
    N = Vseq_mat.shape[0]
    
    # Initialize the phase voltage matrix with zeros.
    Vph_mat = np.zeros((N, 3), dtype=complex)
    
    # Loop over each bus and perform the transformation.
    for i in range(N):
        V0 = Vseq_mat[i, 0]
        V1 = Vseq_mat[i, 1]
        V2 = Vseq_mat[i, 2]
        
        # Compute phase voltages using inverse symmetrical component transformation.
        Va = V0 + V1 + V2
        Vb = V0 + a2 * V1 + a * V2
        Vc = V0 + a * V1 + a2 * V2
        
        # Store the computed phase voltages.
        Vph_mat[i, :] = [Va, Vb, Vc]
        
    return Vph_mat


# ####################################################
# #  Displaying the results in the terminal window   #
# ####################################################
# 2. the DisplayFaultAnalysisResults() function
def DisplayFaultAnalysisResults(Iph, Vph_mat, fault_bus, fault_type, Zf, Vf,bus_to_ind):
    """
    Displays phase fault currents (Iph) at the faulted bus and the phase line-to-ground
    voltages (Vph_mat) at all buses.
    
    Parameters:
        Iph : list or array-like
            A list [I_a, I_b, I_c] containing the phase currents at the faulted bus.
        Vph_mat : np.array
            An N x 3 matrix of phase line-to-ground voltages (each row corresponds to a bus,
            with columns representing phases a, b, and c).
        fault_bus : int
            Bus number where the fault occurs.
        fault_type : int
            Fault type code (0: Three-phase fault, 1: Single line-to-ground fault,
            2: Line-to-line fault, 3: Double line-to-ground fault).
        Zf : complex
            Fault impedance (in per unit).
        Vf : complex
            Prefault voltage (in per unit).
    """
    # Dictionary to translate fault type code to a descriptive string.
    fault_dict = {
        0: "Three-phase fault",
        1: "Single Line-to-Ground fault",
        2: "Line-to-Line fault",
        3: "Double Line-to-Ground fault"
    }
    fault_str = fault_dict.get(fault_type, "Unknown Fault Type")
    
    print("==============================================================")
    print("|                  Fault Analysis Results                    |")
    print("==============================================================")
    print(f" {fault_str} at Bus {fault_bus}")
    print(f" Prefault Voltage: Vf = {Vf:.3f} (pu)")
    print(f" Fault Impedance: Zf = {Zf:.3f} (pu)")
    print("")
    
    # Helper function: returns a string with magnitude and angle (in degrees) using numpy.
    def cpx_mag_ang_str(c):
        mag = np.abs(c)
        ang = np.angle(c, deg=True)  # angle in degrees using numpy
        return f"{mag:.3f} ∠ {ang:6.2f}°"
    
    # --- Display Phase Currents at the Faulted Bus ---
    print("--------------------------------------------------------------")
    print("                   Phase Currents at Fault                   ")
    print("--------------------------------------------------------------")
    Ia_str = cpx_mag_ang_str(Iph[0])
    Ib_str = cpx_mag_ang_str(Iph[1])
    Ic_str = cpx_mag_ang_str(Iph[2])
    
    print(f"  Phase a: {Ia_str}")
    print(f"  Phase b: {Ib_str}")
    print(f"  Phase c: {Ic_str}")
    print("")
    
    # --- Display Phase Voltages at All Buses ---
    print("--------------------------------------------------------------")
    print("           Phase Line-to-Ground Voltages (all buses)         ")
    print("--------------------------------------------------------------")
    print(f"{'Bus':>4s}  {'Phase a':>20s}  {'Phase b':>20s}  {'Phase c':>20s}")
    print("--------------------------------------------------------------")
    
    num_buses = Vph_mat.shape[0]
    for i in range(num_buses):
        bus_number = list(bus_to_ind.keys())[list(bus_to_ind.values()).index(i)]
        Va_str = cpx_mag_ang_str(Vph_mat[i, 0])
        Vb_str = cpx_mag_ang_str(Vph_mat[i, 1])
        Vc_str = cpx_mag_ang_str(Vph_mat[i, 2])
        print(f"{bus_number:4d}  {Va_str:>20s}  {Vb_str:>20s}  {Vc_str:>20s}")
    
    print("==============================================================")
