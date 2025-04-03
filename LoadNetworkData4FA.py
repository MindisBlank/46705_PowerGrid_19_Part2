import numpy as np
import ReadNetworkData as rd

def LoadNetworkData4FA(filename):
    # Make these global so that other functions can use them.
    global Ybus, Ybus0, Ybus2, Zbus0, Zbus1, Zbus2, bus_to_ind

    # Read the network data from the file using the provided library function.
    bus_data, load_data, gen_data, line_data, tran_data, mva_base, bus_to_ind, ind_to_bus = \
        rd.read_network_data_from_file(filename)

    # Number of buses in the system.
    N = len(bus_data)
    
    # Initialize the Ybus matrices (complex-valued) for:
    # Positive-sequence (Ybus), zero-sequence (Ybus0) and negative-sequence (Ybus2)
    Ybus   = np.zeros((N, N), dtype=complex)
    Ybus0  = np.zeros((N, N), dtype=complex)
    Ybus2  = np.zeros((N, N), dtype=complex)

    #############################################################
    # Process line data
    # Each line in line_data is assumed to be formatted as:
    # [fr, to, br_id, R, X, B, MVA_rate, X2, X0]
    # For fault analysis we ignore shunt admittances (B) and use only series impedances.
    #############################################################
    for branch in line_data:
        fr, to, br_id, R, X, B, MVA_rate, X2, X0 = branch
        i = bus_to_ind[fr]
        j = bus_to_ind[to]
        
        # Positive-sequence impedance (Z1)
        Z1 = complex(R, X)
        Y1 = 1.0 / Z1 if Z1 != 0 else 0
        
        # Negative-sequence impedance (Z2)
        Z2 = complex(R, X2)
        Y2 = 1.0 / Z2 if Z2 != 0 else 0
        
        # Zero-sequence impedance (Z0)
        Z0 = complex(R, X0)
        Y0 = 1.0 / Z0 if Z0 != 0 else 0
        
        # Update Ybus matrices using standard Ybus construction
        Ybus[i, i]  += Y1
        Ybus[j, j]  += Y1
        Ybus[i, j]  -= Y1
        Ybus[j, i]  -= Y1
        
        Ybus2[i, i] += Y2
        Ybus2[j, j] += Y2
        Ybus2[i, j] -= Y2
        Ybus2[j, i] -= Y2
        
        Ybus0[i, i] += Y0
        Ybus0[j, j] += Y0
        Ybus0[i, j] -= Y0
        Ybus0[j, i] -= Y0

    #############################################################
    # Process transformer data
    # Each transformer in tran_data is assumed to be formatted as:
    # [fr, to, br_id, R, X, n, ang1, MVA_rate, fr_co, to_co, X2, X0]
    # where n is the tap ratio and ang1 is the phase shift (in degrees).
    # The connection codes (fr_co, to_co) are assumed to be:
    #    1 or 2  --> wye (or wye-grounded),
    #    3       --> delta.
    # For positive and negative sequence, we include the branch (with tap adjustments).
    # For zero sequence, if either winding is delta (code == 3), then we do not include it.
    #############################################################
    for tra in tran_data:
        fr, to, br_id, R, X, n, ang1, MVA_rate, fr_co, to_co, X2, X0 = tra
        i = bus_to_ind[fr]
        j = bus_to_ind[to]

        # Compute the complex tap ratio:
        t = n * np.exp(1j * ang1 * np.pi / 180)
        
        # Positive-sequence branch impedance and admittance.
        Z1 = complex(R, X)
        Y1 = 1.0 / Z1 if Z1 != 0 else 0
        
        # Negative-sequence branch impedance and admittance.
        Z2 = complex(R, X2)
        Y2 = 1.0 / Z2 if Z2 != 0 else 0
        
        # Zero-sequence: check for delta connection.
        if int(fr_co) == 3 or int(to_co) == 3:
            Y0 = 0
        else:
            Z0 = complex(R, X0)
            Y0 = 1.0 / Z0 if Z0 != 0 else 0

        # Update positive-sequence Ybus.
        Ybus[i, i]  += Y1 / (t * np.conjugate(t))
        Ybus[j, j]  += Y1
        Ybus[i, j]  -= Y1 / np.conjugate(t)
        Ybus[j, i]  -= Y1 / t
        
        # Update negative-sequence Ybus2.
        Ybus2[i, i] += Y2 / (t * np.conjugate(t))
        Ybus2[j, j] += Y2
        Ybus2[i, j] -= Y2 / np.conjugate(t)
        Ybus2[j, i] -= Y2 / t
        
        # Update zero-sequence Ybus0 if allowed.
        if Y0 != 0:
            Ybus0[i, i] += Y0 / (t * np.conjugate(t))
            Ybus0[j, j] += Y0
            Ybus0[i, j] -= Y0 / np.conjugate(t)
            Ybus0[j, i] -= Y0 / t

    #############################################################
    # Process generator data: Add generator impedances as shunt admittances.
    # Each generator in gen_data is formatted as:
    # [bus_nr, mva_size, p_gen, p_max, q_max, q_min, X1, X2, X0, Xn, grnd]
    #############################################################
    for gen in gen_data:
        bus_nr, mva_size, p_gen, p_max, q_max, q_min, X1, X2, X0, Xn, grnd = gen
        i = bus_to_ind[bus_nr]
        # Positive-sequence: add shunt admittance 1/(j*X1)
        if X1 != 0:
            Ybus[i, i] += 1.0 / (1j * X1)
        # Negative-sequence: add shunt admittance 1/(j*X2)
        if X2 != 0:
            Ybus2[i, i] += 1.0 / (1j * X2)
        # Zero-sequence: add shunt admittance 1/(j*X0) only if generator is grounded.
        if grnd and X0 != 0:
            Ybus0[i, i] += 1.0 / (1j * X0)

    #############################################################
    # Finally, compute the bus impedance matrices by inverting (using the pseudo-inverse)
    # the Ybus matrices.
    #############################################################
    Zbus0 = np.linalg.inv(Ybus0)
    Zbus1 = np.linalg.inv(Ybus)
    Zbus2 = np.linalg.inv(Ybus2)

    return
    # Alternatively, you can return the matrices if needed:
    # return Zbus0, Zbus1, Zbus2, Ybus, Ybus0, Ybus2, bus_to_ind, ind_to_bus
