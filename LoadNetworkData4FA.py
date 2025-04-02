import numpy as np
import ReadNetworkData as rd


#TODO  Make sure that you dont need the generator sequneces in the matrices

def LoadNetworkData4FA(filename):
    # Make these global so that other functions can use them.
    global Ybus, Ybus0, Ybus2, Zbus0, Zbus1, Zbus2,bus_to_ind

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
        # For off-diagonals, the contribution is -Y, and for diagonals, add Y.
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
        # t = n * exp(j*ang1*pi/180)
        t = n * np.exp(1j * ang1 * np.pi / 180)
        
        # Positive-sequence branch impedance and admittance.
        Z1 = complex(R, X)
        Y1 = 1.0 / Z1 if Z1 != 0 else 0
        
        # Negative-sequence branch impedance and admittance.
        Z2 = complex(R, X2)
        Y2 = 1.0 / Z2 if Z2 != 0 else 0
        
        # Zero-sequence: check for delta connection.
        # If either side is delta (connection code 3), then the transformer branch is open-circuited
        # in the zero-sequence network.
        if int(fr_co) == 3 or int(to_co) == 3:
            Y0 = 0
        else:
            Z0 = complex(R, X0)
            Y0 = 1.0 / Z0 if Z0 != 0 else 0

        # In Ybus construction for transformers with an off-nominal tap,
        # we assume the tap is on the "from" side. The standard formulation is:
        #   Y11 += Y / |t|^2,  Y22 += Y, Y12 = -Y / t*, Y21 = -Y / t.
        # This is applied for both positive and negative sequence networks.
        # For the zero-sequence network, only add if Y0 is nonzero.
        
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
        
        # Update zero-sequence Ybus0 only if transformer winding connections allow it.
        if Y0 != 0:
            Ybus0[i, i] += Y0 / (t * np.conjugate(t))
            Ybus0[j, j] += Y0
            Ybus0[i, j] -= Y0 / np.conjugate(t)
            Ybus0[j, i] -= Y0 / t

    #############################################################
    # Finally, compute the bus impedance matrices by inverting the Ybus matrices.
    # (Note: in practical applications, care must be taken when the Ybus matrix is singular.)
    #############################################################
    Zbus0 = np.linalg.inv(Ybus0)
    Zbus1 = np.linalg.inv(Ybus)
    Zbus2 = np.linalg.inv(Ybus2)

    # (Optionally, you might return these matrices or simply use them as globals.)
    return
    #return Zbus0, Zbus1, Zbus2, Ybus, Ybus0, Ybus2, bus_to_ind, ind_to_bus