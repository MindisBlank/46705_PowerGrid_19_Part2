"""
LoadNetworkData4FA.py

Builds positive-, negative-, and zero-sequence Ybus matrices
and computes the corresponding Zbus matrices (bus impedances)
for fault analysis on a given network data file.
"""

import numpy as np
import ReadNetworkData as rd


def LoadNetworkData4FA(filename):
    """
    Load network data and construct sequence Ybus/Zbus matrices.

    Parameters:
        filename (str): Path to the network data file.

    Globals:
        Ybus (np.ndarray): Positive-sequence bus admittance matrix.
        Ybus0 (np.ndarray): Zero-sequence bus admittance matrix.
        Ybus2 (np.ndarray): Negative-sequence bus admittance matrix.
        Zbus0 (np.ndarray): Zero-sequence bus impedance matrix.
        Zbus1 (np.ndarray): Positive-sequence bus impedance matrix.
        Zbus2 (np.ndarray): Negative-sequence bus impedance matrix.
        bus_to_ind (dict): Mapping from bus numbers to indices.
    """
    global Ybus, Ybus0, Ybus2, Zbus0, Zbus1, Zbus2, bus_to_ind

    # Read raw data
    (bus_data, load_data, gen_data,
     line_data, tran_data,
     mva_base, bus_to_ind, ind_to_bus) = rd.read_network_data_from_file(
        filename
    )

    num_buses = len(bus_data)

    # Initialize zeroed Ybus matrices (complex)
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    Ybus0 = np.zeros((num_buses, num_buses), dtype=complex)
    Ybus2 = np.zeros((num_buses, num_buses), dtype=complex)

    # ---------------------------------------------------------------------
    # 1) Stamp line data 
    # ---------------------------------------------------------------------
    for (fr, to, _br_id,
         r, x, _b_half,
         _mva_rate, x2, x0) in line_data:
        i = bus_to_ind[fr]
        j = bus_to_ind[to]

        z1 = complex(r, x)
        y1 = 1.0 / z1 if z1 != 0 else 0

        z2 = complex(r, x2)
        y2 = 1.0 / z2 if z2 != 0 else 0

        z0 = complex(r, x0)
        y0 = 1.0 / z0 if z0 != 0 else 0

        # Positive-sequence
        Ybus[i, i] += y1
        Ybus[j, j] += y1
        Ybus[i, j] -= y1
        Ybus[j, i] -= y1

        # Negative-sequence
        Ybus2[i, i] += y2
        Ybus2[j, j] += y2
        Ybus2[i, j] -= y2
        Ybus2[j, i] -= y2

        # Zero-sequence
        Ybus0[i, i] += y0
        Ybus0[j, j] += y0
        Ybus0[i, j] -= y0
        Ybus0[j, i] -= y0

        y_sh = 1j * _b_half       # half-shunt at each end
        Ybus0[i,i] += y_sh
        Ybus0[j,j] += y_sh
    # ---------------------------------------------------------------------
    # 2) Stamp transformer data
    # ---------------------------------------------------------------------
    for (fr, to, _br_id,
         r, x, tap, ang,
         _mva_rate, fr_co, to_co,
         x2, x0) in tran_data:
        i = bus_to_ind[fr]
        j = bus_to_ind[to]

        tap_ratio = tap * np.exp(1j * np.deg2rad(ang))

        z1 = complex(r, x)
        y1 = 1.0 / z1 if z1 != 0 else 0

        z2 = complex(r, x2)
        y2 = 1.0 / z2 if z2 != 0 else 0

        if int(fr_co) == 3 or int(to_co) == 3:
            y0 = 0
        else:
            z0 = complex(r, x0)
            y0 = 1.0 / z0 if z0 != 0 else 0

        # Positive-sequence stamping with tap
        Ybus[i, i] += y1 / (tap_ratio * np.conj(tap_ratio))
        Ybus[j, j] += y1
        Ybus[i, j] -= y1 / np.conj(tap_ratio)
        Ybus[j, i] -= y1 / tap_ratio

        # Negative-sequence stamping with tap
        Ybus2[i, i] += y2 / (tap_ratio * np.conj(tap_ratio))
        Ybus2[j, j] += y2
        Ybus2[i, j] -= y2 / np.conj(tap_ratio)
        Ybus2[j, i] -= y2 / tap_ratio

        # Zero-sequence stamping if wye-wye
        if y0:
            Ybus0[i, i] += y0 / (tap_ratio * np.conj(tap_ratio))
            Ybus0[j, j] += y0
            Ybus0[i, j] -= y0 / np.conj(tap_ratio)
            Ybus0[j, i] -= y0 / tap_ratio

    # ---------------------------------------------------------------------
    # 3) Stamp generator shunt admittances (with stator r=0.01 pu)
    # ---------------------------------------------------------------------
    for (bus_nr, mva_gen,
         _p_gen, _p_max,
         _q_max, _q_min,
         x1, x2, x0,
         _x_n, grounded) in gen_data:

        idx = bus_to_ind[bus_nr]

        # 0.01 pu stator resistance on machine base
        r_stator_base = 0.01

        # conversion factor to system base
        scale = mva_base / mva_gen

        # machine impedances on system base
        r_stator = r_stator_base * scale
        x1_sys   = x1 * scale
        x2_sys   = x2 * scale
        x0_sys   = x0 * scale

        # full sequence impedances (r + jx)
        Z1_gen = r_stator + 1j*x1_sys
        Z2_gen = r_stator + 1j*x2_sys
        Z0_gen = (r_stator + 1j*x0_sys) if grounded else None

        # stamp into each Ybus
        Ybus[idx, idx]  += 1.0 / Z1_gen
        Ybus2[idx, idx] += 1.0 / Z2_gen
        if grounded and Z0_gen is not None:
            Ybus0[idx, idx] += 1.0 / Z0_gen

    # ---------------------------------------------------------------------
    # 4) Invert Ybus to Zbus
    # ---------------------------------------------------------------------
    try:
        Zbus0 = np.linalg.inv(Ybus0)
    except np.linalg.LinAlgError:
        print("Error: Ybus0 matrix is singular; cannot invert.")
        Zbus0 = None
    
    try:
        Zbus1 = np.linalg.inv(Ybus)
    except np.linalg.LinAlgError:
        print("Error: Ybus1 matrix is singular; cannot invert.")
        Zbus1 = None
    Zth = Zbus1[bus_to_ind[4011], bus_to_ind[4011]]
    print(f"Zth: {Zth}")
    try:
        Zbus2 = np.linalg.inv(Ybus2)
    except np.linalg.LinAlgError:
        print("Error: Ybus2 matrix is singular; cannot invert.")
        Zbus2 = None

    return
