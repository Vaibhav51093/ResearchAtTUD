#!/usr/bin/python
"""
Vibrational Density of States (VDOS) Analysis

This script processes LAMMPS MD trajectories with velocity information to calculate VDOS. 

The expected dump file format is: ID TYPE x y z vx vy vz. 

It utilizes numpy.correlate and numpy.fft for efficiency. Caution is advised for large datasets 
due to memory usage. For mass-weighted VACF, a mapping file of atomic IDs to masses is required.

"""

__author__ = "Stefan Bringuier"
__affiliation__ = "University of Arizona, Dept. of MSE"
__creation_date__ = "Feb. 12, 2014"
__license__ = "MIT"
__modification_history__ = [
    (
        "Sept. 16, 2019",
        "Added id:mass map reading, mass-weighted VACF, migrated to Python3.",
    ),
    (
        "Jan. 4, 2024",
        "Added robust commandline options, refactor functions and plotting.",
    ),
]
__units__ = {
    "Time": "seconds",
    "Frequency": "THz",
}
__pkg_version_tested__ = {
    "numpy": "1.26.2",
    "scipy": "1.11.4",
    "matplotlib": "3.8.2",
}

import argparse
import glob
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process LAMMPS dump files to calculate VDOS.",
        epilog="NOTES: Program outputs raw files VACF.dat and VDOS.dat as well as figures VDOS.png and VACF.png",
        usage="Example, dump2VDOS.py --dump_file argon.dump --timestep 10.0e-15 --correlation_length 150",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dump_file", type=str, help="Path to the LAMMPS dump file.")
    group.add_argument(
        "--dump_pattern", type=str, help="Pattern for multiple LAMMPS dump files."
    )

    parser.add_argument(
        "--num_columns",
        type=int,
        default=8,
        help="Number of columns in ITEM: ATOMS. Note that vx, vy, and vz must be columns 6-8 (5-7 py)",
    )
    parser.add_argument(
        "--timestep",
        type=float,
        default=1.0e-15,
        help="Timestep in seconds (e.g., 1 fs = 1e-15 s ) between dump snapshots. Note it is the timestep times dump output freq.",
    )
    parser.add_argument(
        "--mass_map_file",
        type=str,
        default=None,
        help="Path to the file mapping atom IDs to masses. Useful for mass weighted VDOS.",
    )
    parser.add_argument(
        "--num_atoms",
        type=int,
        default=None,
        help="Number of atoms. If not provided, it is read from the dump file.",
    )
    parser.add_argument(
        "--num_snapshots",
        type=int,
        default=None,
        help="Number of snapshots. If not provided, it is inferred from the dump file.",
    )
    parser.add_argument(
        "--correlation_length",
        type=int,
        default=500,
        help="Correlation length. This is the number of snapshots to use.",
    )

    # New arguments for windowing function and factors
    parser.add_argument(
        "--window_function",
        type=str,
        default="gaussian",
        choices=["gaussian", "hanning", "hamming", "blackman", "bartlett"],
        help="Windowing function to use for VDOS calculation.",
    )
    parser.add_argument(
        "--std_factor",
        type=float,
        default=15.5,
        help="Standard deviation factor for Gaussian windowing.",
    )
    parser.add_argument(
        "--pad_factor",
        type=int,
        default=15,
        help="Padding factor for FFT calculation in VDOS.",
    )
    parser.add_argument(
        "--vacf_image_file",
        type=str,
        default="VACF.png",
        help="VACF result visualized",
    )
    parser.add_argument(
        "--vdos_image_file",
        type=str,
        default="VDOS.png",
        help="VDOS result visualized",
    )
    parser.add_argument(
        "--vacf_data_file",
        type=str,
        default="VACF.dat",
        help="VACF result data",
    )
    parser.add_argument(
        "--vdos_data_file",
        type=str,
        default="VDOS.dat",
        help="VDOS result data",
    )
    args = parser.parse_args()

    # Check for conditional requirements
    if args.dump_pattern:
        if args.num_atoms is None or args.num_snapshots is None:
            parser.error(
                "--num_atoms and --num_snapshots are required when --dump_pattern is provided"
            )

    return args


def read_lammps_dump(filename, num_atoms, num_fields, num_snapshots):
    """
    Read data from a LAMMPS dump file.

    Parameters:
    - filename (str): Path to the LAMMPS dump file.
    - num_atoms (int): Number of atoms to read from each snapshot.
    - num_fields (int): Number of data fields for each atom in the dump file.
    - num_snapshots (int): Number of snapshots present in the dump file.

    Returns:
    - data (ndarray): A numpy array containing the data from the dump file.
                      The array has a shape of (num_atoms, num_fields, num_snapshots).
    """
    data = np.ndarray((num_atoms, num_fields, num_snapshots), dtype=float)

    with open(filename, "r") as file:
        for snapshot_index in range(num_snapshots):
            # Skip header lines
            for _ in range(9):  # Assuming there are 8 lines in the header
                file.readline()

            # Read data for each atom in the snapshot
            for atom_index in range(num_atoms):
                line = file.readline().strip().split()
                data[atom_index, :, snapshot_index] = [float(field) for field in line]

    return data


def infer_from_dump(file_name):
    """
    Infer the number of atoms, timestep difference, and total snapshots from the LAMMPS dump file.

    Parameters:
    - file_name: str, path to the LAMMPS dump file.

    Returns:
    - num_atoms: int, inferred number of atoms.
    - timestep_diff: int, difference in timesteps between two consecutive snapshots.
    - total_snapshots: int, inferred total number of snapshots.
    """
    num_atoms = 0
    timestep_diff = 0
    total_snapshots = 0
    first_timestep = 0
    second_timestep = 0
    last_timestep = 0

    with open(file_name, "r") as file:
        lines = file.readlines()
        snapshot_count = 0

        for i, line in enumerate(lines):
            if "TIMESTEP" in line:
                snapshot_count += 1
                timestep = int(lines[i + 1].strip())

                if snapshot_count == 1:
                    first_timestep = timestep
                elif snapshot_count == 2:
                    second_timestep = timestep
                last_timestep = timestep

            if "NUMBER OF ATOMS" in line and num_atoms == 0:
                num_atoms = int(lines[i + 1].strip())

    timestep_diff = second_timestep - first_timestep
    total_snapshots = (last_timestep - first_timestep) // timestep_diff + 1

    return num_atoms, total_snapshots


def parsemapfile(Filename, Ntypes):
    """Parse the file containing atom IDs and associated masses. Dump file
    atoms ID need to stay the same, i.e. LAMMP dump_modify sort id"""
    File = open(Filename, "r")
    massmap = dict()

    header = File.readline()
    for n in range(Ntypes):
        line = File.readline()
        try:
            # Remove comments
            id, mass = line.strip("#").split()
            massmap[int(id)] = float(mass)
        except:
            print("Error parsing id:mass file")
            File.close()

    File.close()
    return massmap


def autocorr(X):
    """the convolution is actually being done here
    meaning from -inf to inf so we only want half the
    array"""

    result = np.correlate(X, X, mode="full")
    return result[result.size // 2 :]


def fft_autocorr(AutoCorr, dt):
    """FFT of autocorrelation function"""
    # fft_arry = fftpack.dct(AutoCorr)*dt
    fft_arry = np.fft.rfft(AutoCorr) * dt
    return fft_arry


def calculate_vacf_ensemble_average(
    file_pattern,
    correlation_length,
    num_atoms,
    num_fields,
    total_snapshots,
    mass_map=None,
):
    """
    Calculate the ensemble average of the velocity auto-correlation function (VACF).

    Parameters:
    - file_pattern (str): Either the file name or the pattern or a series of dump files.
    - correlation_length (int): The correlation length. The last 3 fields must be vx, vy, vz.
    - num_atoms (int): Number of atoms in the dump file.
    - num_fields (int): Number of fields (columns) in each record of the dump file.
    - total_snapshots (int): Total number of output configurations in the dump file.
    - mass_map (dict, optional): A dictionary mapping atom IDs to their masses. If provided,
      the mass-weighted VACF is calculated.

    Returns:
    - times (ndarray): Array of time indices.
    - ensemble_avg_vacf (ndarray): Ensemble-averaged VACF.
    """
    ensemble_avg_vacf = np.zeros(correlation_length, dtype=float)
    vacf_components = {
        "vx": np.zeros(correlation_length, dtype=float),
        "vy": np.zeros(correlation_length, dtype=float),
        "vz": np.zeros(correlation_length, dtype=float),
    }

    # Process each file matching the pattern
    for filename in glob.glob(file_pattern):
        data = read_lammps_dump(filename, num_atoms, num_fields, total_snapshots)
        print(f"Processing file: {filename}")

        for component in vacf_components:
            vacf_components[component].fill(0.0)

        blocks = range(correlation_length, total_snapshots, correlation_length)
        for time_index in blocks:
            mass_weighted = mass_map is not None
            mass_sum = (
                sum(mass_map.get(data[i, 1, time_index], 1.0) for i in range(num_atoms))
                if mass_weighted
                else 1.0
            )

            for i in range(num_atoms):
                mass_i = (
                    mass_map.get(data[i, 1, time_index], 1.0) if mass_weighted else 1.0
                )
                for component in ["vx", "vy", "vz"]:
                    velocity_index = 5 + ["vx", "vy", "vz"].index(component)
                    vacf_components[component] += autocorr(
                        mass_i
                        * data[
                            i,
                            velocity_index,
                            time_index - correlation_length : time_index,
                        ]
                    )

        # Averaging over components and time blocks
        for component in vacf_components:
            ensemble_avg_vacf += vacf_components[component]
        ensemble_avg_vacf /= 3 * len(blocks)

    # Average over all files and normalize
    ensemble_avg_vacf /= len(glob.glob("dump.*"))
    ensemble_avg_vacf /= ensemble_avg_vacf[0] if ensemble_avg_vacf[0] != 0 else 1

    times = np.arange(correlation_length)
    return times, ensemble_avg_vacf


def write_vacf(time, vacf, filename="VACF.dat"):
    f = open(filename, "w")
    f.write("# Time(ps) Norm. Int.(No Units) \n")
    for i in range(time.size):
        f.write("%f %f \n" % (time[i], vacf[i]))
    f.close()
    return None


def plot_vacf(time, vacf, image_file="VACF.png"):
    """
    Plot the Velocity Auto-Correlation Function (VACF).

    Parameters:
    - time: ndarray, array of time values.
    - vacf: ndarray, array of VACF values.
    """
    plt.plot(time, vacf, "k")
    plt.xlim((0, None))
    plt.ylabel(r"$\frac{v(t)\cdot v(t)}{v(0)\cdot v(0)}$", fontsize=20)
    plt.xlabel("Time [ps]")
    plt.savefig(image_file)
    plt.close()


def write_vdos(freq, vdos, filename="VDOS.dat"):
    """
    Write the VDOS data (frequency and FFT of VACF) to a file.

    Parameters:
    - freq: ndarray, array of frequency values.
    - fft_v: ndarray, FFT of the VACF.
    - filename: str (optional), name of the file to write the data to.
    """
    with open(filename, "w") as file:
        file.write("# Frequency (THz), VDOS (THz^-1)\n")
        for f, fft_val in zip(freq, vdos):
            file.write(f"{f:.5f}, {fft_val.real:.5e}\n")


def get_vdos(
    vacf,
    corlen,
    dt,
    image_file="VDOS.png",
    data_file="VDOS.dat",
    window="gaussian",
    std_factor=15.5,
    pad_factor=15,
    freq_max=15.0,
):
    """
    Vibrational Density of States (VDOS) using FFT of the autocorrelation data.

     Parameters:
     - vacf: ndarray, Velocity Auto-Correlation Function data.
     - corlen: int, Correlation length.
     - dt: float, Timestep duration in femtoseconds.
    """
    # Zero-padding and windowing
    std = corlen / std_factor
    vsize = vacf.size
    np.lib.pad(vacf, (0, vsize * pad_factor), "constant", constant_values=(0))
    windowed_vacf_v = vacf * signal.get_window((window, std), vacf.size)

    # FFT
    freq = np.fft.rfftfreq(windowed_vacf_v.size, d=dt) / 1.0e12  # Convert to THz
    fft_vacf = fft_autocorr(windowed_vacf_v, dt) * 1.0e12  # Convert to THz^-1

    write_vdos(freq, fft_vacf, filename=data_file)

    # Plotting
    plt.plot(freq, np.abs(fft_vacf))
    plt.xlim((0, freq_max))
    plt.ylim((0, None))
    xx, locs = plt.xticks()
    ll = ["%2.0f" % a for a in xx]
    plt.xticks(xx, ll)
    plt.xlabel("THz")
    plt.ylabel(r"VDOS [THz $^{-1}$]")
    plt.savefig(image_file)
    plt.close()


def main():
    args = parse_arguments()

    # Extract variables from args
    if args.dump_file is not None:
        fil = args.dump_file
    else:
        fil = args.dump_pattern

    ncols = args.num_columns
    dt = args.timestep
    mass_map_file = args.mass_map_file
    num_atoms = args.num_atoms
    num_snapshots = args.num_snapshots
    corlen = args.correlation_length
    window_function = args.window_function
    std_factor = args.std_factor
    pad_factor = args.pad_factor
    vacf_image_file = args.vacf_image_file
    vdos_image_file = args.vdos_image_file
    vacf_data_file = args.vacf_data_file
    vdos_data_file = args.vdos_data_file

    if num_atoms is None or num_snapshots is None:
        if args.dump_file:
            num_atoms, num_snapshots = infer_from_dump(fil)

    x, vacf = calculate_vacf_ensemble_average(
        fil, corlen, num_atoms, ncols, num_snapshots, mass_map=mass_map_file
    )

    time_ps = x * dt * 1.0e12

    write_vacf(time_ps, vacf, filename=vacf_data_file)

    plot_vacf(time_ps, vacf, image_file=vacf_image_file)

    get_vdos(
        vacf,
        corlen,
        dt,
        window=window_function,
        std_factor=std_factor,
        pad_factor=pad_factor,
        image_file=vdos_image_file,
        data_file=vdos_data_file,
    )

    return None


if __name__ == "__main__":
    main()

