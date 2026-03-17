"""
A file to process a LAMMPS trajectory file and generate a
load directory structure for the pyretis simulation
"""

import math
import os
import shutil
import sys

# ── User settings ────────────────────────────────────────────────────────────
TRAJ_FILE     = "biased.lammpstrj"   # input LAMMPS trajectory
ATOM1_ID      = 1                    # id of first diatomic atom
ATOM2_ID      = 2                    # id of second diatomic atom
N_ENSEMBLES   = 8                    # number of ensembles (000 … 007)
LOAD_DIR      = "load"               # output root

# Interface positions from the original TIS article (as shown in Fig. 43)
# lambda_0 is the dividing surface (state A boundary)
# Adjust these to match your retis.rst settings!
INTERFACES = [1.13, 1.3, 1.4, 1.5, 1.55, 1.6, 1.65, 1.7]
# ─────────────────────────────────────────────────────────────────────────────


def parse_lammpstrj(filename):
    """
    Parse a LAMMPS custom dump file.
    Returns a list of frames, each frame being a dict:
        {
          'timestep': int,
          'atoms': {id: {'type': int, 'x': float, 'y': float, 'z': float,
                         'vx': float, 'vy': float, 'vz': float}}
        }
    """
    frames = []
    with open(filename) as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "ITEM: TIMESTEP":
            timestep = int(lines[i + 1].strip())
            # skip NUMBER OF ATOMS and BOX BOUNDS
            natoms = int(lines[i + 3].strip())
            # find ITEM: ATOMS header
            i += 4
            while not lines[i].startswith("ITEM: ATOMS"):
                i += 1
            col_names = lines[i].strip().split()[2:]  # drop "ITEM:" and "ATOMS"
            atoms = {}
            for j in range(natoms):
                vals = lines[i + 1 + j].strip().split()
                d = {col_names[k]: vals[k] for k in range(len(col_names))}
                aid = int(d["id"])
                atoms[aid] = {
                    "type": int(d["type"]),
                    "x":  float(d["x"]),
                    "y":  float(d["y"]),
                    "z":  float(d["z"]),
                    "vx": float(d.get("vx", 0.0)),
                    "vy": float(d.get("vy", 0.0)),
                    "vz": float(d.get("vz", 0.0)),
                }
            frames.append({"timestep": timestep, "atoms": atoms})
            i += 1 + natoms
        else:
            i += 1
    return frames

test_frame = parse_lammpstrj(TRAJ_FILE)
print(test_frame[0]["atoms"])
#sys.exit()


def order_parameter(frame, id1=ATOM1_ID, id2=ATOM2_ID):
    """Distance between the two diatomic atoms (no PBC correction needed for
    the order parameter because LAMMPS unwraps or the distance is small)."""
    a1 = frame["atoms"][id1]
    a2 = frame["atoms"][id2]
    dx = a1["x"] - a2["x"]
    dy = a1["y"] - a2["y"]
    dz = a1["z"] - a2["z"]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

test_dist = order_parameter(frame=test_frame[283], id1=ATOM1_ID, id2=ATOM2_ID)
print(test_dist)
sys.exit()


def find_zero_minus_trajectory(op_values, lambda0):
    """
    Find a [0-] trajectory: starts above lambda_0, dips below, comes back above.
    Returns (start_idx, end_idx) into op_values, or None if not found.
    """
    n = len(op_values)
    for start in range(n):
        if op_values[start] > lambda0:
            # look for a dip below lambda0 then return above
            went_below = False
            for mid in range(start + 1, n):
                if op_values[mid] < lambda0:
                    went_below = True
                if went_below and op_values[mid] > lambda0:
                    # valid [0-] segment found
                    return (start, mid)
    return None


def find_reactive_trajectory(op_values, lambda_first, lambda_last):
    """
    Find a reactive trajectory: starts below lambda_first, ends above lambda_last,
    with intermediate points all between lambda_first and lambda_last.
    Returns (start_idx, end_idx) or None.
    """
    n = len(op_values)
    for start in range(n):
        if op_values[start] < lambda_first:
            for end in range(start + 2, n):
                segment = op_values[start:end + 1]
                # first point below lambda_first, last above lambda_last
                if segment[-1] > lambda_last:
                    interior = segment[1:-1]
                    if all(lambda_first <= v <= lambda_last for v in interior):
                        return (start, end)
    return None


def write_load_dir(ensemble_idx, frames, frame_indices, traj_filename):
    """
    Write order.txt, traj.txt, and accepted/ for one ensemble.
    ensemble_idx : int (0 … N_ENSEMBLES-1)
    frames       : full list of parsed frames
    frame_indices: list of indices into `frames` that make up this path
    traj_filename: basename of the trajectory file placed in accepted/
    """
    ens_dir      = os.path.join(LOAD_DIR, f"{ensemble_idx:03d}")
    accepted_dir = os.path.join(ens_dir, "accepted")
    os.makedirs(accepted_dir, exist_ok=True)

    # ── order.txt ────────────────────────────────────────────────────────────
    with open(os.path.join(ens_dir, "order.txt"), "w") as fh:
        fh.write("# Order parameter values\n")
        for fi in frame_indices:
            op = order_parameter(frames[fi])
            fh.write(f"  {op:.6f}\n")

    # ── traj.txt ─────────────────────────────────────────────────────────────
    with open(os.path.join(ens_dir, "traj.txt"), "w") as fh:
        fh.write("# Step Filename index vel\n")
        for step, fi in enumerate(frame_indices):
            fh.write(f"  {step} {traj_filename} {fi} -1\n")

    # ── accepted/ ────────────────────────────────────────────────────────────
    dest = os.path.join(accepted_dir, traj_filename)
    if not os.path.exists(dest):
        if os.path.exists(TRAJ_FILE):
            shutil.copy2(TRAJ_FILE, dest)
            print(f"  Copied {TRAJ_FILE} → {dest}")
        else:
            print(f"  WARNING: {TRAJ_FILE} not found; copy it manually to {dest}")

    print(f"  Written ensemble {ensemble_idx:03d}  "
          f"(frames {frame_indices[0]}–{frame_indices[-1]}, "
          f"len={len(frame_indices)})")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print(f"Reading {TRAJ_FILE} …")
    frames = parse_lammpstrj(TRAJ_FILE)
    print(f"  {len(frames)} frames found.")

    op_values = [order_parameter(f) for f in frames]
    print(f"  OP range: {min(op_values):.4f} – {max(op_values):.4f}")
    print(f"  Interfaces: {INTERFACES}")

    lambda0    = INTERFACES[0]
    lambda_last = INTERFACES[-1]

    traj_basename = os.path.basename(TRAJ_FILE)
    os.makedirs(LOAD_DIR, exist_ok=True)

    # ── Ensemble 000 : [0-] ──────────────────────────────────────────────────
    result = find_zero_minus_trajectory(op_values, lambda0)
    if result is None:
        print("WARNING: Could not find a [0-] trajectory. "
              "Try running a longer unbiased simulation near state A.")
        zero_minus_indices = list(range(min(50, len(frames))))
    else:
        s, e = result
        zero_minus_indices = list(range(s, e + 1))
    write_load_dir(0, frames, zero_minus_indices, traj_basename)

    # ── Ensembles 001–007 : reactive trajectories ────────────────────────────
    result = find_reactive_trajectory(op_values, lambda0, lambda_last)
    if result is None:
        print("\nWARNING: No reactive trajectory found in the biased run.\n"
              "  The spring pull may need more steps, a stronger force constant,\n"
              "  or a shorter cut distance. Try increasing 'run' to 20000–50000\n"
              "  or increase the spring constant.\n"
              "  Using a fallback: the highest-OP segment of the trajectory.\n")
        # Fallback: just use the segment that reaches the furthest
        peak = op_values.index(max(op_values))
        # walk back to first point below lambda0
        start = peak
        while start > 0 and op_values[start] > lambda0:
            start -= 1
        reactive_indices = list(range(start, peak + 1))
    else:
        s, e = result
        reactive_indices = list(range(s, e + 1))
        print(f"\nReactive trajectory found: frames {s}–{e} "
              f"(OP {op_values[s]:.4f} → {op_values[e]:.4f})")

    for ens in range(1, N_ENSEMBLES):
        write_load_dir(ens, frames, reactive_indices, traj_basename)

    print(f"\nDone. Load directory written to '{LOAD_DIR}/'")
    print("Next steps:")
    print("  1. Verify interface positions in INTERFACES match your retis.rst")
    print("  2. Make sure the double-well potential is active in your lammps.in")
    print("     (remove the comment from  'include \"dw-wca.in\"')")
    print("  3. Run:  pyretisrun -i retis.rst -p")


if __name__ == "__main__":
    main()

"""
# Parameter settings-------------------------------------------------------------------
traj_file   = "/home/jonas/pyretis_venv/examples/11_LAMMPS/init_traj/biased.lammpstrj"
atom_id1    = 1
atom_id2    = 2
target_dir  = "load"
interfaces  = [1.13, 1.3, 1.4, 1.5, 1.55, 1.6, 1.65, 1.7]
n_ensembles = len(interfaces)
#--------------------------------------------------------------------------------------

def parse_trajfile(filename):
    
    Function to parse a LAMMPS trajectory file and create a dictionaries of the timesteps and
    atoms with the respective positions

    Parameters
    ----------
    filename: string
        name/path of the trajectory file

    Returns
    -------
    traj_dict: dictionary
    

    traj_dict = {}
    with open(filename, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("ITEM: TIMESTEP"):
            timestep = int(lines[i+1].strip())
            n_atoms = int(lines[i+3].strip())
            i += 4
    """