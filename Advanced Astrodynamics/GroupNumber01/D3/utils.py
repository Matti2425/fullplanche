import re
import numpy as np
from datetime import datetime

def read_glonass_rinex(filepath):
    """
    Parse a RINEX GLONASS Navigation file (v2.x or v3.x) into structured data.
    Handles missing spaces and compact scientific notation.
    """
    eph_data = {}

    with open(filepath, "r") as f:
        content = f.read()

    # Split into lines and skip header
    lines = content.splitlines()
    start_idx = 0
    for i, line in enumerate(lines):
        if "END OF HEADER" in line:
            start_idx = i + 1
            break

    # Pattern for start of GLONASS block (PRN + date/time)
    header_pattern = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)"
    )

    # Pattern for floats including missing spaces
    float_pattern = re.compile(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?")

    i = start_idx
    while i < len(lines):
        line = lines[i]
        m = header_pattern.match(line)
        if not m:
            i += 1
            continue

        try:
            prn = int(m.group(1))
            year = int(m.group(2))
            month = int(m.group(3))
            day = int(m.group(4))
            hour = int(m.group(5))
            minute = int(m.group(6))
            sec = float(m.group(7))

            if year < 80:
                year += 2000
            else:
                year += 1900

            epoch = datetime(year, month, day, hour, minute, int(sec))

            # Get next 3 lines (position, velocity, acceleration)
            Xline = lines[i + 1] if i + 1 < len(lines) else ""
            Yline = lines[i + 2] if i + 2 < len(lines) else ""
            Zline = lines[i + 3] if i + 3 < len(lines) else ""

            def parse_line(line):
                vals = [float(v.replace("D", "E")) for v in float_pattern.findall(line)]
                return vals if len(vals) >= 3 else [0.0, 0.0, 0.0]

            X = parse_line(Xline)
            Y = parse_line(Yline)
            Z = parse_line(Zline)

            eph_data.setdefault(prn, []).append({
                "epoch": epoch,
                "pos": np.array([X[0], Y[0], Z[0]]),
                "vel": np.array([X[1], Y[1], Z[1]]),
                "acc": np.array([X[2], Y[2], Z[2]])
            })

        except Exception as e:
            print(f"⚠️ Parse issue at line {i}: {e}")

        i += 4  # Move to next satellite record

    return eph_data




def get_ephemeris(eph_data, prn, date_time):
    """
    Retrieve the ephemeris entry for a given PRN and UTC datetime.
    Finds the record with the closest matching epoch.
    """
    if prn not in eph_data:
        raise ValueError(f"PRN {prn} not found in data.")

    recs = eph_data[prn]
    # Find nearest epoch
    rec = min(recs, key=lambda e: abs((e["epoch"] - date_time).total_seconds()))
    return rec
