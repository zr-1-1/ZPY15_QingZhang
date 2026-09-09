"""Export Debris satellites' initial Cartesian states from an ATK scenario.

Uses only the Python standard library. Values retain the XML's precision and
units: position in m, velocity in m/s, StepSize in s, and timestamps in UTC.
"""

import argparse
import csv
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parent.parent
TIME_FIELDS = ("StartUTC", "StopUTC", "OrbEpoch")
NUMERIC_FIELDS = (
    "StepSize", "PositionX", "PositionY", "PositionZ",
    "VelocityX", "VelocityY", "VelocityZ", "GravityModel",
    "MaxDegree", "MaxOrder", "UseDrag", "UseFluxGeoFile", "DragCoefficient",
)
FIELDS = ("DebrisID", "Name") + TIME_FIELDS + NUMERIC_FIELDS


def read_debris(xml_path, encoding="gb2312"):
    """Read direct Orbit fields only (exclude nested third-body settings)."""
    # Decode explicitly: ElementTree's byte parser does not support GB2312.
    root = ET.fromstring(Path(xml_path).read_text(encoding=encoding))
    rows = []
    seen = set()
    for satellite in root.findall("Satellite"):
        name = satellite.get("Name", "")
        match = re.fullmatch(r"Debris([0-9]+)", name)
        if not match:
            continue
        debris_id = int(match.group(1))
        if debris_id in seen:
            raise ValueError(f"Duplicate debris ID: {debris_id}")
        seen.add(debris_id)
        orbit = satellite.find("Orbit")
        if orbit is None:
            raise ValueError(f"{name}: missing Orbit")
        row = {"DebrisID": debris_id, "Name": name}
        for field in TIME_FIELDS + NUMERIC_FIELDS:
            value = orbit.findtext(field, "").strip()
            if not value:
                raise ValueError(f"{name}: missing Orbit/{field}")
            if field in NUMERIC_FIELDS:
                try:
                    valid = math.isfinite(float(value))
                except ValueError:
                    valid = False
                if not valid:
                    raise ValueError(f"{name}: invalid Orbit/{field}: {value!r}")
            row[field] = value
        rows.append(row)
    if not rows:
        raise ValueError("No Satellite elements named Debris<number> found")
    return sorted(rows, key=lambda row: row["DebrisID"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xml", nargs="?", type=Path,
                        default=PROJECT_ROOT / "ZPY15.xml")
    parser.add_argument("-o", "--output", type=Path,
                        default=PROJECT_ROOT / "data/debris_orbits.csv")
    parser.add_argument("--encoding", default="gb2312",
                        help="XML encoding (default: gb2312, as in ZPY15.xml)")
    args = parser.parse_args()
    try:
        if args.xml.resolve() == args.output.resolve():
            raise ValueError("Output must differ from input XML")
        # Validate every row before opening the output file.
        rows = read_debris(args.xml, args.encoding)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
    except (OSError, ValueError, LookupError, ET.ParseError) as exc:
        parser.exit(1, f"Export failed: {exc}\n")
    print(f"Exported {len(rows)} debris satellites to {args.output}")


if __name__ == "__main__":
    main()
