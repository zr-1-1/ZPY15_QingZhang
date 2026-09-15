"""Extract saved ATK segment states without changing the scenario (SI units)."""
from pathlib import Path
import csv
import json
import sys
import xml.etree.ElementTree as ET


def extract(path):
    root = ET.fromstring(Path(path).read_text(encoding="gb18030"))
    rows = []
    for satellite in root.iter("Satellite"):
        name = satellite.get("Name", "")
        if name not in ("Mother1", "Mother2", "Mother3"):
            continue
        for segment in satellite.find("Planning").iter("Segment"):
            kind = segment.get("ComponentType")
            if kind not in ("CMCSInitialState", "CMCSPropagate", "CMCSManeuver"):
                continue
            row = {"mother": int(name[-1]), "kind": kind}
            for tag in ("InitialState", "FinalState"):
                state = segment.find(tag)
                row[tag] = dict(
                    utc=state.findtext("UTC"),
                    frame=state.find("CoordSystem").get("Name"),
                    rv=[float(state.find(v).get(axis)) for v in ("Pos", "Vel") for axis in "XYZ"],
                )
            burn = segment.find("Burn")
            if burn is not None:
                row["burn"] = [float(burn.get(axis)) for axis in "XYZ"]
                row["burn_config"] = segment.find("OrbBurn").attrib
            rows.append(row)
    return rows


if __name__ == "__main__":
    out = Path("data/atk_correction")
    out.mkdir(exist_ok=True)
    rows = extract(sys.argv[1] if len(sys.argv)>1 else "ZPY15.xml")
    target = Path(sys.argv[2]) if len(sys.argv)>2 else out / "original_segments.json"
    target.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    if len(sys.argv)>1:
        sys.exit(0)
    sets = []
    reports = []
    for mother in range(1, 4):
        lines = Path(f'ATK_result/Mother{mother}.txt').read_text(encoding='gb18030').splitlines()[1:]
        records = [line.split() for line in lines if line.strip()]
        eligible = sorted({int(r[1].removeprefix('Debris')) for r in records if float(r[19]) < .15})
        sets.append(set(eligible))
        reports.append(dict(Mother=mother, EncounterRows=len(records),
                            UniqueEncounterDebris=len({r[1] for r in records}),
                            EligibleDebris=eligible))
    summary = dict(Reports=reports, Union12=sorted(sets[0] | sets[1]), Union123=sorted(set.union(*sets)))
    (out / 'original_report_counts.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    for row in rows:
        if row["kind"] == "CMCSManeuver":
            delta = [row["FinalState"]["rv"][i+3]-row["InitialState"]["rv"][i+3] for i in range(3)]
            print(row["mother"], row["InitialState"]["utc"], row["burn_config"],
                  "actual-minus-J2000:", [delta[i]-row["burn"][i] for i in range(3)])
