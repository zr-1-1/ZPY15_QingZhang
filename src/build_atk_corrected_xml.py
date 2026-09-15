"""Build a separate corrected scenario, preserving all unrelated XML bytes."""
from pathlib import Path
import csv
import re


def build():
    original = Path('ZPY15.xml').read_bytes()
    text = original.decode('gb18030')
    assert text.encode('gb18030') == original
    with Path('data/atk_correction/maneuvers.csv').open(encoding='utf-8-sig') as f:
        pulses = list(csv.DictReader(f))
    for mother in range(1, 4):
        pattern = rf'<Satellite Name="Mother{mother}"[^>]*>.*?</Satellite>'
        matches = list(re.finditer(pattern, text, re.S))
        assert len(matches) == 1
        match = matches[0]; block = match.group()
        rows = [p for p in pulses if int(p['Mother']) == mother]
        assert len(rows) == 6
        idx = 0
        def replace_burn(m):
            nonlocal idx
            p = rows[idx]; idx += 1
            return '<Burn ' + ' '.join(f'{axis}="{float(p[field]):.17g}"' for axis, field in
                zip('XYZ', ('DeltaVx_mps','DeltaVy_mps','DeltaVz_mps'))) + ' />'
        block = re.sub(r'<Burn X="[^"]*" Y="[^"]*" Z="[^"]*"\s*/>', replace_burn, block)
        assert idx == 6
        # Stored segment endpoints are from the failed run; require a fresh run.
        block = re.sub(r'<LastRunCode>[^<]*</LastRunCode>', '<LastRunCode>-1</LastRunCode>', block)
        text = text[:match.start()] + block + text[match.end():]
    Path('ZPY15_corrected.xml').write_bytes(text.encode('gb18030'))


if __name__ == '__main__':
    build()
