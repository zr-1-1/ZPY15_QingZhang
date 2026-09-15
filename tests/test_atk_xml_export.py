"""The correction must alter only 18 burns and invalidate mother run caches."""
from pathlib import Path
import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from audit_atk_xml import extract

original = Path('ZPY15.xml').read_bytes().decode('gb18030')
corrected = Path('ZPY15_corrected.xml').read_bytes().decode('gb18030')
def neutralize(text):
    for ship in range(1, 4):
        pattern = rf'<Satellite Name="Mother{ship}"[^>]*>.*?</Satellite>'
        m = re.search(pattern, text, re.S)
        block = re.sub(r'<Burn X="[^"]*" Y="[^"]*" Z="[^"]*"\s*/>', '<Burn />', m.group())
        block = re.sub(r'<LastRunCode>[^<]*</LastRunCode>', '<LastRunCode />', block)
        text = text[:m.start()] + block + text[m.end():]
    return text
assert neutralize(original) == neutralize(corrected), 'Unrelated scenario data changed'
old, new = extract('ZPY15.xml'), extract('ZPY15_corrected.xml')
assert len(old) == len(new) == 42
for a, b in zip(old, new):
    assert a['InitialState'] == b['InitialState'] and a['FinalState'] == b['FinalState']
    if a['kind'] == 'CMCSManeuver':
        assert a['burn_config'] == b['burn_config'] and a['burn'] != b['burn']
print('ATK XML export: only 18 burns and mother cache status changed; other bytes preserved.')
