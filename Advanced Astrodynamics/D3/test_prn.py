import re

# Test the regex pattern
header_pattern = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)")

test_line = " 2 09  8 20  0 15  0.0-0.197580084205E-04-0.272848410532E-11 0.300000000000E+02"

m = header_pattern.match(test_line)
if m:
    print('Match found:')
    for i in range(1, len(m.groups()) + 1):
        print(f'Group {i}: "{m.group(i)}"')
    try:
        sec = float(m.group(7))
        print(f'Seconds: {sec}')
    except Exception as e:
        print(f'Error converting seconds: {e}')
else:
    print('No match')

# Now test the actual parsing
import sys
sys.path.insert(0, '.')
from utils import read_glonass_rinex

try:
    data = read_glonass_rinex('brdc2320.09g.txt')
    print('Available PRNs:', sorted(data.keys()))
    if 2 in data:
        print('PRN 2 data:', len(data[2]), 'records')
        print('First record epoch:', data[2][0]['epoch'])
    else:
        print('PRN 2 not found')
except Exception as e:
    print('Error:', e)
    import traceback
    traceback.print_exc()