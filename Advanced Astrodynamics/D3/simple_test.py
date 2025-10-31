with open('debug_output.txt', 'w') as f:
    f.write("Starting test\n")
    from utils import read_glonass_rinex
    f.write("Import successful\n")
    data = read_glonass_rinex('brdc2320.09g.txt')
    f.write(f'Available PRNs: {sorted(data.keys())}\n')
    f.write("Test complete\n")