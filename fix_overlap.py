import sys
import re

def truncate_log(path, max_state):
    lines = []
    with open(path, 'r') as f:
        for line in f:
            m = re.match(r'^([0-9]+)\t', line)
            if m and int(m.group(1)) >= max_state:
                break
            lines.append(line)
    with open(path, 'w') as f:
        f.writelines(lines)

def truncate_trees(path, max_state):
    lines = []
    with open(path, 'r') as f:
        for line in f:
            m = re.match(r'^tree STATE_([0-9]+)', line.strip())
            if m and int(m.group(1)) >= max_state:
                break
            lines.append(line)
    lines.append('End;\n')
    with open(path, 'w') as f:
        f.writelines(lines)

truncate_log('results/beast/GSS/ucln_exp/H5N1_HA_panel_postQC.log.prior', 75000000)
truncate_trees('results/beast/GSS/ucln_exp/H5N1_HA_panel_postQC.trees.prior', 75000000)
print("Truncated successfully")
