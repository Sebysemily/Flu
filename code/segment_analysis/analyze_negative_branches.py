import sys, re
import xml.etree.ElementTree as ET

def analyze_tree(path):
    with open(path, 'r') as f:
        text = f.read()

    tree_line = None
    for line in text.split('\n'):
        if line.startswith('tree '):
            tree_line = line
            break

    if not tree_line:
        print("No tree line found")
        sys.exit(1)
        
    matches = re.findall(r'(\w+)?\[&([^\]]+)\]:(-[0-9]*\.?[0-9]+(?:E-?[0-9]+)?)', tree_line)
    
    print(f"Found {len(matches)} negative branches in {path}.")
    for m in matches:
        tip = m[0] if m[0] else "Internal Node"
        annots = m[1]
        length = m[2]
        
        # extract height_median or height_mean
        h_mean = re.search(r'height_mean=([0-9.]+)', annots)
        h_median = re.search(r'height_median=([0-9.]+)', annots)
        
        h_mean_val = h_mean.group(1) if h_mean else "N/A"
        h_median_val = h_median.group(1) if h_median else "N/A"
        
        print(f"Node: {tip}, Length: {length}, Height_Mean: {h_mean_val}, Height_Median: {h_median_val}")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        analyze_tree(sys.argv[1])
    else:
        analyze_tree('results/beast/final_run/H5N1_HA_panel_postQC.mcc.tree')
