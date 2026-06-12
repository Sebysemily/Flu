import argparse
import os
import csv
import re

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    parser.add_argument('logs', nargs='+')
    args = parser.parse_args()

    mle_dict = {}
    for log_path in args.logs:
        scenario = os.path.basename(os.path.dirname(log_path))
        with open(log_path, 'r') as f:
            lines = f.readlines()
            for line in reversed(lines):
                if line.startswith("log marginal likelihood"):
                    match = re.search(r'=\s*([-\d\.]+)', line)
                    if match:
                        mle_dict[scenario] = float(match.group(1))
                    break

    if not mle_dict:
        print("No MLE values found.")
        return

    max_mle = max(mle_dict.values())
    max_model = max(mle_dict, key=mle_dict.get)
    
    diff_dict = {model: 2 * (max_mle - val) for model, val in mle_dict.items()}

    models = list(mle_dict.keys())
    
    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Metric"] + models)
        
        row_mle = ["MLE"] + [mle_dict[m] for m in models]
        writer.writerow(row_mle)
        
        row_diff = [f"2 * (L({max_model}) - L)"] + [diff_dict[m] for m in models]
        writer.writerow(row_diff)

if __name__ == "__main__":
    main()
