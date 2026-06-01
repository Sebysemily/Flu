import sys
import os

def main():
    out_file = sys.argv[1]
    in_files = sys.argv[2:]
    
    with open(out_file, "w") as out:
        for f in in_files:
            seg = os.path.basename(f).split("_")[-1].replace(".fasta", "")
            with open(f, "r") as inf:
                for line in inf:
                    if line.startswith(">"):
                        out.write(line.strip() + f"_{seg}\n")
                    else:
                        out.write(line)

if __name__ == "__main__":
    main()
