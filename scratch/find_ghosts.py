import os
import sys

def main():
    summary_file = "scratch/snakemake_summary.tsv"
    if not os.path.exists(summary_file):
        print(f"Summary file not found: {summary_file}")
        return

    tracked_files = set()
    with open(summary_file, 'r') as f:
        # Skip header
        next(f, None)
        for line in f:
            parts = line.split('\t')
            if parts:
                path = parts[0].strip()
                tracked_files.add(os.path.abspath(path))

    # We also want to consider parent directories of tracked files as "tracked directories"
    tracked_dirs = set()
    for tf in tracked_files:
        d = tf
        while d != '/' and d != '':
            tracked_dirs.add(d)
            d = os.path.dirname(d)

    search_roots = ['data/phylogeny', 'data/processed_alignments', 'data/beast', 'results']
    
    ghost_dirs = []
    ghost_files = []

    # Known auxiliary extensions that we should ignore if they are untracked
    # but sit next to tracked files.
    ignore_exts = {
        '.iqtree', '.log', '.bionj', '.mldist', '.contree', '.ckp.gz', '.splits.nex',
        '.state', '.trees', '.ops', '.part_', '.xml.state', '.nex', '.model.gz'
    }

    for root_dir in search_roots:
        if not os.path.exists(root_dir):
            continue
        for dirpath, dirnames, filenames in os.walk(root_dir):
            abs_dirpath = os.path.abspath(dirpath)
            
            # If the directory itself is not tracked (meaning no tracked files inside it or its subdirs)
            if abs_dirpath not in tracked_dirs:
                # This whole directory is a ghost!
                ghost_dirs.append(dirpath)
                # Clear dirnames so we don't recurse into it and list its subdirs too
                dirnames[:] = []
                continue
            
            # The directory is tracked. Are there any ghost files inside it?
            for fname in filenames:
                fpath = os.path.join(dirpath, fname)
                abs_fpath = os.path.abspath(fpath)
                
                if abs_fpath not in tracked_files:
                    # It's an untracked file. Is it an expected auxiliary file?
                    is_aux = False
                    for ext in ignore_exts:
                        if fname.endswith(ext) or '.part_' in fname:
                            is_aux = True
                            break
                    if not is_aux:
                        ghost_files.append(fpath)

    print("=== GHOST DIRECTORIES ===")
    for gd in sorted(ghost_dirs):
        print(gd)
        
    print("\n=== GHOST FILES ===")
    for gf in sorted(ghost_files):
        print(gf)

if __name__ == '__main__':
    main()
