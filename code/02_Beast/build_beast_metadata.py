#!/usr/bin/env python3
import argparse
import csv
import os
import sys
import time
import requests

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build BEAST metadata TSV from H5N1_context.csv with latitude and longitude"
    )
    parser.add_argument("--metadata", required=True, help="Path to metadata/H5N1_context.csv")
    parser.add_argument("--out", required=True, help="Output metadata TSV")
    parser.add_argument("--aln", required=False, help="Path to FASTA alignment to filter taxa")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    valid_taxa = None
    if args.aln:
        valid_taxa = set()
        with open(args.aln, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    valid_taxa.add(line[1:].strip())

    seen = set()
    rows = []

    # Cache for coordinates to avoid redundant API calls
    coord_cache = {}
    
    BASE_URL = "https://nominatim.openstreetmap.org/search"
    HEADERS = {"User-Agent": "FluCentroidApp/1.0 (sebascalvas@outlook.com)"}

    print("Reading metadata and querying OpenStreetMap Nominatim for coordinates...", file=sys.stderr)

    with open(args.metadata, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for index, row in enumerate(reader):
            file_name = (row.get("file_name") or "").strip()
            date_val = (row.get("collection_date") or "").strip()
            country = (row.get("country") or "").strip()
            expected_role = (row.get("expected_role") or "").strip()
            host_group = (row.get("host_type") or "?").strip()

            if file_name and file_name not in seen:
                if valid_taxa is not None and file_name not in valid_taxa:
                    continue

                if country == "FalklandIslands":
                    country = "Falkland Islands"

                query_string = country
                resolution = "Country"
                
                latitude = ""
                longitude = ""

                if query_string:
                    if query_string in coord_cache:
                        latitude, longitude, resolution = coord_cache[query_string]
                    else:
                        print(f"Querying -> {query_string} (Resolution: {resolution})", file=sys.stderr)
                        params = {"q": query_string, "format": "json", "limit": 1}
                        try:
                            response = requests.get(BASE_URL, params=params, headers=HEADERS)
                            if response.status_code == 200:
                                data = response.json()
                                if data:
                                    latitude = data[0]['lat']
                                    longitude = data[0]['lon']
                                    coord_cache[query_string] = (latitude, longitude, resolution)
                                else:
                                    print(f"   ⚠️ No coordinates found for: {query_string}", file=sys.stderr)
                                    coord_cache[query_string] = ("", "", resolution)
                            else:
                                print(f"   ❌ API Error (Status {response.status_code})", file=sys.stderr)
                                coord_cache[query_string] = ("", "", resolution)
                        except Exception as e:
                            print(f"   ❌ Request failed: {e}", file=sys.stderr)
                            coord_cache[query_string] = ("", "", resolution)
                        
                        # Nominatim's usage policy requires a 1-second sleep between requests
                        time.sleep(1)

                location = country
                if country.lower() == "ecuador":
                    if expected_role == "flu_costa":
                        location = "ecuador_coastal"
                    elif expected_role == "flu_sierra":
                        location = "ecuador_andine"
                    elif expected_role == "flu_amazonia":
                        location = "ecuador_amazonia"
                    else:
                        location = "Ecuador"
                
                if not location:
                    location = "?"

                rows.append((file_name, date_val, latitude, longitude, location, host_group))
                seen.add(file_name)

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("Taxon\tDate\tLatitude\tLongitude\tLocation\tHost\n")
        for file_name, date_val, lat, lon, loc, host in sorted(rows):
            handle.write(f"{file_name}\t{date_val}\t{lat}\t{lon}\t{loc}\t{host}\n")

    print(f"Wrote {len(rows)} entries to {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
