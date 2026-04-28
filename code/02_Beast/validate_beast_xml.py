#!/usr/bin/env python3
import argparse
import os
import xml.etree.ElementTree as ET


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate that a BEAST XML is well formed.")
    parser.add_argument("--xml", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ET.parse(args.xml)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(f"xml\t{os.path.abspath(args.xml)}\n")
        handle.write("status\tvalidated\n")


if __name__ == "__main__":
    main()
