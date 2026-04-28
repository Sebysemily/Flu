#!/usr/bin/env python3
import argparse
import csv
import os
import xml.etree.ElementTree as ET


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare a BEAST run XML from a versioned base template.")
    parser.add_argument("--template-xml", required=True)
    parser.add_argument("--scenario-name", required=True)
    parser.add_argument("--output-xml", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--panel-taxa", default=None)
    parser.add_argument("--chain-length", type=int, required=True)
    parser.add_argument("--log-every", type=int, required=True)
    parser.add_argument("--tree-every", type=int, required=True)
    parser.add_argument("--echo-every", type=int, required=True)
    parser.add_argument("--checkpoint-every", type=int, default=None)
    return parser.parse_args()


def require_attr(node: ET.Element | None, description: str) -> ET.Element:
    if node is None:
        raise ValueError(f"No se encontro {description} en el XML template.")
    return node


def find_with_id(root: ET.Element, tag: str, expected_id: str) -> ET.Element:
    for node in root.iter(tag):
        if node.get("id") == expected_id:
            return node
    raise ValueError(f"No se encontro <{tag} id=\"{expected_id}\"> en el XML template.")


def prefixed_log_name(output_prefix: str, original_name: str) -> str:
    suffix = os.path.basename(original_name)
    if suffix.startswith("H5N1_HA."):
        suffix = suffix[len("H5N1_HA.") :]
    return f"{output_prefix}.{suffix}"


def read_panel_taxa(path: str) -> set[str]:
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {
            str(row.get("taxon", "")).strip()
            for row in reader
            if str(row.get("taxon", "")).strip()
        }


def validate_template_taxa(root: ET.Element, panel_taxa_path: str) -> None:
    panel_taxa = read_panel_taxa(panel_taxa_path)
    xml_taxa = {
        node.get("id")
        for node in root.iter("taxon")
        if node.get("id")
    }

    if panel_taxa != xml_taxa:
        only_panel = sorted(panel_taxa - xml_taxa)
        only_xml = sorted(xml_taxa - panel_taxa)
        raise ValueError(
            "El XML template no coincide con el panel final. "
            f"Solo en panel: {only_panel[:5]}; solo en XML: {only_xml[:5]}"
        )

    for alignment in root.iter("alignment"):
        if not alignment.get("id"):
            continue
        sequence_count = sum(1 for child in alignment if child.tag == "sequence")
        if sequence_count == 0:
            continue
        if sequence_count != len(panel_taxa):
            raise ValueError(
                f"El alignment {alignment.get('id')} tiene {sequence_count} secuencias, "
                f"pero el panel final tiene {len(panel_taxa)} taxa."
            )


def main() -> None:
    args = parse_args()

    tree = ET.parse(args.template_xml)
    root = tree.getroot()

    if args.panel_taxa:
        validate_template_taxa(root, args.panel_taxa)

    mcmc = require_attr(root.find("mcmc"), "<mcmc>")
    screen_log = find_with_id(root, "log", "screenLog")
    file_log = find_with_id(root, "log", "fileLog")
    tree_log = find_with_id(root, "logTree", "treeFileLog")
    checkpoint_log = find_with_id(root, "logCheckpoint", "checkpointFileLog")

    mcmc.set("chainLength", str(args.chain_length))
    mcmc.set("operatorAnalysis", f"{args.output_prefix}.ops")

    screen_log.set("logEvery", str(args.echo_every))
    file_log.set("logEvery", str(args.log_every))
    file_log.set("fileName", f"{args.output_prefix}.log")
    tree_log.set("logEvery", str(args.tree_every))
    tree_log.set("fileName", f"{args.output_prefix}.trees")

    managed_logs = {"screenLog", "fileLog"}
    for log_node in root.iter("log"):
        log_id = log_node.get("id")
        log_file = log_node.get("fileName")
        if log_id in managed_logs or not log_file:
            continue
        log_node.set("fileName", prefixed_log_name(args.output_prefix, log_file))

    checkpoint_log.set("fileName", f"{args.output_prefix}.chkpt")
    checkpoint_log.set("checkpointFinal", str(args.chain_length))
    if args.checkpoint_every is not None:
        checkpoint_log.set("checkpointEvery", str(args.checkpoint_every))

    output_dir = os.path.dirname(args.output_xml)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    ET.indent(tree, space="  ")
    tree.write(args.output_xml, encoding="utf-8", xml_declaration=True)


if __name__ == "__main__":
    main()
