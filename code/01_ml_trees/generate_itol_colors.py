#!/usr/bin/env python3
import argparse
import re
import os

def main():
    parser = argparse.ArgumentParser(description="Genera plantillas de colores iTOL a partir de un árbol Newick.")
    parser.add_argument("-i", "--tree", required=True, help="Ruta al archivo de árbol Newick de entrada.")
    parser.add_argument("--tree-colors", required=True, help="Ruta de salida para TREE_COLORS.txt.")
    parser.add_argument("--colorstrip", required=True, help="Ruta de salida para DATASET_COLORSTRIP.txt.")
    args = parser.parse_args()

    if not os.path.exists(args.tree):
        raise FileNotFoundError(f"No se encontró el árbol: {args.tree}")

    # Leer el archivo de árbol
    with open(args.tree, "r") as f:
        content = f.read()

    # Extraer taxa/tips usando regex robusto
    taxa = re.findall(r"(?<=[(,])([^():,;]+):", content)
    taxa = [t.strip() for t in taxa if t.strip()]

    # Crear directorio padre si no existe
    os.makedirs(os.path.dirname(os.path.abspath(args.tree_colors)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.colorstrip)), exist_ok=True)

    # Escribir TREE_COLORS
    with open(args.tree_colors, "w") as out:
        out.write("TREE_COLORS\n")
        out.write("SEPARATOR COMMA\n")
        out.write("DATA\n")
        for t in taxa:
            t_lower = t.lower()
            if "flu-epi_isl" in t_lower or "flu-0406" in t_lower:
                color, style = "#FF0000", "bold"
            elif "flu-" in t_lower:
                color, style = "#00008B", "bold"
            elif "american_anchor" in t_lower:
                color, style = "#800080", "normal"
            else:
                color, style = "#008000", "normal"
            out.write(f"{t},label,{color},{style}\n")
            out.write(f"{t},branch,{color},normal\n")

    # Escribir DATASET_COLORSTRIP
    with open(args.colorstrip, "w") as out:
        out.write("DATASET_COLORSTRIP\n")
        out.write("SEPARATOR COMMA\n")
        out.write("DATASET_LABEL,clades_colors\n")
        out.write("COLOR,#ff0000\n")
        out.write("LEGEND_TITLE,Clades\n")
        out.write("LEGEND_SHAPES,1,1,1,1\n")
        out.write("LEGEND_COLORS,#FF0000,#00008B,#800080,#008000\n")
        out.write("LEGEND_LABELS,costa,sierra,american_anchor,regional_context\n")
        out.write("DATA\n")
        for t in taxa:
            t_lower = t.lower()
            if "flu-epi_isl" in t_lower or "flu-0406" in t_lower:
                color, label = "#FF0000", "costa"
            elif "flu-" in t_lower:
                color, label = "#00008B", "sierra"
            elif "american_anchor" in t_lower:
                color, label = "#800080", "american_anchor"
            else:
                color, label = "#008000", "regional_context"
            out.write(f"{t},{color},{label}\n")

    print(f"Plantillas de colores iTOL generadas con éxito para {len(taxa)} taxa.")

if __name__ == "__main__":
    main()
