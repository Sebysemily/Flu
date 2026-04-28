from __future__ import annotations

import copy
import os
import shutil
import tempfile
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


P_NS = "http://schemas.openxmlformats.org/presentationml/2006/main"
A_NS = "http://schemas.openxmlformats.org/drawingml/2006/main"
R_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_NS = "http://schemas.openxmlformats.org/package/2006/relationships"

NS = {"p": P_NS, "a": A_NS, "r": R_NS, "rel": PKG_NS}

ET.register_namespace("a", A_NS)
ET.register_namespace("p", P_NS)
ET.register_namespace("r", R_NS)


ROOT = Path("/home/sebas/MIRA_NGS/flu_aviar")
PPTX_PATH = ROOT / "Ecuador_H5N1_Mutation_Characterization_Results.pptx"
RTT_PNG = Path("/tmp/flu_aviar_slide_assets/root_to_tip_regression.png")
RULEGRAPH_PNG = Path("/tmp/flu_aviar_slide_assets/pipeline_rulegraph.png")
FOOTER = "Avances en panel temporal y BEAST exploratorio | H5N1 Ecuador"


def qname(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"


def slide_path(num: int) -> Path:
    return Path("ppt/slides") / f"slide{num}.xml"


def rels_path(num: int) -> Path:
    return Path("ppt/slides/_rels") / f"slide{num}.xml.rels"


def parse_xml(path: Path) -> ET.ElementTree:
    return ET.parse(path)


def find_sp(root: ET.Element, shape_id: int) -> ET.Element:
    for sp in root.findall(".//p:sp", NS):
        c_nv = sp.find("./p:nvSpPr/p:cNvPr", NS)
        if c_nv is not None and int(c_nv.get("id")) == shape_id:
            return sp
    raise ValueError(f"Shape id {shape_id} not found")


def find_graphic_frame(root: ET.Element, frame_id: int) -> ET.Element:
    for gf in root.findall(".//p:graphicFrame", NS):
        c_nv = gf.find("./p:nvGraphicFramePr/p:cNvPr", NS)
        if c_nv is not None and int(c_nv.get("id")) == frame_id:
            return gf
    raise ValueError(f"GraphicFrame id {frame_id} not found")


def find_pic(root: ET.Element, pic_id: int) -> ET.Element:
    for pic in root.findall(".//p:pic", NS):
        c_nv = pic.find("./p:nvPicPr/p:cNvPr", NS)
        if c_nv is not None and int(c_nv.get("id")) == pic_id:
            return pic
    raise ValueError(f"Picture id {pic_id} not found")


def find_sp_by_name(root: ET.Element, name: str) -> ET.Element | None:
    for sp in root.findall(".//p:sp", NS):
        c_nv = sp.find("./p:nvSpPr/p:cNvPr", NS)
        if c_nv is not None and c_nv.get("name") == name:
            return sp
    return None


def find_pic_by_name(root: ET.Element, name: str) -> ET.Element | None:
    for pic in root.findall(".//p:pic", NS):
        c_nv = pic.find("./p:nvPicPr/p:cNvPr", NS)
        if c_nv is not None and c_nv.get("name") == name:
            return pic
    return None


def set_txbody_lines(txbody: ET.Element, lines: list[str]) -> None:
    body_pr = txbody.find("a:bodyPr", NS)
    lst_style = txbody.find("a:lstStyle", NS)
    template_p = txbody.find("a:p", NS)
    if body_pr is None:
        body_pr = ET.Element(qname(A_NS, "bodyPr"))
        txbody.insert(0, body_pr)
    if lst_style is None:
        lst_style = ET.Element(qname(A_NS, "lstStyle"))
        if len(txbody) == 1:
            txbody.append(lst_style)
        else:
            txbody.insert(1, lst_style)

    for p in list(txbody.findall("a:p", NS)):
        txbody.remove(p)

    if template_p is None:
        template_p = ET.Element(qname(A_NS, "p"))
        template_p.append(ET.Element(qname(A_NS, "endParaRPr")))

    for line in lines:
        p = copy.deepcopy(template_p)
        ppr = p.find("a:pPr", NS)
        end_para = p.find("a:endParaRPr", NS)
        rpr_template = None
        first_r = p.find("a:r", NS)
        if first_r is not None:
            rpr_template = first_r.find("a:rPr", NS)
        for child in list(p):
            if child.tag in {qname(A_NS, "r"), qname(A_NS, "br"), qname(A_NS, "fld")}:
                p.remove(child)
        r = ET.Element(qname(A_NS, "r"))
        if rpr_template is not None:
            r.append(copy.deepcopy(rpr_template))
        else:
            r.append(ET.Element(qname(A_NS, "rPr")))
        t = ET.Element(qname(A_NS, "t"))
        if line.startswith(" ") or line.endswith(" ") or "  " in line:
            t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        t.text = line
        r.append(t)
        insert_at = 1 if ppr is not None else 0
        p.insert(insert_at, r)
        if end_para is None:
            p.append(ET.Element(qname(A_NS, "endParaRPr")))
        txbody.append(p)


def set_shape_text(root: ET.Element, shape_id: int, lines: list[str]) -> None:
    sp = find_sp(root, shape_id)
    txbody = sp.find("p:txBody", NS)
    if txbody is None:
        raise ValueError(f"Shape {shape_id} has no text body")
    set_txbody_lines(txbody, lines)


def set_table_text(root: ET.Element, frame_id: int, rows: list[list[str]]) -> None:
    gf = find_graphic_frame(root, frame_id)
    tr_list = gf.findall(".//a:tbl/a:tr", NS)
    if len(tr_list) != len(rows):
        raise ValueError(f"Frame {frame_id} expected {len(tr_list)} rows, got {len(rows)}")
    for tr, values in zip(tr_list, rows, strict=True):
        tc_list = tr.findall("a:tc", NS)
        if len(tc_list) != len(values):
            raise ValueError(f"Frame {frame_id} row expected {len(tc_list)} cols, got {len(values)}")
        for tc, value in zip(tc_list, values, strict=True):
            txbody = tc.find("a:txBody", NS)
            if txbody is None:
                raise ValueError("Table cell missing txBody")
            set_txbody_lines(txbody, [value])


def update_footer(root: ET.Element, shape_id: int) -> None:
    set_shape_text(root, shape_id, [FOOTER])


def max_shape_id(root: ET.Element) -> int:
    max_id = 1
    for tag in ("p:sp", "p:pic", "p:graphicFrame"):
        for el in root.findall(f".//{tag}", NS):
            c_nv = (
                el.find("./p:nvSpPr/p:cNvPr", NS)
                or el.find("./p:nvPicPr/p:cNvPr", NS)
                or el.find("./p:nvGraphicFramePr/p:cNvPr", NS)
            )
            if c_nv is not None:
                max_id = max(max_id, int(c_nv.get("id")))
    return max_id


def add_textbox(
    root: ET.Element,
    shape_id: int,
    name: str,
    x: int,
    y: int,
    cx: int,
    cy: int,
    lines: list[str],
    font_size: str = "1800",
) -> None:
    existing = find_sp_by_name(root, name)
    if existing is not None:
        xfrm = existing.find("./p:spPr/a:xfrm", NS)
        if xfrm is not None:
            off = xfrm.find("a:off", NS)
            ext = xfrm.find("a:ext", NS)
            off.set("x", str(x))
            off.set("y", str(y))
            ext.set("cx", str(cx))
            ext.set("cy", str(cy))
        txbody = existing.find("p:txBody", NS)
        if txbody is not None:
            set_txbody_lines(txbody, lines)
        return

    sp = ET.Element(qname(P_NS, "sp"))
    nv = ET.SubElement(sp, qname(P_NS, "nvSpPr"))
    ET.SubElement(nv, qname(P_NS, "cNvPr"), {"id": str(shape_id), "name": name})
    ET.SubElement(nv, qname(P_NS, "cNvSpPr"), {"txBox": "1"})
    ET.SubElement(nv, qname(P_NS, "nvPr"))

    sp_pr = ET.SubElement(sp, qname(P_NS, "spPr"))
    xfrm = ET.SubElement(sp_pr, qname(A_NS, "xfrm"))
    ET.SubElement(xfrm, qname(A_NS, "off"), {"x": str(x), "y": str(y)})
    ET.SubElement(xfrm, qname(A_NS, "ext"), {"cx": str(cx), "cy": str(cy)})
    geom = ET.SubElement(sp_pr, qname(A_NS, "prstGeom"), {"prst": "rect"})
    ET.SubElement(geom, qname(A_NS, "avLst"))
    ET.SubElement(sp_pr, qname(A_NS, "noFill"))

    txbody = ET.SubElement(sp, qname(P_NS, "txBody"))
    body_pr = ET.SubElement(txbody, qname(A_NS, "bodyPr"), {"wrap": "square"})
    ET.SubElement(body_pr, qname(A_NS, "spAutoFit"))
    ET.SubElement(txbody, qname(A_NS, "lstStyle"))
    template_p = ET.Element(qname(A_NS, "p"))
    p_pr = ET.SubElement(template_p, qname(A_NS, "pPr"))
    ET.SubElement(
        p_pr,
        qname(A_NS, "defRPr"),
        {"sz": font_size},
    )
    ET.SubElement(template_p, qname(A_NS, "endParaRPr"))
    txbody.append(template_p)
    set_txbody_lines(txbody, lines)
    root.find("p:cSld/p:spTree", NS).append(sp)


def load_rel_root(path: Path) -> ET.ElementTree:
    return ET.parse(path)


def add_image_rel(rel_root: ET.Element, rel_id: str, target: str) -> None:
    ET.SubElement(
        rel_root.getroot(),
        qname(PKG_NS, "Relationship"),
        {
            "Id": rel_id,
            "Type": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image",
            "Target": target,
        },
    )


def next_rel_id(rel_root: ET.Element) -> str:
    used = []
    for rel in rel_root.getroot().findall("rel:Relationship", NS):
        rid = rel.get("Id", "")
        if rid.startswith("rId"):
            try:
                used.append(int(rid[3:]))
            except ValueError:
                pass
    return f"rId{max(used, default=0) + 1}"


def set_picture_rel(pic: ET.Element, rel_id: str, descr: str) -> None:
    blip = pic.find(".//a:blip", NS)
    if blip is None:
        raise ValueError("Picture missing blip")
    blip.set(qname(R_NS, "embed"), rel_id)
    c_nv = pic.find("./p:nvPicPr/p:cNvPr", NS)
    if c_nv is not None:
        c_nv.set("descr", descr)


def set_picture_geometry(pic: ET.Element, x: int, y: int, cx: int, cy: int) -> None:
    xfrm = pic.find("./p:spPr/a:xfrm", NS)
    if xfrm is None:
        raise ValueError("Picture missing transform")
    off = xfrm.find("a:off", NS)
    ext = xfrm.find("a:ext", NS)
    off.set("x", str(x))
    off.set("y", str(y))
    ext.set("cx", str(cx))
    ext.set("cy", str(cy))


def remove_graphic_frame(root: ET.Element, frame_id: int) -> None:
    sp_tree = root.find("p:cSld/p:spTree", NS)
    try:
        gf = find_graphic_frame(root, frame_id)
    except ValueError:
        return
    sp_tree.remove(gf)


def add_picture_from_template(
    root: ET.Element,
    template_pic: ET.Element,
    new_shape_id: int,
    name: str,
    rel_id: str,
    descr: str,
    x: int,
    y: int,
    cx: int,
    cy: int,
) -> None:
    existing = find_pic_by_name(root, name)
    if existing is not None:
        c_nv = existing.find("./p:nvPicPr/p:cNvPr", NS)
        c_nv.set("descr", descr)
        set_picture_rel(existing, rel_id, descr)
        set_picture_geometry(existing, x, y, cx, cy)
        return

    pic = copy.deepcopy(template_pic)
    c_nv = pic.find("./p:nvPicPr/p:cNvPr", NS)
    c_nv.set("id", str(new_shape_id))
    c_nv.set("name", name)
    c_nv.set("descr", descr)
    set_picture_rel(pic, rel_id, descr)
    set_picture_geometry(pic, x, y, cx, cy)
    root.find("p:cSld/p:spTree", NS).append(pic)


def write_tree(tree: ET.ElementTree, path: Path) -> None:
    tree.write(path, encoding="UTF-8", xml_declaration=True)


def update_slides(workdir: Path) -> None:
    slide_roots: dict[int, ET.ElementTree] = {}
    for num in range(1, 14):
        slide_roots[num] = parse_xml(workdir / slide_path(num))

    # Slide 1
    root = slide_roots[1].getroot()
    set_shape_text(root, 4, ["Avances en panel temporal y BEAST exploratorio"])
    set_shape_text(
        root,
        5,
        [
            "Pregunta principal: fechar la insercion del clado de Ecuador",
            "Objetivo: construir un panel temporalmente informativo y correr BEAST exploratorio",
        ],
    )
    set_shape_text(root, 6, ["74", "panel relajado"])
    set_shape_text(root, 7, ["R^2 0.88", "RTT pre-filter"])
    set_shape_text(root, 8, ["1", "outlier temporal"])
    set_shape_text(root, 9, ["Baseline", "strict + constant"])
    update_footer(root, 11)

    # Slide 2
    root = slide_roots[2].getroot()
    set_shape_text(root, 2, ["Flujo general"])
    set_shape_text(
        root,
        3,
        ["Preparacion del panel y corrida bayesiana separadas para iterar sin rehacer BEAST"],
    )
    set_shape_text(root, 5, ["Build subset", "QC", "Root-to-tip", "Filtro de outliers", "Export final"])
    set_shape_text(
        root,
        6,
        ["02_pre_beast prepara el panel", "03_beast corre strict + constant y UCLN + constant"],
    )
    remove_graphic_frame(root, 10)
    update_footer(root, 8)

    # Slide 3 (old slide 3 table slide now panel decisions)
    root = slide_roots[3].getroot()
    set_shape_text(root, 2, ["Decisiones del panel en esta iteracion"])
    set_table_text(
        root,
        3,
        [
            ["Parametro", "Antes", "Ahora"],
            ["regional_context", "35", "40"],
            ["american_anchor", "2", "6"],
            ["panel total", "65", "74"],
        ],
    )
    set_shape_text(
        root,
        4,
        [
            "Se mantuvo solo el perfil relajado.",
            "Mas contexto regional y mas anchors americanos mejoraron la senal temporal.",
        ],
    )
    set_table_text(
        root,
        5,
        [
            ["Decision", "Estado", "Nota"],
            ["perfil", "relajado", "switch eliminado"],
            ["RTT outliers", "filtrado", "sin rerun RTT"],
            ["export final", "data/beast", "panel limpio"],
        ],
    )
    set_shape_text(root, 6, ["El panel relajado fue el mas util para esta iteracion."])
    update_footer(root, 8)

    # Slide 4 (old slide 4 flow)
    root = slide_roots[4].getroot()
    set_shape_text(root, 39, ["Como se construye el subset"])
    set_shape_text(
        root,
        40,
        ["Parte de un nucleo de Ecuador y suma contexto filogenetico cercano con control de redundancia"],
    )
    set_shape_text(root, 41, ["1", "Nucleo Ecuador"])
    set_shape_text(root, 42, ["seed del clado"])
    set_shape_text(root, 44, ["2", "Contexto regional"])
    set_shape_text(root, 45, ["cercania filogenetica"])
    set_shape_text(root, 47, ["3", "Control redundancia"])
    set_shape_text(root, 48, ["pais y mes"])
    set_shape_text(root, 50, ["4", "Anchors Americas"])
    set_shape_text(root, 51, ["mejorar senal temporal"])
    set_shape_text(root, 53, ["5", "Panel final"])
    set_shape_text(root, 54, ["diverso y defendible"])
    set_shape_text(root, 55, ["No se arma al azar: parte de Ecuador y anade contexto cercano."])
    set_shape_text(root, 56, ["Se fuerza diversidad temporal y geografica para evitar paneles redundantes."])
    set_shape_text(root, 57, ["Nucleo ecuatoriano + contexto regional + anchors americanos enmarcan mejor el evento de interes."])
    update_footer(root, 59)

    # Slide 5
    root = slide_roots[5].getroot()
    set_shape_text(root, 2, ["QC antes de BEAST"])
    set_shape_text(root, 3, ["1", "QC fuente"])
    set_shape_text(root, 4, ["secuencias origen"])
    set_shape_text(root, 6, ["2", "Filtro extremo"])
    set_shape_text(root, 7, ["problemas graves"])
    set_shape_text(root, 9, ["3", "Root-to-tip"])
    set_shape_text(root, 10, ["primer barrido"])
    set_shape_text(root, 12, ["4", "Outliers RTT"])
    set_shape_text(root, 13, ["discordancia temporal"])
    set_shape_text(root, 15, ["5", "Excluir"])
    set_shape_text(root, 16, ["antes del panel"])
    set_shape_text(root, 18, ["6", "Exportar"])
    set_shape_text(root, 19, ["a data/beast"])
    set_shape_text(
        root,
        20,
        [
            "El pipeline combina QC de secuencia y TreeTime para detectar outliers en un RTT pre-filter y excluirlos antes de exportar el panel final, sin recalcular RTT.",
        ],
    )
    update_footer(root, 22)

    # Slide 6 (old slide 6 now next runs)
    root = slide_roots[6].getroot()
    set_shape_text(root, 2, ["Proxima tanda de corridas"])
    set_table_text(
        root,
        3,
        [
            ["Escenario", "Plan"],
            ["strict + constant", "100M"],
            ["UCLN + constant", "120M"],
            ["log/tree/echo", "10000"],
            ["tiempo por replica", "~5-6 h"],
            ["prioridad", "profundizar constant"],
            ["siguiente salto", "luego exponential"],
        ],
    )
    set_shape_text(root, 4, ["100M", "strict_constant"])
    set_shape_text(root, 5, ["120M", "ucln_constant"])
    set_shape_text(
        root,
        6,
        [
            "Antes de pasar a exponential, conviene estabilizar reloj y datacion bajo constant.",
        ],
    )
    update_footer(root, 8)

    # Slide 7 (old slide 7 results table)
    root = slide_roots[7].getroot()
    set_shape_text(root, 2, ["Resultado del primer exploratorio"])
    set_shape_text(
        root,
        3,
        ["strict + constant mezcla mejor; UCLN ajusta mejor likelihood pero aun no madura"],
    )
    set_table_text(
        root,
        4,
        [
            ["Escenario", "Metrica", "Valor", "Lectura", "Decision"],
            ["strict + constant", "treeLikelihood", "~ -28606", "baseline estable", "referencia"],
            ["strict + constant", "rootHeight ESS", "~296 combinado", "temporal usable", "ok"],
            ["strict + constant", "clock.rate ESS", "~385 combinado", "buena mezcla", "ok"],
            ["strict + constant", "constant.popSize ESS", "~354 combinado", "estable", "ok"],
            ["UCLN + constant", "treeLikelihood", "~ -28558", "likelihood mejor", "no basta"],
            ["UCLN + constant", "rootHeight ESS", "~50 combinado", "temporal flojo", "seguir"],
            ["UCLN + constant", "meanRate ESS", "~55 combinado", "mezcla pobre", "seguir"],
            ["UCLN + constant", "ucld.mean/stdev", "~135 / ~172", "mas flexible", "mas cadena"],
        ],
    )
    set_shape_text(root, 5, ["Baseline actual: strict + constant"])
    update_footer(root, 7)

    # Slide 8 (old slide 8 image slide now RTT)
    root = slide_roots[8].getroot()
    set_shape_text(root, 2, ["Root-to-tip regression"])
    update_footer(root, 5)
    new_id = max_shape_id(root) + 1
    add_textbox(
        root,
        new_id,
        "TextBox RTT",
        320040,
        1188720,
        3200400,
        2743200,
        [
            "- Pendiente positiva",
            "- R^2 = 0.88",
            "- Root date = 2018.05",
            "- 1 outlier temporal marcado en rojo",
        ],
        font_size="1800",
    )

    # Slide 9 (old slide 9 conclusion layout)
    root = slide_roots[9].getroot()
    set_shape_text(root, 2, ["Conclusion"])
    set_shape_text(root, 3, ["Senal temporal"])
    set_shape_text(root, 4, ["R^2 = 0.88 y root date 2018.05 en RTT pre-filter."])
    set_shape_text(root, 5, ["Filtro RTT"])
    set_shape_text(root, 6, ["1 outlier temporal detectado en RTT y excluido antes del panel final."])
    set_shape_text(root, 7, ["Baseline BEAST"])
    set_shape_text(root, 8, ["strict + constant mezcla mejor y sirve como referencia temporal inicial."])
    set_shape_text(
        root,
        9,
        ["El panel ahora llega mas limpio y defendible a BEAST; la siguiente decision es si UCLN mejora con cadenas mas largas."],
    )
    update_footer(root, 11)

    # Slide 10
    root = slide_roots[10].getroot()
    set_shape_text(root, 2, ["XMLs y particiones"])
    set_table_text(
        root,
        3,
        [
            ["Aspecto", "Decision"],
            ["Templates base", "versionados"],
            ["Regeneracion", "no desde BEAUti cada iteracion"],
            ["Particiones de sitio", "se preservan del template"],
        ],
    )
    set_table_text(
        root,
        4,
        [
            ["Solo se reescribe", "Uso"],
            ["chainLength", "profundidad de cadena"],
            ["logEvery / treeEvery", "muestreo y salida"],
            ["echoEvery / outputs", "seguimiento reproducible"],
            ["comparabilidad", "mismo particionado entre escenarios"],
        ],
    )
    update_footer(root, 6)

    # Slide 11
    root = slide_roots[11].getroot()
    set_shape_text(root, 2, ["Outlier temporal"])
    set_shape_text(
        root,
        3,
        [
            "- Contexto regional, no Ecuador core",
            "- Residual = 4.51",
            "- Detectado en RTT pre-filter",
            "- Excluido antes del panel final",
        ],
    )
    set_shape_text(
        root,
        4,
        ["1775-N_OR878092.1__regional_context/", "Brazil/2023-06-23"],
    )
    update_footer(root, 6)

    # Slide 12
    root = slide_roots[12].getroot()
    set_shape_text(root, 2, ["Export final para BEAST"])
    update_footer(root, 5)

    # Slide 13
    root = slide_roots[13].getroot()
    set_shape_text(root, 2, ["Primer BEAST exploratorio"])
    update_footer(root, 5)

    # Add text boxes to slides 12 and 13
    for num, lines, font_size in [
        (
            12,
            [
                "- Intermedios -> data/pre_beast",
                "- Final limpio -> data/beast",
                "- panel_main_taxa.final.tsv",
                "- panel_main_concat.final.fasta",
                "- final_panel_segment/H5N1_*.fasta",
                "- 72 taxa kept tras excluir 1 outlier RTT",
            ],
            "1800",
        ),
        (
            13,
            [
                "- Escenarios: strict + constant y UCLN + constant",
                "- Dos replicas por escenario",
                "- Objetivo: evaluar mezcla y estabilidad",
                "- Exploratorio: no cerrar aun la pregunta final",
            ],
            "1900",
        ),
    ]:
        root = slide_roots[num].getroot()
        new_id = max_shape_id(root) + 1
        add_textbox(
            root,
            new_id,
            f"TextBox Added {num}",
            731520,
            1097280 if num == 12 else 1629070,
            10789920,
            5017562 if num == 12 else 3355017,
            lines,
            font_size=font_size,
        )

    for num, tree in slide_roots.items():
        write_tree(tree, workdir / slide_path(num))


def add_images_and_relationships(workdir: Path) -> None:
    media_dir = workdir / "ppt/media"
    existing = sorted(media_dir.glob("image*.*"))
    next_idx = 1
    for path in existing:
        stem = path.stem.replace("image", "")
        if stem.isdigit():
            next_idx = max(next_idx, int(stem) + 1)

    rtt_name = f"image{next_idx}.png"
    rule_name = f"image{next_idx + 1}.png"
    shutil.copy2(RTT_PNG, media_dir / rtt_name)
    shutil.copy2(RULEGRAPH_PNG, media_dir / rule_name)

    # Slide 8: update existing picture to RTT
    slide8 = parse_xml(workdir / slide_path(8))
    rel8 = load_rel_root(workdir / rels_path(8))
    pic = find_pic(slide8.getroot(), 3)
    rid = "rId2"
    for rel in rel8.getroot().findall("rel:Relationship", NS):
        if rel.get("Id") == rid:
            rel.set("Target", f"../media/{rtt_name}")
            break
    set_picture_rel(pic, rid, "root_to_tip_regression.png")
    write_tree(slide8, workdir / slide_path(8))
    write_tree(rel8, workdir / rels_path(8))

    # Slide 2: add picture with rulegraph
    slide2 = parse_xml(workdir / slide_path(2))
    rel2 = load_rel_root(workdir / rels_path(2))
    template_pic = find_pic(parse_xml(workdir / slide_path(8)).getroot(), 3)
    new_rid = next_rel_id(rel2)
    add_image_rel(rel2, new_rid, f"../media/{rule_name}")
    root2 = slide2.getroot()
    add_picture_from_template(
        root2,
        template_pic,
        max_shape_id(root2) + 1,
        "Picture Rulegraph",
        new_rid,
        "pipeline_rulegraph.png",
        249382,
        1179576,
        6567054,
        4818888,
    )
    write_tree(slide2, workdir / slide_path(2))
    write_tree(rel2, workdir / rels_path(2))


def rezip(src_dir: Path, dest_pptx: Path) -> None:
    with zipfile.ZipFile(dest_pptx, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(src_dir.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(src_dir))


def reorder_slides(workdir: Path) -> None:
    presentation = parse_xml(workdir / "ppt/presentation.xml")
    root = presentation.getroot()
    sld_id_lst = root.find("p:sldIdLst", NS)
    current = list(sld_id_lst)

    rels = parse_xml(workdir / "ppt/_rels/presentation.xml.rels").getroot()
    relmap = {
        rel.get("Id"): rel.get("Target")
        for rel in rels.findall("rel:Relationship", NS)
    }

    slide_by_title: dict[str, ET.Element] = {}
    extras: list[ET.Element] = []
    desired_titles = [
        "Avances en panel temporal y BEAST exploratorio",
        "Flujo general",
        "Como se construye el subset",
        "Decisiones del panel en esta iteracion",
        "QC antes de BEAST",
        "Root-to-tip regression",
        "Outlier temporal",
        "Export final para BEAST",
        "XMLs y particiones",
        "Primer BEAST exploratorio",
        "Resultado del primer exploratorio",
        "Proxima tanda de corridas",
        "Conclusion",
    ]

    for sld in current:
        rel_id = sld.get(qname(R_NS, "id"))
        target = relmap.get(rel_id)
        if not target:
            extras.append(sld)
            continue
        slide_xml = parse_xml(workdir / "ppt" / target).getroot()
        texts = [t.text for t in slide_xml.findall(".//a:t", NS) if t.text]
        title = texts[0] if texts else ""
        if title in desired_titles and title not in slide_by_title:
            slide_by_title[title] = sld
        else:
            extras.append(sld)

    reordered = [slide_by_title[title] for title in desired_titles if title in slide_by_title] + extras
    sld_id_lst[:] = reordered
    write_tree(presentation, workdir / "ppt/presentation.xml")


def main() -> None:
    if not RTT_PNG.exists():
        raise FileNotFoundError(RTT_PNG)
    if not RULEGRAPH_PNG.exists():
        raise FileNotFoundError(RULEGRAPH_PNG)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        extract_dir = tmpdir_path / "pptx"
        with zipfile.ZipFile(PPTX_PATH) as zf:
            zf.extractall(extract_dir)

        update_slides(extract_dir)
        add_images_and_relationships(extract_dir)
        reorder_slides(extract_dir)

        out_path = tmpdir_path / PPTX_PATH.name
        rezip(extract_dir, out_path)
        shutil.copy2(out_path, PPTX_PATH)


if __name__ == "__main__":
    main()
