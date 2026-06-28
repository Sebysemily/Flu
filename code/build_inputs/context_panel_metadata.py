"""Build GISAID context panel rows from the raw context FASTA."""

from __future__ import annotations

import re

from build_inputs.process_raw_to_segments import (
    dedupe_metadata_rows,
    filter_complete_context_isolates,
    parse_context_isolates,
)

PANEL_COLUMNS = [
    "file_name", "host", "host_type", "collection_date",
    "country", "province", "expected_role",
]

# ---------------------------------------------------------------------------
# Host-type classification
# ---------------------------------------------------------------------------

# Known non-avian mammal host compact keys (lowercase alphanumeric only)
_MAMMAL_HOST_KEYS: frozenset[str] = frozenset({
    # Pinnipeds
    "sealion", "southamericansealion", "southamericanfurseal", "southernsealion",
    "otariaflavescens", "arctocephalusaustralis",
    "elephantseal", "southernelephantseal", "pinniped",
    # Cetaceans
    "dolphin", "porpoise", "burmeistersporpoise",
    # Mustelids / otters
    "marineotter", "riverotter", "chungungo",
    # Felids / other carnivores / mammals
    "feline", "lion", "nasuanasua", "mink",
})

# Substrings that, when found inside the compact host key, indicate avian
_AVIAN_KEYWORDS: frozenset[str] = frozenset({
    "bird", "avian", "chicken", "duck", "turkey", "goose", "gull", "eagle",
    "hawk", "owl", "tern", "pelican", "vulture", "crow", "raven", "swan",
    "teal", "mallard", "wigeon", "merganser", "pheasant", "pigeon", "heron",
    "cormorant", "gannet", "albatross", "penguin", "skua", "petrel", "booby",
    "frigatebird", "condor", "falcon", "parrot", "ibis", "oystercatcher",
    "skimmer", "grebe", "grabe", "curlew", "whimbrel", "sanderling", "plover",
    "stork", "caracara", "eider", "brant", "shoveler", "gadwall", "scaup",
    "goldeneye", "flamingo", "quail", "dove", "cuckoo", "kingfisher",
    "swift", "swallow", "finch", "sparrow", "warbler", "thrush", "starling",
    "magpie", "blackbird", "coot", "crane", "kite", "harrier", "osprey",
    "kestrel", "merlin", "lark", "pipit", "wagtail", "murre", "puffin",
    "razorbill", "guillemot", "auklet", "murrelet", "kittiwake", "jaeger",
    "noddy", "tropicbird", "shearwater", "fulmar", "prion", "poultry",
    "backyard", "wildbird", "seabird", "emu", "spoonbill", "hornero",
    "guanay", "egret",
})

# First word of Latin binomials known to be avian genera (compact lowercase)
_AVIAN_GENERA: frozenset[str] = frozenset({
    # Cormorants / gannets / boobies / frigatebirds
    "phalacrocorax", "morus", "sula", "fregata", "phaethon",
    # Terns / gulls / skimmers
    "thalasseus", "sterna", "onychoprion", "anous", "hydroprogne",
    "larus", "leucophaeus", "chroicocephalus", "ichthyaetus", "rissa",
    "rynchops",
    # Pelicans / herons / ibis
    "pelecanus", "ardea", "egretta", "nycticorax", "bubulcus", "threskiornis",
    "platalea",
    # Shorebirds / waders
    "pluvialis", "calidris", "tringa", "numenius", "limosa", "arenaria",
    "phalaropus", "charadrius",
    # Tubenoses (petrels / shearwaters / albatrosses)
    "procellaria", "ardenna", "puffinus", "fulmarus", "daption", "pterodroma",
    "oceanites", "diomedea", "thalassarche", "macronectes",
    # Ducks / geese / swans
    "coscoroba", "cygnus", "dendrocygna", "anser", "branta", "chen",
    "aythya", "bucephala", "mergus", "lophodytes", "oxyura",
    # Galliformes / ratites
    "gallus", "numida", "meleagris", "dromaius",
    # Passerines / raptors / owls / others
    "furnarius", "megascops", "buteogallus", "coragyps", "cathartes",
    "buteo", "accipiter", "falco", "haliaeetus", "pandion",
    "strix", "tyto", "otus",
    "columba", "streptopelia", "zenaida",
    "turdus", "sialia", "mimus", "sturnus",
    "passer", "spizella", "zonotrichia", "setophaga",
})


def classify_host_type(host: str) -> str:
    """Classify a host string into specific host groups."""
    if not host or str(host).lower() in {"", "unknown", "na", "nan", "none", "avian", "bird", "environment"}:
        return "?"

    host_lower = host.lower()
    
    # Domesticated bird
    if any(kw in host_lower for kw in ["chicken", "turkey", "gallus", "numida", "guineafowl", "pheasant", "poultry"]):
        return "domesticated bird"
        
    # Wild bird
    elif any(kw in host_lower for kw in [
        "duck", "goose", "eider", "wigeon", "teal", "swan", "merganser", "shoveler", "scaup", "mallard", 
        "gadwall", "goldeneye", "brant", "cygnus", "dendrocygna", "mergus", "coscoroba", "pelican", 
        "booby", "gull", "tern", "pelecanus", "skimmer", "cormorant", "phalacrocorax", "frigatebird", 
        "sterna", "thalasseus", "sula", "fregata", "penguin", "fulmar", "shearwater", "sanderling", 
        "plover", "seabird", "puffinus", "calidris", "seagull", "vulture", "eagle", "owl", "hawk", 
        "falcon", "caracara", "condor", "coragyps", "buteogallus", "skua", "megascops", "crow", "raven", 
        "furnarius", "dromaius", "emu", "backyard_bird", "backyard bird", "wild_bird", "wild bird",
        "strix", "tyto", "otus", "columba", "streptopelia", "zenaida", "turdus", "sialia", "mimus", "sturnus",
        "passer", "spizella", "zonotrichia", "setophaga"
    ]):
        return "wild bird"
        
    # Domesticated mammal
    elif any(kw in host_lower for kw in ["cow", "cattle", "bovine", "bos", "cat", "feline", "dog", "canine", "mink", "swine", "pig", "sus"]):
        return "domesticated mammal"
        
    # Wild mammal
    elif any(kw in host_lower for kw in [
        "sea_lion", "sea lion", "pinniped", "dolphin", "porpoise", "elephant_seal", "elephant seal", 
        "fur_seal", "fur seal", "marine_otter", "otter", "otaria", "chungungo", "seal", "arctocephalus", 
        "lion", "nasua", "coati", "fox", "vulpes", "bear", "ursus", "skunk", "racoon", "raccoon", "procyon", 
        "mustela", "martes", "puma", "leopard", "tiger", "panthera"
    ]):
        return "wild mammal"
        
    # If "human", map to "?"
    elif any(kw in host_lower for kw in ["human", "homo", "person", "patient"]):
        return "?"

    return "?"


def read_epi_isls_from_fasta(path: str) -> frozenset[str]:
    """Return all EPI_ISL identifiers found in FASTA headers at *path*."""
    if not path:
        return frozenset()
    try:
        with open(path, encoding="utf-8") as fh:
            return frozenset(
                m.group()
                for line in fh
                if line.startswith(">")
                for m in [re.search(r"EPI_ISL_\d+", line)]
                if m
            )
    except FileNotFoundError:
        return frozenset()


# ---------------------------------------------------------------------------
# Ecuador / GISAID role helpers
# ---------------------------------------------------------------------------

def is_ecuador_country(country: str) -> bool:
    key = str(country or "").strip().lower().replace("_", "").replace(" ", "")
    return key == "ecuador"


def context_expected_role(country: str, default_role: str) -> str:
    """Ecuadorian GISAID context isolates are coastal Ecuador panel samples."""
    if is_ecuador_country(country):
        return "flu_costa"
    return default_role


# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------

def build_gisaid_context_rows(
    context_fasta: str,
    local_epi_isls: set[str],
    human_epi_isls: frozenset[str] = frozenset(),
    avian_epi_isls: frozenset[str] = frozenset(),
) -> list[dict[str, str]]:
    """One metadata row per GISAID EPI_ISL in the context FASTA (non-local isolates).

    host_type priority for each EPI_ISL:
        1. In human_epi_isls  → 'human'
        2. In avian_epi_isls  → 'avian'
        3. classify_host_type(host) fallback
    """
    isolates_data = parse_context_isolates(context_fasta)
    complete_context, context_dates, context_places, context_types, context_provinces = (
        filter_complete_context_isolates(isolates_data, local_epi_isls)
    )

    rows: list[dict[str, str]] = []
    for isolate, segs in sorted(complete_context.items()):
        country = context_places[isolate]
        province = context_provinces[isolate]
        date_value = context_dates[isolate]
        role = context_expected_role(country, context_types[isolate])
        parts = isolate.split("/")
        host = parts[1] if (len(parts) > 1 and parts[0] == "A") else "unknown"

        for _seg, (epi_isl, _seq, _hdr) in segs.items():
            if epi_isl in human_epi_isls:
                host_type = "Human"
            elif epi_isl in avian_epi_isls:
                host_type = classify_host_type(host)
                if host_type == "?":
                    host_type = "Terrestrial Birds" # Fallback if forced avian but unknown group
            else:
                host_type = classify_host_type(host)

            rows.append(
                {
                    "file_name": epi_isl,
                    "host": host,
                    "host_type": host_type,
                    "collection_date": date_value,
                    "country": country,
                    "province": province,
                    "expected_role": role,
                }
            )

    return dedupe_metadata_rows(rows)
