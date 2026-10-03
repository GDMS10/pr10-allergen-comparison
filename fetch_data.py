#!/usr/bin/env python3
"""Download every PR-10 structure from the RCSB and write the entry metadata.

The PR-10 family is InterPro IPR000916 (Bet v I / major latex protein) plus
IPR024949 (Bet v I type allergen), which is how the 133 entries were found.
Structures go to pdb_all/, which is gitignored: 140 MB of coordinates does not
belong in the repository, and fetch_data.py rebuilds it in about a minute.

Run:  python3 fetch_data.py
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import os
import urllib.error
import urllib.request

PDB_DIR = "pdb_all"
META = "meta.json"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA = "https://data.rcsb.org/rest/v1/core"
DOWNLOAD = "https://files.rcsb.org/download/{}.pdb"

INTERPRO = ["IPR000916", "IPR024949"]


def _get_json(url, data=None, timeout=60):
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    try:
        return json.load(urllib.request.urlopen(req, timeout=timeout))
    except urllib.error.HTTPError as e:
        return {"__error__": e.code, "body": e.read().decode()[:300]}


def find_entries():
    """Query the RCSB search API for every entry carrying a PR-10 InterPro hit."""
    entries = set()
    for ipr in INTERPRO:
        body = {
            "query": {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                "operator": "exact_match", "value": ipr}},
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 1000},
                                "results_verbosity": "compact"},
        }
        res = _get_json(SEARCH, json.dumps(body).encode())
        if "__error__" in res:
            raise RuntimeError(f"RCSB search failed for {ipr}: {res}")
        entries |= set(res.get("result_set", []))
    return sorted(entries)


def entry_meta(entry):
    """Title, method, resolution, and the PR-10-annotated polymer entity for one entry."""
    d = _get_json(f"{DATA}/entry/{entry}")
    if "__error__" in d:
        return entry, {"error": "entry fetch failed"}

    entities = []
    for pe in d.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []):
        e = _get_json(f"{DATA}/polymer_entity/{entry}/{pe}")
        if "__error__" in e:
            continue
        annotations = {a["annotation_id"] for a in e.get("rcsb_polymer_entity_annotation", [])}
        if not annotations & set(INTERPRO):
            continue  # the Fab, a ligand, or something else in the same entry
        src = (e.get("rcsb_entity_source_organism") or [{}])[0]
        entities.append(dict(
            entity=pe,
            name=e["rcsb_polymer_entity"].get("pdbx_description"),
            seq=e["entity_poly"].get("pdbx_seq_one_letter_code_can"),
            org=src.get("ncbi_scientific_name") or src.get("scientific_name"),
            uniprot=e.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids"),
        ))
    return entry, dict(
        title=d["struct"]["title"],
        method=[m["method"] for m in d.get("exptl", [])],
        res=d.get("rcsb_entry_info", {}).get("resolution_combined"),
        entities=entities,
    )


def fetch(entry):
    os.makedirs(PDB_DIR, exist_ok=True)
    out = os.path.join(PDB_DIR, f"{entry}.pdb")
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        return entry, os.path.getsize(out)
    try:
        data = urllib.request.urlopen(urllib.request.Request(DOWNLOAD.format(entry)),
                                      timeout=60).read()
        open(out, "wb").write(data)
        return entry, len(data)
    except Exception as e:
        return entry, f"ERR {str(e)[:60]}"


def main():
    entries = find_entries()
    print(f"{len(entries)} PR-10 entries")

    meta = {}
    with cf.ThreadPoolExecutor(10) as ex:
        for entry, value in ex.map(entry_meta, entries):
            meta[entry] = value
    kept = {k: v for k, v in meta.items() if v.get("entities")}
    print(f"{len(kept)} have a PR-10 entity annotated")
    json.dump(kept, open(META, "w"), indent=1)

    failed = []
    with cf.ThreadPoolExecutor(8) as ex:
        for entry, size in ex.map(fetch, entries):
            if isinstance(size, str):
                failed.append((entry, size))
    print(f"structures: {len(entries) - len(failed)} downloaded, {len(failed)} failed {failed}")
    print(f"wrote {META} and {PDB_DIR}/")


if __name__ == "__main__":
    main()