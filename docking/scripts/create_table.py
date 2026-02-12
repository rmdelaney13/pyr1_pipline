import argparse
import csv
import fnmatch
import os
import re
import sys
from configparser import ConfigParser


def _debug(enabled, message):
    if enabled:
        print(f"[create_table][DEBUG] {message}")


def _safe_get(config, section, key, fallback=None):
    if section in config and key in config[section]:
        return config[section][key]
    if key in config["DEFAULT"]:
        return config["DEFAULT"][key]
    return fallback


def _safe_getbool(config, section, key, fallback=False):
    raw = _safe_get(config, section, key, None)
    if raw is None:
        return fallback
    return str(raw).strip().lower() in {"1", "true", "yes", "on"}


def _split_patterns(raw_patterns):
    return [item for item in raw_patterns.split(" ") if item.strip()]


def _parse_target_atom_triplets(raw):
    if raw is None:
        return []
    items = re.split(r"[;\n]+", raw.strip())
    out = []
    for item in items:
        token = item.strip()
        if not token:
            continue
        atoms = tuple([x.strip() for x in token.split("-") if x.strip()])
        if len(atoms) != 3:
            raise ValueError(f"Target atom triplet must have 3 atoms, got: '{token}'")
        out.append(atoms)
    return out


def _parse_acceptor_modes(raw):
    if raw is None:
        return ["auto"]
    tokens = [token.strip().lower() for token in re.split(r"[,\s]+", raw) if token.strip()]
    if not tokens:
        return ["auto"]
    return tokens


def _read_params_atom_names(params_path):
    atom_names = []
    with open(params_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            # Rosetta atom names used by downstream alignment.
            if line.startswith("ATOM "):
                fields = line.split()
                if len(fields) >= 2:
                    atom_names.append(fields[1].strip())
    return atom_names


def _read_sdf_atom_symbols(molecule_sdf):
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(molecule_sdf, removeHs=False)
    for entry in supplier:
        if entry is not None:
            return [atom.GetSymbol() for atom in entry.GetAtoms()]
    return []


def _read_sdf_atom_records(molecule_sdf):
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(molecule_sdf, removeHs=False)
    for entry in supplier:
        if entry is None:
            continue
        conf = entry.GetConformer()
        records = []
        for atom in entry.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            records.append(
                {
                    "idx": idx,
                    "symbol": atom.GetSymbol().upper(),
                    "coord": (float(pos.x), float(pos.y), float(pos.z)),
                }
            )
        return records
    return []


def _read_pdb_atom_records(pdb_path):
    records = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            name = line[12:16].strip()
            elem = line[76:78].strip().upper()
            if not elem:
                letters = "".join([c for c in name if c.isalpha()])
                elem = letters[:1].upper() if letters else ""
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                continue
            records.append(
                {
                    "name": name,
                    "symbol": elem,
                    "coord": (x, y, z),
                }
            )
    return records


def _coord_centroid(records):
    if len(records) == 0:
        return (0.0, 0.0, 0.0)
    sx = sum([r["coord"][0] for r in records])
    sy = sum([r["coord"][1] for r in records])
    sz = sum([r["coord"][2] for r in records])
    n = float(len(records))
    return (sx / n, sy / n, sz / n)


def _sq_dist(a, b):
    return (
        (a[0] - b[0]) * (a[0] - b[0])
        + (a[1] - b[1]) * (a[1] - b[1])
        + (a[2] - b[2]) * (a[2] - b[2])
    )


def _build_sdf_idx_to_pdb_name_map(molecule_sdf, pdb_path, debug=False):
    sdf_records = _read_sdf_atom_records(molecule_sdf)
    pdb_records = _read_pdb_atom_records(pdb_path)
    if len(sdf_records) == 0 or len(pdb_records) == 0:
        _debug(
            debug,
            f"coord_map unavailable sdf_records={len(sdf_records)} pdb_records={len(pdb_records)}",
        )
        return {}

    sdf_ctr = _coord_centroid(sdf_records)
    pdb_ctr = _coord_centroid(pdb_records)

    def centered(coord, ctr):
        return (coord[0] - ctr[0], coord[1] - ctr[1], coord[2] - ctr[2])

    candidates_by_sdf = []
    for s in sdf_records:
        cand = []
        s_coord = centered(s["coord"], sdf_ctr)
        for p_idx, p in enumerate(pdb_records):
            if p["symbol"] != s["symbol"]:
                continue
            p_coord = centered(p["coord"], pdb_ctr)
            cand.append((_sq_dist(s_coord, p_coord), p_idx))
        cand.sort(key=lambda x: x[0])
        candidates_by_sdf.append((s["idx"], s["symbol"], cand))

    # Deterministic greedy assignment: hardest atoms first.
    candidates_by_sdf.sort(key=lambda item: (len(item[2]), item[0]))
    used_pdb = set()
    idx_to_name = {}
    assigned_dists = []
    for sdf_idx, symbol, cand in candidates_by_sdf:
        chosen = None
        for dist2, p_idx in cand:
            if p_idx in used_pdb:
                continue
            chosen = (dist2, p_idx)
            break
        if chosen is None:
            continue
        dist2, p_idx = chosen
        used_pdb.add(p_idx)
        idx_to_name[sdf_idx] = pdb_records[p_idx]["name"]
        assigned_dists.append(dist2)

    if assigned_dists:
        max_dist = max(assigned_dists) ** 0.5
        mean_dist = (sum(assigned_dists) / len(assigned_dists)) ** 0.5
    else:
        max_dist = 0.0
        mean_dist = 0.0
    _debug(
        debug,
        f"coord_map assigned={len(idx_to_name)}/{len(sdf_records)} "
        f"max_dist={max_dist:.4f} mean_dist={mean_dist:.4f}",
    )
    return idx_to_name


def _get_rdkit_acceptor_triplets(
    molecule_sdf,
    include_reverse=True,
    acceptor_smarts=None,
    acceptor_mode="auto",
    debug=False,
):
    from rdkit import Chem, RDConfig
    from rdkit.Chem import ChemicalFeatures

    supplier = Chem.SDMolSupplier(molecule_sdf, removeHs=False)
    mol = None
    for entry in supplier:
        if entry is not None:
            mol = entry
            break
    if mol is None:
        return []

    mode = str(acceptor_mode or "generic").strip().lower()
    if mode not in {
        "auto",
        "generic",
        "carbonyl",
        "sulfonyl",
        "nitro",
        "hydroxyl_ether",
        "ring_n",
    }:
        raise ValueError(
            "AcceptorMode must be one of: auto, generic, carbonyl, sulfonyl, nitro, hydroxyl_ether, ring_n"
        )

    heavy_idxs = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1}
    acceptor_idxs = set()

    # "auto" delegates acceptor detection to RDKit's built-in feature definitions.
    if mode == "auto" and not acceptor_smarts:
        fdef = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
        factory = ChemicalFeatures.BuildFeatureFactory(fdef)
        for feat in factory.GetFeaturesForMol(mol):
            if feat.GetFamily() != "Acceptor":
                continue
            for atom_idx in feat.GetAtomIds():
                if atom_idx in heavy_idxs:
                    acceptor_idxs.add(atom_idx)
        _debug(
            debug,
            f"triplet_search mode={mode} acceptor_atoms={len(acceptor_idxs)} sample={sorted(acceptor_idxs)[:10]}",
        )
    else:
        # Restrict by mode only when no explicit SMARTS is provided.
        if not acceptor_smarts:
            mode_smarts = {
                "generic": "[#7,#8;!$([N+]);!$([O+])]",
                # Includes both carbonyl oxygens (O=C) and carboxylate oxygens (O- C(=O)).
                "carbonyl": "[#8;!$([O+]);$([O]=[#6]),$([O-]-[#6](~[#8])),$([O]-[#6](~[#8]))]",
                "sulfonyl": "[#8;!$([O+]);$([O]-[#16](~[#8]))]",
                "nitro": "[#8;!$([O+]);$([O]-[#7+](~[#8]))]",
                "hydroxyl_ether": "[#8X2;!$([O]-[#6](~[#8]))]",
                # Ring nitrogens that can accept, aromatic or kekulized.
                "ring_n": "[#7;r;H0;!+]",
            }
            acceptor_smarts = mode_smarts[mode]
        _debug(debug, f"triplet_search mode={mode} smarts='{acceptor_smarts}'")

        acceptor_pattern = Chem.MolFromSmarts(acceptor_smarts)
        if acceptor_pattern is None:
            raise ValueError(f"Invalid AcceptorSMARTS pattern: '{acceptor_smarts}'")

        matches = mol.GetSubstructMatches(acceptor_pattern)
        _debug(
            debug,
            f"mode={mode} atoms={mol.GetNumAtoms()} heavy_atoms={len(heavy_idxs)} matches={len(matches)}",
        )
        for match in matches:
            if len(match) < 1:
                continue
            acc_idx = match[0]
            if acc_idx in heavy_idxs:
                acceptor_idxs.add(acc_idx)

    triplets = set()
    for acc_idx in sorted(acceptor_idxs):
        acc_atom = mol.GetAtomWithIdx(acc_idx)
        hub_idxs = [nbr.GetIdx() for nbr in acc_atom.GetNeighbors() if nbr.GetIdx() in heavy_idxs]
        for hub_idx in hub_idxs:
            hub_atom = mol.GetAtomWithIdx(hub_idx)
            plane_idxs = [
                nbr.GetIdx()
                for nbr in hub_atom.GetNeighbors()
                if nbr.GetIdx() in heavy_idxs and nbr.GetIdx() != acc_idx
            ]
            for plane_idx in plane_idxs:
                triplets.add((acc_idx, hub_idx, plane_idx))
                if include_reverse:
                    triplets.add((acc_idx, plane_idx, hub_idx))
    result = sorted(triplets)
    _debug(
        debug,
        f"mode={mode} atoms={mol.GetNumAtoms()} heavy_atoms={len(heavy_idxs)} "
        f"acceptor_atoms={len(acceptor_idxs)} triplets={len(result)}",
    )
    if debug and acceptor_idxs:
        _debug(debug, f"mode={mode} sample_acceptors={sorted(acceptor_idxs)[:10]}")
    if debug and result:
        _debug(debug, f"mode={mode} sample_triplets={result[:5]}")
    return result


def _map_idx_triplets_to_atom_names(idx_triplets, atom_names, debug=False, idx_to_name=None):
    mapped = []
    skipped_out_of_range = 0
    skipped_unmapped = 0
    skipped_hydrogen_name = 0
    for a, b, c in idx_triplets:
        if idx_to_name is not None and len(idx_to_name) > 0:
            if a not in idx_to_name or b not in idx_to_name or c not in idx_to_name:
                skipped_unmapped += 1
                continue
            names = (idx_to_name[a], idx_to_name[b], idx_to_name[c])
        else:
            if max(a, b, c) >= len(atom_names):
                skipped_out_of_range += 1
                continue
            names = (atom_names[a], atom_names[b], atom_names[c])
        if any(name.upper().startswith("H") for name in names):
            skipped_hydrogen_name += 1
            continue
        mapped.append(names)
    _debug(
        debug,
        "mapping "
        f"input={len(idx_triplets)} mapped={len(mapped)} "
        f"skipped_out_of_range={skipped_out_of_range} skipped_unmapped={skipped_unmapped} "
        f"skipped_hydrogen_name={skipped_hydrogen_name}",
    )
    return mapped


def generate_params_pdb_and_table(
    mtp,
    csv_file_name,
    path_to_conformers,
    molecule_sdfs,
    use_mol_id=False,
    no_name=False,
    dynamic_acceptor_alignment=False,
    target_atom_triplets=None,
    max_dynamic_alignments=0,
    include_reverse_neighbors=True,
    acceptor_smarts=None,
    acceptor_modes=None,
    dynamic_alignment_debug=False,
):
    """
    Used internally to generate csv files for reading in by yield_ligand_pose, params files,
    and pdb objects for conformers, will create directories inside of path_to_conformers 
    which are either mol_id/ or mol_name/

    Arguments:
    csv_file: File to create or append new entries to
    path_to_confomers: directory where ligand conformer and params file directories should be held
    molecule_sdfs: list of strings of molec Moule sdfs to read in and generate params file/pdbs for
    use_mol_id: boolean, if a mol id is present on the third line of each molecule in the sdf
        uses that to name directories and files
    no_name: boolean, if the sdf file has no discernable name then just name them sequentially

    Creates:
    csv_file with the following columns:
        Molecule Name, Molecule ID, Conformer Range, Molecule Atoms, Residue Atoms
        Molecule name: Name of the molecule as determined by first line of SDF
        Molecule ID: ID of the molecule as determiend by the third line of SDF (optional)
        Conformer Range: How many conformers were used in generating the params file
            e.g. 1-1 for 1 conformer, 1-100 for 100 conformers etc
            (Note this does assume that your conformers are in separate files)
        Molecule Atoms: MANUAL ENTRY atom labels for atoms on the conformer pdbs that correspond to target atoms
            e.g. C1-C4-C6 denotes three atoms C1, C4, C6 in that order
        Target Atoms: MANUAL ENTRY atom labels for target atoms on a Residue Object
            e.g. CD2-CZ2-CZ3 denote three atoms CD2, CZ2, CZ3 in that order on the target residue

    """
    try:
        os.mkdir(path_to_conformers)
    except:
        print(f"Directory {path_to_conformers} already made")
        pass


    if target_atom_triplets is None:
        target_atom_triplets = []

    if acceptor_modes is None:
        acceptor_modes = ["generic"]

    with open(csv_file_name, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Molecule Name", "Molecule ID", "Conformer Range", "Molecule Atoms", "Target Atoms"])
        for i, molecule_sdf in enumerate(molecule_sdfs):
            _debug(dynamic_alignment_debug, f"processing_sdf idx={i} path='{molecule_sdf}'")
            with open(molecule_sdf, "r", encoding='utf-8-sig') as sdf:
                lines = list(sdf)
                if len(lines) == 0:
                    continue
                mol_name = lines[0].strip()
                

                if use_mol_id:
                    mol_id = lines[2].split(" ")[0].strip()
                    file_stem = f"{mol_id}/{mol_id}"
                    dir_name = f"{path_to_conformers}/{mol_id}"
                elif not no_name:
                    mol_id = mol_name
                    file_stem =f"{mol_name}/{mol_name}"
                    dir_name = f"{path_to_conformers}/{mol_name}"

                else:
                    mol_name = i
                    mol_id = i
                    file_stem =f"{i}/{i}"
                    dir_name = f"{path_to_conformers}/{i}"

                try:
                    os.mkdir(dir_name)
                except:
                    print(f"Directory {dir_name} already made")
                    pass
                _debug(
                    dynamic_alignment_debug,
                    f"mol_name='{mol_name}' mol_id='{mol_id}' file_stem='{file_stem}'",
                )
                mtp.main([f"{molecule_sdf}", "-n", f"{path_to_conformers}/{file_stem}"])
                count = str(lines.count("$$$$\n"))
                conformer_range = f"1_{count}_"
                _debug(
                    dynamic_alignment_debug,
                    f"conformer_range='{conformer_range}' params_base='{path_to_conformers}/{file_stem}'",
                )

                if dynamic_acceptor_alignment and target_atom_triplets:
                    params_file = f"{path_to_conformers}/{file_stem}.params"
                    first_conformer_pdb = f"{path_to_conformers}/{file_stem}_0001.pdb"
                    atom_names = _read_params_atom_names(params_file)
                    _debug(
                        dynamic_alignment_debug,
                        f"params_file='{params_file}' atom_count={len(atom_names)} sample={atom_names[:10]}",
                    )
                    idx_to_name = _build_sdf_idx_to_pdb_name_map(
                        molecule_sdf,
                        first_conformer_pdb,
                        debug=dynamic_alignment_debug,
                    )
                    if len(idx_to_name) == 0:
                        _debug(
                            dynamic_alignment_debug,
                            "coord_map unavailable; falling back to params index mapping",
                        )
                    idx_triplets = []
                    for acceptor_mode in acceptor_modes:
                        mode_triplets = _get_rdkit_acceptor_triplets(
                            molecule_sdf,
                            include_reverse=include_reverse_neighbors,
                            acceptor_smarts=acceptor_smarts,
                            acceptor_mode=acceptor_mode,
                            debug=dynamic_alignment_debug,
                        )
                        _debug(
                            dynamic_alignment_debug,
                            f"mode='{acceptor_mode}' raw_triplets={len(mode_triplets)}",
                        )
                        idx_triplets.extend(mode_triplets)
                    idx_triplets = list(dict.fromkeys(idx_triplets))
                    _debug(
                        dynamic_alignment_debug,
                        f"deduped_idx_triplets={len(idx_triplets)} include_reverse={include_reverse_neighbors}",
                    )
                    mol_triplets = _map_idx_triplets_to_atom_names(
                        idx_triplets,
                        atom_names,
                        debug=dynamic_alignment_debug,
                        idx_to_name=idx_to_name,
                    )

                    if dynamic_alignment_debug:
                        atom_symbols = _read_sdf_atom_symbols(molecule_sdf)
                        print(
                            f"DEBUG {mol_id}: params_atoms={len(atom_names)} "
                            f"modes={','.join(acceptor_modes)} "
                            f"idx_triplets={len(idx_triplets)} mapped_triplets={len(mol_triplets)}"
                        )
                        if mol_triplets:
                            for idx_triplet, triplet in zip(idx_triplets[:10], mol_triplets[:10]):
                                idx_symbol_triplet = "-".join(
                                    [
                                        atom_symbols[idx] if idx < len(atom_symbols) else "?"
                                        for idx in idx_triplet
                                    ]
                                )
                                print(
                                    f"DEBUG {mol_id}: idx_triplet={idx_triplet} "
                                    f"idx_symbols={idx_symbol_triplet} "
                                    f"atom_triplet={'-'.join(triplet)}"
                                )

                    if max_dynamic_alignments > 0:
                        _debug(
                            dynamic_alignment_debug,
                            f"apply_max_dynamic_alignments={max_dynamic_alignments}",
                        )
                        mol_triplets = mol_triplets[:max_dynamic_alignments]

                    if len(mol_triplets) == 0:
                        _debug(dynamic_alignment_debug, "no_mapped_triplets writing blank row")
                        writer.writerow([mol_name, mol_id, conformer_range, "", ""])
                        continue

                    for mol_atoms in mol_triplets:
                        for tgt_atoms in target_atom_triplets:
                            writer.writerow(
                                [
                                    mol_name,
                                    mol_id,
                                    conformer_range,
                                    "-".join(mol_atoms),
                                    "-".join(tgt_atoms),
                                ]
                            )
                else:
                    writer.writerow([mol_name, mol_id, conformer_range, "", ""])


def main(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("config_file",
                        help = "your config file",
                        default = "my_conf.txt")
    
    if len(argv) == 0:
        print(parser.print_help())
        return

    args = parser.parse_args(argv)

    config = ConfigParser()
    with open(args.config_file, "r", encoding="utf-8-sig") as handle:
        config.read_file(handle)
    spec = config["create_table"]

    script_dir = os.path.dirname(os.path.abspath(__file__))
    legacy_dir = os.path.normpath(os.path.join(script_dir, "..", "legacy"))
    if legacy_dir not in sys.path:
        sys.path.append(legacy_dir)
    import molfile_to_params as mtp

    csv_file_name = _safe_get(config, "create_table", "CSVFileName")
    path_to_conformers = _safe_get(config, "create_table", "PathToConformers")
    if not csv_file_name:
        raise KeyError("Missing CSVFileName in [create_table] or [DEFAULT]")
    if not path_to_conformers:
        raise KeyError("Missing PathToConformers in [create_table] or [DEFAULT]")
    input_molecule_sdfs = spec["MoleculeSDFs"].split(" ")


    # Grabbing all of the Molecule SDFs
    molecule_sdfs = []
    for inp in input_molecule_sdfs:
        if "*" in inp:
            if "/" in inp:
                directory = "/".join([e for e in inp.split("/")[:-1]])
            else:
                directory = "."

            for file in os.listdir(directory):
                if fnmatch.fnmatch(file.lower(), inp.split("/")[-1]):
                    molecule_sdfs.append(f"{directory}/{file}")
        else:
            molecule_sdfs.append(inp)

    use_mol_id = _safe_getbool(config, "create_table", "UseMoleculeID", fallback=False)
    no_name = _safe_getbool(config, "create_table", "NoName", fallback=False)
    dynamic_acceptor_alignment = _safe_getbool(
        config, "create_table", "DynamicAcceptorAlignment", fallback=False
    )
    include_reverse_neighbors = _safe_getbool(
        config, "create_table", "IncludeReverseNeighborOrder", fallback=True
    )
    max_dynamic_alignments = int(
        _safe_get(config, "create_table", "MaxDynamicAlignments", fallback="0")
    )
    target_atom_triplets = _parse_target_atom_triplets(
        _safe_get(config, "create_table", "TargetAtomTriplets", fallback="")
    )
    acceptor_smarts = _safe_get(config, "create_table", "AcceptorSMARTS", fallback="")
    acceptor_modes = _parse_acceptor_modes(
        _safe_get(config, "create_table", "AcceptorMode", fallback="auto")
    )
    dynamic_alignment_debug = _safe_getbool(
        config, "create_table", "DynamicAlignmentDebug", fallback=False
    )

    if dynamic_acceptor_alignment and len(target_atom_triplets) == 0:
        raise ValueError(
            "DynamicAcceptorAlignment=True but no TargetAtomTriplets were provided."
        )

    if dynamic_acceptor_alignment:
        try:
            import rdkit  # noqa: F401
        except ImportError as exc:
            raise ImportError(
                "RDKit is required when DynamicAcceptorAlignment=True"
            ) from exc
    _debug(
        dynamic_alignment_debug,
        "config "
        f"csv='{csv_file_name}' conformers='{path_to_conformers}' "
        f"dynamic_acceptor_alignment={dynamic_acceptor_alignment} "
        f"acceptor_modes={acceptor_modes} acceptor_smarts='{acceptor_smarts}' "
        f"include_reverse_neighbors={include_reverse_neighbors} "
        f"max_dynamic_alignments={max_dynamic_alignments} "
        f"target_atom_triplets={target_atom_triplets}",
    )
    _debug(
        dynamic_alignment_debug,
        f"input_sdfs_count={len(molecule_sdfs)} sample={molecule_sdfs[:5]}",
    )

    generate_params_pdb_and_table(
        mtp,
        csv_file_name,
        path_to_conformers,
        molecule_sdfs,
        use_mol_id,
        no_name,
        dynamic_acceptor_alignment=dynamic_acceptor_alignment,
        target_atom_triplets=target_atom_triplets,
        max_dynamic_alignments=max_dynamic_alignments,
        include_reverse_neighbors=include_reverse_neighbors,
        acceptor_smarts=acceptor_smarts,
        acceptor_modes=acceptor_modes,
        dynamic_alignment_debug=dynamic_alignment_debug,
    )
    print(f"Succesfully generated table at {csv_file_name} and conformers at {path_to_conformers}")
    
    
if __name__ == '__main__':
    main(sys.argv[1:])
