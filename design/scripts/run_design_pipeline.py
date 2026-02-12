#!/usr/bin/env python3
"""
Design Pipeline Orchestrator: Docking → MPNN → Rosetta → Filter → AF3 Prep

This script automates the design workflow that was previously manual:
1. LigandMPNN design on docked structures
2. Rosetta relax of MPNN sequences
3. Score aggregation
4. Filtering by Rosetta metrics
5. (Optional) Iteration (repeat 1-4)
6. FASTA generation
7. AF3 JSON preparation with SMILES extraction

Usage:
    # Full pipeline (docking must already be complete)
    python run_design_pipeline.py config.txt

    # Skip to specific stage
    python run_design_pipeline.py config.txt --skip-mpnn
    python run_design_pipeline.py config.txt --skip-rosetta

    # Run only AF3 prep
    python run_design_pipeline.py config.txt --af3-prep-only

    # Specify iteration round
    python run_design_pipeline.py config.txt --iteration 2
"""

import argparse
import os
import sys
import subprocess
import json
import shutil
from configparser import ConfigParser
from pathlib import Path

# Add docking scripts to path for utility imports
SCRIPT_DIR = Path(__file__).parent.absolute()
PIPELINE_ROOT = SCRIPT_DIR.parent.parent
DOCKING_SCRIPTS = PIPELINE_ROOT / "docking" / "scripts"
sys.path.insert(0, str(DOCKING_SCRIPTS))

try:
    import docking_pipeline_utils as dpu
except ImportError:
    print("Warning: Could not import docking_pipeline_utils")
    dpu = None


class DesignPipelineConfig:
    """Configuration holder for design pipeline"""

    def __init__(self, config_file):
        self.config = ConfigParser()
        with open(config_file, 'r', encoding='utf-8-sig') as f:
            self.config.read_file(f)
        self.config_file = config_file

        # Core paths
        self.pipe_root = Path(self._get('DEFAULT', 'PIPE_ROOT', PIPELINE_ROOT))
        self.campaign_root = Path(self._get('DEFAULT', 'CAMPAIGN_ROOT'))
        self.scratch_root = Path(self._get('DEFAULT', 'SCRATCH_ROOT'))

        # Docking output (input to design)
        self.docked_dir = Path(self._get('grade_conformers', 'OutputDir'))
        self.clustered_dir = Path(self._get('grade_conformers', 'ClusteredOutputDir',
                                           self.docked_dir / 'clustered_final'))

        # Design settings
        self.design_root = self.scratch_root / self._get('design', 'DesignRoot', 'design')
        self.iterations = int(self._get('design', 'DesignIterationRounds', '1'))

        # MPNN settings
        self.mpnn_script = self.pipe_root / self._get('design', 'MPNNScript',
            'design/instructions/ligand_alignment_mpnni_grouped.sh')
        self.mpnn_batch_size = int(self._get('design', 'MPNNBatchSize', '40'))
        self.mpnn_temperature = float(self._get('design', 'MPNNTemperature', '0.3'))
        self.mpnn_omit_file = self._get('design', 'MPNNOmitDesignFile', '')
        self.mpnn_bias_file = self._get('design', 'MPNNBiasFile', '')
        self.design_residues = self._get('design', 'DesignResidues',
            '59 79 81 90 92 106 108 115 118 120 139 157 158 161 162 165')

        # Rosetta settings
        self.rosetta_script = self.pipe_root / self._get('design', 'RosettaScript',
            'design/instructions/submit_pyrosetta_general_threading_relax.sh')
        self.rosetta_relax_py = self.pipe_root / self._get('design', 'RosettaRelaxPy',
            'design/instructions/general_relax.py')
        self.ligand_params = Path(self._get('design', 'LigandParams'))

        # Filtering settings
        self.filter_script = self.pipe_root / self._get('design', 'FilterScript',
            'design/instructions/relax_2_filter__allpolar_unsats.py')
        self.filter_target_n = int(self._get('design', 'FilterTargetN', '1000'))
        self.filter_max_unsat = int(self._get('design', 'FilterMaxUnsats', '1'))
        self.filter_max_per_parent = int(self._get('design', 'FilterMaxPerParent', '20'))
        self.filter_ignore_o1 = self._get('design', 'FilterIgnoreO1', 'False').lower() == 'true'
        self.filter_ignore_o2 = self._get('design', 'FilterIgnoreO2', 'False').lower() == 'true'
        self.filter_ignore_charge = self._get('design', 'FilterIgnoreCharge', 'False').lower() == 'true'

        # AF3 settings
        self.af3_binary_template = Path(self._get('design', 'AF3BinaryTemplate',
            self.pipe_root / 'design/templates/pyr1_binary_template.json'))
        self.af3_ternary_template = Path(self._get('design', 'AF3TernaryTemplate',
            self.pipe_root / 'design/templates/pyr1_ternary_template.json'))
        self.ligand_sdf = Path(self._get('design', 'LigandSDF', ''))

        # Helper scripts
        self.aggregate_scores = self.pipe_root / self._get('design', 'AggregateScoresScript',
            'design/instructions/aggregate_scores.py')
        self.split_to_fasta = self.pipe_root / self._get('design', 'SplitToFastaScript',
            'design/instructions/split_and_mutate_to_fasta.py')
        self.make_af3_jsons = self.pipe_root / self._get('design', 'MakeAF3JSONsScript',
            'design/instructions/make_af3_jsons.py')

    def _get(self, section, key, default=None):
        """Safely get config value"""
        if dpu:
            return dpu.cfg_get(self.config, section, key, default)
        else:
            try:
                return self.config.get(section, key)
            except:
                if default is not None:
                    return default
                raise

    def get_iteration_dir(self, iteration_num):
        """Get directory for specific iteration"""
        return self.design_root / f"iteration_{iteration_num}"

    def get_mpnn_output_dir(self, iteration_num):
        """Get MPNN output directory for iteration"""
        return self.get_iteration_dir(iteration_num) / "mpnn_output"

    def get_rosetta_output_dir(self, iteration_num):
        """Get Rosetta output directory for iteration"""
        return self.get_iteration_dir(iteration_num) / "rosetta_output"

    def get_scores_dir(self, iteration_num):
        """Get scores directory for iteration"""
        return self.get_iteration_dir(iteration_num) / "scores"

    def get_filtered_dir(self, iteration_num):
        """Get filtered results directory for iteration"""
        return self.get_iteration_dir(iteration_num) / "filtered"

    def get_input_pdbs_dir(self, iteration_num):
        """Get input PDB directory for iteration (docked or previous filtered)"""
        if iteration_num == 1:
            return self.clustered_dir
        else:
            return self.get_filtered_dir(iteration_num - 1)


def extract_smiles_from_sdf(sdf_path):
    """Extract SMILES from SDF file using RDKit"""
    try:
        from rdkit import Chem
        supplier = Chem.SDMolSupplier(str(sdf_path))
        mol = next(iter(supplier))
        if mol:
            return Chem.MolToSmiles(mol)
    except ImportError:
        print("Warning: RDKit not available, cannot extract SMILES from SDF")
    except Exception as e:
        print(f"Warning: Could not extract SMILES from {sdf_path}: {e}")
    return None


def update_json_smiles(template_path, output_path, smiles):
    """Update SMILES in AF3 JSON template"""
    with open(template_path, 'r') as f:
        data = json.load(f)

    # Find and update ligand SMILES
    # Typically at sequences[N]["ligand"]["smiles"]
    updated = False
    if 'sequences' in data:
        for seq in data['sequences']:
            if 'ligand' in seq and 'smiles' in seq['ligand']:
                seq['ligand']['smiles'] = smiles
                updated = True
                break

    if not updated:
        print(f"Warning: Could not find SMILES field in {template_path}")

    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)

    return updated


def run_mpnn_design(cfg, iteration_num, submit_slurm=True):
    """Run LigandMPNN design on input PDBs"""
    print("\n" + "=" * 80)
    print(f"STAGE 1: LigandMPNN Design (Iteration {iteration_num})")
    print("=" * 80)

    input_dir = cfg.get_input_pdbs_dir(iteration_num)
    output_dir = cfg.get_mpnn_output_dir(iteration_num)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input PDBs: {input_dir}")
    print(f"Output dir: {output_dir}")

    # Count input PDBs
    pdb_count = len(list(input_dir.glob("*.pdb")))
    print(f"Found {pdb_count} input PDBs")

    if pdb_count == 0:
        print("ERROR: No PDB files found in input directory")
        return False

    # Generate MPNN submit script with updated paths
    mpnn_script_custom = output_dir / "submit_mpnn.sh"

    with open(cfg.mpnn_script, 'r') as f:
        script_content = f.read()

    # Replace directory paths
    script_content = script_content.replace(
        'PDB_DIR="/scratch/alpine/ryde3462/',
        f'PDB_DIR="{input_dir}"\n# ORIGINAL: PDB_DIR="/scratch/alpine/ryde3462/'
    )
    script_content = script_content.replace(
        'OUTPUT_BASE="/scratch/alpine/ryde3462/',
        f'OUTPUT_BASE="{output_dir}"\n# ORIGINAL: OUTPUT_BASE="/scratch/alpine/ryde3462/'
    )

    # Update array count based on number of PDBs
    script_content = script_content.replace(
        '#SBATCH --array=1-',
        f'#SBATCH --array=1-{pdb_count}'
    )

    with open(mpnn_script_custom, 'w') as f:
        f.write(script_content)

    print(f"Generated custom MPNN script: {mpnn_script_custom}")

    if submit_slurm:
        # Submit to SLURM
        cmd = ['sbatch', str(mpnn_script_custom)]
        print(f"Submitting: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.returncode != 0:
            print(f"ERROR: {result.stderr}")
            return False

        # Extract job ID
        job_id = result.stdout.strip().split()[-1]
        print(f"✓ MPNN job submitted: {job_id}")
        return job_id
    else:
        print(f"SLURM submission skipped (dry-run mode)")
        print(f"To submit manually: sbatch {mpnn_script_custom}")
        return True


def run_rosetta_relax(cfg, iteration_num, submit_slurm=True):
    """Run Rosetta relax on MPNN designs"""
    print("\n" + "=" * 80)
    print(f"STAGE 2: Rosetta Relax (Iteration {iteration_num})")
    print("=" * 80)

    template_dir = cfg.get_input_pdbs_dir(iteration_num)
    mpnn_output = cfg.get_mpnn_output_dir(iteration_num)
    rosetta_output = cfg.get_rosetta_output_dir(iteration_num)
    rosetta_output.mkdir(parents=True, exist_ok=True)

    print(f"Template PDBs: {template_dir}")
    print(f"MPNN outputs: {mpnn_output}")
    print(f"Rosetta output: {rosetta_output}")

    # Count FASTA files
    fasta_count = len(list(mpnn_output.rglob("*.fa")))
    print(f"Found {fasta_count} FASTA files")

    if fasta_count == 0:
        print("ERROR: No FASTA files found in MPNN output")
        return False

    # Generate Rosetta submit script with updated paths
    rosetta_script_custom = rosetta_output / "submit_rosetta.sh"

    with open(cfg.rosetta_script, 'r') as f:
        script_content = f.read()

    # Replace directory paths
    script_content = script_content.replace(
        'TEMPLATE_DIR="/scratch/alpine/ryde3462/',
        f'TEMPLATE_DIR="{template_dir}"\n# ORIGINAL: TEMPLATE_DIR="/scratch/alpine/ryde3462/'
    )
    script_content = script_content.replace(
        'MPNN_OUTPUT_BASE="/scratch/alpine/ryde3462/',
        f'MPNN_OUTPUT_BASE="{mpnn_output}"\n# ORIGINAL: MPNN_OUTPUT_BASE="/scratch/alpine/ryde3462/'
    )
    script_content = script_content.replace(
        'OUTPUT_DIR="/scratch/alpine/ryde3462/',
        f'OUTPUT_DIR="{rosetta_output}"\n# ORIGINAL: OUTPUT_DIR="/scratch/alpine/ryde3462/'
    )

    # Update ligand params path
    script_content = script_content.replace(
        'LIGAND_PARAMS="/projects/ryde3462/',
        f'LIGAND_PARAMS="{cfg.ligand_params}"\n# ORIGINAL: LIGAND_PARAMS="/projects/ryde3462/'
    )

    # Update array count
    script_content = script_content.replace(
        '#SBATCH --array=1-500',
        f'#SBATCH --array=1-{fasta_count}'
    )

    with open(rosetta_script_custom, 'w') as f:
        f.write(script_content)

    print(f"Generated custom Rosetta script: {rosetta_script_custom}")

    if submit_slurm:
        # Submit to SLURM
        cmd = ['sbatch', str(rosetta_script_custom)]
        print(f"Submitting: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.returncode != 0:
            print(f"ERROR: {result.stderr}")
            return False

        job_id = result.stdout.strip().split()[-1]
        print(f"✓ Rosetta job submitted: {job_id}")
        return job_id
    else:
        print(f"SLURM submission skipped (dry-run mode)")
        print(f"To submit manually: sbatch {rosetta_script_custom}")
        return True


def aggregate_rosetta_scores(cfg, iteration_num):
    """Aggregate Rosetta .sc files into CSV"""
    print("\n" + "=" * 80)
    print(f"STAGE 3: Aggregate Scores (Iteration {iteration_num})")
    print("=" * 80)

    rosetta_output = cfg.get_rosetta_output_dir(iteration_num)
    scores_dir = cfg.get_scores_dir(iteration_num)
    scores_dir.mkdir(parents=True, exist_ok=True)

    output_csv = scores_dir / f"iteration_{iteration_num}_scores.csv"

    print(f"Scanning: {rosetta_output}")
    print(f"Output CSV: {output_csv}")

    # Run aggregate script
    cmd = [
        sys.executable,
        str(cfg.aggregate_scores),
        str(rosetta_output),
        '--output', str(output_csv)
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True)

    print(f"✓ Scores aggregated to {output_csv}")
    return output_csv


def filter_designs(cfg, iteration_num, scores_csv):
    """Filter designs by Rosetta metrics"""
    print("\n" + "=" * 80)
    print(f"STAGE 4: Filter Designs (Iteration {iteration_num})")
    print("=" * 80)

    rosetta_output = cfg.get_rosetta_output_dir(iteration_num)
    filtered_dir = cfg.get_filtered_dir(iteration_num)
    filtered_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input CSV: {scores_csv}")
    print(f"Rosetta PDBs: {rosetta_output}")
    print(f"Filtered output: {filtered_dir}")

    cmd = [
        sys.executable,
        str(cfg.filter_script),
        str(scores_csv),
        str(rosetta_output),
        str(filtered_dir),
        '--target_n', str(cfg.filter_target_n),
        '--max_unsat', str(cfg.filter_max_unsat),
        '--max_per_parent', str(cfg.filter_max_per_parent),
        '--output_csv_name', 'filtered.csv'
    ]

    if cfg.filter_ignore_o1:
        cmd.append('--ignore_o1')
    if cfg.filter_ignore_o2:
        cmd.append('--ignore_o2')
    if cfg.filter_ignore_charge:
        cmd.append('--ignore_charge')

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True)

    filtered_csv = filtered_dir / 'filtered.csv'
    print(f"✓ Filtered designs saved to {filtered_dir}")

    # Count filtered PDBs
    pdb_count = len(list(filtered_dir.glob("*.pdb")))
    print(f"  {pdb_count} PDBs passed filters")

    return filtered_csv


def generate_fasta(cfg, final_filtered_dir):
    """Generate FASTA from filtered PDBs"""
    print("\n" + "=" * 80)
    print("STAGE 5: Generate FASTA")
    print("=" * 80)

    output_fasta = final_filtered_dir / "filtered.fasta"

    print(f"Input dir: {final_filtered_dir}")
    print(f"Output FASTA: {output_fasta}")

    # Check if script exists
    if not cfg.split_to_fasta.exists():
        print(f"Warning: Script not found: {cfg.split_to_fasta}")
        print("Skipping FASTA generation")
        return None

    cmd = [
        sys.executable,
        str(cfg.split_to_fasta),
        str(final_filtered_dir),
        str(output_fasta)
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True)

    print(f"✓ FASTA generated: {output_fasta}")
    return output_fasta


def prepare_af3_inputs(cfg, fasta_file, smiles=None):
    """Prepare AF3 JSON inputs (binary and ternary)"""
    print("\n" + "=" * 80)
    print("STAGE 6: Prepare AF3 Inputs")
    print("=" * 80)

    af3_root = cfg.design_root / "af3_inputs"
    binary_dir = af3_root / "binary"
    ternary_dir = af3_root / "ternary"

    binary_dir.mkdir(parents=True, exist_ok=True)
    ternary_dir.mkdir(parents=True, exist_ok=True)

    # Extract SMILES if not provided
    if smiles is None and cfg.ligand_sdf.exists():
        print(f"Extracting SMILES from {cfg.ligand_sdf}")
        smiles = extract_smiles_from_sdf(cfg.ligand_sdf)
        if smiles:
            print(f"  SMILES: {smiles}")

    # Prepare templates with SMILES
    binary_template = cfg.af3_binary_template
    ternary_template = cfg.af3_ternary_template

    if smiles:
        # Create updated templates
        binary_template_updated = binary_dir / "template_with_smiles.json"
        ternary_template_updated = ternary_dir / "template_with_smiles.json"

        update_json_smiles(binary_template, binary_template_updated, smiles)
        update_json_smiles(ternary_template, ternary_template_updated, smiles)

        binary_template = binary_template_updated
        ternary_template = ternary_template_updated
        print(f"✓ Updated templates with SMILES")

    # Generate binary JSONs
    print(f"\nGenerating binary JSONs...")
    cmd = [
        sys.executable,
        str(cfg.make_af3_jsons),
        '--template', str(binary_template),
        '--fasta', str(fasta_file),
        '--outdir', str(binary_dir)
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    binary_count = len(list(binary_dir.glob("*.json"))) - 1  # Exclude template
    print(f"✓ Generated {binary_count} binary JSONs in {binary_dir}")

    # Generate ternary JSONs
    print(f"\nGenerating ternary JSONs...")
    cmd = [
        sys.executable,
        str(cfg.make_af3_jsons),
        '--template', str(ternary_template),
        '--fasta', str(fasta_file),
        '--outdir', str(ternary_dir)
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    ternary_count = len(list(ternary_dir.glob("*.json"))) - 1
    print(f"✓ Generated {ternary_count} ternary JSONs in {ternary_dir}")

    return binary_dir, ternary_dir


def run_full_pipeline(cfg, args):
    """Run complete design pipeline"""

    # Determine iteration range
    start_iter = args.iteration if args.iteration else 1
    end_iter = start_iter if args.iteration else cfg.iterations

    print("\n" + "=" * 80)
    print("DESIGN PIPELINE ORCHESTRATOR")
    print("=" * 80)
    print(f"Config: {cfg.config_file}")
    print(f"Design root: {cfg.design_root}")
    print(f"Iterations: {start_iter} to {end_iter}")
    print("=" * 80)

    # Track outputs across iterations
    final_filtered_dir = None

    for iter_num in range(start_iter, end_iter + 1):
        print(f"\n{'=' * 80}")
        print(f"ITERATION {iter_num} / {end_iter}")
        print(f"{'=' * 80}")

        # Stage 1: MPNN Design
        if not args.skip_mpnn and not args.af3_prep_only:
            job_id = run_mpnn_design(cfg, iter_num, submit_slurm=not args.dry_run)
            if not job_id:
                print("ERROR: MPNN stage failed")
                return False

            if not args.dry_run and args.wait:
                print("\nWaiting for MPNN job to complete...")
                print("(Press Ctrl+C to continue without waiting)")
                try:
                    subprocess.run(['squeue', '-j', str(job_id)], check=False)
                except KeyboardInterrupt:
                    print("\nContinuing without waiting...")

        # Stage 2: Rosetta Relax
        if not args.skip_rosetta and not args.af3_prep_only:
            job_id = run_rosetta_relax(cfg, iter_num, submit_slurm=not args.dry_run)
            if not job_id:
                print("ERROR: Rosetta stage failed")
                return False

            if not args.dry_run and args.wait:
                print("\nWaiting for Rosetta job to complete...")
                try:
                    subprocess.run(['squeue', '-j', str(job_id)], check=False)
                except KeyboardInterrupt:
                    print("\nContinuing without waiting...")

        # Stage 3: Aggregate Scores
        if not args.skip_aggregate and not args.af3_prep_only:
            if args.dry_run:
                print("\nStage 3: Aggregate Scores (skipped in dry-run)")
            else:
                scores_csv = aggregate_rosetta_scores(cfg, iter_num)

        # Stage 4: Filter
        if not args.skip_filter and not args.af3_prep_only:
            if args.dry_run:
                print("\nStage 4: Filter Designs (skipped in dry-run)")
                final_filtered_dir = cfg.get_filtered_dir(iter_num)
            else:
                filter_designs(cfg, iter_num, scores_csv)
                final_filtered_dir = cfg.get_filtered_dir(iter_num)

    # Stage 5 & 6: FASTA and AF3 prep (after all iterations)
    if not args.skip_af3_prep and final_filtered_dir:
        if args.dry_run:
            print("\nStage 5-6: FASTA & AF3 Prep (skipped in dry-run)")
        else:
            fasta_file = generate_fasta(cfg, final_filtered_dir)
            if fasta_file:
                prepare_af3_inputs(cfg, fasta_file)

    print("\n" + "=" * 80)
    print("PIPELINE COMPLETE")
    print("=" * 80)

    if final_filtered_dir:
        print(f"\nFinal designs: {final_filtered_dir}")
        if (final_filtered_dir / "filtered.fasta").exists():
            print(f"FASTA: {final_filtered_dir / 'filtered.fasta'}")
        print(f"AF3 inputs: {cfg.design_root / 'af3_inputs'}")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Design Pipeline Orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('config', help='Config file path')
    parser.add_argument('--iteration', type=int, help='Run specific iteration only')
    parser.add_argument('--skip-mpnn', action='store_true', help='Skip MPNN stage')
    parser.add_argument('--skip-rosetta', action='store_true', help='Skip Rosetta stage')
    parser.add_argument('--skip-aggregate', action='store_true', help='Skip score aggregation')
    parser.add_argument('--skip-filter', action='store_true', help='Skip filtering')
    parser.add_argument('--skip-af3-prep', action='store_true', help='Skip AF3 preparation')
    parser.add_argument('--af3-prep-only', action='store_true',
                       help='Only prepare AF3 inputs (skip MPNN/Rosetta)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Generate scripts but do not submit SLURM jobs')
    parser.add_argument('--wait', action='store_true',
                       help='Wait for SLURM jobs to complete before continuing')

    args = parser.parse_args()

    # Load config
    cfg = DesignPipelineConfig(args.config)

    # Run pipeline
    success = run_full_pipeline(cfg, args)

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
