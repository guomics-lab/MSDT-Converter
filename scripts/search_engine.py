import logging
import os
import shutil
import subprocess
from pathlib import Path
import json
from dataclasses import dataclass

# Configure logger (accessible by caller)
logger = logging.getLogger(__name__)
sage_script = './linux_sage'


def generate_sage_search_result_fn(param):
    """
    Generate a rawspectrum file.
    Return values:
        0: Successfully generated
        1: Already exists, no need to generate
        2: Input file does not exist
        -1: Generation failed
    """
    workdir = param.get('workdir')
    fasta = param.get('fasta')
    data_path = param.get('data_path')
    config_path = param.get('config_path')
    if data_path.endswith('.mzML'):
        fn = data_path.split('/')[-1][:-5]
    elif data_path.endswith('.d'):
        fn = data_path.split('/')[-1][:-2]
    output_path = os.path.join(workdir, fn + '_search_result.tsv')

    if not os.path.exists(config_path):
        logger.error(f"Input file does not exist: {config_path}")
        return 2

    if os.path.exists(output_path):
        logger.info(f"Output already exists, skipping: {output_path}")
        return 1

    # change sage config
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except json.JSONDecodeError:
        logger.error(f"ERROR: {config_path} is not json format")
        return -1
    except Exception as e:
        logger.error(f"ERROR: wrong occurs when reading {e}")
        return -1

    data['database']['fasta'] = fasta
    data['mzml_paths'] = [data_path]
    data['output_directory'] = workdir

    with open(config_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)

    # Create the output directory if it does not exist
    new_run_exe_path = os.path.join(workdir, 'linux_sage')
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(sage_script, new_run_exe_path)

    cmd = [new_run_exe_path, config_path]
    try:
        logger.info(f"Generating: {output_path}")
        logger.info(f"Running command: {' '.join(cmd)}")
        env = os.environ.copy()
        env['RUST_MIN_STACK'] = '8388608'
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
        result_sage_file_path = os.path.join(workdir, 'results.sage.tsv')
        shutil.move(result_sage_file_path, output_path)
        logger.info(f"Successfully generated: {output_path}")
        os.remove(new_run_exe_path)
        return 0
    except subprocess.CalledProcessError as e:
        logger.error(f"Generation failed: {e.stderr}")
        return -1
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        return -1


fragpipe_exe_path = os.path.join(os.getcwd(), 'FragPipe-21.1', 'bin', 'fragpipe')


@dataclass(frozen=True)
class ManifestEntry:
    path: Path
    experiment: str = "exp"
    bioreplicate: str = ""
    data_type: str = "DDA"


def read_file_list(file_list_path):
    """Read one-path-per-line or official four-column FragPipe input files."""
    source = Path(file_list_path)
    if not source.is_file():
        raise FileNotFoundError(f"FragPipe file_list not found: {source}")
    entries = []
    seen_paths = set()
    for line_number, raw_line in enumerate(
        source.read_text(encoding="utf-8-sig").splitlines(), start=1
    ):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = raw_line.rstrip("\r\n").split("\t")
        if fields[0].strip().lower() in {"path", "file_path", "data_path"}:
            continue
        if len(fields) == 1:
            entry = ManifestEntry(Path(fields[0].strip()))
        elif len(fields) == 4:
            entry = ManifestEntry(
                Path(fields[0].strip()),
                fields[1].strip(),
                fields[2].strip(),
                fields[3].strip(),
            )
        else:
            raise ValueError(
                f"file_list line {line_number} must have 1 or 4 tab-separated columns"
            )
        if not entry.path.exists():
            raise FileNotFoundError(
                f"Input from file_list does not exist: {entry.path}"
            )
        resolved = entry.path.resolve()
        if resolved in seen_paths:
            raise ValueError(f"Duplicate input in file_list: {entry.path}")
        seen_paths.add(resolved)
        entries.append(entry)
    if not entries:
        raise ValueError(f"file_list contains no inputs: {source}")
    return entries


def build_manifest(file_path, fragpipe_output_path):
    if isinstance(file_path, ManifestEntry):
        entries = [file_path]
    elif isinstance(file_path, (str, os.PathLike)):
        entries = [ManifestEntry(Path(file_path))]
    else:
        entries = list(file_path)
    if not entries:
        raise ValueError("Cannot build an empty FragPipe manifest")

    Path(fragpipe_output_path).mkdir(parents=True, exist_ok=True)
    manifest_path = os.path.join(fragpipe_output_path, 'fragpipe-files.fp-manifest')
    logger.info(f'Processing build fragpipe manifest, path is {manifest_path}')
    with open(manifest_path, 'w', encoding='utf-8') as f:
        for entry in entries:
            f.write(
                f'{entry.path}\t{entry.experiment}\t{entry.bioreplicate}\t{entry.data_type}\n'
            )
    logger.info(f'Finished build fragpipe manifest')
    return Path(manifest_path)


def prepare_workflow(workflow_path, destination_path, database_path):
    """Copy a workflow and force the v2 database/Percolator settings."""
    source = Path(workflow_path)
    if not source.is_file():
        raise FileNotFoundError(f"FragPipe workflow not found: {source}")
    destination = Path(destination_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    settings = {
        "database.db-path": str(database_path),
        "percolator.keep-tsv-files": "true",
        "percolator.run-percolator": "true",
    }
    output_lines = []
    found = set()
    for line in source.read_text(encoding="utf-8").splitlines():
        key = line.split("=", 1)[0] if "=" in line else None
        if key in settings:
            output_lines.append(f"{key}={settings[key]}")
            found.add(key)
        else:
            output_lines.append(line)
    for key, value in settings.items():
        if key not in found:
            output_lines.append(f"{key}={value}")
    destination.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return destination


def find_fragpipe_outputs(entry, workdir):
    """Find one complete FragPipe PSM result set for a manifest entry."""
    result_group = entry.experiment
    if entry.bioreplicate:
        result_group = f"{result_group}_{entry.bioreplicate}"
    experiment_root = (
        Path(workdir) / result_group if result_group else Path(workdir)
    )
    if not experiment_root.is_dir():
        return None

    edited_pins = list(experiment_root.rglob(f"{entry.path.stem}_edited.pin"))
    pin_candidates = edited_pins or list(
        experiment_root.rglob(f"{entry.path.stem}.pin")
    )
    result_sets = []
    for pin in pin_candidates:
        target = pin.parent / f"{entry.path.stem}_percolator_target_psms.tsv"
        decoy = pin.parent / f"{entry.path.stem}_percolator_decoy_psms.tsv"
        # FragPipe/Percolator creates output files before it has finished
        # writing them.  A killed run can therefore leave a PIN and two
        # zero-byte TSVs that must not be treated as a reusable result set.
        outputs = (pin, target, decoy)
        if all(path.is_file() and path.stat().st_size > 0 for path in outputs):
            result_sets.append({"pin": pin, "target": target, "decoy": decoy})
    if len(result_sets) > 1:
        raise RuntimeError(
            f"Multiple complete FragPipe result sets found for {entry.path}: "
            + ", ".join(str(result["pin"].parent) for result in result_sets)
        )
    return result_sets[0] if result_sets else None


def run_cmd(cmd, cwd=None):
    env = os.environ.copy()
    java_bin_path = os.path.join(os.getcwd(), 'jdk-11.0.26', 'bin')
    env['PATH'] = java_bin_path + os.pathsep + env['PATH']
    java_home = os.path.join(os.getcwd(), 'jdk-11.0.26')
    env['JAVA_HOME'] = java_home
    cmd_str = ' '.join(cmd)
    logger.info(f'Run cmd: {cmd_str}')
    process = subprocess.Popen(
        cmd,
        env=env,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding='utf-8',
        errors='replace',
    )
    if process.stdout is not None:
        for line in process.stdout:
            message = line.rstrip()
            if message:
                logger.info(message)
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, cmd)


def run_fragpipe(manifest_path, workflow_path, fragpipe_output_path, exe_abs_path, thread_num, fasta_path):
    frag_base_dir = os.path.dirname(os.path.dirname(exe_abs_path))

    ion_quant_exe_path = os.path.join(frag_base_dir, 'IonQuant-1.10.27', 'IonQuant-1.10.27.jar')
    msfrag_exe_path = os.path.join(frag_base_dir, 'MSFragger-4.0', 'MSFragger-4.0.jar')
    philosopher_exe_path = os.path.join(frag_base_dir, 'philosopher-v5.1.1')

    # Philosopher
    fasta_dir = os.path.dirname(fasta_path)
    if not os.path.exists(fasta_dir):
        logger.error(f"Fasta file does not exist: {fasta_path}")

    decoy_pattern = '*.fasta.fas'
    previous_decoys = {
        path.resolve(): (path.stat().st_mtime_ns, path.stat().st_size)
        for path in Path(fasta_dir).glob(decoy_pattern)
    }

    philosopher_cmd1 = [philosopher_exe_path, 'workspace', '--clean', '--nocheck']
    subprocess.run(philosopher_cmd1, cwd=fasta_dir, check=True)
    
    philosopher_cmd2 = [philosopher_exe_path, 'workspace', '--init', '--nocheck']
    subprocess.run(philosopher_cmd2, cwd=fasta_dir, check=True)
    
    philosopher_cmd3 = [philosopher_exe_path, 'database', '--custom', fasta_path]
    subprocess.run(philosopher_cmd3, cwd=fasta_dir, check=True)
    
    philosopher_cmd4 = [philosopher_exe_path, 'workspace', '--clean', '--nocheck']
    subprocess.run(philosopher_cmd4, cwd=fasta_dir, check=True)

    # get decoyfasta
    current_decoys = list(Path(fasta_dir).glob(decoy_pattern))
    generated_decoys = [
        path
        for path in current_decoys
        if path.resolve() not in previous_decoys
        or (path.stat().st_mtime_ns, path.stat().st_size)
        != previous_decoys[path.resolve()]
    ]
    if not generated_decoys:
        raise FileNotFoundError(
            f"Philosopher did not create or update a decoy FASTA in {fasta_dir}"
        )
    if len(generated_decoys) != 1:
        raise RuntimeError(
            "Philosopher created or updated multiple decoy FASTA files; "
            "cannot choose safely: "
            + ", ".join(str(path) for path in generated_decoys)
        )
    decoyfasta = generated_decoys[0]
    
    workflow_dest = os.path.join(fragpipe_output_path, os.path.basename(workflow_path))
    prepare_workflow(workflow_path, workflow_dest, decoyfasta)
    logger.info(f'new workflow file has been generated: {workflow_dest}')
    logger.info(f'fasta with decoy has been generated: {decoyfasta}')
    
    cmd = [exe_abs_path, '--headless', '--workflow', workflow_dest, '--manifest',
           str(manifest_path), '--workdir', fragpipe_output_path,
           '--config-ionQuant', ion_quant_exe_path, '--config-msfragger', msfrag_exe_path,
           '--config-philosopher', philosopher_exe_path, '--threads', str(thread_num)]
    logger.info(f'Processing run fragpipe, command is {cmd}')
    run_cmd(cmd)
    logger.info(f'Finished run fragpipe, result path is {fragpipe_output_path}')


def generate_fp_search_result_fn(param):
    """
    Return values:
        0: Successfully generated
        1: Already exists, no need to generate
        2: Input file does not exist
        -1: Generation failed
    """
    try:
        file_list = param.get('file_list')
        data_path = param.get('data_path')
        if file_list:
            entries = read_file_list(file_list)
        elif data_path:
            input_path = Path(data_path)
            if not input_path.exists():
                logger.error(f"Input file does not exist: {input_path}")
                return 2
            entries = [ManifestEntry(input_path)]
        else:
            raise ValueError("FragPipe search requires data_path or file_list")

        thread_num = int(param.get('thread_num', 1))
        fasta_path = param.get('fasta_path')
        workflow_path = param.get('workflow_path')
        workdir = Path(param.get('workdir'))
        workdir.mkdir(parents=True, exist_ok=True)

        if all(find_fragpipe_outputs(entry, workdir) is not None for entry in entries):
            logger.info("All FragPipe PSM outputs already exist, skipping search")
            return 1

        manifest_path = build_manifest(entries, workdir)
        run_fragpipe(
            manifest_path,
            workflow_path,
            str(workdir),
            fragpipe_exe_path,
            thread_num,
            fasta_path,
        )
        missing = [
            str(entry.path)
            for entry in entries
            if find_fragpipe_outputs(entry, workdir) is None
        ]
        if missing:
            logger.error(
                "FragPipe completed without complete PIN/Percolator outputs for: "
                + ", ".join(missing)
            )
            return -1
        logger.info("generate FragPipe search results success")
        return 0
    except Exception as error:
        logger.exception(f"FragPipe search failed: {error}")
        return -1
