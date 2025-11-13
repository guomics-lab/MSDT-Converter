import logging
import os
import shutil
import subprocess
from pathlib import Path
import json

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
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
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


def build_manifest(file_path, fragpipe_output_path):
    manifest_path = os.path.join(fragpipe_output_path, 'fragpipe-files.fp-manifest')
    logger.info(f'Processing build fragpipe manifest, path is {manifest_path}')
    #
    with open(manifest_path, 'w+') as f:
        f.write(f'{file_path}\texp\t\tDDA')
    logger.info(f'Finished build fragpipe manifest')
    return manifest_path


def run_cmd(cmd, cwd=None):
    env = os.environ.copy()
    java_bin_path = os.path.join(os.getcwd(), 'jdk-11.0.26', 'bin')
    env['PATH'] = java_bin_path + os.pathsep + env['PATH']
    java_home = os.path.join(os.getcwd(), 'jdk-11.0.26')
    env['JAVA_HOME'] = java_home
    cmd_str = ' '.join(cmd)
    logger.info(f'Run cmd: {cmd_str}')
    p = subprocess.Popen(cmd, env=env, cwd=cwd, stdout=subprocess.PIPE)
    stdout = p.stdout
    while True:
        output = stdout.readline()
        if output == b'' or (output == '' and p.poll() is not None):
            break
        logger.info(output)
        if output:
            #
            info_msg = output.decode('utf-8')
            info_msg = info_msg.rstrip()
            if len(info_msg) == 0:
                continue


def run_fragpipe(manifest_path, workflow_path, fragpipe_output_path, exe_abs_path, thread_num):
    frag_base_dir = os.path.dirname(os.path.dirname(exe_abs_path))

    ion_quant_exe_path = os.path.join(frag_base_dir, 'IonQuant-1.10.27', 'IonQuant-1.10.27.jar')
    msfrag_exe_path = os.path.join(frag_base_dir, 'MSFragger-4.0', 'MSFragger-4.0.jar')
    philosopher_exe_path = os.path.join(frag_base_dir, 'philosopher-v5.1.1')

    cmd = [exe_abs_path, '--headless', '--workflow', workflow_path, '--manifest',
           manifest_path, '--workdir', fragpipe_output_path,
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
    file_path = param.get('data_path')
    thread_num = param.get('thread_num')
    workflow_path = param.get('workflow_path')
    workdir = param.get('workdir')
    os.makedirs(workdir, exist_ok=True)
    manifest_path = build_manifest(file_path, workdir)
    run_fragpipe(manifest_path, workflow_path, workdir, fragpipe_exe_path, thread_num)
    # check bin file
    base_file_name = os.path.basename(file_path)
    # find _edited.pin file
    pin_file_list = list(Path(workdir).rglob(f'{base_file_name}_edited.pin'))
    if len(pin_file_list) > 0:
        pin_file_path = pin_file_list[0]
    else:
        return -1
    if os.path.exists(pin_file_path):
        logger.info(f"generate fragpipe success:")
        return 0
    else:
        logger.error(f"generate fragpipe fail")
        return -1
