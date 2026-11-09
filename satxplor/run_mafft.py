import os
import subprocess
import multiprocessing
import utils.paths as paths
from utils.constant_loader import constants as constants
from logging_config import logger
from concurrent.futures import ThreadPoolExecutor, as_completed


def run_mafft(input_file, output_file, threads=None):
    # Get the number of available threads
    available_threads = multiprocessing.cpu_count()

    # Calculate half of the available threads
    half_threads = max(1, available_threads // 2)

    if threads is None:
        threads = half_threads

    # Build the MAFFT command
    mafft_command = [
        'mafft',
        '--adjustdirection',
        '--reorder',
        '--quiet',
        '--thread',
        str(threads),
        input_file
    ]
    try:
        result = subprocess.run(mafft_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, check=True)
        # Redirect output to the specified file
        if result.stderr:
            print(f"MAFFT STDERR for {input_file}:\n{result.stderr}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT for {input_file}. Error message: {e}")
        raise e

    # Save output to the specified file
    with open(output_file, 'w') as output:
        output.write(result.stdout)


def align_one(input_path, output_path, threads):
    try:
        run_mafft(input_path, output_path, threads=threads)
        logger.info(f"Alignment completed for {input_path}. Output saved to {output_path}")
        return True
    except Exception as e:
        logger.warning(f"Error in mafft for {input_path}, error message: {e}")
        return False


def align_directory(directory, suffix, workers, threads):
    input_files = [f for f in os.listdir(directory) if f.endswith(suffix)]
    if not input_files:
        return

    logger.info(f"Aligning {len(input_files)} files ending with {suffix} "
                f"using {workers} worker(s) x {threads} thread(s)")

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for input_file in input_files:
            output_file = os.path.splitext(input_file)[0] + '_aligned.fasta'
            input_path = directory + input_file
            output_path = directory + output_file
            futures[executor.submit(align_one, input_path, output_path, threads)] = input_path
        for future in as_completed(futures):
            future.result()


def main():
    cpus = max(1, int(os.environ.get("SATXPLOR_CPUS", multiprocessing.cpu_count())))
    workers = max(1, int(os.environ.get("SATXPLOR_MAFFT_WORKERS", "4")))
    threads = max(1, cpus // workers)

    # Monomer alignment
    align_directory(paths.SEQ_SAVE_PATH, "_monomers.fasta", workers, threads)

    # Microhomology alignment
    align_directory(paths.FLANKS_SAVE_ROOT, "_microhomology.fasta", workers, threads)


if __name__ == "__main__":
    main()
