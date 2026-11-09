import os
import shutil
import tempfile
import subprocess
from logging_config import logger
import utils.paths as paths
from utils.constant_loader import constants as constants
import utils.utils as utils

def run_pca_script(alignment):

    logger.info(f"Running PCA umap for {alignment}")
    command = ['Rscript', './satxplor/r/pca_umap_plots.R',
                    "-a", alignment,
                    "-d", constants.DIMENSION_RED_MODE,
                    "-o", paths.PCA_UMAP_SAVE_ROOT]
    
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode == 0:
        logger.info(f"R script executed successfully for {alignment}.")
        return True

    print(f"Error executing R script for {alignment}. Return code:", result.returncode)
    print("STDERR:\n", result.stderr)
    return False

    
def run_networks_script(distance_matrix):

    logger.info(f"Running networks {distance_matrix}")
    command = ['Rscript', './satxplor/r/network_plots.R',
                    "-a", paths.ARRAYS_OUT_PATH,
                    "-m", distance_matrix,
                    "-k", "5",
                    "-o", paths.NETWORKS_SAVE_ROOT]
    
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode == 0:
        logger.info(f"R script executed successfully for {distance_matrix}.")
        return True

    print(f"Error executing R script for {distance_matrix}. Return code:", result.returncode)
    print("STDERR:\n", result.stderr)
    return False

        
def run_flank_kmers(flank_fasta):
    """Compute Mash distances directly from the raw flank FASTA and plot them."""
    logger.info(f"Running k-mer flank distances for {flank_fasta}")

    stem = os.path.splitext(os.path.basename(flank_fasta))[0]
    workdir = tempfile.mkdtemp(prefix="satxplor_flankdist_")
    matrix_path = os.path.join(workdir, stem + ".kmer_dist.bin")
    cpus = max(1, int(os.environ.get("SATXPLOR_CPUS", os.cpu_count() or 1)))

    command = ['python3', './satxplor/kmer_distances.py',
               '-a', flank_fasta,
               '-o', matrix_path,
               '--threads', str(cpus)]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error computing k-mer distances for {flank_fasta}. Return code:", result.returncode)
        print("STDERR:\n", result.stderr)
        shutil.rmtree(workdir, ignore_errors=True)
        return False

    command = ['Rscript', './satxplor/r/distance_maps.R',
               '-a', matrix_path,
               '-o', paths.DISTANCE_SAVE_ROOT]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode == 0:
        logger.info(f"K-mer distance map executed successfully for {flank_fasta}.")
        ok = True
    else:
        print(f"Error executing R script for {flank_fasta}. Return code:", result.returncode)
        print("STDERR:\n", result.stderr)
        ok = False

    shutil.rmtree(workdir, ignore_errors=True)
    return ok


def run_microhomology():
    search_string = 'microhomology_aligned'
    matching_files = utils.get_files_containing_string(paths.FLANKS_SAVE_ROOT, search_string)
    print(matching_files)
    failed = []
    for mhl in matching_files:
        logger.info(f"Running flank alignment {mhl}")
        command = ['Rscript', './satxplor/r/microhomology_analysis.R',
                    "-a", mhl,
                    "-o", paths.MICROHOMOLOGY_SAVE_ROOT]
    
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode == 0:
            logger.info(f"R script executed successfully for {mhl}.")
        else:
            print(f"Error executing R script for {mhl}. Return code:", result.returncode)
            print("STDERR:\n", result.stderr)
            failed.append(mhl)

    if failed:
        raise RuntimeError(
            f"Microhomology analysis failed for {len(failed)} file(s): "
            + ", ".join(failed)
        )
    
