import os
import shutil
import subprocess
from utils.constant_loader import constants as constants
import utils.utils as utils
import utils.paths as paths
import argparse
import sys
from logging_config import logger
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from r_script_runners import (run_pca_script,
                                  run_flank_kmers,
                                  run_microhomology,
                                  run_networks_script )
import json
from preprocess import sanitize_and_filter_sequences


def _available_memory_gb():
    try:
        with open('/proc/meminfo') as handle:
            for line in handle:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) / (1024.0 * 1024.0)
    except OSError:
        pass
    return None


def _parallel_workers(per_worker_gb, budget_fraction=0.5):
    """Bound worker count by free memory as well as by CPU count.

    The R stages build dense N x N matrices and use roughly 10-14 GB per
    process.  One worker per CPU overcommits memory and the kernel kills R
    with SIGKILL (return code -9), so also budget against MemAvailable.
    """
    cpus = max(1, int(os.environ.get('SATXPLOR_CPUS', os.cpu_count() or 1)))
    available = _available_memory_gb()
    if available is None:
        return cpus
    workers = int((available * budget_fraction) / per_worker_gb)
    return max(1, min(cpus, workers))


def _clear_directory(path):
    """Remove the contents of a directory without removing the directory.

    Needed when ./results is a bind mount: rmtree() would fail with EBUSY
    when it tries to remove the mount point itself.
    """
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
        return
    for entry in os.listdir(path):
        item = os.path.join(path, entry)
        if os.path.isdir(item) and not os.path.islink(item):
            shutil.rmtree(item)
        else:
            os.remove(item)


def save_constants(data):
    with open('satxplor/utils/constants.json', 'w') as f:
        json.dump(data, f, indent=4)

def update_constants(args):
    # Update constants with new values from command-line arguments (if provided)
    constants = {
        'MAX_PLOT_LEN': args.max_plot_len,
        'HISTOGRAM_BIN_WIDTH': args.histogram_bin_width,
        'NKERNEL_BINS': args.nkernel_bins,
        'FLANK_SIZE': args.flank_size,
        'ARRAY_HOR_PERC': args.array_hor_perc,
        'MONOMER_NUMBER': args.monomer_number,
        'CONTIG_FILTER': args.contig_filter,
        'DIMENSION_RED_MODE': args.dimension_red_mode,
        'PERC_ID_FILTER': args.perc_id_filter,
        'QCOVHSP_FILTER': args.qcovhsp_filter,
        'SQUISH': args.squish,
        'SQUISH_LEN': args.squish_len
    }

    # Save updated constants back to constants.json
    save_constants(constants)

    # Print updated constants
    print(f"Updated Constants:\n{json.dumps(constants, indent=4)}")


with open("run_config.json", "r") as f:
    config = json.load(f)
checkpoints = ["nothing",
               "blast",
                   "arrays",
                   "density",
                   "extract",
                   "filenames",
                   "kmercount",
                   "output",
                   "alignment",
                   "PCA"]

def audit_checkpoints() -> str: 
    checkpoint_file = "./checkpoints.tsv"
    # Check if checkpoint file exists or is empty
    if not os.path.exists(checkpoint_file):
        with open(checkpoint_file, mode='w') as file:
            file.write("Checkpoint\tDatetime\n")
        return("nothing")
    else:
        # Read the content of the last row
        with open(checkpoint_file, mode='r') as file:
            lines = file.readlines()
            last_row = lines[-1].strip()

            # Extract the last checkpoint name and datetime
            last_checkpoint, last_datetime = last_row.split('\t')
            
            if last_checkpoint=="Checkpoint":
                logger.info("Checkpoint file was created, but nothing was in it, starting from scratch.")
                return("nothing")
            else: 
                logger.info(f"Last checkpoint is: {last_checkpoint} at {last_datetime}")
                return(last_checkpoint)
    

def create_checkpoint(checkpoint_name):
    checkpoint_file = "./checkpoints.tsv"
    
    # Append checkpoint and datetime to the checkpoint file
    with open(checkpoint_file, mode='a') as file:
        file.write(f"{checkpoint_name}\t{datetime.now()}\n")


#new blast will always reinitialize folders
        
def initialize_preprocess_blast():
    # Create or clear the 'results' folder
    logger.info('Removing data in existing results and tmp directories.')
    results_folder = './results'
    _clear_directory(results_folder)

    #data saving roots
    os.makedirs(paths.DATA_SAVE_ROOT)
    os.makedirs(paths.TABLE_SAVE_ROOT)

    #sequence saving roots
    os.makedirs(paths.SEQ_SAVE_PATH)
    os.makedirs(paths.FLANKS_SAVE_ROOT)

    #picture saving roots
    os.makedirs(paths.PIC_SAVE_ROOT)
    os.makedirs(paths.PCA_UMAP_SAVE_ROOT)
    os.makedirs(paths.DISTANCE_SAVE_ROOT)
    os.makedirs(paths.NETWORKS_SAVE_ROOT)
    os.makedirs(paths.MICROHOMOLOGY_SAVE_ROOT)
    
    logger.info('Data deleted and folders created.')

    logger.info(
        f"Running preprocessing script for the input.\n  Removing all sequences shorter than {constants.CONTIG_FILTER}.")
    
    input_file = config["INPUT_GENOME_PATH"]
    output_file = config["GENOME_PATH"]
    
    sanitize_and_filter_sequences(input_file, output_file,min_length=constants.CONTIG_FILTER)

    input_file = config["SAT_RAW"]
    output_file = config["SAT_FASTA_PATH"]
    
    
    logger.info("Sanitizing input names")
    
    sanitize_and_filter_sequences(input_file, output_file,min_length=0)


    logger.info( f"Output sequences written to { config['GENOME_PATH'] }")

    command = ['python3', './satxplor/blast.py',
               config["SAT_FASTA_PATH"],
               config["GENOME_PATH"]
               ]
    
    print(f"Blast command:{command}")
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        logger.warning(f"Error. Return code:{ result.returncode}")
        print("STDERR:\n", result.stderr)
        quit()


def run_create_arrays_script():
    # Run the 'create_arrays.py' script

    command = ['python3', './satxplor/create_arrays.py']
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()
    

def run_density_plots():
    # Run the 'create_arrays.py' script
    command = ['python3', './satxplor/density_plots.py']
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        #this is an optional output, if error, continue normaly
        logger.warning(f"Error with density plots. Return code:{ result.returncode}")
        
    

def run_extract_naive():
    # Run the 'create_arrays.py' script
    command = ['python3', './satxplor/extract_naive.py',
               config["GENOME_PATH"]]
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()

def run_filename_mapping():
    # Run the 'create_arrays.py' script

    command = ['python3', './satxplor/filename_mapping.py']
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()

def run_kmer_edge_finder():
    # Run the 'create_arrays.py' script
    command = ['./satxplor/executables/kmer_edge_finder']
    
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()

def run_rust_output_processing():
    # Run the 'create_arrays.py' script
    command = ['python', './satxplor/rust_output_processing.py',
                              config["GENOME_PATH"]]

    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        print("STDERR:\n", result.stderr)
        quit()

def run_mafft():
    # Run the 'create_arrays.py' script
    command = ['python3', './satxplor/run_mafft.py']
    result = subprocess.run(command)
    if result.returncode == 0:
        pass
    else:
        logger.warning(f"Error. Return code:{ result.returncode}. Error in running MAFFT for some sequences. ")

def run_pca():
    search_string = 'monomers_aligned'
    matching_files = utils.get_files_containing_string(paths.SEQ_SAVE_PATH, search_string)

    workers = _parallel_workers(per_worker_gb=16.0, budget_fraction=0.5)
    logger.info(f"Running PCA/UMAP with {workers} parallel worker(s)")
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(run_pca_script, matching_files))
    failed = [f for f, ok in zip(matching_files, results) if not ok]
    if failed:
        raise RuntimeError(
            f"PCA/UMAP failed for {len(failed)} file(s): " + ", ".join(failed)
        )

def run_networks():
    search_string = 'matrix.csv.gz'
    matching_files = utils.get_files_containing_string(paths.PCA_UMAP_SAVE_ROOT, search_string)

    
    workers = _parallel_workers(per_worker_gb=20.0, budget_fraction=0.4)
    logger.info(f"Running network plots with {workers} parallel worker(s)")
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(run_networks_script, matching_files))
    failed = [f for f, ok in zip(matching_files, results) if not ok]
    if failed:
        raise RuntimeError(
            f"Network plots failed for {len(failed)} file(s): " + ", ".join(failed)
        )

def run_flanks():
    # Compute Mash distances directly on the raw 500 bp flanks.  No MSA is
    # needed, so run_mafft.py no longer aligns the flank files.
    matching_files = [
        path for path in utils.get_files_containing_string(
            paths.FLANKS_SAVE_ROOT, '_flanks.fasta')
        if path.endswith('_flanks.fasta')
    ]
    matching_files.sort()

    workers = min(3, _parallel_workers(per_worker_gb=25.0, budget_fraction=0.4))
    logger.info(f"Running k-mer flank distances with {workers} parallel worker(s)")
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(run_flank_kmers, matching_files))
    failed = [f for f, ok in zip(matching_files, results) if not ok]
    if failed:
        raise RuntimeError(
            f"Flank distance maps failed for {len(failed)} file(s): "
            + ", ".join(failed)
        )


def main():
    logger.info("Python Version:")
    logger.info(sys.version)

    logger.info(config)
    
    if os.path.exists(config["FINAL_RESULTS_DIR"]):
        logger.info("The output folder already exists")
        if not config["OVERWRITE"]:

            logger.error(f"Destination folder '{config[ 'FINAL_RESULTS_DIR' ] }' already exists. but OVERWRITE is off. Exiting the application.")
        else: 
            logger.info("Since OVERWRITE is on, trying to delete the folder!")
            _clear_directory(config["FINAL_RESULTS_DIR"])
            logger.info("Folder successfully deleted")



    check_point = audit_checkpoints()

    starting_id = checkpoints.index(check_point)

    parser = argparse.ArgumentParser(description='Update constants in constants.json')

    # Hard-coded default values in argparse
    parser.add_argument('--max_plot_len', type=int, default=5000, help='Maximum plot length')
    parser.add_argument('--histogram_bin_width', type=int, default=50, help='Histogram bin width')
    parser.add_argument('--nkernel_bins', type=int, default=100, help='Number of kernel bins')
    parser.add_argument('--flank_size', type=int, default=500, help='Flank size')
    parser.add_argument('--array_hor_perc', type=float, default=0.01, help='Array horizontal percentage')
    parser.add_argument('--monomer_number', type=int, default=4, help='Monomer number')
    parser.add_argument('--contig_filter', type=int, default=100000, help='Contig filter')
    parser.add_argument('--dimension_red_mode', type=str, default='both', help='Dimension reduction mode')
    parser.add_argument('--perc_id_filter', type=float, default=70.0, help='Percentage identity filter')
    parser.add_argument('--qcovhsp_filter', type=float, default=70.0, help='Query coverage filter')
    parser.add_argument('--squish', type=bool, default=True, help='Squish flag')
    parser.add_argument('--squish_len', type=int, default=2500, help='Squish length')

    args = parser.parse_args()

    update_constants(args)


    logger.info("The output folder already exists")
    fun_list = [initialize_preprocess_blast,
                run_create_arrays_script,
                run_density_plots,
                run_extract_naive,
                run_filename_mapping,
                run_kmer_edge_finder,
                run_rust_output_processing,
                run_mafft,
                run_flanks,
                run_pca,
                run_networks,
                run_microhomology]
    
    for i in range(starting_id,len(checkpoints)+2):
        fun_list[i]()
        try:
            create_checkpoint(checkpoints[i])
        except IndexError:
            pass
        except Exception as e:
            print("Stopping execution due to:", e)

    
    
    logger.info(f"Running finished at {datetime.now()}, deleting the checkpoints file. Copying results file to output directory {config['FINAL_RESULTS_DIR']}")
    os.remove("checkpoints.tsv")


        

# Copy the entire folder from source to destination.  When ./results is
    # already a bind mount of FINAL_RESULTS_DIR, write in place and skip the
    # expensive (and redundant) copy.
    same_results = False
    try:
        same_results = os.path.samefile("./results", config["FINAL_RESULTS_DIR"])
    except OSError:
        same_results = (os.path.abspath("./results") ==
                        os.path.abspath(config["FINAL_RESULTS_DIR"]))
    if not same_results:
        shutil.copytree("./results", config["FINAL_RESULTS_DIR"])
    else:
        logger.info("Results already written directly to FINAL_RESULTS_DIR; skipping copy")
    utils.move_file(config["GENOME_PATH"],config['FINAL_RESULTS_DIR'])
    utils.move_file(config["SAT_FASTA_PATH"],config['FINAL_RESULTS_DIR'])
    

if __name__ == '__main__':

    main()