import json
import argparse
import os


def update_config(input_genome_path, sat_raw, final_results_dir, overwrite):
    # Load the existing configuration
    with open("run_config.json", 'r') as config_file:
        config = json.load(config_file)

    def resolve(path):
        # Paths are relative to the /mnt/data mount unless already absolute.
        if os.path.isabs(path):
            return path
        return "/mnt/data/" + path

    # Update the configuration with new paths
    config['INPUT_GENOME_PATH'] = resolve(input_genome_path)
    config['SAT_RAW'] = resolve(sat_raw)
    # The outer launcher bind-mounts the staging directory at ./results and
    # asks the controller to write there directly.  In that case keep the
    # relative path so it resolves to the mount, not to /mnt/data/results.
    if os.path.isabs(final_results_dir) or final_results_dir == "./results":
        config['FINAL_RESULTS_DIR'] = final_results_dir
    else:
        config['FINAL_RESULTS_DIR'] = "/mnt/data/" + final_results_dir
    config['OVERWRITE'] = overwrite

    # Save the updated configuration
    with open("run_config.json", 'w') as config_file:
        json.dump(config, config_file, indent=4)


def main():
    parser = argparse.ArgumentParser(description="Update paths in the run_config.json file.")
    parser.add_argument('--input_genome_path', type=str, required=True, help='Path to the input genome file')
    parser.add_argument('--sat_raw', type=str, required=True, help='Path to the satellite raw file')
    parser.add_argument('--final_results_dir', type=str, required=True, help='Path to the final results directory')
    parser.add_argument('--overwrite', type=bool, default=True, help='Overwrite existing files')

    args = parser.parse_args()

    update_config(
        args.input_genome_path,
        args.sat_raw,
        args.final_results_dir,
        args.overwrite
    )


if __name__ == '__main__':
    print("Updating config file for your settings")
    main()
    print("Config file updated successfully")
