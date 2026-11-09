import os
import polars as pl
import utils.paths as paths
from utils.constant_loader import constants as constants
from utils.utils import read_gff_output, extract_subsequences
from Bio import SeqIO
import argparse

from logging_config import logger


def extract_array_flanks(df: pl.DataFrame, fasta_records: list,
                         flank_size: int = constants.FLANK_SIZE,
                         flank_type: str = "flanks") -> None:
    grouped = df.group_by(["feature"])

    for group_key, group_df in grouped:

        # extract left flanks
        left_flanks = extract_subsequences(fasta_records, group_df, flanks="left", flank_size=flank_size)
        # extract right flanks
        right_flanks = extract_subsequences(fasta_records, group_df, flanks="right", flank_size=flank_size)
        # merge them together
        left_flanks.extend(right_flanks)

        # Write the subsequences to a new FASTA file
        output_file = paths.FLANKS_SAVE_ROOT + group_key[0] + "_" + flank_type + ".fasta"
        SeqIO.write(left_flanks, output_file, "fasta")


def kmer_file_name(path: str, squish: bool) -> str:
    name = path.split("/")[-1]
    name = name[:-4]
    name = name.split("_")[0:4] if squish else name.split("_")[0:3]
    return "_".join(name)


def summarize_kmer_file(path: str, squish: bool):
    """Return min/max index where roll_mean < 5 for one array's kmer file.

    Streaming one file at a time keeps memory at O(one array) instead of
    O(all kmer rows).  The original code concatenated every file and grouped
    afterwards, which used ~17 GB at 48 M rows and would need ~100 GB at the
    279 M rows of the largest assembly.
    """
    name = kmer_file_name(path, squish)
    tmp = pl.read_csv(path, dtypes={"index": int, "actual": float,
                                    "roll_mean": float}, separator="\t")
    tmp = tmp.filter(pl.col("roll_mean") < 5.0)
    if tmp.height == 0:
        return None
    return {
        "name": name,
        "min_value": int(tmp["index"].min()),
        "max_value": int(tmp["index"].max()),
    }


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('genome_path', help='Path to the subject file in FASTA format')

    args = parser.parse_args()

    fasta_records = list(SeqIO.parse(args.genome_path, 'fasta'))

    logger.info("\n Loaded the fasta")

    file_list = os.listdir(paths.KMER_ANALYSIS_OUT + "data")
    file_list = [paths.KMER_ANALYSIS_OUT + "data/" + string for string in file_list]

    logger.info(f"Streaming {len(file_list)} kmer files")

    rows = []
    for i, path in enumerate(file_list, 1):
        row = summarize_kmer_file(path, constants.SQUISH)
        if row is not None:
            rows.append(row)
        if i % 10000 == 0:
            logger.info(f"  processed {i}/{len(file_list)} kmer files")

    if not rows:
        logger.warning("No k-mer regions passed the roll_mean filter; writing empty table")
        df_tot = pl.DataFrame(schema={
            "name": pl.String,
            "min_value": pl.Int64,
            "max_value": pl.Int64,
        })
    else:
        df_tot = pl.DataFrame(rows)

    df_tot = df_tot.with_columns(
        query=pl.col("name").map_elements(lambda x: x.split("_")[1], return_dtype=pl.String)
    )
    df_tot = df_tot.rename({
        "min_value": "start",
        "max_value": "end",
        "query": "feature",
        "name": "seqnames"
    })

    df_tot = df_tot.with_columns(
        [
            (pl.lit("feature")).alias("group"),
            (pl.lit("*")).alias("strand"),
            (pl.lit("EuSatXplor")).alias("source"),
            (pl.lit(".")).alias("frame"),
            (pl.lit("100")).alias("score")
        ]
    )

    df_tot = df_tot.select(
        ["seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "group"]
    )

    df_tot.write_csv(paths.KMER_ANALYSIS_OUT + "final_array_table_array_maps.tsv",
                     include_header=False, separator="\t")

    logger.info("Created the real edges annotations")

    # read in the original arrays
    df = read_gff_output(paths.ARRAYS_OUT_PATH, headers=False)

    df_real_edges = df_tot
    if constants.SQUISH:
        df = df.with_columns(
            name=pl.col("seqnames") + "_" + pl.col("group") + "_" + pl.col("start").cast(pl.String) + "_" + pl.col("end").cast(pl.String)
        ).with_columns(
            len=pl.col("end") - pl.col("start")
        )

        df = df.join(df_real_edges, left_on="name", right_on="seqnames")

        df = df.with_columns(
            real_start=pl.col("start") - (500 - pl.col("start_right"))).with_columns(
            real_end=pl.when(pl.col("len") < 5000)
                        .then(pl.col("end_right") - pl.col("start_right") + pl.col("real_start"))
                        .otherwise(pl.col("end") - (4500 - pl.col("end_right")))
        ).select(
            ["seqnames", "source", "feature", "real_start", "real_end", "score", "strand", "frame", "group"]
        ).rename(
            {
                "real_start": "start",
                "real_end": "end"
            }
        )

        df.write_csv(paths.KMER_ANALYSIS_OUT + "real_edges_in_genome.gff", include_header=False, separator="\t")

        logger.info("Extracted the real flanks")

        extract_array_flanks(df=df, fasta_records=fasta_records)

        logger.info("Extracted the microhomology regions")

        extract_array_flanks(df=df, fasta_records=fasta_records,
                             flank_size=20, flank_type="microhomology")

    else:
        df = df.with_columns(
            name=pl.col("seqnames") + "_" + pl.col("group") + "_" + pl.col("start").cast(pl.String)
        )

        df = (df.join(df_real_edges, left_on="name", right_on="seqnames")
              .with_columns(
                  real_start=pl.col("start") - (500 - pl.col("start_right")))
              .with_columns(
                  real_end=pl.col("end_right") - pl.col("start_right") + pl.col("real_start")
              )
              .select(
                  ["seqnames", "source", "feature", "real_start", "real_end", "score", "strand", "frame", "group"]
              )
              .rename(
                  {
                      "real_start": "start",
                      "real_end": "end"
                  }
              )
              )

        df.write_csv(paths.KMER_ANALYSIS_OUT + "real_edges_in_genome.gff", include_header=False, separator="\t")

        logger.info("Extracted the real flanks")

        extract_array_flanks(df=df, fasta_records=fasta_records)

        logger.info("Extracted the microhomology regions")

        extract_array_flanks(df=df, fasta_records=fasta_records,
                             flank_size=20, flank_type="microhomology")
