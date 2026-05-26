import argparse
import logging
import re
from pathlib import Path

import pandas as pd

BASE_HEADER = ["feature_id", "feature_type"]
ACTIVITY_MAP = {
    "EVIDENCE": "ACTIVE",
    "NO EVIDENCE": "INACTIVE",
}
DEFAULT_CHUNK_SIZE = 100_000
OPEN_CHROMATIN_PEAKS_COLUMN = "base_epigenome_open_chromatin_peaks"
REGULATORY_ACTIVITY_MARKER = ".Regulatory_activity"
SPECIES_COLUMN = "species_name"
LABEL_COLUMNS = {
    SPECIES_COLUMN,
    OPEN_CHROMATIN_PEAKS_COLUMN,
    "term",
    "term_for_label",
    "sex_for_label",
    "age_for_label",
}
LOGGER = logging.getLogger(__name__)


def parse_gff_attributes(value: str) -> dict[str, str]:
    attributes = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Malformed GFF attribute: {item!r}")
        key, attribute_value = item.split("=", 1)
        attributes[key] = attribute_value

    return attributes


def extract_id(value: str) -> str:
    values = parse_gff_attributes(value)
    feature_id = values.get("ID")
    if not feature_id:
        raise ValueError(f"GFF attributes do not contain an ID: {value!r}")

    id_match = re.search(r"ENS\D+[0-9]{11}", feature_id)
    if not id_match:
        raise ValueError(f"GFF ID does not look like an Ensembl ID: {feature_id}")

    return feature_id


def read_gff_feature_types(gff: Path) -> dict[str, str]:
    feature_types = {}
    with gff.open() as gff_file:
        for line_number, line in enumerate(gff_file, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{gff}:{line_number} has {len(fields)} columns; expected 9."
                )

            feature_id = extract_id(fields[8])
            if feature_id in feature_types:
                raise ValueError(f"Duplicate GFF feature ID: {feature_id}")
            feature_types[feature_id] = fields[2]

    if not feature_types:
        raise ValueError(f"No feature records found in {gff}.")

    return feature_types


def clean_label_value(value) -> str:
    if pd.isna(value):
        return ""
    return str(value).replace("\xa0", " ").strip()


def activity_epigenome_key(activity_file: str | Path) -> str:
    filename = Path(activity_file).name
    prefix, marker, _ = filename.partition(REGULATORY_ACTIVITY_MARKER)
    if not marker:
        raise ValueError(
            f"Could not extract epigenome key from activity file: {filename}"
        )

    key = prefix.split(".")[-1]
    if not key:
        raise ValueError(
            f"Could not extract epigenome key from activity file: {filename}"
        )

    return key


def open_chromatin_peak_key(peaks_file: str) -> str:
    peaks_file = clean_label_value(peaks_file)
    if not peaks_file:
        raise ValueError("Empty open chromatin peaks file value.")

    first_peaks_file = peaks_file.strip("{}").split(",")[0]
    filename = Path(first_peaks_file).name
    stem = filename.removesuffix(".bed.gz").removesuffix(".bed")
    stem = re.sub(r"_[pq]_\d+_\d+$", "", stem)

    for marker in ("-dnase-seq-", "-atac-seq-"):
        if marker in stem:
            key = stem.rsplit(marker, 1)[1]
            if key:
                return key

    raise ValueError(
        f"Could not extract epigenome key from open chromatin peaks file: {peaks_file}"
    )


def epigenome_label(row: pd.Series) -> str:
    label = clean_label_value(row.get("term_for_label")) or clean_label_value(
        row.get("term")
    )
    if not label:
        raise ValueError("Could not build epigenome label without a term value.")

    details = [
        clean_label_value(row.get("sex_for_label")),
        clean_label_value(row.get("age_for_label")),
    ]
    details = [detail for detail in details if detail]

    if details:
        return f"{label} ({', '.join(details)})"
    return label


def list_activity_files(reg_activity_folder: Path) -> list[Path]:
    activity_files = sorted(
        path for path in reg_activity_folder.iterdir() if path.suffix.lower() == ".csv"
    )
    if not activity_files:
        raise ValueError(f"No activity CSV files found in {reg_activity_folder}.")

    return activity_files


def build_epigenome_labels(
    epigenome_csv: Path, species_name: str, activity_keys: set[str]
) -> dict[str, str]:
    species_name = clean_label_value(species_name)
    if not species_name:
        raise ValueError("species_name must not be empty.")

    LOGGER.info(
        "Reading epigenome metadata for %s from %s",
        species_name,
        epigenome_csv,
    )
    epigenomes = pd.read_csv(epigenome_csv, usecols=sorted(LABEL_COLUMNS))

    missing_columns = LABEL_COLUMNS - set(epigenomes.columns)
    if missing_columns:
        raise ValueError(
            f"{epigenome_csv} is missing required columns: "
            f"{', '.join(sorted(missing_columns))}"
        )

    epigenomes = epigenomes[
        epigenomes[SPECIES_COLUMN].map(clean_label_value) == species_name
    ]
    if epigenomes.empty:
        raise ValueError(
            f"{epigenome_csv} does not contain metadata for species {species_name!r}."
        )

    labels: dict[str, str] = {}
    for _, row in epigenomes.iterrows():
        peaks_file = clean_label_value(row[OPEN_CHROMATIN_PEAKS_COLUMN])
        if not peaks_file:
            continue

        try:
            key = open_chromatin_peak_key(peaks_file)
        except ValueError:
            continue
        if key not in activity_keys:
            continue

        label = epigenome_label(row)
        existing_label = labels.get(key)
        if existing_label and existing_label != label:
            raise ValueError(
                f"Conflicting labels for epigenome key {key!r}: "
                f"{existing_label!r} and {label!r}"
            )
        labels[key] = label

    missing_keys = activity_keys - set(labels.keys())
    if missing_keys:
        raise ValueError(
            f"{epigenome_csv} does not contain labels for species "
            f"{species_name!r} and activity keys: "
            f"{', '.join(sorted(missing_keys))}"
        )

    return labels


def epigenome_header(activity_file: str | Path, labels: dict[str, str]) -> str:
    key = activity_epigenome_key(activity_file)
    return labels[key]


def read_activity_file(
    activity_file: Path,
    column_index: int,
    column_count: int,
    data: dict[str, list[str]],
    chunk_size: int,
) -> int:
    seen_feature_ids = set()
    row_count = 0

    LOGGER.info("Reading activity file %s", activity_file)
    chunks = pd.read_csv(
        activity_file,
        names=["regulatory_feature_id", "activity"],
        skiprows=1,
        dtype=str,
        chunksize=chunk_size,
    )

    for chunk in chunks:
        chunk["regulatory_feature_id"] = chunk["regulatory_feature_id"].map(
            clean_label_value
        )
        chunk["activity"] = chunk["activity"].map(clean_label_value)

        empty_ids = chunk["regulatory_feature_id"] == ""
        if empty_ids.any():
            raise ValueError(f"{activity_file} contains empty feature IDs.")

        duplicated_in_chunk = chunk["regulatory_feature_id"].duplicated()
        if duplicated_in_chunk.any():
            duplicate = chunk.loc[duplicated_in_chunk, "regulatory_feature_id"].iloc[0]
            raise ValueError(
                f"{activity_file} contains duplicate feature ID {duplicate!r}."
            )

        already_seen = set(chunk["regulatory_feature_id"]) & seen_feature_ids
        if already_seen:
            duplicate = sorted(already_seen)[0]
            raise ValueError(
                f"{activity_file} contains duplicate feature ID {duplicate!r}."
            )
        seen_feature_ids.update(chunk["regulatory_feature_id"])

        unknown_activities = set(chunk["activity"]) - set(ACTIVITY_MAP)
        if unknown_activities:
            raise ValueError(
                f"{activity_file} contains unknown activity values: "
                f"{', '.join(sorted(unknown_activities))}"
            )

        chunk["activity"] = chunk["activity"].map(ACTIVITY_MAP)
        for feature_id, activity in zip(
            chunk["regulatory_feature_id"], chunk["activity"]
        ):
            activities = data.setdefault(feature_id, [""] * column_count)
            activities[column_index] = activity

        row_count += len(chunk)

    return row_count


def read_activity_data(
    activity_files: list[Path], chunk_size: int
) -> dict[str, list[str]]:
    if chunk_size < 1:
        raise ValueError("chunk_size must be greater than zero.")

    data: dict[str, list[str]] = {}
    column_count = len(activity_files)
    for column_index, activity_file in enumerate(activity_files):
        row_count = read_activity_file(
            activity_file,
            column_index,
            column_count,
            data,
            chunk_size,
        )
        LOGGER.info("Read %s activity rows from %s", row_count, activity_file)

    return data


def validate_row_width(row: list[str], expected_width: int) -> None:
    if len(row) != expected_width:
        raise ValueError(
            f"Output row has {len(row)} columns; expected {expected_width}."
        )


def write_activity_output(
    output_filename: Path,
    header: list[str],
    feature_types: dict[str, str],
    data: dict[str, list[str]],
) -> None:
    missing_features = set(data) - set(feature_types)
    if missing_features:
        raise ValueError(
            "Activity files contain feature IDs missing from the GFF: "
            f"{', '.join(sorted(missing_features)[:20])}"
        )

    expected_width = len(header)
    missing_activity_cells = sum(row.count("") for row in data.values())
    if missing_activity_cells:
        LOGGER.warning(
            "Writing %s empty activity cells for IDs absent from some files.",
            missing_activity_cells,
        )

    LOGGER.info("Writing merged activity file to %s", output_filename)
    with output_filename.open("w") as output:
        output.write("\t".join(header) + "\n")
        for feature_id, feature_type in feature_types.items():
            if feature_id not in data:
                continue

            row = [feature_id, feature_type, *data[feature_id]]
            validate_row_width(row, expected_width)
            output.write("\t".join(row) + "\n")


def validate_inputs(
    reg_activity_folder: Path, gff3_filename: Path, epigenome_csv: Path
) -> None:
    if not reg_activity_folder.is_dir():
        raise NotADirectoryError(f"{reg_activity_folder} is not a directory.")
    if not gff3_filename.is_file():
        raise FileNotFoundError(f"{gff3_filename} is not a file.")
    if not epigenome_csv.is_file():
        raise FileNotFoundError(f"{epigenome_csv} is not a file.")


def merge_activities(
    reg_activity_folder: Path,
    gff3_filename: Path,
    epigenome_csv: Path,
    species_name: str,
    output_filename: Path,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> None:
    validate_inputs(reg_activity_folder, gff3_filename, epigenome_csv)

    activity_files = list_activity_files(reg_activity_folder)
    activity_keys = {activity_epigenome_key(af) for af in activity_files}
    labels = build_epigenome_labels(epigenome_csv, species_name, activity_keys)
    header = BASE_HEADER + [epigenome_header(af, labels) for af in activity_files]

    LOGGER.info("Reading GFF feature types from %s", gff3_filename)
    feature_types = read_gff_feature_types(gff3_filename)
    data = read_activity_data(activity_files, chunk_size)
    write_activity_output(output_filename, header, feature_types, data)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge regulatory activity CSV files into a TSV."
    )
    parser.add_argument(
        "reg_activity_folder",
        type=Path,
        help="Path to a folder containing all the CSV activity files",
    )
    parser.add_argument(
        "gff3",
        type=Path,
        help="Regulatory build GFF3 file",
    )
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "-m",
        "--epigenome-metadata",
        required=True,
        type=Path,
        help="CSV containing base epigenome metadata",
    )
    parser.add_argument(
        "-s",
        "--species-name",
        required=True,
        help=("Species name to select from the metadata CSV, e.g. 'Homo sapiens'"),
    )
    parser.add_argument(
        "--chunk-size",
        default=DEFAULT_CHUNK_SIZE,
        type=int,
        help="Rows to read at a time from each activity CSV",
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.INFO,
    )
    args = parse_args(argv)
    merge_activities(
        args.reg_activity_folder,
        args.gff3,
        args.epigenome_metadata,
        args.species_name,
        args.output,
        chunk_size=args.chunk_size,
    )


if __name__ == "__main__":
    main()
