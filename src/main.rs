use clap::{Parser, Subcommand};
use ensreg::{
    feature::{FeatureRegistry, RegulatoryFeature, RegulatoryFeatureType},
    stable_id::{EnsemblSerialId, SpeciesPrefix},
};
use gannot::{format::Bed6Row, genome::GenomicRange};
use std::io::Write;
use std::{fs::File, io::BufWriter, path::PathBuf};

#[derive(Parser)]
#[command(version)]
struct Args {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(name = "regbuild")]
    RegBuild {
        #[arg(value_enum)]
        species: SpeciesPrefix,
        #[arg(short, long)]
        ctcf_bed: Option<PathBuf>,
        #[clap(long, short, default_value_t = false)]
        rewrite_id: bool,
        features_bed: PathBuf,
    },
}

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();
    let args = Args::parse();
    match args.command {
        Commands::RegBuild {
            species,
            ctcf_bed,
            rewrite_id,
            features_bed,
        } => {
            let mut registry = FeatureRegistry::new(species);

            if let Some(ctcf_bed) = ctcf_bed {
                // read in CTCF features
                let ctcf_file = File::open(ctcf_bed)?;
                let mut ctcf_csv = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .comment(Some(b'#'))
                    .has_headers(false)
                    .from_reader(ctcf_file);

                for row in ctcf_csv.deserialize() {
                    let row: Bed6Row = row?;
                    let location = GenomicRange::from_bed_row(&row);
                    let ctcf_feature = RegulatoryFeature::with_genomic_range(
                        location,
                        row.strand(),
                        RegulatoryFeatureType::CTCF,
                    );
                    registry.register_stable_id(ctcf_feature);
                }
            }

            let features_file = File::open(features_bed)?;
            let mut features_csv = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .comment(Some(b'#'))
                .has_headers(false)
                .from_reader(features_file);

            // read in other features
            for row in features_csv.deserialize() {
                if row.is_err() {
                    log::error!("could not deserialise feature: {:?}", row.err());
                    continue;
                }
                let row = row.unwrap();
                let feature = RegulatoryFeature::from_gff_row(&row);
                let gff_id = row.attributes.get("ID");
                if rewrite_id && feature.feature_type() != &RegulatoryFeatureType::EMAR {
                    registry.register_stable_id(feature);
                } else if let Some(gff_id) = gff_id {
                    let serial_id = EnsemblSerialId::try_from_str(gff_id);
                    if let Ok(serial_id) = serial_id {
                        registry.register_serial_id(serial_id, feature);
                    } else {
                        log::error!("could not parse serial id {gff_id} in row {row:?}");
                    }
                } else {
                    log::error!(
                        "could not register a feature with a serial ID because it has no ID defined {row:?}"
                    );
                }
            }

            let stdout = std::io::stdout();
            let lock = stdout.lock();
            let mut writer = BufWriter::new(lock);

            writeln!(writer, "##gff-version 3")?;
            for (feature_id, feature) in registry {
                feature.to_gff_row(&mut writer, &*feature_id)?;
            }

            Ok(())
        }
    }
}
