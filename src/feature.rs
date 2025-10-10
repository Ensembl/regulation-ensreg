use crate::stable_id::{EnsemblSerialId, FeatureId, SpeciesPrefix};
use gannot::{
    format::{Gff3Row, Strand},
    genome::GenomicRange,
};
use indexmap::IndexMap;
use serde::Deserialize;
use std::{
    collections::{BTreeMap, btree_map},
    fmt::{self, Display},
};

use super::stable_id::StableId;

#[derive(Debug, PartialEq, Deserialize, Clone, Copy)]
#[serde(rename_all = "lowercase")]
pub enum RegulatoryFeatureType {
    Promoter,
    Enhancer,
    #[serde(rename = "open_chromatin_region")]
    OpenChromatin,
    #[serde(rename = "CTCF_binding_site")]
    CTCF,
    #[serde(rename = "EMAR")]
    EMAR,
}

impl RegulatoryFeatureType {
    pub fn rgb(&self) -> String {
        match self {
            RegulatoryFeatureType::Promoter => "#ff0000",
            RegulatoryFeatureType::Enhancer => "#faca00",
            RegulatoryFeatureType::OpenChromatin => "#d9d9d9",
            RegulatoryFeatureType::CTCF => "#40e0d0",
            RegulatoryFeatureType::EMAR => "#004d40",
        }
        .to_owned()
    }
}

impl Display for RegulatoryFeatureType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                RegulatoryFeatureType::Promoter => "promoter",
                RegulatoryFeatureType::Enhancer => "enhancer",
                RegulatoryFeatureType::OpenChromatin => "open_chromatin_region",
                RegulatoryFeatureType::CTCF => "CTCF_binding_site",
                RegulatoryFeatureType::EMAR => "EMAR",
            }
        )
    }
}

#[derive(Debug)]
pub struct RegulatoryFeature {
    feature_type: RegulatoryFeatureType,
    location: GenomicRange,
    extended_location: Option<GenomicRange>, // extended_location here should have the same seqid as location
    strand: Strand,
    other_attributes: Option<IndexMap<String, String>>,
}

impl RegulatoryFeature {
    pub fn with_genomic_range(
        location: GenomicRange,
        strand: Strand,
        feature_type: RegulatoryFeatureType,
    ) -> RegulatoryFeature {
        RegulatoryFeature {
            feature_type,
            location,
            extended_location: None,
            strand,
            other_attributes: None,
        }
    }

    pub(crate) fn original_id(&self) -> Option<&str> {
        self.other_attributes
            .as_ref()
            .and_then(|attr| attr.get("ID"))
            .map(|id| id.as_str())
    }

    pub fn from_gff_row(row: &Gff3Row<RegulatoryFeatureType>) -> RegulatoryFeature {
        let location = GenomicRange::from_gff_row(row);
        let mut other_attributes = row.attributes.clone();
        let mut extended_location = if row.attributes.contains_key("extended_start")
            || row.attributes.contains_key("extended_end")
        {
            if let (Some(start), Some(end)) = (
                other_attributes.shift_remove("extended_start"),
                other_attributes.shift_remove("extended_end"),
            ) {
                let start: Result<u64, _> = start.parse();
                let end: Result<u64, _> = end.parse();
                if let (Ok(start), Ok(end)) = (start, end) {
                    let grange = GenomicRange::from_1closed(row.seqid.clone(), start..=end);
                    if let Ok(grange) = grange {
                        if row.feature_type == RegulatoryFeatureType::Promoter {
                            Some(grange)
                        } else if grange == location {
                            None
                        } else {
                            log::warn!(
                                "extended_start and extended_end defined for an unexpected feature type in {row:?}"
                            );
                            Some(grange)
                        }
                    } else {
                        log::error!("could not parse extended_start or extended_end in {row:?}");
                        None
                    }
                } else {
                    log::error!("could not parse extended_start or extended_end in {row:?}");
                    None
                }
            } else {
                log::error!("only one of extended_start or extended_end is defined in {row:?}");
                None
            }
        } else {
            None
        };
        if row.feature_type == RegulatoryFeatureType::Promoter && extended_location.is_none() {
            extended_location = Some(location.clone());
        }
        if let Some(ref extended_location) = extended_location {
            let core_range = location.range_1closed();
            let extended_range = extended_location.range_1closed();
            // Check that core_range is fully contained within extended_range (inclusive)
            if core_range.start() < extended_range.start()
                || core_range.end() > extended_range.end()
            {
                log::error!(
                    "core range ({core_range:?}) is not fully contained within the extended range {extended_range:?}"
                );
            }
        }

        RegulatoryFeature {
            feature_type: row.feature_type,
            location,
            extended_location,
            strand: row.strand,
            other_attributes: if other_attributes.is_empty() {
                None
            } else {
                Some(other_attributes)
            },
        }
    }

    pub fn to_gff_row<T: std::io::Write>(
        &self,
        writer: &mut T,
        id: &dyn FeatureId,
    ) -> Result<(), anyhow::Error> {
        let mut attributes = vec![format!("ID={}", id)];
        let gff_range = self.location.range_1closed();
        if let Some(ref extended_location) = self.extended_location {
            let range = extended_location.range_1closed();
            attributes.push(format!("extended_start={}", range.start()));
            attributes.push(format!("extended_end={}", range.end()));
        };
        if let Some(ref other_attributes) = self.other_attributes {
            for (key, value) in other_attributes {
                match key.as_str() {
                    "ID" => {
                        if value != &id.to_string() {
                            log::info!("changed ID from {value} to {id}");
                        }
                    }
                    "color" => {}
                    _ => attributes.push(format!("{key}={value}")),
                }
            }
        }
        attributes.push(format!("color={}", self.feature_type.rgb()));
        writeln!(
            writer,
            "{}\tEnsembl\t{}\t{}\t{}\t.\t{}\t.\t{}",
            self.location.seqid(),
            self.feature_type,
            gff_range.start(),
            gff_range.end(),
            self.strand,
            attributes.join(";")
        )?;
        Ok(())
    }

    pub fn location(&self) -> &GenomicRange {
        &self.location
    }

    pub fn extended_location(&self) -> &Option<GenomicRange> {
        &self.extended_location
    }

    pub fn strand(&self) -> &Strand {
        &self.strand
    }

    pub fn feature_type(&self) -> &RegulatoryFeatureType {
        &self.feature_type
    }
}

/// A registry for storing and managing regulatory features by their stable and serial IDs for a given species.
pub struct FeatureRegistry {
    species: SpeciesPrefix,
    stable_ids: BTreeMap<StableId, RegulatoryFeature>,
    serial_ids: BTreeMap<EnsemblSerialId, RegulatoryFeature>,
}

pub struct FeatureRegistryIntoIter {
    stable_iter: btree_map::IntoIter<StableId, RegulatoryFeature>,
    serial_iter: btree_map::IntoIter<EnsemblSerialId, RegulatoryFeature>,
}

impl Iterator for FeatureRegistryIntoIter {
    type Item = (Box<dyn FeatureId>, RegulatoryFeature);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((key, value)) = self.stable_iter.next() {
            Some((Box::new(key), value))
        } else if let Some((key, value)) = self.serial_iter.next() {
            Some((Box::new(key), value))
        } else {
            None
        }
    }
}

impl IntoIterator for FeatureRegistry {
    type Item = (Box<dyn FeatureId>, RegulatoryFeature);
    type IntoIter = FeatureRegistryIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        FeatureRegistryIntoIter {
            stable_iter: self.stable_ids.into_iter(),
            serial_iter: self.serial_ids.into_iter(),
        }
    }
}

impl FeatureRegistry {
    pub fn new(species: SpeciesPrefix) -> Self {
        Self {
            species,
            stable_ids: BTreeMap::new(),
            serial_ids: BTreeMap::new(),
        }
    }

    pub fn register_serial_id(&mut self, id: EnsemblSerialId, feature: RegulatoryFeature) {
        assert_eq!(id.species(), self.species);
        let existing_feature = self.serial_ids.remove(&id);
        match existing_feature {
            None => {
                self.serial_ids.insert(id, feature);
            }
            Some(existing_feature) => {
                log::error!(
                    "two features with the same serial id.\
                {feature:?} and {existing_feature:?}"
                )
            }
        }
    }

    pub fn register_stable_id(&mut self, feature: RegulatoryFeature) {
        let stable_id = StableId::from_genomic_range(self.species, feature.location());
        let existing_feature = self.stable_ids.remove(&stable_id);
        match existing_feature {
            None => {
                self.stable_ids.insert(stable_id, feature);
            }
            Some(mut existing_feature) => {
                if existing_feature.feature_type == feature.feature_type
                    && feature.feature_type == RegulatoryFeatureType::CTCF
                {
                    assert!(feature.strand != Strand::None);
                    let new_location = feature.location.combine(existing_feature.location());
                    let new_strand = if existing_feature.strand == feature.strand {
                        existing_feature.strand
                    } else {
                        Strand::None
                    };
                    match new_location {
                        Err(err) => log::error!(
                            "error combining two features with the same stable id: {stable_id}.\
                                    Features are {existing_feature:?} and {feature:?}.\
                                    Error: {err}"
                        ),
                        Ok(new_location) => {
                            let new_stable_id =
                                StableId::from_genomic_range(self.species, &new_location);
                            existing_feature.location = new_location;
                            existing_feature.strand = new_strand;
                            log::warn!(
                                "combined two features with the same stable id: {stable_id}. New stable id is {new_stable_id}."
                            );
                            self.stable_ids.insert(new_stable_id, existing_feature);
                        }
                    }
                } else if existing_feature.feature_type == RegulatoryFeatureType::CTCF
                    && (feature.feature_type == RegulatoryFeatureType::Enhancer
                        || feature.feature_type == RegulatoryFeatureType::OpenChromatin)
                {
                    log::warn!(
                        "CTCF feature with same stable id ({}) as another feature. Discarding {} with original ID {}.",
                        stable_id,
                        feature.feature_type,
                        feature.original_id().unwrap_or_default()
                    );
                    self.stable_ids.insert(stable_id, existing_feature);
                } else if (existing_feature.feature_type == RegulatoryFeatureType::Enhancer
                    || existing_feature.feature_type == RegulatoryFeatureType::OpenChromatin)
                    && feature.feature_type == RegulatoryFeatureType::CTCF
                {
                    log::warn!(
                        "CTCF feature with same stable id ({}) as another feature. Discarding {} with original ID {}.",
                        stable_id,
                        existing_feature.feature_type,
                        feature.original_id().unwrap_or_default()
                    );
                    self.stable_ids.insert(stable_id, feature);
                } else {
                    log::error!(
                        "two features with the same stable id ({stable_id}) that cannot be combined.\
                    {feature:?} and {existing_feature:?}"
                    )
                }
            }
        }
    }
}
