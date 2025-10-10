use anyhow::anyhow;
use log::warn;
use regex::Regex;
use std::hash::Hash;
use std::{cmp, ops::Range};

use gannot::genome::{GenomicRange, SeqId};

const ID_RADIX: u64 = 27;
const ID_DIGITS: [char; ID_RADIX as usize] = [
    '2', '3', '4', '5', '6', '7', '8', '9', 'B', 'C', 'D', 'F', 'G', 'H', 'J', 'K', 'L', 'M', 'N',
    'P', 'Q', 'R', 'S', 'T', 'W', 'X', 'Z',
];
const MAX_CHROMOSOME_LENGTH: u64 = 1_000_000_000;
const MAX_BIN_ID_CHARS: usize = (MAX_CHROMOSOME_LENGTH.ilog(ID_RADIX) + 1) as usize;

pub trait FeatureId: std::fmt::Display {}
impl FeatureId for EnsemblSerialId {}
impl FeatureId for StableId {}

#[derive(clap::ValueEnum, Debug, Hash, Clone, Copy, PartialOrd, Ord, Eq, PartialEq)]
pub enum SpeciesPrefix {
    Human,
    Mouse,
    Pig,
    Cow,
    Chicken,
    Salmon,
    Trout,
    Turbot,
    Carp,
    Seabass,
}

impl std::fmt::Display for SpeciesPrefix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let prefix = match self {
            Self::Human => "ENSR",
            Self::Mouse => "ENSMUSR",
            Self::Pig => "ENSSSCR",
            Self::Cow => "ENSBTAR",
            Self::Chicken => "ENSGALR",
            Self::Salmon => "ENSSSAR",
            Self::Trout => "ENSOMYR",
            Self::Turbot => "ENSSMAR",
            Self::Carp => "ENSCCRR",
            Self::Seabass => "ENSDLAR",
        };
        write!(f, "{prefix}")
    }
}

#[derive(PartialOrd, Ord, Eq, PartialEq)]
pub struct EnsemblSerialId {
    prefix: SpeciesPrefix,
    serial_num: u32,
}

impl EnsemblSerialId {
    pub fn try_from_str(value: &str) -> Result<Self, anyhow::Error> {
        let re = Regex::new(r"^(ENS[A-Z]*R)([0-9]+)$").unwrap();
        let parts = re
            .captures(value)
            .ok_or_else(|| anyhow!("could not parse id {value}"))?;

        let first = parts.get(1).expect("regex should have parsed id");
        let second = parts.get(2).expect("regex should have parsed id");
        let serial_num: u32 = second.as_str().parse()?;
        let prefix = match first.as_str() {
            "ENSR" => SpeciesPrefix::Human,
            "ENSMUSR" => SpeciesPrefix::Mouse,
            "ENSSSCR" => SpeciesPrefix::Pig,
            "ENSBTAR" => SpeciesPrefix::Cow,
            "ENSGALR" => SpeciesPrefix::Chicken,
            "ENSSSAR" => SpeciesPrefix::Salmon,
            "ENSOMYR" => SpeciesPrefix::Trout,
            "ENSSMAR" => SpeciesPrefix::Turbot,
            "ENSCCRR" => SpeciesPrefix::Carp,
            "ENSDLAR" => SpeciesPrefix::Seabass,
            _ => return Err(anyhow!("unknown prefix for id {value}")),
        };
        Ok(EnsemblSerialId { prefix, serial_num })
    }

    pub fn species(&self) -> SpeciesPrefix {
        self.prefix
    }
}

impl std::fmt::Display for EnsemblSerialId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{:011}", self.prefix, self.serial_num)
    }
}

#[derive(Debug, Hash, PartialEq, Eq)]
pub struct StableId {
    species: SpeciesPrefix,
    seqid: SeqId,
    bin_exp: IdInt,
    bin_num: IdInt,
}

impl Ord for StableId {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.seqid.cmp(&other.seqid) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        // order is reversed so that higher bin_exp come first
        match other.bin_exp.cmp(&self.bin_exp) {
            cmp::Ordering::Equal => {}
            ord => return ord,
        }
        self.bin_num.cmp(&other.bin_num)
    }
}

impl PartialOrd for StableId {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl StableId {
    pub fn from_genomic_range(species: SpeciesPrefix, grange: &GenomicRange) -> StableId {
        let (bin_exp_u32, bin_num_u64) = assign_bin(grange.range_0halfopen());
        let bin_num = IdInt::from(bin_num_u64);
        if bin_num.string_len() > MAX_BIN_ID_CHARS {
            log::error!("bin id with more chars than expected");
        }
        let bin_exp = IdInt::from(bin_exp_u32 as u64);
        let seqid = grange.seqid().to_owned();

        StableId {
            species,
            seqid,
            bin_exp,
            bin_num,
        }
    }
    pub fn seqid(&self) -> &SeqId {
        &self.seqid
    }
}

impl std::fmt::Display for StableId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}_{}{}",
            self.species, self.seqid, self.bin_exp, self.bin_num
        )
    }
}

fn assign_bin(range: Range<u64>) -> (u32, u64) {
    // quarter bin size is determined by 2^bin_exp
    // bin centre is closest to start of one quarter size bin
    // final bin is the two quarter bins on both sides of the bin centre
    let length = range.end - range.start;
    let bin_exp = u32::saturating_sub(length.ilog2(), 1); // might not work for small features
    let quarter_binsize = 2_u64.pow(bin_exp); // there is probably a maximum size before overflow
    let half_binsize = quarter_binsize * 2;
    let binsize = 2 * half_binsize;
    let bin_centre = (range.start + length / 2) / quarter_binsize;
    let id_coord_start = u64::saturating_sub(bin_centre * quarter_binsize, half_binsize);
    let id_coord_end = id_coord_start + binsize;

    // info!("Quarter bin size: {}. ID coverage: {}", quarter_binsize, binsize);
    // info!("Goal is:   {}\t{}", range.start, range.end);
    // info!("ID covers: {}\t{}", id_coord_start, id_coord_end);
    let left_overhang = u64::saturating_sub(range.start, id_coord_start);
    let right_overhang = u64::saturating_sub(id_coord_end, range.end);
    let bases_covered = if id_coord_start < range.start {
        std::cmp::min(range.start + length, id_coord_end) - range.start
    } else if id_coord_end > range.start + length {
        range.start + length - std::cmp::max(range.start, id_coord_start)
    } else {
        length
    };
    let quarter_length = length / 4;
    let core_left = range.start + quarter_length;
    let core_right = range.end - quarter_length;
    let left = std::cmp::max(core_left, id_coord_start);
    let right = std::cmp::min(core_right, id_coord_end);
    let core_bases_covered = u64::saturating_sub(right, left);
    let core_bases_covered_perc = 50.0 * core_bases_covered as f64 / quarter_length as f64;

    let bases_covered_perc = 100.0 * bases_covered as f64 / length as f64;
    let left_overhang_perc = 100.0 * left_overhang as f64 / binsize as f64;
    let right_overhang_perc = 100.0 * right_overhang as f64 / binsize as f64;
    let total_overhang_perc = left_overhang_perc + right_overhang_perc;
    if core_bases_covered_perc < 100.0 || total_overhang_perc > 50.0 {
        warn!(
            "Range: {range:?}. Core bases covered: {core_bases_covered_perc:.1}%. Bases covered: {bases_covered_perc:.1}%. Left overhang {left_overhang_perc:.1}%. Right overhang {right_overhang_perc:.1}%."
        );
    }
    assert_eq!(left_overhang + right_overhang + bases_covered, binsize);
    (bin_exp, bin_centre)
}

#[derive(Debug)]
pub struct IdInt {
    string: String,
    value: u64,
}

impl IdInt {
    pub fn string_len(&self) -> usize {
        self.string.len()
    }

    pub fn value(&self) -> u64 {
        self.value
    }
}

impl Hash for IdInt {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.value.hash(state);
    }
}

impl PartialEq for IdInt {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl Eq for IdInt {}

impl Ord for IdInt {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.value.cmp(&other.value)
    }
}

impl PartialOrd for IdInt {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::fmt::Display for IdInt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.string.fmt(f)
    }
}

impl TryFrom<&str> for IdInt {
    type Error = anyhow::Error;

    fn try_from(str: &str) -> Result<Self, Self::Error> {
        if str.is_empty() {
            return Err(anyhow!("input string for IdInt cannot be empty"));
        }
        let mut num = 0;
        for (i, digit) in str.chars().rev().enumerate() {
            let value = ID_DIGITS
                .iter()
                .position(|&x| x == digit)
                .ok_or_else(|| anyhow!("digit not found"))?;
            num += value as u64 * ID_RADIX.pow(i as u32);
        }
        Ok(IdInt {
            value: num,
            string: str.to_owned(),
        })
    }
}

impl From<u64> for IdInt {
    fn from(value: u64) -> Self {
        let mut result = vec![];
        let mut start = value;
        loop {
            let digit = start % ID_RADIX;
            start /= ID_RADIX;
            result.push(ID_DIGITS[digit as usize]);
            if start == 0 {
                break;
            }
        }
        let string: String = result.into_iter().rev().collect();
        IdInt { string, value }
    }
}

impl From<IdInt> for u64 {
    fn from(id_int: IdInt) -> Self {
        id_int.value
    }
}
