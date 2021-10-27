use norad::{Font, Glyph};

use crate::Opt;

const AREA_KEY: &str = "com.ht.spacer.area";
const DEPTH_KEY: &str = "com.ht.spacer.depth";
const OVERSHOOT_KEY: &str = "com.ht.spacer.overshoot";

pub struct SpacingParameters {
    pub area: f64,
    pub depth: f64,
    pub overshoot: f64,
}

#[derive(Debug)]
pub enum SpacingParametersError {
    ExpectedRealNumberGlyph(String),
    ExpectedRealNumberFont(String),
}

impl SpacingParameters {
    pub(crate) fn try_new_determine(
        cli: &Opt,
        font: &Font,
        glyph: &Glyph,
    ) -> Result<Self, SpacingParametersError> {
        Ok(Self {
            area: Self::cascading_lookup(AREA_KEY, cli.area, glyph, font)?,
            depth: Self::cascading_lookup(DEPTH_KEY, cli.depth, glyph, font)?,
            overshoot: Self::cascading_lookup(OVERSHOOT_KEY, cli.overshoot, glyph, font)?,
        })
    }

    #[cfg(test)]
    pub(crate) fn try_new_determine_or_default(
        font: &Font,
        glyph: &Glyph,
    ) -> Result<Self, SpacingParametersError> {
        Ok(Self {
            area: Self::cascading_lookup(AREA_KEY, 400.0, glyph, font)?,
            depth: Self::cascading_lookup(DEPTH_KEY, 15.0, glyph, font)?,
            overshoot: Self::cascading_lookup(OVERSHOOT_KEY, 0.0, glyph, font)?,
        })
    }

    fn cascading_lookup(
        key: &str,
        default: f64,
        glyph: &Glyph,
        font: &Font,
    ) -> Result<f64, SpacingParametersError> {
        // Look for a key first in the glyph lib, then the font lib, then fall back to the CLI argument.
        match glyph.lib.get(key) {
            Some(v) => match v.as_real() {
                Some(r) => Ok(r),
                None => Err(SpacingParametersError::ExpectedRealNumberGlyph(key.into())),
            },
            None => match font.lib.get(key) {
                Some(v) => match v.as_real() {
                    Some(r) => Ok(r),
                    None => Err(SpacingParametersError::ExpectedRealNumberFont(key.into())),
                },
                None => Ok(default),
            },
        }
    }
}

impl std::error::Error for SpacingParametersError {}

impl std::fmt::Display for SpacingParametersError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpacingParametersError::ExpectedRealNumberGlyph(key) => {
                write!(f, "Expected glyph lib key '{}' to be a number.", key)
            }
            SpacingParametersError::ExpectedRealNumberFont(key) => {
                write!(f, "Expected font lib key '{}' to be a number.", key)
            }
        }
    }
}
