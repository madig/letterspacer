pub const AREA_KEY: &str = "com.ht.spacer.area";
pub const DEPTH_KEY: &str = "com.ht.spacer.depth";
pub const OVERSHOOT_KEY: &str = "com.ht.spacer.overshoot";

pub struct SpacingParameters {
    pub area: f64,
    pub depth: f64,
    pub overshoot: f64,
    pub sample_frequency: usize,
}

#[derive(Debug)]
pub enum SpacingParametersError {
    ExpectedRealNumberGlyph(String),
    ExpectedRealNumberFont(String),
}

impl Default for SpacingParameters {
    fn default() -> Self {
        Self {
            area: 400.0,
            depth: 15.0,
            overshoot: 0.0,
            sample_frequency: 5,
        }
    }
}

impl SpacingParameters {
    pub fn try_new_from_glyph(
        lib: &norad::Plist,
        fallback: &Self,
    ) -> Result<Self, SpacingParametersError> {
        Ok(Self {
            area: glyph_lookup_f64(lib, AREA_KEY, fallback.area)?,
            depth: glyph_lookup_f64(lib, DEPTH_KEY, fallback.depth)?,
            overshoot: glyph_lookup_f64(lib, OVERSHOOT_KEY, fallback.overshoot)?,
            sample_frequency: fallback.sample_frequency,
        })
    }
}

pub fn font_lookup_f64(
    lib: &norad::Plist,
    key: &str,
    fallback: f64,
) -> Result<f64, SpacingParametersError> {
    if let Some(v) = lib.get(key) {
        if let Some(v) = v.as_real() {
            Ok(v)
        } else if let Some(i) = v.as_signed_integer() {
            Ok(i as f64)
        } else {
            Err(SpacingParametersError::ExpectedRealNumberFont(key.into()))
        }
    } else {
        Ok(fallback)
    }
}

pub fn glyph_lookup_f64(
    lib: &norad::Plist,
    key: &str,
    fallback: f64,
) -> Result<f64, SpacingParametersError> {
    if let Some(v) = lib.get(key) {
        if let Some(v) = v.as_real() {
            Ok(v)
        } else if let Some(i) = v.as_signed_integer() {
            Ok(i as f64)
        } else {
            Err(SpacingParametersError::ExpectedRealNumberGlyph(key.into()))
        }
    } else {
        Ok(fallback)
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
