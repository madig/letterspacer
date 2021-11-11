use std::collections::{HashMap, HashSet};

use kurbo::{BezPath, Line, ParamCurve, Point, Rect, Shape};
use log::{error, warn};
use norad::{Font, Glyph, GlyphName, Layer};
use structopt::StructOpt;

use parameters::SpacingParameters;

pub mod config;
pub mod drawing;
pub mod parameters;

// TODO:
// - Write deslanter for italic sidebearings
// - Make spacing polygons BezPaths for free `area` fn?

#[derive(Debug, StructOpt)]
#[structopt(
    name = "letterspacer",
    about = "A configurable letter spacer for UFO font sources."
)]
pub(crate) struct Opt {
    /// Set area (overrides font lib key).
    #[structopt(long)]
    area: Option<f64>,

    /// Set depth (overrides font lib key).
    #[structopt(long)]
    depth: Option<f64>,

    /// Set overshoot (overrides font lib key).
    #[structopt(long)]
    overshoot: Option<f64>,

    /// Set sample frequency.
    #[structopt(long)]
    sample_frequency: Option<usize>,

    /// The path to the UFOs to space.
    #[structopt(parse(from_os_str))]
    input: std::path::PathBuf,
}

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn"))
        .format_timestamp(None)
        .init();
    let args = Opt::from_args();

    let mut font = norad::Font::load(&args.input)?;

    // Lookup order: CLI > font lib > default value.
    let parameters = SpacingParameters::default();
    let parameters = SpacingParameters {
        area: match args.area {
            Some(v) => v,
            None => parameters::font_lookup_f64(&font.lib, parameters::AREA_KEY, parameters.area)?,
        },
        depth: match args.depth {
            Some(v) => v,
            None => {
                parameters::font_lookup_f64(&font.lib, parameters::DEPTH_KEY, parameters.depth)?
            }
        },
        overshoot: match args.overshoot {
            Some(v) => v,
            None => parameters::font_lookup_f64(
                &font.lib,
                parameters::OVERSHOOT_KEY,
                parameters.overshoot,
            )?,
        },
        sample_frequency: parameters.sample_frequency,
    };

    space_default_layer(&mut font, &parameters);

    font.meta.creator = Some("org.linebender.norad".into());
    font.save(std::path::PathBuf::from("/tmp").join(args.input.file_name().unwrap()))?;

    Ok(())
}

fn space_default_layer(font: &mut norad::Font, parameters: &SpacingParameters) {
    let font_metrics = FontMetrics::from_font(font);
    let default_layer = font.default_layer_mut();
    let glyphs = prepare_glyphs(default_layer, parameters);
    let composites = composite_graph(default_layer);

    let new_side_bearings = calculate_sidebearings(&glyphs, &font_metrics);

    for (name, (delta_left, advance_width)) in new_side_bearings {
        let glyph = default_layer.get_glyph_mut(&name).unwrap();
        move_glyph_x(glyph, delta_left);
        glyph.width = advance_width;

        // Go though all composites using the glyph and counter-move the component of it to
        // keep it in place.
        if let Some(glyphs_to_countermove) = composites.get(&name) {
            for glyph_name in glyphs_to_countermove {
                let dependent_glyph = default_layer.get_glyph_mut(glyph_name).unwrap();
                for component in dependent_glyph
                    .components
                    .iter_mut()
                    .filter(|c| c.base == name)
                {
                    component.transform.x_offset -= delta_left;
                }
            }
        }
    }
}

/// Returns which glyphs are used as a component in which other glyphs.
///
/// E.g. `{"a": {"aacute", "acircumflex"}, ...}`.
fn composite_graph(layer: &Layer) -> HashMap<GlyphName, HashSet<GlyphName>> {
    let mut composite_graph: HashMap<GlyphName, HashSet<GlyphName>> = HashMap::new();

    for glyph in layer.iter() {
        for component in &glyph.components {
            composite_graph
                .entry(component.base.clone())
                .or_default()
                .insert(glyph.name.clone());
        }
    }

    composite_graph
}

/// Prepares the BezPath, bounding box and glyph-specific spacing parameters for all glyphs in a layer.
fn prepare_glyphs(
    layer: &norad::Layer,
    global_parameters: &SpacingParameters,
) -> HashMap<GlyphName, (BezPath, Rect, SpacingParameters)> {
    let mut paths = HashMap::with_capacity(layer.len());

    for glyph in layer.iter() {
        match drawing::path_for_glyph(glyph, layer) {
            Ok(path) => {
                let bbox = path.bounding_box();
                let parameters =
                    match SpacingParameters::try_new_from_glyph(&glyph.lib, global_parameters) {
                        Ok(p) => p,
                        Err(e) => {
                            error!(
                                "Error while extracting spacing parameters from {}: {:?}",
                                glyph.name, e
                            );
                            continue;
                        }
                    };
                paths.insert(glyph.name.clone(), (path, bbox, parameters));
            }
            Err(e) => {
                error!("Error while drawing {}: {:?}", glyph.name, e);
            }
        };
    }

    paths
}

// TODO: refactor to work per-glyph?
/// Returns a map of glyph name to the delta to move the glyph on the left side by and the advance width.
///
/// First move the glyph by the left delta, then set the advance width.
fn calculate_sidebearings(
    glyphs: &HashMap<GlyphName, (BezPath, Rect, SpacingParameters)>,
    font_metrics: &FontMetrics,
) -> HashMap<GlyphName, (f64, f64)> {
    let mut new_side_bearings = HashMap::new();

    for (glyph_name, (paths, bounds, parameters)) in glyphs.iter() {
        let (reference_name, factor) = config::config_for_glyph(glyph_name);
        let (_, bounds_reference, _) = match glyphs.get(reference_name) {
            Some(path) => path,
            None => {
                warn!(
                    "Skipping {} because reference glyph {} wasn't found.",
                    glyph_name, reference_name
                );
                continue;
            }
        };

        let overshoot = font_metrics.xheight * parameters.overshoot / 100.0;

        let bounds_reference_lower = (bounds_reference.min_y() - overshoot).round();
        let bounds_reference_upper = (bounds_reference.max_y() + overshoot).round();

        // Stash new metrics away. We have to do a second iteration so we can get a
        // mut ref to glyphs to modify them. Doing it in this loop is complicated by
        // misc. functions needing a normal ref to layer while we hold a mut ref...
        // TODO: Still correct explanation when using `GlyphName`?
        if let Some((new_left, new_right)) = calculate_spacing(
            paths,
            bounds,
            (bounds_reference_lower, bounds_reference_upper),
            font_metrics,
            parameters,
            factor,
        ) {
            // Discard bogus computation results (typically from bogus geometry) here
            // rather than in `calculate_spacing` because we know which glyph is affected
            // here.
            use std::num::FpCategory::{Infinite, Nan, Subnormal};
            if matches!(new_left.classify(), Infinite | Nan | Subnormal)
                || matches!(new_right.classify(), Infinite | Nan | Subnormal)
            {
                warn!(
                    "Glyph '{}': failed to compute spacing, skipping.",
                    glyph_name
                );
                continue;
            }

            let delta_left = new_left - bounds.min_x();
            let advance_width = bounds.max_x() + delta_left + new_right;
            new_side_bearings.insert(glyph_name.clone(), (delta_left, advance_width));
        }
    }

    new_side_bearings
}

struct FontMetrics {
    angle: f64,
    units_per_em: f64,
    xheight: f64,
}

impl FontMetrics {
    fn from_font(font: &Font) -> Self {
        Self {
            angle: font
                .font_info
                .italic_angle
                .map(|v| -v.get())
                .unwrap_or_default(),
            units_per_em: font
                .font_info
                .units_per_em
                .map(|v| v.get())
                .unwrap_or(1000.0),
            xheight: font.font_info.x_height.map(|v| v.get()).unwrap_or_default(),
        }
    }
}

/// Shift anchors, contours and components of a glyph horizontally by `delta`.
fn move_glyph_x(glyph: &mut Glyph, delta: f64) {
    for contour in &mut glyph.contours {
        for point in &mut contour.points {
            point.x += delta;
        }
    }
    for component in &mut glyph.components {
        component.transform.x_offset += delta;
    }
    for anchor in &mut glyph.anchors {
        anchor.x += delta;
    }
}

fn calculate_spacing(
    paths: &BezPath,
    bounds: &Rect,
    (bounds_reference_lower, bounds_reference_upper): (f64, f64),
    font_metrics: &FontMetrics,
    parameters: &SpacingParameters,
    factor: f64,
) -> Option<(f64, f64)> {
    if paths.is_empty() {
        return None;
    }

    let (left, extreme_left_full, extreme_left, right, extreme_right_full, extreme_right) =
        spacing_polygons(
            paths,
            bounds,
            (bounds_reference_lower, bounds_reference_upper),
            font_metrics,
            parameters.sample_frequency,
            parameters.depth,
        );

    // Difference between extreme points full and in zone.
    let distance_left = (extreme_left.x - extreme_left_full.x).ceil();
    let distance_right = (extreme_right_full.x - extreme_right.x).ceil();

    let new_left = (-distance_left
        + calculate_sidebearing_value(
            factor,
            (bounds_reference_lower, bounds_reference_upper),
            parameters.area,
            &left,
            font_metrics,
        ))
    .ceil();
    let new_right = (-distance_right
        + calculate_sidebearing_value(
            factor,
            (bounds_reference_lower, bounds_reference_upper),
            parameters.area,
            &right,
            font_metrics,
        ))
    .ceil();

    Some((new_left, new_right))
}

fn calculate_sidebearing_value(
    factor: f64,
    (bounds_reference_lower, bounds_reference_upper): (f64, f64),
    area: f64,
    polygon: &[Point],
    font_metrics: &FontMetrics,
) -> f64 {
    let amplitude_y = bounds_reference_upper - bounds_reference_lower;
    let area_upm = area * (font_metrics.units_per_em / 1000.0).powi(2);
    let white_area = area_upm * factor * 100.0;
    let prop_area = (amplitude_y * white_area) / font_metrics.xheight;
    let valor = prop_area - calculate_area(polygon);
    valor / amplitude_y
}

fn calculate_area(points: &[Point]) -> f64 {
    // https://mathopenref.com/coordpolygonarea2.html
    points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .fold(0.0, |sum, (prev, next)| {
            sum + (prev.x * next.y - next.x * prev.y)
        })
        .abs()
        / 2.0
}

fn spacing_polygons(
    paths: &BezPath,
    bounds: &Rect,
    (bounds_reference_lower, bounds_reference_upper): (f64, f64),
    font_metrics: &FontMetrics,
    sample_frequency: usize,
    depth: f64,
) -> (Vec<Point>, Point, Point, Vec<Point>, Point, Point) {
    // For deskewing angled glyphs. Makes subsequent processing easier.
    let skew_offset = font_metrics.xheight / 2.0;
    let tan_angle = font_metrics.angle.to_radians().tan();

    // First pass: Collect the outer intersections of a horizontal line with the glyph on both sides, going bottom
    // to top. The spacing polygon is vertically limited to lower_bound_reference..=upper_bound_reference,
    // but we need to collect the extreme points on both sides for the full stretch for spacing later.

    // A glyph can over- or undershoot its reference bounds. Measure the tallest stretch.
    let bounds_sampling_lower = bounds.min_y().round().min(bounds_reference_lower) as isize;
    let bounds_sampling_upper = bounds.max_y().round().max(bounds_reference_upper) as isize;

    let mut left = Vec::new();
    let left_bounds = bounds.min_x();
    let mut extreme_left_full: Option<Point> = None;
    let mut extreme_left: Option<Point> = None;
    let mut right = Vec::new();
    let right_bounds = bounds.max_x();
    let mut extreme_right_full: Option<Point> = None;
    let mut extreme_right: Option<Point> = None;
    for y in (bounds_sampling_lower..=bounds_sampling_upper)
        .step_by(sample_frequency)
        .map(|v| v as f64)
    {
        let line = Line::new((left_bounds, y), (right_bounds, y));
        let in_reference_zone = bounds_reference_lower <= y && y <= bounds_reference_upper;

        let mut hits = intersections_for_line(paths, line);
        if hits.is_empty() {
            if in_reference_zone {
                // Treat no hits as hits deep off the other side.
                left.push(Point::new(f64::INFINITY, y));
                right.push(Point::new(-f64::INFINITY, y));
            }
        } else {
            hits.sort_by_key(|k| k.x.round() as i32);
            let mut first = *hits.first().unwrap();
            let mut last = *hits.last().unwrap();
            if font_metrics.angle != 0.0 {
                // De-slant both points to simulate an upright glyph. Makes things easier.
                first = Point::new(first.x - (y - skew_offset) * tan_angle, first.y);
                last = Point::new(last.x - (y - skew_offset) * tan_angle, last.y);
            }
            if in_reference_zone {
                left.push(first);
                right.push(last);

                extreme_left = extreme_left
                    .map(|l| if l.x < first.x { l } else { first })
                    .or(Some(first));
                extreme_right = extreme_right
                    .map(|r| if r.x > last.x { r } else { last })
                    .or(Some(last));
            }

            extreme_left_full = extreme_left_full
                .map(|l| if l.x < first.x { l } else { first })
                .or(Some(first));
            extreme_right_full = extreme_right_full
                .map(|r| if r.x > last.x { r } else { last })
                .or(Some(last));
        }
    }

    // A glyph may have bogus geometry, skip it later.
    let extreme_left_full = extreme_left_full.unwrap_or_default();
    let extreme_left = extreme_left.unwrap_or_default();
    let extreme_right_full = extreme_right_full.unwrap_or_default();
    let extreme_right = extreme_right.unwrap_or_default();

    // Second pass: Cap the margin samples to a maximum depth from the outermost point in to get our depth cut-in.
    let depth = font_metrics.xheight * depth / 100.0;
    let max_depth = extreme_left.x + depth;
    let min_depth = extreme_right.x - depth;
    left.iter_mut().for_each(|s| s.x = s.x.min(max_depth));
    right.iter_mut().for_each(|s| s.x = s.x.max(min_depth));

    // Third pass: Close open counterforms at 45 degrees.
    let dx_max = sample_frequency as f64;

    for i in 0..left.len() - 1 {
        if left[i + 1].x - left[i].x > dx_max {
            left[i + 1].x = left[i].x + dx_max;
        }
        if right[i + 1].x - right[i].x < -dx_max {
            right[i + 1].x = right[i].x - dx_max;
        }
    }
    for i in (0..left.len() - 1).rev() {
        if left[i].x - left[i + 1].x > dx_max {
            left[i].x = left[i + 1].x + dx_max;
        }
        if right[i].x - right[i + 1].x < -dx_max {
            right[i].x = right[i + 1].x - dx_max;
        }
    }

    left.insert(
        0,
        Point {
            x: extreme_left.x,
            y: bounds_reference_lower as f64,
        },
    );
    left.push(Point {
        x: extreme_left.x,
        y: bounds_reference_upper as f64,
    });
    right.insert(
        0,
        Point {
            x: extreme_right.x,
            y: bounds_reference_lower as f64,
        },
    );
    right.push(Point {
        x: extreme_right.x,
        y: bounds_reference_upper as f64,
    });

    (
        left,
        extreme_left_full,
        extreme_left,
        right,
        extreme_right_full,
        extreme_right,
    )
}

fn intersections_for_line(paths: &BezPath, line: Line) -> Vec<Point> {
    paths
        .segments()
        .flat_map(|s| {
            s.intersect_line(line)
                .into_iter()
                .map(move |h| s.eval(h.segment_t).round())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn space_bogus() {
        // Glyphs with bogus geometry should be skipped.

        let zero_point = || {
            norad::ContourPoint::new(
                0.0,
                0.0,
                norad::PointType::OffCurve,
                false,
                None,
                None,
                None,
            )
        };

        let mut font = Font::new();
        let mut glyph = Glyph::new_named("bogus");
        glyph.width = 1000.0;
        glyph.contours.push(norad::Contour::new(
            vec![zero_point(), zero_point()],
            None,
            None,
        ));
        let glyph_clone = glyph.clone();
        font.default_layer_mut().insert_glyph(glyph);

        let parameters = SpacingParameters::default();
        space_default_layer(&mut font, &parameters);

        assert_eq!(font.get_glyph("bogus").unwrap().as_ref(), &glyph_clone);
    }

    // XXX: rewrite https://gist.github.com/madig/9553f66257c4a1abe793b0e82ada84d4
    #[test]
    fn space_mutatorsans() {
        let font = norad::Font::load("testdata/mutatorSans/MutatorSansLightWide.ufo").unwrap();
        let parameters = SpacingParameters::default();
        let expected = [
            ("A", Some((23.0, 23.0))),
            ("Aacute", Some((23.0, 23.0))),
            ("Adieresis", Some((23.0, 23.0))),
            ("B", Some((92.0, 43.0))),
            ("C", Some((49.0, 43.0))),
            ("D", Some((92.0, 49.0))),
            ("E", Some((92.0, 33.0))),
            ("F", Some((92.0, 32.0))),
            ("G", Some((49.0, 66.0))),
            ("H", Some((92.0, 92.0))),
            ("I", Some((33.0, 33.0))),
            ("I.narrow", Some((92.0, 92.0))),
            ("IJ", Some((71.0, 72.0))),
            ("J", Some((41.0, 75.0))),
            ("J.narrow", Some((24.0, 72.0))),
            ("K", Some((92.0, 25.0))),
            ("L", Some((92.0, 25.0))),
            ("M", Some((92.0, 92.0))),
            ("N", Some((92.0, 92.0))),
            ("O", Some((49.0, 49.0))),
            ("P", Some((92.0, 46.0))),
            ("Q", Some((49.0, 49.0))),
            ("R", Some((92.0, 49.0))),
            ("S", Some((38.0, 44.0))),
            ("S.closed", Some((43.0, 42.0))),
            ("T", Some((25.0, 25.0))),
            ("U", Some((72.0, 72.0))),
            ("V", Some((23.0, 23.0))),
            ("W", Some((26.0, 26.0))),
            ("X", Some((19.0, 25.0))),
            ("Y", Some((22.0, 22.0))),
            ("Z", Some((33.0, 33.0))),
            ("acute", Some((79.0, 79.0))),
            ("arrowdown", Some((89.0, 91.0))),
            ("arrowleft", Some((95.0, 111.0))),
            ("arrowright", Some((110.0, 96.0))),
            ("arrowup", Some((91.0, 89.0))),
            ("colon", Some((88.0, 88.0))),
            ("comma", Some((94.0, 91.0))),
            ("dieresis", Some((80.0, 80.0))),
            ("dot", Some((80.0, 80.0))),
            ("period", Some((96.0, 96.0))),
            ("quotedblbase", Some((94.0, 91.0))),
            ("quotedblleft", Some((91.0, 94.0))),
            ("quotedblright", Some((94.0, 91.0))),
            ("quotesinglbase", Some((94.0, 91.0))),
            ("semicolon", Some((88.0, 86.0))),
            ("space", None),
        ];
        let glyphs = prepare_glyphs(font.default_layer(), &parameters);

        check_expectations(&font, &glyphs, &expected);
    }

    fn check_expectations(
        font: &norad::Font,
        glyphs: &HashMap<GlyphName, (BezPath, Rect, SpacingParameters)>,
        expected: &[(&str, Option<(f64, f64)>)],
    ) {
        let font_metrics = FontMetrics::from_font(font);
        let sidebearing_deltas = calculate_sidebearings(glyphs, &font_metrics);

        for (name, margins) in expected {
            match margins {
                Some((left, right)) => {
                    let glyph = font.get_glyph(*name).unwrap();
                    let drawing = drawing::path_for_glyph(glyph, font.default_layer()).unwrap();
                    let bbox = drawing.bounding_box();

                    let (glyph_reference, factor) = config::config_for_glyph(name);
                    let (delta_left, advance_width) = sidebearing_deltas[*name];
                    let (new_left, new_right) = (
                        bbox.min_x() + delta_left,
                        advance_width - bbox.max_x() - delta_left,
                    );

                    assert!(
                        (left - new_left).abs() <= 1.0,
                        "Glyph {}: expected left {} but got {} (glyph reference {}, factor {})",
                        name,
                        left,
                        new_left,
                        glyph_reference,
                        factor
                    );

                    assert!(
                        (right - new_right).abs() <= 1.0,
                        "Glyph {}: expected right {} but got {} (glyph reference {}, factor {})",
                        name,
                        right,
                        new_right,
                        glyph_reference,
                        factor
                    );
                }
                None => assert!(!sidebearing_deltas.contains_key(*name)),
            }
        }
    }
}
