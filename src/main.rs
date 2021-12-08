use std::collections::{HashMap, HashSet};

use kurbo::{BezPath, Rect, Shape};
use log::{error, warn};
use norad::{Font, Glyph, GlyphName, Layer};
use structopt::StructOpt;

use parameters::SpacingParameters;

pub mod computation;
pub mod config;
pub mod drawing;
pub mod parameters;

// TODO:
// - Make spacing polygons BezPaths instead of Vec<Point> for free `area` fn?
// - Make SpacingParameters also hold FontMetrics?

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
                    // Apply affine counter-translation instead of changing the x_offset
                    // directly. Otherwise, a comma flipped horizontally and vertically
                    // would move in the opposite direction of the delta.
                    let transform: kurbo::Affine = component.transform.into();
                    component.transform =
                        (transform * kurbo::Affine::translate((-delta_left, 0.0))).into();
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
        let bounds_reference = match glyphs.get(reference_name) {
            Some(path) => &path.1,
            None => {
                warn!(
                    "Reference glyph '{}' does not exist, spacing '{}' with own bounds.",
                    reference_name, glyph_name
                );
                &bounds
            }
        };

        let overshoot = font_metrics.xheight * parameters.overshoot / 100.0;

        let bounds_reference_lower = (bounds_reference.min_y() - overshoot).round();
        let bounds_reference_upper = (bounds_reference.max_y() + overshoot).round();

        // Stash new metrics away. We have to do a second iteration so we can get a
        // mut ref to glyphs to modify them. Doing it in this loop is complicated by
        // misc. functions needing a normal ref to layer while we hold a mut ref...
        // TODO: Still correct explanation when using `GlyphName`?
        if let Some((new_left, new_right)) = computation::compute_sidebearings(
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
            use std::num::FpCategory::{Normal, Zero};
            if !(matches!(new_left.classify(), Zero | Normal)
                && matches!(new_right.classify(), Zero | Normal))
            {
                warn!(
                    "Glyph '{}': failed to compute spacing, skipping.",
                    glyph_name
                );
                continue;
            }

            let (delta_left, advance_width) = if font_metrics.angle == 0.0 {
                let delta_left = new_left - bounds.min_x();
                (delta_left, bounds.max_x() + delta_left + new_right)
            } else {
                deslant_sidebearings(
                    &paths,
                    new_left,
                    new_right,
                    font_metrics.angle,
                    font_metrics.xheight,
                )
            };
            new_side_bearings.insert(glyph_name.clone(), (delta_left, advance_width));
        }
    }

    new_side_bearings
}

fn deslant_sidebearings(
    paths: &BezPath,
    left: f64,
    right: f64,
    angle: f64,
    xheight: f64,
) -> (f64, f64) {
    // let bbox_slanted = paths.bounding_box();
    let mut paths2 = paths.clone();
    paths2.apply_affine(slant_paths(-angle, (left, xheight / 2.0)));
    let bbox_deslanted = paths2.bounding_box();

    let left_delta = left - bbox_deslanted.min_x();
    let new_advance_width = bbox_deslanted.max_x() + left_delta + right;

    (left_delta.round(), new_advance_width.round())
}

fn slant_paths(angle: f64, offset: (f64, f64)) -> kurbo::Affine {
    let (dx, dy) = offset;
    let t1 = kurbo::Affine::translate((dx, dy));
    let t2 = kurbo::Affine::new([1.0, 0.0, angle.to_radians().tan(), 1.0, 0.0, 0.0]);
    let t3 = kurbo::Affine::translate((-dx, -dy));
    t1 * t2 * t3
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

#[cfg(test)]
mod tests {
    use norad::{Component, Contour};

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
        let mut glyph = Glyph::new_named("a");
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

        assert_eq!(font.get_glyph("a").unwrap().as_ref(), &glyph_clone);
    }

    // Sanity check: empty glyphs should have a zero bounding box.
    #[test]
    fn empty_glyph_has_zero_bbox() {
        let layer = norad::Layer::new("public.default".into(), None);
        let mut glyph = norad::Glyph::new_named("space");
        glyph.width = 250.0;

        let bbox = drawing::path_for_glyph(&glyph, &layer)
            .unwrap()
            .bounding_box();

        assert_eq!(bbox, Rect::ZERO);
    }

    fn glyph_bounds(glyph: &Glyph, glyph_set: &Layer) -> Rect {
        drawing::path_for_glyph(glyph, glyph_set)
            .unwrap()
            .bounding_box()
    }

    fn contour_bounds(contour: &Contour) -> Rect {
        let mut path = BezPath::new();
        for element in drawing::contour_segments(contour).unwrap() {
            path.push(element);
        }
        path.bounding_box()
    }

    fn component_bounds(component: &Component, glyph_set: &Layer) -> Rect {
        let glyph = glyph_set.get_glyph(&component.base).unwrap();
        let mut path = drawing::path_for_glyph(glyph, glyph_set).unwrap();
        let transform: kurbo::Affine = component.transform.into();
        path.apply_affine(transform);
        path.bounding_box()
    }

    fn assert_approx_eq(a: f64, b: f64) {
        assert!((a - b).abs() <= f64::EPSILON, "{:?} != {:?}", a, b);
    }

    #[test]
    fn test_components_pinned() {
        let font = norad::Font::load("testdata/NestedComponents.ufo").unwrap();
        let glyph_set = font.default_layer();

        let mut font_rs = font.clone();
        let parameters = SpacingParameters::default();
        space_default_layer(&mut font_rs, &parameters);
        let glyph_set_rs = font_rs.default_layer();

        // C: test offset of components from each other stays the same.
        let glyph_c = font.get_glyph("C").unwrap();
        let bbox_c_1 = component_bounds(&glyph_c.components[0], glyph_set);
        let bbox_c_2 = component_bounds(&glyph_c.components[1], glyph_set);
        let bbox_c_x_offset = bbox_c_1.min_x() - bbox_c_2.min_x();

        let glyph_c_rs = font_rs.get_glyph("C").unwrap();
        let bbox_c_1_rs = component_bounds(&glyph_c_rs.components[0], glyph_set_rs);
        let bbox_c_2_rs = component_bounds(&glyph_c_rs.components[1], glyph_set_rs);
        let bbox_c_x_offset_new = bbox_c_1_rs.min_x() - bbox_c_2_rs.min_x();

        assert_approx_eq(bbox_c_x_offset, bbox_c_x_offset_new);

        // D: test glyph bbox x-length stays the same to test x-moving of flipped
        // component.
        let bbox_d = glyph_bounds(font.get_glyph("D").unwrap(), glyph_set);
        let bbox_d_len = bbox_d.max_x() - bbox_d.min_x();

        let bbox_d_rs = glyph_bounds(font_rs.get_glyph("D").unwrap(), glyph_set_rs);
        let bbox_d_len_new = bbox_d_rs.max_x() - bbox_d_rs.min_x();

        assert_approx_eq(bbox_d_len, bbox_d_len_new);

        // E: test offset of components from each other stays the same.
        let glyph_e = font.get_glyph("E").unwrap();
        let bbox_e_1 = component_bounds(&glyph_e.components[0], glyph_set);
        let bbox_e_2 = component_bounds(&glyph_e.components[1], glyph_set);
        let bbox_e_3 = component_bounds(&glyph_e.components[2], glyph_set);
        let bbox_e_4 = component_bounds(&glyph_e.components[3], glyph_set);
        let bbox_e_x_offset1 = bbox_e_1.min_x() - bbox_e_2.min_x();
        let bbox_e_x_offset2 = bbox_e_1.min_x() - bbox_e_3.min_x();
        let bbox_e_x_offset3 = bbox_e_1.min_x() - bbox_e_4.min_x();

        let glyph_e_rs = font_rs.get_glyph("E").unwrap();
        let bbox_e_1_rs = component_bounds(&glyph_e_rs.components[0], glyph_set_rs);
        let bbox_e_2_rs = component_bounds(&glyph_e_rs.components[1], glyph_set_rs);
        let bbox_e_3_rs = component_bounds(&glyph_e_rs.components[2], glyph_set_rs);
        let bbox_e_4_rs = component_bounds(&glyph_e_rs.components[3], glyph_set_rs);
        let bbox_e_x_offset1_new = bbox_e_1_rs.min_x() - bbox_e_2_rs.min_x();
        let bbox_e_x_offset2_new = bbox_e_1_rs.min_x() - bbox_e_3_rs.min_x();
        let bbox_e_x_offset3_new = bbox_e_1_rs.min_x() - bbox_e_4_rs.min_x();

        assert_approx_eq(bbox_e_x_offset1, bbox_e_x_offset1_new);
        assert_approx_eq(bbox_e_x_offset2, bbox_e_x_offset2_new);
        assert_approx_eq(bbox_e_x_offset3, bbox_e_x_offset3_new);

        // F: test offset of component and outline from each other stays the same.
        let glyph_f = font.get_glyph("F").unwrap();
        let bbox_f_1 = contour_bounds(&glyph_f.contours[0]);
        let bbox_f_2 = component_bounds(&glyph_f.components[0], glyph_set);
        let bbox_f_x_offset = bbox_f_1.min_x() - bbox_f_2.min_x();

        let glyph_f_rs = font_rs.get_glyph("F").unwrap();
        let bbox_f_1_rs = contour_bounds(&glyph_f_rs.contours[0]);
        let bbox_f_2_rs = component_bounds(&glyph_f_rs.components[0], glyph_set_rs);
        let bbox_f_x_offset_new = bbox_f_1_rs.min_x() - bbox_f_2_rs.min_x();

        assert_approx_eq(bbox_f_x_offset, bbox_f_x_offset_new);
    }

    #[test]
    fn space_mutatorsans() {
        let mut font = norad::Font::load("testdata/mutatorSans/MutatorSansLightWide.ufo").unwrap();
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

        check_expectations(&mut font, &parameters, &expected);
    }

    fn check_expectations(
        font: &mut norad::Font,
        parameters: &SpacingParameters,
        expected: &[(&str, Option<(f64, f64)>)],
    ) {
        space_default_layer(font, parameters);

        for (name, margins) in expected {
            let glyph = font.get_glyph(*name).unwrap();
            let (glyph_reference, factor) = config::config_for_glyph(name);
            let bbox = drawing::path_for_glyph(glyph, font.default_layer())
                .unwrap()
                .bounding_box();

            match margins {
                Some((left, right)) => {
                    let (new_left, new_right) = (bbox.min_x(), glyph.width - bbox.max_x());

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
                None => assert_eq!(bbox, Rect::ZERO),
            }
        }
    }

    #[test]
    fn space_merriweather() {
        let mut font = norad::Font::load("testdata/Merriweather-LightItalic.ufo").unwrap();
        let parameters = SpacingParameters::default();
        space_default_layer(&mut font, &parameters);

        let font_cmp = norad::Font::load("testdata/Merriweather-LightItalic-Respaced.ufo").unwrap();
        for glyph in font.default_layer().iter() {
            dbg!(&glyph.name);

            let glyph_cmp = font_cmp.get_glyph(&glyph.name).unwrap();

            let bbox = drawing::path_for_glyph(glyph, font.default_layer())
                .unwrap()
                .bounding_box();
            let bbox_cmp = drawing::path_for_glyph(glyph_cmp, font_cmp.default_layer())
                .unwrap()
                .bounding_box();

            // if &*glyph.name == "DZcaron" {
            //     // TODO: x1 is += 43???
            //     assert_eq!(bbox_cmp, Rect { x0: -25.0, y0: -16.0, x1: 2692.0, y1: 1988.0 })
            // }

            if bbox == Rect::ZERO {
                assert_eq!(bbox_cmp, Rect::ZERO)
            } else {
                let (left, right) = (bbox.min_x(), glyph.width - bbox.max_x());
                let (new_left, new_right) = (bbox_cmp.min_x(), glyph_cmp.width - bbox_cmp.max_x());

                assert!(
                    (left - new_left).abs() <= 1.0,
                    "Glyph {}: expected left {} but got {}",
                    glyph.name,
                    left,
                    new_left
                );

                // Glyph E: expected right 76 but got 78
                // Glyph acutecomb: expected left 146 but got 113
                assert!(
                    (right - new_right).abs() <= 2.0,
                    "Glyph {}: expected right {} but got {}",
                    glyph.name,
                    right,
                    new_right
                );
            }
        }
    }
}
