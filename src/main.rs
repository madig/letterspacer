use log::error;

use std::collections::HashMap;

use kurbo::{Affine, BezPath, Line, ParamCurve, PathEl, Point, Rect, Shape};
use norad::{Component, Contour, ContourPoint, Glyph, Layer, PointType};
use structopt::StructOpt;

pub mod config;

// TODO:
// - Write `set_(left|right)_margin` plus `move`
// - Write deslanter for italic sidebearings
// - Write new_layer
// - Make spacing polygons BezPaths for free `area` fn?

#[derive(Debug, StructOpt)]
#[structopt(
    name = "letterspacer",
    about = "A configurable letter spacer for UFO font sources."
)]
struct Opt {
    /// Set area.
    #[structopt(long, default_value = "400.0")]
    area: f64,

    /// Set depth.
    #[structopt(long, default_value = "15.0")]
    depth: f64,

    /// Set overshoot.
    #[structopt(long, default_value = "0.0")]
    overshoot: f64,

    /// Set sample frequency.
    #[structopt(long, default_value = "5")]
    sample_frequency: usize,

    /// The path to the UFOs to space.
    #[structopt(parse(from_os_str))]
    input: std::path::PathBuf,
}

fn main() {
    env_logger::init();
    let args = Opt::from_args();

    let mut ufo = match norad::Font::load(&args.input) {
        Ok(v) => v,
        Err(e) => {
            error!("Loading UFO failed: {}", e);
            std::process::exit(1);
        }
    };

    // TODO: Read from glyph lib, fall back to UFO lib, fall back to args
    let param_area = args.area;
    let param_depth = args.depth;
    let param_overshoot = args.overshoot;
    let param_sample_frequency = args.sample_frequency;

    let units_per_em: f64 = ufo
        .font_info
        .units_per_em
        .map(|v| v.get())
        .unwrap_or(1000.0);
    let angle: f64 = ufo
        .font_info
        .italic_angle
        .map(|v| -v.get())
        .unwrap_or_default();
    let xheight: f64 = ufo.font_info.x_height.map(|v| v.get()).unwrap_or_default();

    let overshoot = xheight * param_overshoot / 100.0;

    let default_layer = ufo.default_layer();
    // let mut background_glyphs: Vec<Glyph> = Vec::new();

    let mut new_side_bearings = HashMap::new();
    for glyph in default_layer.iter() {
        let (glyph_reference, factor) = config::config_for_glyph(&glyph.name);
        let glyph_reference = default_layer.get_glyph(glyph_reference).unwrap();

        let paths = match path_for_glyph(glyph, default_layer) {
            Ok(path) => path,
            Err(e) => {
                error!("Error while drawing {}: {:?}", glyph.name, e);
                continue;
            }
        };
        let bounds = paths.bounding_box();
        let paths_reference = match path_for_glyph(glyph_reference, default_layer) {
            Ok(path) => path,
            Err(e) => {
                error!("Error while drawing {}: {:?}", glyph_reference.name, e);
                continue;
            }
        };
        let bounds_reference = paths_reference.bounding_box();
        let bounds_reference_lower = (bounds_reference.min_y() - overshoot).round();
        let bounds_reference_upper = (bounds_reference.max_y() + overshoot).round();

        let (new_left, new_right) = calculate_spacing(
            paths,
            bounds,
            (bounds_reference_lower, bounds_reference_upper),
            angle,
            xheight,
            param_sample_frequency,
            param_depth,
            // glyph.name.clone(),
            // &mut background_glyphs,
            factor,
            param_area,
            units_per_em,
        );

        // Stash new metrics away. We have to do a second iteration so we can get a
        // mut ref to glyphs to modify them. Doing it in this loop is complicated by
        // misc. functions needing a normal ref to default_layer while we hold a mut ref...
        if let (Some(new_left), Some(new_right)) = (new_left, new_right) {
            let delta_left = new_left - bounds.min_x();
            let delta_right = bounds.max_x() + delta_left + new_right;
            new_side_bearings.insert(glyph.name.clone(), (delta_left, delta_right));
        }
    }

    let default_layer = ufo.default_layer_mut();
    for (name, (left, right)) in new_side_bearings {
        let glyph = default_layer.get_glyph_mut(&name).unwrap();
        move_glyph_x(glyph, left);
        glyph.width = right;
        // TODO: go though all composites using `glyph` and counter-move them to
        // keep them in place. Need some sort of topological sorter to move
        // component glyphs in the right order, once.
    }

    // Write out background layer.
    // let mut background_layer = ufo.layers.get_or_create("public.background");
    // for glyph in background_glyphs {
    //     background_layer.insert_glyph(glyph)
    // }

    ufo.meta.creator = Some("org.linebender.norad".into());
    ufo.save(std::path::PathBuf::from("/tmp").join(args.input.file_name().unwrap()))
        .unwrap();
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
    paths: BezPath,
    bounds: Rect,
    (bounds_reference_lower, bounds_reference_upper): (f64, f64),
    angle: f64,
    xheight: f64,
    param_sample_frequency: usize,
    param_depth: f64,
    // glyph_name: impl Into<GlyphName>,
    // background_glyphs: &mut Vec<Glyph>,
    factor: f64,
    param_area: f64,
    units_per_em: f64,
) -> (Option<f64>, Option<f64>) {
    if paths.is_empty() {
        return (None, None);
    }

    let (left, extreme_left_full, extreme_left, right, extreme_right_full, extreme_right) =
        spacing_polygons(
            &paths,
            &bounds,
            (bounds_reference_lower, bounds_reference_upper),
            angle,
            xheight,
            param_sample_frequency,
            param_depth,
        );

    // let background_glyph = draw_glyph_outer_outline_into_glyph(glyph_name, (&left, &right));
    // background_glyphs.push(background_glyph);

    // Difference between extreme points full and in zone.
    let distance_left = (extreme_left.x - extreme_left_full.x).ceil();
    let distance_right = (extreme_right_full.x - extreme_right.x).ceil();

    let new_left = (-distance_left
        + calculate_sidebearing_value(
            factor,
            (bounds_reference_lower, bounds_reference_upper),
            param_area,
            &left,
            units_per_em,
            xheight,
        ))
    .ceil();
    let new_right = (-distance_right
        + calculate_sidebearing_value(
            factor,
            (bounds_reference_lower, bounds_reference_upper),
            param_area,
            &right,
            units_per_em,
            xheight,
        ))
    .ceil();

    (Some(new_left), Some(new_right))
}

fn calculate_sidebearing_value(
    factor: f64,
    (bounds_reference_lower, bounds_reference_upper): (f64, f64),
    param_area: f64,
    polygon: &[Point],
    units_per_em: f64,
    xheight: f64,
) -> f64 {
    let amplitude_y = bounds_reference_upper - bounds_reference_lower;
    let area_upm = param_area * (units_per_em / 1000.0).powi(2);
    let white_area = area_upm * factor * 100.0;
    let prop_area = (amplitude_y * white_area) / xheight;
    let valor = prop_area - area(polygon);
    valor / amplitude_y
}

fn area(points: &[Point]) -> f64 {
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
    angle: f64,
    xheight: f64,
    scan_frequency: usize,
    depth_cut: f64,
) -> (Vec<Point>, Point, Point, Vec<Point>, Point, Point) {
    // For deskewing angled glyphs. Makes subsequent processing easier.
    let skew_offset = xheight / 2.0;
    let tan_angle = angle.to_radians().tan();

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
        .step_by(scan_frequency)
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
            if angle != 0.0 {
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

    let extreme_left_full = extreme_left_full.unwrap();
    let extreme_left = extreme_left.unwrap();
    let extreme_right_full = extreme_right_full.unwrap();
    let extreme_right = extreme_right.unwrap();

    // Second pass: Cap the margin samples to a maximum depth from the outermost point in to get our depth cut-in.
    let depth = xheight * depth_cut / 100.0;
    let max_depth = extreme_left.x + depth;
    let min_depth = extreme_right.x - depth;
    left.iter_mut().for_each(|s| s.x = s.x.min(max_depth));
    right.iter_mut().for_each(|s| s.x = s.x.max(min_depth));

    // Third pass: Close open counterforms at 45 degrees.
    let dx_max = scan_frequency as f64;

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

// fn draw_glyph_outer_outline_into_glyph(
//     glyph_name: impl Into<GlyphName>,
//     outlines: (&Vec<Point>, &Vec<Point>),
// ) -> Glyph {
//     let mut builder = GlyphBuilder::new(glyph_name, GlifVersion::V2);
//     let mut outline_builder = OutlineBuilder::new();
//     outline_builder.begin_path(None).unwrap();
//     for left in outlines.0 {
//         outline_builder
//             .add_point(
//                 (left.x.round() as f32, left.y.round() as f32),
//                 PointType::Line,
//                 false,
//                 None,
//                 None,
//             )
//             .unwrap();
//     }
//     outline_builder.end_path().unwrap();
//     outline_builder.begin_path(None).unwrap();
//     for right in outlines.1 {
//         outline_builder
//             .add_point(
//                 (right.x.round() as f32, right.y.round() as f32),
//                 PointType::Line,
//                 false,
//                 None,
//                 None,
//             )
//             .unwrap();
//     }
//     outline_builder.end_path().unwrap();
//     let (outline, identifiers) = outline_builder.finish().unwrap();
//     builder.outline(outline, identifiers).unwrap();
//     builder.finish().unwrap()
// }

/// Returns a Vec of decomposed components of a composite. Ignores incoming identifiers and libs
/// and dangling components; contours are in no particular order.
// XXX: deal with cycles?!
fn decomposed_components(glyph: &Glyph, glyphset: &Layer) -> Vec<Contour> {
    let mut contours = Vec::new();
    let mut stack: Vec<(&Component, Affine)> = Vec::new();

    for component in &glyph.components {
        stack.push((component, component.transform.into()));

        while let Some((component, transform)) = stack.pop() {
            let new_glyph = match glyphset.get_glyph(&component.base) {
                Some(g) => g,

                None => continue,
            };

            for contour in &new_glyph.contours {
                let mut decomposed_contour = Contour::default();
                for point in &contour.points {
                    let new_point = transform * Point::new(point.x, point.y);
                    decomposed_contour.points.push(ContourPoint::new(
                        new_point.x,
                        new_point.y,
                        point.typ.clone(),
                        point.smooth,
                        point.name.clone(),
                        None,
                        None,
                    ))
                }
                contours.push(decomposed_contour);
            }

            for new_component in new_glyph.components.iter().rev() {
                let new_transform: Affine = new_component.transform.into();
                stack.push((new_component, transform * new_transform));
            }
        }
    }

    contours
}

fn path_for_glyph(glyph: &Glyph, glyphset: &Layer) -> Result<BezPath, ContourDrawingError> {
    let mut path = BezPath::new();
    for contour in glyph
        .contours
        .iter()
        .chain(decomposed_components(glyph, glyphset).iter())
    {
        for element in contour_segments(contour)? {
            path.push(element);
        }
    }
    Ok(path)
}

#[derive(Debug)]
enum ContourDrawingError {
    IllegalPointCount(PointType, usize),
    IllegalMove,
    TrailingOffCurves,
}

fn contour_segments(contour: &Contour) -> Result<Vec<PathEl>, ContourDrawingError> {
    let mut points: Vec<&ContourPoint> = contour.points.iter().collect();
    let mut segments = Vec::new();

    let closed;
    let start: &ContourPoint;
    let implied_oncurve: ContourPoint;

    // Phase 1: Preparation
    match points.len() {
        // Empty contours cannot be represented by segments.
        0 => return Ok(segments),
        // Single points are converted to open MoveTos because closed single points of any
        // PointType make no sense.
        1 => {
            segments.push(PathEl::MoveTo(Point::new(
                points[0].x as f64,
                points[0].y as f64,
            )));
            return Ok(segments);
        }
        // Contours with two or more points come in three flavors...:
        _ => {
            // 1. ... Open contours begin with a Move. Start the segment on the first point
            // and don't close it. Note: Trailing off-curves are an error.
            if let PointType::Move = points[0].typ {
                closed = false;
                // Pop off the Move here so the segmentation loop below can just error out on
                // encountering any other Move.
                start = points.remove(0);
            } else {
                closed = true;
                // 2. ... Closed contours begin with anything else. Locate the first on-curve
                // point and rotate the point list so that it _ends_ with that point. The first
                // point could be a curve with its off-curves at the end; moving the point
                // makes always makes all associated off-curves reachable in a single pass
                // without wrapping around. Start the segment on the last point.
                if let Some(first_oncurve) =
                    points.iter().position(|e| e.typ != PointType::OffCurve)
                {
                    points.rotate_left(first_oncurve + 1);
                    start = points.last().unwrap();
                // 3. ... Closed all-offcurve quadratic contours: Rare special case of
                // TrueType's “implied on-curve points” principle. Compute the last implied
                // on-curve point and append it, so we can handle this normally in the loop
                // below. Start the segment on the last, computed point.
                } else {
                    let first = points.first().unwrap();
                    let last = points.last().unwrap();
                    implied_oncurve = ContourPoint::new(
                        0.5 * (last.x + first.x),
                        0.5 * (last.y + first.y),
                        PointType::QCurve,
                        false,
                        None,
                        None,
                        None,
                    );
                    points.push(&implied_oncurve);
                    start = &implied_oncurve;
                }
            }
        }
    }

    // Phase 1.5: Always need a MoveTo as the first element.
    segments.push(PathEl::MoveTo(Point::new(start.x as f64, start.y as f64)));

    // Phase 2: Conversion
    let mut controls: Vec<Point> = Vec::new();
    for point in points {
        let p = Point::new(point.x as f64, point.y as f64);
        match point.typ {
            PointType::OffCurve => controls.push(p),
            // The first Move is removed from the points above, any other Move we encounter is illegal.
            PointType::Move => return Err(ContourDrawingError::IllegalMove),
            // A line must have 0 off-curves preceeding it.
            PointType::Line => match controls.len() {
                0 => segments.push(PathEl::LineTo(p)),
                _ => {
                    return Err(ContourDrawingError::IllegalPointCount(
                        PointType::Line,
                        controls.len(),
                    ))
                }
            },
            // A quadratic curve can have any number of off-curves preceeding it. Zero means it's
            // a line, numbers > 1 mean we must expand “implied on-curve points”.
            PointType::QCurve => match controls.len() {
                0 => segments.push(PathEl::LineTo(p)),
                1 => {
                    segments.push(PathEl::QuadTo(controls[0], p));
                    controls.clear()
                }
                _ => {
                    // TODO: make iterator? controls.iter().zip(controls.iter().cycle().skip(1))
                    for i in 0..=controls.len() - 2 {
                        let c = controls[i];
                        let cn = controls[i + 1];
                        let pi = Point::new(0.5 * (c.x + cn.x), 0.5 * (c.y + cn.y));
                        segments.push(PathEl::QuadTo(c, pi));
                    }
                    segments.push(PathEl::QuadTo(controls[controls.len() - 1], p));
                    controls.clear()
                }
            },
            // A curve can have 0, 1 or 2 off-curves preceeding it according to the UFO specification.
            // Zero means it's a line, one means it's a quadratic curve, two means it's a cubic curve.
            PointType::Curve => match controls.len() {
                0 => segments.push(PathEl::LineTo(p)),
                1 => {
                    segments.push(PathEl::QuadTo(controls[0], p));
                    controls.clear()
                }
                2 => {
                    segments.push(PathEl::CurveTo(controls[0], controls[1], p));
                    controls.clear()
                }
                _ => {
                    return Err(ContourDrawingError::IllegalPointCount(
                        PointType::Curve,
                        controls.len(),
                    ))
                }
            },
        }
    }
    // If we have control points left at this point, we are an open contour, which must end on
    // an on-curve point.
    if !controls.is_empty() {
        debug_assert!(!closed);
        return Err(ContourDrawingError::TrailingOffCurves);
    }
    if closed {
        segments.push(PathEl::ClosePath);
    }

    Ok(segments)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn space_mutatorsans() {
        let ufo = norad::Font::load("testdata/mutatorSans/MutatorSansLightWide.ufo").unwrap();
        let default_layer = ufo.default_layer();

        let units_per_em: f64 = ufo
            .font_info
            .units_per_em
            .map(|v| v.get())
            .unwrap_or(1000.0);
        let angle: f64 = ufo
            .font_info
            .italic_angle
            .map(|v| -v.get())
            .unwrap_or_default();
        let xheight: f64 = ufo.font_info.x_height.map(|v| v.get()).unwrap_or_default();
        let param_area: f64 = 400.0;
        let param_depth: f64 = 15.0;
        let param_overshoot: f64 = 0.0;
        let overshoot = xheight * param_overshoot / 100.0;
        let param_sample_frequency: usize = 5;

        // let mut background_glyphs = Vec::new();

        // TODO: make own config instead of polluting global one
        for (name, left, right) in &[
            ("A", Some(23.0), Some(23.0)),
            ("Aacute", Some(23.0), Some(23.0)),
            ("Adieresis", Some(23.0), Some(23.0)),
            ("B", Some(92.0), Some(43.0)),
            ("C", Some(49.0), Some(43.0)),
            ("D", Some(92.0), Some(49.0)),
            ("E", Some(92.0), Some(33.0)),
            ("F", Some(92.0), Some(32.0)),
            ("G", Some(49.0), Some(66.0)),
            ("H", Some(92.0), Some(92.0)),
            ("I", Some(33.0), Some(33.0)),
            ("I.narrow", Some(92.0), Some(92.0)),
            ("IJ", Some(71.0), Some(72.0)),
            ("J", Some(41.0), Some(75.0)),
            ("J.narrow", Some(24.0), Some(72.0)),
            ("K", Some(92.0), Some(25.0)),
            ("L", Some(92.0), Some(25.0)),
            ("M", Some(92.0), Some(92.0)),
            ("N", Some(92.0), Some(92.0)),
            ("O", Some(49.0), Some(49.0)),
            ("P", Some(92.0), Some(46.0)),
            ("Q", Some(49.0), Some(49.0)),
            ("R", Some(92.0), Some(49.0)),
            ("S", Some(38.0), Some(44.0)),
            ("S.closed", Some(43.0), Some(42.0)),
            ("T", Some(25.0), Some(25.0)),
            ("U", Some(72.0), Some(72.0)),
            ("V", Some(23.0), Some(23.0)),
            ("W", Some(26.0), Some(26.0)),
            ("X", Some(19.0), Some(25.0)),
            ("Y", Some(22.0), Some(22.0)),
            ("Z", Some(33.0), Some(33.0)),
            ("acute", Some(79.0), Some(79.0)),
            ("arrowdown", Some(89.0), Some(91.0)),
            ("arrowleft", Some(95.0), Some(111.0)),
            ("arrowright", Some(110.0), Some(96.0)),
            ("arrowup", Some(91.0), Some(89.0)),
            ("colon", Some(88.0), Some(88.0)),
            ("comma", Some(94.0), Some(91.0)),
            ("dieresis", Some(80.0), Some(80.0)),
            ("dot", Some(80.0), Some(80.0)),
            ("period", Some(96.0), Some(96.0)),
            ("quotedblbase", Some(94.0), Some(91.0)),
            ("quotedblleft", Some(91.0), Some(94.0)),
            ("quotedblright", Some(94.0), Some(91.0)),
            ("quotesinglbase", Some(94.0), Some(91.0)),
            ("semicolon", Some(88.0), Some(86.0)),
            ("space", None, None),
        ] {
            let glyph = ufo.get_glyph(*name).unwrap();

            let (glyph_ref_name, factor) = config::config_for_glyph(&glyph.name);
            let glyph_ref = ufo.get_glyph(glyph_ref_name).unwrap();

            let paths = path_for_glyph(glyph, default_layer).unwrap();
            let bounds = paths.bounding_box();
            let paths_reference = path_for_glyph(glyph_ref, default_layer).unwrap();
            let bounds_reference = paths_reference.bounding_box();
            let bounds_reference_lower = (bounds_reference.min_y() - overshoot).round();
            let bounds_reference_upper = (bounds_reference.max_y() + overshoot).round();

            let (new_left, new_right) = calculate_spacing(
                paths,
                bounds,
                (bounds_reference_lower, bounds_reference_upper),
                angle,
                xheight,
                param_sample_frequency,
                param_depth,
                // name.clone(),
                // &mut background_glyphs,
                factor,
                param_area,
                units_per_em,
            );

            match (left, new_left) {
                (Some(v), Some(new_v)) => assert!(
                    (*v - new_v).abs() <= 1.0,
                    "Glyph {}: expected left {} but got {} (factor {})",
                    *name,
                    v,
                    new_v,
                    factor
                ),
                (None, None) => (),
                _ => panic!(
                    "Glyph {}, left side: expected {:?}, got {:?} (factor {})",
                    *name, left, new_left, factor
                ),
            }
            match (right, new_right) {
                (Some(v), Some(new_v)) => assert!(
                    (*v - new_v).abs() <= 1.0,
                    "Glyph {}: expected right {} but got {} (factor {})",
                    *name,
                    v,
                    new_v,
                    factor
                ),
                (None, None) => (),
                _ => panic!(
                    "Glyph {}, right side: expected {:?}, got {:?} (factor {})",
                    *name, right, new_right, factor
                ),
            }
        }
    }
}
