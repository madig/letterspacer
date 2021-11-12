use kurbo::{BezPath, Line, ParamCurve, Point, Rect};

use crate::parameters::SpacingParameters;
use crate::FontMetrics;

pub(crate) fn compute_sidebearings(
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
        + compute_sidebearing_value(
            factor,
            (bounds_reference_lower, bounds_reference_upper),
            parameters.area,
            &left,
            font_metrics,
        ))
    .ceil();
    let new_right = (-distance_right
        + compute_sidebearing_value(
            factor,
            (bounds_reference_lower, bounds_reference_upper),
            parameters.area,
            &right,
            font_metrics,
        ))
    .ceil();

    Some((new_left, new_right))
}

fn compute_sidebearing_value(
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
    let valor = prop_area - compute_area(polygon);
    valor / amplitude_y
}

fn compute_area(points: &[Point]) -> f64 {
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
