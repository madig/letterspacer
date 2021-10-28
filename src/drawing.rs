//! Code for drawing norad Glyphs into kurbo BezPaths.

use kurbo::{Affine, BezPath, PathEl, Point};
use norad::{Component, Contour, ContourPoint, Glyph, Layer, PointType};

pub fn path_for_glyph(glyph: &Glyph, glyphset: &Layer) -> Result<BezPath, ContourDrawingError> {
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

/// Returns a Vec of decomposed components of a composite. Ignores incoming identifiers and libs
/// and dangling components; contours are in no particular order.
// TODO: deal with cycles, extend ContourDrawingError.
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

#[derive(Debug)]
pub enum ContourDrawingError {
    IllegalPointCount(PointType, usize),
    IllegalMove,
    TrailingOffCurves,
}

impl std::error::Error for ContourDrawingError {}

impl std::fmt::Display for ContourDrawingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ContourDrawingError::IllegalPointCount(typ, count) => match typ {
                PointType::Line => write!(
                    f,
                    "A line must have zero off-curves preceeding it, but found {}.",
                    count
                ),
                PointType::Curve => write!(
                    f,
                    "A curve must have zero to two off-curves preceeding it, but found {}.",
                    count
                ),
                _ => unreachable!(),
            },
            ContourDrawingError::IllegalMove => {
                write!(f, "A move must only occur at the start of a contour.")
            }
            ContourDrawingError::TrailingOffCurves => {
                write!(f, "Open contours must not have trailing off-curve points.")
            }
        }
    }
}
