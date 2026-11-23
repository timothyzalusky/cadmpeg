// SPDX-License-Identifier: Apache-2.0
//! Regression tests over the controlled specimen pairs in `corpus/sldprt/`.

use std::fs;
use std::path::PathBuf;

use cadmpeg_core::decode::{DecodeArena, DecodeContext, DecodePolicy};
use cadmpeg_ir::codec::{Codec, DecodeOptions};
use cadmpeg_ir::features::{AngularTermination, ExtrudeDirection, FeatureOperation, RevolveExtent};
use cadmpeg_ir::sketches::SketchGeometryDefinition;
use cadmpeg_measure::{MeasureInput, MeasureOptions, Measurer, Selection};

use crate::SldprtCodec;

fn measure_first_body_bbox(ir: &cadmpeg_ir::CadIr) -> cadmpeg_measure::value::BoxValue {
    let arena = DecodeArena::new();
    let policy = DecodePolicy::service();
    let ctx = DecodeContext::for_loaded_input(&arena, &policy, 0).expect("measure context");
    let input = MeasureInput {
        ir,
        fidelity: None,
        geometry_transferred: true,
    };
    let report = Measurer::new(input, MeasureOptions::default())
        .measure(&ctx, &Selection::default())
        .expect("measured report");
    ctx.finish_session().expect("measure session");
    *report.bodies()[0]
        .bbox()
        .measured()
        .expect("measured bbox")
        .value()
}

fn rodrigues_rotate(
    point: [f64; 3],
    origin: [f64; 3],
    unit_axis: [f64; 3],
    angle: f64,
) -> [f64; 3] {
    let rel_point = [
        point[0] - origin[0],
        point[1] - origin[1],
        point[2] - origin[2],
    ];
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let axis_dot =
        unit_axis[0] * rel_point[0] + unit_axis[1] * rel_point[1] + unit_axis[2] * rel_point[2];
    let axis_cross = [
        unit_axis[1] * rel_point[2] - unit_axis[2] * rel_point[1],
        unit_axis[2] * rel_point[0] - unit_axis[0] * rel_point[2],
        unit_axis[0] * rel_point[1] - unit_axis[1] * rel_point[0],
    ];
    let one_minus_cos = 1.0 - cos_a;
    [
        origin[0]
            + rel_point[0] * cos_a
            + axis_cross[0] * sin_a
            + unit_axis[0] * axis_dot * one_minus_cos,
        origin[1]
            + rel_point[1] * cos_a
            + axis_cross[1] * sin_a
            + unit_axis[1] * axis_dot * one_minus_cos,
        origin[2]
            + rel_point[2] * cos_a
            + axis_cross[2] * sin_a
            + unit_axis[2] * axis_dot * one_minus_cos,
    ]
}

fn revolve_swept_bounds(ir: &cadmpeg_ir::CadIr) -> ([f64; 3], [f64; 3]) {
    const STEPS: u32 = 360;

    let sketch = &ir.model.sketches[0];
    let (origin, normal, u_axis) = sketch.resolved_placement().expect("resolved placement");
    let normal_vec = normal.get();
    let u_vec = u_axis.get();
    let v_vec = [
        normal_vec.y * u_vec.z - normal_vec.z * u_vec.y,
        normal_vec.z * u_vec.x - normal_vec.x * u_vec.z,
        normal_vec.x * u_vec.y - normal_vec.y * u_vec.x,
    ];
    let origin_pt = origin.get();
    let to_3d = |pt: cadmpeg_ir::math::Point2| {
        [
            origin_pt.x + pt.u * u_vec.x + pt.v * v_vec[0],
            origin_pt.y + pt.u * u_vec.y + pt.v * v_vec[1],
            origin_pt.z + pt.u * u_vec.z + pt.v * v_vec[2],
        ]
    };

    let mut sketch_points = Vec::new();
    for entity in &ir.model.sketch_entities {
        match entity.geometry.definition() {
            SketchGeometryDefinition::Point { position } => sketch_points.push(to_3d(**position)),
            SketchGeometryDefinition::Line { start, end } => {
                sketch_points.push(to_3d(**start));
                sketch_points.push(to_3d(**end));
            }
            _ => {}
        }
    }
    assert!(!sketch_points.is_empty(), "expected sketch profile points");

    let revolve = ir
        .model
        .features
        .iter()
        .find(|feature| feature.name.as_deref() == Some("Revolve1"))
        .expect("Revolve1");
    let FeatureOperation::Revolve { construction, .. } =
        revolve.evaluation.definition().operation()
    else {
        panic!("expected Revolve");
    };
    let cadmpeg_ir::features::RevolveConstruction::Resolved { axis, extent, .. } = construction
    else {
        panic!("expected resolved revolve");
    };
    let RevolveExtent::OneSided {
        termination: AngularTermination::Angle { angle },
    } = extent
    else {
        panic!("expected one-sided angular revolve");
    };
    let axis_origin = [
        axis.origin.get().x,
        axis.origin.get().y,
        axis.origin.get().z,
    ];
    let axis_dir = [
        axis.direction.get().x,
        axis.direction.get().y,
        axis.direction.get().z,
    ];
    let total_angle = angle.get();

    let distance_to_axis = |pt: [f64; 3]| {
        let delta = [
            pt[0] - axis_origin[0],
            pt[1] - axis_origin[1],
            pt[2] - axis_origin[2],
        ];
        let proj = delta[0] * axis_dir[0] + delta[1] * axis_dir[1] + delta[2] * axis_dir[2];
        let perp = [
            delta[0] - proj * axis_dir[0],
            delta[1] - proj * axis_dir[1],
            delta[2] - proj * axis_dir[2],
        ];
        (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt()
    };
    let off_axis_points: Vec<[f64; 3]> = sketch_points
        .into_iter()
        .filter(|pt| distance_to_axis(*pt) > 1e-6)
        .collect();
    assert!(
        !off_axis_points.is_empty(),
        "expected off-axis sketch profile points"
    );

    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for step in 0..=STEPS {
        let theta = total_angle * f64::from(step) / f64::from(STEPS);
        for &pt in &off_axis_points {
            let rotated = rodrigues_rotate(pt, axis_origin, axis_dir, theta);
            for idx in 0..3 {
                min[idx] = min[idx].min(rotated[idx]);
                max[idx] = max[idx].max(rotated[idx]);
            }
        }
    }
    (min, max)
}

fn decode_specimen(relative: &str) -> cadmpeg_ir::codec::DecodeResult {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../corpus/sldprt")
        .join(relative);
    let mut file = fs::File::open(&path).expect("specimen file exists");
    SldprtCodec
        .decode(&mut file, &DecodeOptions::default())
        .expect("specimen decodes")
}

#[test]
fn pair1_front_flip_decodes_reversed_profile_normal_di16() {
    let blind = decode_specimen("Pair1/plane_front_blind.SLDPRT");
    let flip = decode_specimen("Pair1/plane_front_flip.SLDPRT");

    for (result, expected_reversed) in [(&blind, false), (&flip, true)] {
        let extrude = result
            .ir()
            .model
            .features
            .iter()
            .find(|feature| feature.name.as_deref() == Some("Boss-Extrude1"))
            .expect("Boss-Extrude1");
        let FeatureOperation::Extrude { direction, .. } =
            extrude.evaluation.definition().operation()
        else {
            panic!("expected Extrude operation");
        };
        if expected_reversed {
            assert!(matches!(
                direction,
                ExtrudeDirection::ReversedProfileNormal {}
            ));
        } else {
            assert!(matches!(direction, ExtrudeDirection::ProfileNormal {}));
        }
        assert!(
            !result.report().losses.iter().any(|loss| {
                loss.code.to_string() == "sldprt/feature.extrude-direction-unresolved"
            }),
            "unexpected extrude-direction-unresolved loss"
        );
    }
}

#[test]
fn pair1_right_plane_retains_front_fallback_frame_sf02() {
    let right = decode_specimen("Pair1/plane_right_blind.SLDPRT");
    let sketch = &right.ir().model.sketches[0];
    let (_, normal, u_axis) = sketch.resolved_placement().expect("resolved placement");
    assert_eq!(normal.get().x, 0.0);
    assert_eq!(normal.get().y, -1.0);
    assert_eq!(normal.get().z, 0.0);
    assert_eq!(u_axis.get().x, 0.0);
    assert_eq!(u_axis.get().y, 0.0);
    assert_eq!(u_axis.get().z, -1.0);

    let bbox = measure_first_body_bbox(right.ir());
    assert!((bbox.min.get().x - 0.0).abs() < 1e-6);
    assert!((bbox.min.get().y - 0.0).abs() < 1e-6);
    assert!((bbox.min.get().z + 20.0).abs() < 1e-6);
    assert!((bbox.max.get().x - 10.0).abs() < 1e-6);
    assert!((bbox.max.get().y - 10.0).abs() < 1e-6);
    assert!((bbox.max.get().z - 0.0).abs() < 1e-6);
}

#[test]
fn pair3_concentric_circles_decode_both_radii_di70_di71() {
    for path in [
        "Pair3/sketch_concentric_one_dim.SLDPRT",
        "Pair3/sketch_concentric_ring.SLDPRT",
    ] {
        let result = decode_specimen(path);
        let mut radii: Vec<f64> = result
            .ir()
            .model
            .sketch_entities
            .iter()
            .filter_map(|entity| match entity.geometry.definition() {
                SketchGeometryDefinition::Circle { radius, .. } => Some(radius.get()),
                _ => None,
            })
            .collect();
        radii.sort_by(f64::total_cmp);
        radii.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
        assert_eq!(radii.len(), 2, "{path}: expected both r=5 and r=10 circles");
        assert!((radii[0] - 5.0).abs() < 1e-6);
        assert!((radii[1] - 10.0).abs() < 1e-6);
    }
}

#[test]
fn pair4_arcs_cw_ccw_decodes_ffff_arc_ev26() {
    let result = decode_specimen("Pair4/sketch_arcs_cw_ccw.SLDPRT");
    let has_r10_arc = result.ir().model.sketch_entities.iter().any(|entity| {
        matches!(
            entity.geometry.definition(),
            SketchGeometryDefinition::Arc { center, radius, .. }
                if (center.u - 5.0).abs() < 1e-4
                    && (center.v - 18.660_254).abs() < 1e-4
                    && (radius.get() - 10.0).abs() < 1e-4
        )
    });
    assert!(
        has_r10_arc,
        "expected 0xffff arc centered at (5, 18.660254) with r=10; got: {:#?}",
        result
            .ir()
            .model
            .sketch_entities
            .iter()
            .map(|e| (e.id().to_string(), e.geometry.definition()))
            .collect::<Vec<_>>()
    );
}

#[test]
fn pair5_revolve_90_default_vs_flipped_direction_ev25() {
    let default = decode_specimen("Pair5/revolve_90_default.SLDPRT");
    let flipped = decode_specimen("Pair5/revolve_90_flipped.SLDPRT");

    let axis_y = |result: &cadmpeg_ir::codec::DecodeResult| {
        let revolve = result
            .ir()
            .model
            .features
            .iter()
            .find(|feature| feature.name.as_deref() == Some("Revolve1"))
            .expect("Revolve1");
        let FeatureOperation::Revolve { construction, .. } =
            revolve.evaluation.definition().operation()
        else {
            panic!("expected Revolve");
        };
        let cadmpeg_ir::features::RevolveConstruction::Resolved { axis, .. } = construction else {
            panic!("expected resolved revolve");
        };
        axis.direction.get().y
    };

    assert!((axis_y(&default) - 1.0).abs() < 1e-9);
    assert!((axis_y(&flipped) + 1.0).abs() < 1e-9);

    for (label, result) in [("default", &default), ("flipped", &flipped)] {
        let (swept_min, swept_max) = revolve_swept_bounds(result.ir());
        let body_bbox = measure_first_body_bbox(result.ir());
        let body_min = [
            body_bbox.min.get().x,
            body_bbox.min.get().y,
            body_bbox.min.get().z,
        ];
        let body_max = [
            body_bbox.max.get().x,
            body_bbox.max.get().y,
            body_bbox.max.get().z,
        ];
        for axis in 0..3 {
            assert!(
                (swept_min[axis] - body_min[axis]).abs() < 1e-3,
                "{label}: swept min[{axis}]={} != body min[{axis}]={}",
                swept_min[axis],
                body_min[axis]
            );
            assert!(
                (swept_max[axis] - body_max[axis]).abs() < 1e-3,
                "{label}: swept max[{axis}]={} != body max[{axis}]={}",
                swept_max[axis],
                body_max[axis]
            );
        }
    }

    assert!(
        !flipped
            .report()
            .losses
            .iter()
            .any(|loss| loss.code.to_string() == "sldprt/feature.input-archive-walk-stopped"),
        "flipped revolve should not stop the CArchive walk"
    );
}
