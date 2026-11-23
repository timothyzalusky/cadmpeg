// SPDX-License-Identifier: Apache-2.0
//! Cross-codec ground-truth verification for Pair5 revolve specimens (`EV-25`).

#![cfg(all(feature = "sldprt", feature = "step"))]

use std::fs;
use std::path::PathBuf;

use cadmpeg_codec_sldprt::SldprtCodec;
use cadmpeg_codec_step::StepCodec;
use cadmpeg_core::decode::{DecodeArena, DecodeContext, DecodePolicy};
use cadmpeg_ir::codec::{Codec, DecodeOptions};
use cadmpeg_ir::features::{
    AngularTermination, FeatureOperation, RevolveConstruction, RevolveExtent,
};
use cadmpeg_ir::sketches::SketchGeometryDefinition;
use cadmpeg_measure::{MeasureInput, MeasureOptions, Measurer, Selection};

fn corpus_path(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../corpus/sldprt")
        .join(relative)
}

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
    let RevolveConstruction::Resolved { axis, extent, .. } = construction else {
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

#[test]
fn pair5_revolve_swept_region_matches_step_ground_truth_ev25() {
    for stem in ["Pair5/revolve_90_default", "Pair5/revolve_90_flipped"] {
        let mut sldprt_file =
            fs::File::open(corpus_path(&format!("{stem}.SLDPRT"))).expect("SLDPRT specimen opens");
        let sldprt = SldprtCodec
            .decode(&mut sldprt_file, &DecodeOptions::default())
            .expect("SLDPRT decodes");

        let mut step_file =
            fs::File::open(corpus_path(&format!("{stem}.STEP"))).expect("STEP specimen opens");
        let step = StepCodec::default()
            .decode(&mut step_file, &DecodeOptions::default())
            .expect("STEP decodes");

        let (swept_min, swept_max) = revolve_swept_bounds(sldprt.ir());
        let step_bbox = measure_first_body_bbox(step.ir());
        let step_min = [
            step_bbox.min.get().x,
            step_bbox.min.get().y,
            step_bbox.min.get().z,
        ];
        let step_max = [
            step_bbox.max.get().x,
            step_bbox.max.get().y,
            step_bbox.max.get().z,
        ];
        for axis in 0..3 {
            assert!(
                (swept_min[axis] - step_min[axis]).abs() < 1e-3,
                "{stem}: swept min[{axis}]={} != step min[{axis}]={}",
                swept_min[axis],
                step_min[axis]
            );
            assert!(
                (swept_max[axis] - step_max[axis]).abs() < 1e-3,
                "{stem}: swept max[{axis}]={} != step max[{axis}]={}",
                swept_max[axis],
                step_max[axis]
            );
        }
    }
}
