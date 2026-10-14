use crate::astrometry::matcher::{AstrometricMatch, AstrometricSolution};
use crate::astrometry::quality::{
    evaluate_astrometric_quality, AstrometricQuality, ReductionStatus,
};
use crate::astrometry::wcs::Wcs;
use crate::reduction::SourceMeasurement;
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct PlateSolveResult {
    pub run_id: Option<String>,
    pub success: bool,
    pub status: ReductionStatus,
    pub failure_code: Option<String>,
    pub wcs: Option<Wcs>,
    pub num_matched: u32,
    pub num_catalog: u32,
    pub residual_rms: Option<f64>,
    pub backend: Option<String>,
    pub message: String,
    pub matches: Vec<AstrometricMatch>,
    pub quality: Option<AstrometricQuality>,
    pub manual_review_confirmed: bool,
}

pub fn missing_hint(message: &str) -> PlateSolveResult {
    PlateSolveResult {
        run_id: None,
        success: false,
        status: ReductionStatus::Rejected,
        failure_code: Some("missing_hint".into()),
        wcs: None,
        num_matched: 0,
        num_catalog: 0,
        residual_rms: None,
        backend: None,
        message: message.to_string(),
        matches: Vec::new(),
        quality: None,
        manual_review_confirmed: false,
    }
}

pub fn match_failed(num_catalog: u32, message: String) -> PlateSolveResult {
    PlateSolveResult {
        run_id: None,
        success: false,
        status: ReductionStatus::Rejected,
        failure_code: Some("match_failed".into()),
        wcs: None,
        num_matched: 0,
        num_catalog,
        residual_rms: None,
        backend: Some("triangle invariants + robust TAN/CD".to_string()),
        message: format!("归算失败：{message}"),
        matches: Vec::new(),
        quality: None,
        manual_review_confirmed: false,
    }
}

pub fn solved(
    num_catalog: u32,
    solution: AstrometricSolution,
    sources: &[SourceMeasurement],
) -> PlateSolveResult {
    let mut quality = evaluate_astrometric_quality(
        &solution.matches,
        sources,
        solution.wcs.image_width,
        solution.wcs.image_height,
    );
    quality.plate_model_requested = solution.plate_model_requested.as_str().to_string();
    quality.plate_model_applied = solution.plate_model_applied.as_str().to_string();
    quality.sip_order = solution.sip_order;
    quality.sip_escalation_attempted = solution.sip_escalation_attempted;
    quality.sip_escalation_succeeded = solution.sip_escalation_succeeded;
    quality.sip_fit_converged = solution.sip_fit_converged;
    quality.plate_model_downgraded = solution.plate_model_downgraded;
    quality.plate_model_downgrade_reason = solution.plate_model_downgrade_reason.clone();
    if let Some(reason) = &solution.sip_rejection_reason {
        quality.status = ReductionStatus::Rejected;
        quality.reasons.push(reason.clone());
    }

    let accepted = quality.status == ReductionStatus::Accepted;
    let sip_label = solution.sip_order.map(|order| match order {
        2 => "二阶",
        3 => "三阶",
        _ => "SIP",
    });
    let downgrade_note = solution
        .plate_model_downgrade_reason
        .as_deref()
        .map(|reason| format!("；板模型已降级：{reason}"))
        .unwrap_or_default();
    let rejection_reason = solution.sip_rejection_reason.clone();
    let failure_code = (!accepted).then(|| {
        match rejection_reason.as_deref() {
            Some("implausible_distortion") => "implausible_distortion",
            Some("sip_fit_not_converged") => "sip_fit_not_converged",
            Some("distortion_suspected") => "distortion_suspected",
            _ if quality.distortion_suspected => "distortion_suspected",
            _ => "quality_gate",
        }
        .to_string()
    });
    let message = match quality.status {
        ReductionStatus::Accepted => match sip_label {
            Some(label) => format!(
                "归算通过（{label} SIP 畸变校正）：匹配 {} 颗参考星，RMS {:.3} arcsec，P95 {:.3} arcsec。{downgrade_note}",
                solution.matches.len(),
                quality.residual_rms_arcsec,
                quality.residual_p95_arcsec
            ),
            None => format!(
                "归算通过：匹配 {} 颗参考星，RMS {:.3} arcsec，P95 {:.3} arcsec。{downgrade_note}",
                solution.matches.len(),
                quality.residual_rms_arcsec,
                quality.residual_p95_arcsec
            ),
        },
        ReductionStatus::ReviewRequired => format!(
            "归算需要复核：匹配 {} 颗参考星，RMS {:.3} arcsec；{}{downgrade_note}",
            solution.matches.len(),
            quality.residual_rms_arcsec,
            quality.reasons.join("；")
        ),
        ReductionStatus::Rejected => match rejection_reason.as_deref() {
            Some("implausible_distortion") => {
                format!("归算被拒绝：拟合出的畸变模型不物理（角点偏移过大或径向修正非单调）{downgrade_note}")
            }
            Some("sip_fit_not_converged") => {
                format!("归算被拒绝：SIP 畸变拟合未收敛，未采用不确定的畸变模型{downgrade_note}")
            }
            _ => format!("归算被拒绝：{}{downgrade_note}", quality.reasons.join("；")),
        },
    };
    let num_matched = solution.matches.len() as u32;
    let residual_rms = Some(solution.rms_arcsec);
    PlateSolveResult {
        run_id: None,
        success: accepted,
        status: quality.status,
        failure_code,
        wcs: Some(solution.wcs),
        num_matched,
        num_catalog,
        residual_rms,
        backend: Some(
            "extended Delaunay / hinted pair voting + iterative robust TAN/CD".to_string(),
        ),
        message,
        matches: solution.matches,
        quality: Some(quality),
        manual_review_confirmed: false,
    }
}

impl PlateSolveResult {
    /// Promote only a reviewable solution after an explicit operator action.
    /// Rejected solutions remain rejected, so weak geometry or excessive
    /// residuals cannot be bypassed through the manual-alignment workflow.
    pub fn confirm_review(&mut self) {
        if self.status != ReductionStatus::ReviewRequired || self.wcs.is_none() {
            return;
        }
        self.success = true;
        self.status = ReductionStatus::Accepted;
        self.failure_code = None;
        self.manual_review_confirmed = true;
        if let Some(quality) = &mut self.quality {
            quality.status = ReductionStatus::Accepted;
        }
        self.message = format!("人工复核通过：{}", self.message);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astrometry::matcher::PlateModel;

    #[test]
    fn sip_rejection_cannot_be_accepted_or_manually_confirmed() {
        let sources: Vec<_> = (0..16)
            .map(|index| SourceMeasurement {
                x: 80.0 + (index % 4) as f64 * 250.0,
                y: 80.0 + (index / 4) as f64 * 250.0,
                peak: 100.0,
                flux: 1_000.0,
                fwhm: 3.0,
                ellipticity: 0.1,
                npix: 20,
                flags: 0,
                saturated: false,
                snr: Some(30.0),
                x_error_px: Some(0.05),
                y_error_px: Some(0.05),
                centroid_refined: true,
            })
            .collect();
        for count in [8, 16] {
            let matches: Vec<_> = (0..count)
                .map(|index| AstrometricMatch {
                    image_source_index: index,
                    catalog_source_index: index,
                    residual_arcsec: 0.12,
                    residual_x_arcsec: 0.08,
                    residual_y_arcsec: 0.09,
                    weight: 1.0,
                    used: true,
                    rejection_reason: None,
                })
                .collect();
            let ordinary = evaluate_astrometric_quality(&matches, &sources, 1_000, 1_000);
            assert_eq!(
                ordinary.status,
                if count == 8 {
                    ReductionStatus::ReviewRequired
                } else {
                    ReductionStatus::Accepted
                }
            );
            for reason in [
                "sip_fit_not_converged",
                "implausible_distortion",
                "distortion_suspected",
            ] {
                let solution = AstrometricSolution {
                    wcs: Wcs {
                        crpix1: 500.0,
                        crpix2: 500.0,
                        crval1: 180.0,
                        crval2: 30.0,
                        cd1_1: -1.0 / 3_600.0,
                        cd1_2: 0.0,
                        cd2_1: 0.0,
                        cd2_2: 1.0 / 3_600.0,
                        image_width: 1_000,
                        image_height: 1_000,
                        sip: None,
                    },
                    matches: matches.clone(),
                    rms_arcsec: 0.12,
                    plate_model_requested: PlateModel::Cubic,
                    plate_model_applied: PlateModel::Linear,
                    sip_order: None,
                    sip_escalation_attempted: true,
                    sip_escalation_succeeded: false,
                    sip_fit_converged: Some(reason != "sip_fit_not_converged"),
                    plate_model_downgraded: false,
                    plate_model_downgrade_reason: None,
                    sip_rejection_reason: Some(reason.into()),
                };
                let mut result = solved(16, solution, &sources);
                result.confirm_review();
                assert!(!result.success);
                assert!(!result.manual_review_confirmed);
                assert_eq!(result.status, ReductionStatus::Rejected);
                assert_eq!(result.failure_code.as_deref(), Some(reason));
                let quality = result.quality.unwrap();
                assert_eq!(quality.status, ReductionStatus::Rejected);
                assert!(quality.reasons.iter().any(|item| item == reason));
            }
        }
    }
}
