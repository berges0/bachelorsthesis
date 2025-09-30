#ifndef METRIC_SIMPLIFICATION_HPP
#define METRIC_SIMPLIFICATION_HPP
#include "utils/pch.hpp"

namespace METRIC_SIMPLIFICATION {

// Kombination für symmetrische Metrik
enum class CombineMode {
    LengthWeightedMean,  // empfohlen
    ArithmeticMean,
    MaxOfDirected
};

// Steuerparameter für Orthogonalenlänge (relativ zur gemeinsamen BBox-Diagonale)
struct RangeParams {
    // Modus A (Standard): r = clamp(alpha_len * |edge|, min_factor*diag, max_factor*diag)
    double alpha_len   = 0.5;   ///< Basis: r = alpha_len * |edge|
    double min_factor  = 1e-6;  ///< min  = min_factor * diag(BBox)
    double max_factor  = 0.05;  ///< max  = max_factor * diag(BBox)

    // Modus B (fix): r = fixed_range_abs  ODER  r = fixed_range_factor * diag(BBox)
    bool   use_fixed_range     = false;  ///< true => Kantenlänge ignorieren
    double fixed_range_abs     = 0.0;    ///< absolute fixe Reichweite (Einheiten deiner Daten)
    double fixed_range_factor  = 0.0;    ///< falls >0: r = fixed_range_factor * diag(BBox)
};

// Visualisierungs-Parameter (prozentual + Mindestovershoot)
struct VizParams {
    double overshoot_factor      = 0.10; ///< Anteil des Trefferabstands (z.B. 0.10 = +10 %)
    double min_overshoot_factor  = 1e-5; ///< Mindestovershoot = faktor * diag(BBox)
    bool   draw_misses           = false;///< bei Miss: Strahl bis r zeichnen
};

// ===================== Metrik ohne Viz =====================

double weighted_angle_deviation_bounded(
    const std::vector<PWH>& Aset,
    const std::vector<PWH>& Bset,
    bool penalize_miss,
    const RangeParams& rp = {}
);

double weighted_angle_deviation_bounded_symmetric(
    const std::vector<PWH>& Aset,
    const std::vector<PWH>& Bset,
    bool penalize_miss,
    const RangeParams& rp = {},
    CombineMode mode = CombineMode::LengthWeightedMean
);

// ===================== Metrik MIT Viz =====================

double weighted_angle_deviation_bounded_with_viz(
    const std::vector<PWH>& Aset,
    const std::vector<PWH>& Bset,
    bool penalize_miss,
    const RangeParams& rp,
    const VizParams& vz,
    std::vector<Segment>& out_viz_AB // Segmente der "gewinnenden" Normalen (A->B)
);

double weighted_angle_deviation_bounded_symmetric_with_viz(
    const std::vector<PWH>& Aset,
    const std::vector<PWH>& Bset,
    bool penalize_miss,
    const RangeParams& rp,
    const VizParams& vz,
    std::vector<Segment>& out_viz_AB, // A->B
    std::vector<Segment>& out_viz_BA, // B->A
    CombineMode mode = CombineMode::LengthWeightedMean
);

} // namespace METRIC_SIMPLIFICATION

#endif