#include "starkware/algebra/fft/fft_with_precompute.h"
#include "starkware/algebra/fft/multiplicative_group_ordering.h"
#include "starkware/algebra/field_operations.h"
#include "starkware/algebra/lde/multiplicative_lde.h"
#include "starkware/algebra/polynomials.h"

template <typename FieldElementT>
auto EcLde<FieldElementT>::AddFromCoefficients(
    std::vector<FieldElementT>&& coefficients) -> EcLde {
  return EcLde(std::move(coefficients));
}

template <typename FieldElementT>
auto EcLde<FieldElementT>::AddFromEvaluation(
    const BasesT& bases, std::vector<FieldElementT>&& evaluation,
    starkware::FftWithPrecomputeBase* fft_precomputed) -> EcLde {
  starkware::MaybeOwnedPtr<EcFftWithPrecompute<FieldElementT>> maybe_precomputed =
      UseOwned(dynamic_cast<EcFftWithPrecompute<FieldElementT>*>(fft_precomputed));
  if (fft_precomputed == nullptr) {
    maybe_precomputed = UseMovedValue(EcFftWithPrecompute<FieldElementT>(bases));
  }

  if (bases.NumLayers() > 0) {

    (*maybe_precomputed).SFft(evaluation, evaluation);

/*    using DualBasesT = decltype(GetDualBases(bases));
    MaybeOwnedPtr<FftWithPrecompute<DualBasesT>> maybe_precomputed =
        UseOwned(dynamic_cast<FftWithPrecompute<DualBasesT>*>(fft_precomputed));
    if (fft_precomputed == nullptr) {
      maybe_precomputed = UseMovedValue(FftWithPrecompute<DualBasesT>(GetDualBases(bases)));
    }
    (*maybe_precomputed).Fft(evaluation, evaluation);

    auto lde_size_inverse = FieldElementT::FromUint(evaluation.size()).Inverse();

    TaskManager& task_manager = TaskManager::GetInstance();
    task_manager.ParallelFor(
        evaluation.size(),
        [&evaluation, lde_size_inverse](const TaskInfo& task_info) {
          for (size_t i = task_info.start_idx; i < task_info.end_idx; i++) {
            evaluation[i] *= lde_size_inverse;
          }
        },
        evaluation.size());
*/
  }
  return AddFromCoefficients(std::move(evaluation));
}

template <typename FieldElementT>
void EcLde<FieldElementT>::EvalAtCoset(
    const EcFftWithPrecompute<FieldElementT>& fft_precompute, gsl::span<FieldElementT> result) const {
  fft_precompute.SIFft(polynomial_, result);
}

/*template <typename FieldElementT>
void EcLde<FieldElementT>::EvalAtPoints(
    gsl::span<FieldElementT> points, gsl::span<FieldElementT> outputs) const {
  switch (Order) {
    case MultiplicativeGroupOrdering::kBitReversedOrder:
      OptimizedBatchHornerEval<FieldElementT>(points, polynomial_, outputs);
      break;
    case MultiplicativeGroupOrdering::kNaturalOrder:
      // Natural Order Lde stores polynomial_ in BitReversed order.
      BatchHornerEvalBitReversed<FieldElementT>(points, polynomial_, outputs);
      break;
  }
}*/

/*template <typename FieldElementT>
int64_t EcLde<FieldElementT>::GetDegree() const {
  // Natural Order Lde stores polynomial_ in BitReversed order.
  bool bit_revese = Order == MultiplicativeGroupOrdering::kNaturalOrder;
  size_t log_n = SafeLog2(polynomial_.size());
  for (int64_t deg = polynomial_.size() - 1; deg >= 0; deg--) {
    if (polynomial_[bit_revese ? BitReverse(deg, log_n) : deg] != FieldElementT::Zero()) {
      return deg;
    }
  }
  return -1;
}*/

template <typename FieldElementT>
EcFftWithPrecompute<FieldElementT>
EcLde<FieldElementT>::FftPrecompute(
    const BasesT& bases, const PointT& offset_compensation,
    const PointT& new_offset) {
  return EcFftWithPrecompute<FieldElementT>(bases.GetShiftedBases(new_offset + offset_compensation));
}

template <typename FieldElementT>
std::unique_ptr<starkware::FftWithPrecomputeBase> EcLde<FieldElementT>::IfftPrecompute(
    const BasesT& bases) {
  return std::make_unique<EcFftWithPrecompute<FieldElementT>>(
      EcFftWithPrecompute<FieldElementT>(bases));
}
