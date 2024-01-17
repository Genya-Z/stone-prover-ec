#ifndef NETHERMIND_EC_LDE_H_
#define NETHERMIND_EC_LDE_H_

#include <memory>
#include <utility>
#include <vector>

#include "starkware/algebra/fft/fft_with_precompute.h"
#include "starkware/algebra/fft/multiplicative_group_ordering.h"
#include "nethermind/ec_fft_bases.h"
#include "nethermind/ec_fft_with_precompute.h"

template <typename FieldElementT>
class EcLde {
 public:
  using BasesT = EcFftBases<FieldElementT>;
  using PointT = starkware::EcPoint<FieldElementT>;
  using PrecomputeType = EcFftWithPrecompute<FieldElementT>;
  using T = FieldElementT;

  /*
    Constructs an LDE from the coefficients of the polynomial (obtained by GetCoefficients()).
  */
  static EcLde AddFromCoefficients(std::vector<FieldElementT>&& coefficients);

  /*
    Constructs an LDE from the evaluation of the polynomial on the domain bases[0].
  */
  static EcLde AddFromEvaluation(
      const BasesT& bases, std::vector<FieldElementT>&& evaluation,
      starkware::FftWithPrecomputeBase* fft_precomputed);

  void EvalAtCoset(
      const EcFftWithPrecompute<FieldElementT>& fft_precompute, gsl::span<FieldElementT> result) const;

  //void EvalAtPoints(gsl::span<FieldElementT> points, gsl::span<FieldElementT> outputs) const;

  int64_t GetDegree() const;

  const std::vector<FieldElementT>& GetCoefficients() const { return polynomial_; }

  static EcFftWithPrecompute<FieldElementT> FftPrecompute(
      const BasesT& bases, const PointT& offset_compensation,
      const PointT& new_offset);
  static std::unique_ptr<starkware::FftWithPrecomputeBase> IfftPrecompute(const BasesT& bases);

 private:
  explicit EcLde(std::vector<FieldElementT>&& polynomial)
      : polynomial_(std::move(polynomial)) {}


  std::vector<FieldElementT> polynomial_;
};

#include "nethermind/ec_lde.inl"

#endif  // NETHERMIND_EC_LDE_H_
