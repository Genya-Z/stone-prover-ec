#ifndef NETHERMIND_EC_FFT_WITH_PRECOMPUTE_H_
#define NETHERMIND_EC_FFT_WITH_PRECOMPUTE_H_

#include <algorithm>
#include <utility>
#include <vector>

#include "starkware/algebra/fft/details.h"

//namespace fft_tuning_params {
//const size_t kPrecomputeDepth = 22;  // Found empirically though benchmarking.
//}  // namespace fft_tuning_params

#if 0
class FftWithPrecomputeBase {
 public:
  virtual ~FftWithPrecomputeBase() = default;
  virtual void ShiftTwiddleFactors(const FieldElement& offset, const FieldElement& prev_offset) = 0;
};

/*
  Make this a mock if FftWithPrecomputeBase becomes more complex.
*/
class DummyFftWithPrecompute : public FftWithPrecomputeBase {
  void ShiftTwiddleFactors(
      const FieldElement& /*offset*/, const FieldElement& /*prev_offset*/) override {}
};
#endif

// Notation based on
// https://arxiv.org/pdf/2107.08473.pdf and
// https://www.math.toronto.edu/swastik/ECFFT2.pdf
template <typename FieldElementT>
class EcFftWithPrecompute : public starkware::FftWithPrecomputeBase {
  using BasesT = EcFftBases<FieldElementT>;
  using PointT = starkware::EcPoint<FieldElementT>;

 public:
  explicit EcFftWithPrecompute(BasesT bases);

  void SFft(gsl::span<const FieldElementT> src, gsl::span<FieldElementT> dst) const;
  void MFft(gsl::span<const FieldElementT> src, gsl::span<FieldElementT> dst, size_t level = 0) const;
  void SIFft(gsl::span<const FieldElementT> src, gsl::span<FieldElementT> dst) const;
  void MIFft(gsl::span<const FieldElementT> src, gsl::span<FieldElementT> dst, size_t level = 0) const;

  const std::vector<FieldElementT>& GetTwiddleFactors() const { return twiddle_factors_; }
  const std::vector<FieldElementT>& GetOmega() const { return omega_; }
  const std::vector<FieldElementT>& GetIOmega() const { return iomega_; }
  const std::vector<FieldElementT>& GetZeta() const { return zeta_; }
  const std::vector<FieldElementT>& GetIZeta() const { return izeta_; }

  // Shifts the twiddle factors by c, to accommodate for evaluation.
  void ShiftTwiddleFactors(const starkware::FieldElement& /*offset*/, const starkware::FieldElement& /*prev_offset*/) override {
    //ASSERT_RELEASE(false,"method not implemented");
/*    if (twiddle_factors_.empty()) {
      return;
    }
    fft::details::ParallelFromOtherTwiddle<FieldElementT, BasesT>(
        BasesT::GroupT::GroupOperation(
            offset.As<FieldElementT>(),
            BasesT::GroupT::GroupOperationInverse(prev_offset.As<FieldElementT>())),
        bases_, twiddle_factors_);*/
  }

 private:
  //void FftNaturalOrder(gsl::span<const FieldElementT> src, gsl::span<FieldElementT> dst) const;
  //void FftReversedOrder(gsl::span<const FieldElementT> src, gsl::span<FieldElementT> dst) const;

  const BasesT bases_;
  std::vector<FieldElementT> twiddle_factors_;
  std::vector<FieldElementT> omega_;
  std::vector<FieldElementT> iomega_;
  std::vector<FieldElementT> zeta_;
  std::vector<FieldElementT> izeta_;
  bool icun_;
};

#include "nethermind/ec_fft_with_precompute.inl"

#endif  // NETHERMIND_EC_FFT_WITH_PRECOMPUTE_H_
