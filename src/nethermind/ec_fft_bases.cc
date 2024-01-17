#include <vector>

#include "nethermind/ec_fft_bases.h"

template <typename FieldElementT>
EcFftBases<FieldElementT>::EcFftBases(const DomainT &base_domain)
{
    this->basis_.push_back(base_domain);

    size_t log_n = base_domain.BasisSize();

    EcPointT current_gen = base_domain.Basis()[0];

    EcPointT current_offset = base_domain.StartOffset();

    EcT curve = base_domain.Curve();

    for (size_t i = 1; i < log_n; ++i)
    {
        DomainT domain = MakeEcFftDomain(
            curve, current_gen, log_n - i, current_offset);
        this->bases_.push_back(domain);
        current_offset = domain.ApplyTwoIsogeny(current_offset);
        current_gen = domain.ApplyTwoIsogeny(current_gen);
        curve = curve.TwoIsogenyCodomain(domain.TwoTor());
    }

    this->bases_.emplace_back(std::vector<FieldElementT>{}, current_offset);

}


template <typename FieldElementT>
EcFftBases<FieldElementT>::EcFftBases(const EcPointT &generator,
           size_t log_n,
           const EcPointT &start_offset,
           const EcT &start_curve)
{
    EcPointT current_gen = generator;

    EcPointT current_offset = start_offset;

    EcT curve = start_curve;

    for (size_t i = 0; i < log_n - 1; ++i)
    {
        DomainT domain = MakeEcFftDomain(
            curve, current_gen, log_n - i, current_offset);
        this->bases_.push_back(domain);
        current_offset = domain.ApplyTwoIsogeny(current_offset);
        current_gen = domain.ApplyTwoIsogeny(current_gen);
        curve = curve.TwoIsogenyCodomain(domain.TwoTor());
    }

    this->bases_.emplace_back(std::vector<EcPointT>{current_gen}, current_offset,
                              curve);
}

template <typename FieldElementT>
std::tuple<std::unique_ptr<starkware::FftBases>, std::vector<starkware::FieldElement>>
EcFftBases<FieldElementT>::SplitToCosets(size_t /*n_log_cosets*/) const {
  ASSERT_RELEASE(false,"Method not implemented");
  /*ASSERT_RELEASE(!bases_.empty(), "Can't split empty bases");
  DomainT domain = bases_[0];
  const std::vector<EcPointT>& basis = domain.Basis();
  ASSERT_RELEASE(basis.size() >= n_log_cosets, "Too many cosets requested");

  auto [coset_bases, offsets_domain] = SplitDomain(domain, n_log_cosets);  // NOLINT
  std::vector<FieldElement> offsets({offsets_domain.begin(), offsets_domain.end()});
  ASSERT_RELEASE(offsets.size() == Pow2(n_log_cosets), "Wrong number of offsets");

  return {std::make_unique<EcFftBases>(std::move(coset_bases)), std::move(offsets)};*/
}
