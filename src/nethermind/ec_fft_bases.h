#ifndef NETHERMIND_EC_FFT_BASES_H_
#define NETHERMIND_EC_FFT_BASES_H_

#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "starkware/algebra/fft/multiplicative_group_ordering.h"
#include "starkware/fft_utils/fft_domain.h"
#include "starkware/fft_utils/fft_bases.h"
#include "starkware/algebra/elliptic_curve/elliptic_curve.h"
#include "nethermind/ec_fft_domain.h"

// This may already exist, but I don't know where.
template <typename T>
inline std::vector<T> subvector(std::vector<T> v, size_t index)
{
    return std::vector<T>{v.begin() + index, v.end()};
}

/*
  Creates a series of EC FFT domains.
  The i'th domain is psi_i of the (i-1)'st one.

  The domains are ordered according to the Order paramater.
*/
template <typename FieldElementT>
class EcFftBases
    : public starkware::FftBases
{

public:
    using DomainT = EcFftDomain<FieldElementT>;
    using EcPointT = starkware::EcPoint<FieldElementT>;
    using EcT = EC<FieldElementT>;

    EcFftBases(
        const EcPointT &generator,
        size_t log_n,
        const EcPointT &start_offset,
        const EcT &start_curve);
    
    EcFftBases(const DomainT &base_domain);

    /*
     Returns the number of layers in the instance, not including the last empty domain at the end.
   */
    size_t NumLayers() const override
    {
        return bases_.size() - 1;
    }

    EcFftBases FromLayer(size_t idx) const
    {
        ASSERT_RELEASE(idx < bases_.size(), "index out of range");
        return EcFftBases(subvector(bases_, idx));
    }

    std::unique_ptr<FftBases> FromLayerAsUniquePtr(size_t idx) const override
    {
        return std::make_unique<EcFftBases>(FromLayer(idx));
    }

    const DomainT &operator[](size_t idx) const { return bases_.at(idx); }
    starkware::Field GetField() const override { return starkware::Field::Create<FieldElementT>(); }

    /*
      Same as operator[]. This is more readable when the object is given as a pointer.
    */
    const starkware::FftDomainBase &At(size_t idx) const override { return (*this)[idx]; }

    /*
      Returns an FftBasesImpl instance that is derived from the original FftBasesImpl
      by changing the offsets in all the layers.

      The offset at layer i is obtained from the offset at layer i-1 using a 2 to 1 mapping.

      The result is independent of the offset in the original FftBasesImpl instance.
    */
    EcFftBases GetShiftedBases(const EcPointT &offset) const
    {
        DomainT b = bases_[0];
        return EcFftBases(b.Basis()[0],b.BasisSize(),offset,b.Curve());
        //return EcFftBases(bases_[0].GetShiftedDomain(offset));
    }

    std::unique_ptr<FftBases> GetShiftedBasesAsUniquePtr(const starkware::FieldElement &offset) const override
    {
        return std::make_unique<EcFftBases>(GetShiftedBases(offset.AsEc<FieldElementT>()));
    }

    FieldElementT ApplyBasisTransformTmpl(const FieldElementT &point, size_t ind) const
    {
        return bases_[ind].ApplyTwoIsogeny(point);
    }

    std::tuple<std::unique_ptr<FftBases>, std::vector<starkware::FieldElement>> SplitToCosets(
        size_t n_log_cosets) const override;

    starkware::FieldElement ApplyBasisTransform(const starkware::FieldElement& point, size_t layer_index) const override {
        ASSERT_RELEASE(layer_index < this->NumLayers(), "Layer index out of range.");
        return starkware::FieldElement(
             ApplyBasisTransformTmpl(point.As<FieldElementT>(), layer_index));
    }


private:
    std::vector<DomainT> bases_;

    /*
      Private so no one will give us a bad domain.
    */
    explicit EcFftBases(std::vector<DomainT> bases)
        : bases_(std::move(bases))
    {
        ASSERT_RELEASE(
            !bases_.empty() && bases_.back().BasisSize() == 0, "bases must end in an empty domain");
    }
};

#include "ec_fft_bases.cc"

#endif // NETHERMIND_EC_FFT_BASES_H_
