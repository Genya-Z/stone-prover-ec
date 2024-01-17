
#ifndef NETHERMIND_EC_FFT_DOMAIN_H_
#define NETHERMIND_EC_FFT_DOMAIN_H_

#include <memory>
#include <stack>
#include <tuple>
#include <utility>
#include <vector>

#include "third_party/gsl/gsl-lite.hpp"

#include "starkware/algebra/polymorphic/field_element.h"
#include "starkware/algebra/elliptic_curve/elliptic_curve.h"
#include "starkware/fft_utils/fft_domain.h"
#include "starkware/math/math.h"
#include "nethermind/ec_group.h"

//using namespace starkware; // NOLINT

/* Grey code, A006068 at oeis */
inline size_t grey(size_t n)
{
    size_t s = 1;
    while (true)
    {
        size_t ns = n >> s;
        if (ns == 0)
            break;
        n ^= ns;
        s <<= 1;
    }
    return n;
}

inline size_t inv_grey(size_t n)
{
    return n ^ n>>1;
}

/*
    An fft domain on the projection of an elliptic curve onto the x-axis.
    This is not the same as the domain of the elliptic curve group,
    which is why we can't just use FftDomain<EC>
*/
template <typename FieldElementT>
class EcFftDomain : public starkware::FftDomainBase
{
public:
    using EcPointT = starkware::EcPoint<FieldElementT>;
    using EcT = EC<FieldElementT>;

    EcFftDomain(std::vector<EcPointT> basis, const EcPointT &start_offset,
                const EcT &curve)
        : basis_(std::move(basis)), start_offset_(start_offset),
        curve_(curve) {}

    const std::vector<EcPointT> &Basis() const { return basis_; }
    const EcPointT &StartOffset() const { return start_offset_; }
    const EcT &Curve() const { return curve_; }
    const FieldElementT &TwoTor() const {return basis_.back().x; }

    EcPointT ApplyTwoIsogeny(const EcPointT &pt) const
    {
        return curve_.TwoIsogeny(TwoTor(),pt);
    }

    FieldElementT ApplyTwoIsogeny(const FieldElementT &pt) const
    {
        return curve_.TwoIsogeny(TwoTor(),pt);
    }

    bool IsClosedUnderNegation() const
    {
        return curve_.Double(start_offset_) == basis_[0];
    }

    /*
        Returns a new instance of FftDomain with the same basis as the original domain,
        but with a different offset.

        The offset in the original domain is ignored.
    */
    EcFftDomain GetShiftedDomain(const EcPointT& offset) const {
      return EcFftDomain(basis_, offset, curve_);
    }


    EcPointT operator[](uint64_t index) const
    {
        EcPointT ret = start_offset_;
        ASSERT_VERIFIER(index < Size(), "Index out of range.");
        for (const auto &b : basis_)
        {
            if ((index & 1) == 1)
            {
                ret += b;
            }
            index >>= 1;
        }
        return ret;
    }

    // precompute?
    FieldElementT GetFieldElementAt(uint64_t idx,bool use_grey) const
    {
        ASSERT_RELEASE(idx < Size(), "Index out of range.");
        // Fix for bit reverse order?
        return (operator[](use_grey ? grey(idx) : idx)).x;
    }

    starkware::FieldElement GetFieldElementAt(uint64_t idx) const override
    {
        return starkware::FieldElement(GetFieldElementAt(idx,IsClosedUnderNegation()));
    }

    uint64_t BasisSize() const override { return basis_.size(); }

    /*
      Returns the number of elements in the domain.
    */
    uint64_t Size() const override { return starkware::Pow2(basis_.size()); }

    std::unique_ptr<FftDomainBase> RemoveFirstBasisElementsAsUniquePtr(size_t n) const override
    {
        return std::make_unique<EcFftDomain>(RemoveFirstBasisElements(n));
    }

    std::unique_ptr<FftDomainBase> RemoveLastBasisElementsAsUniquePtr(size_t ) const override
    {
        ASSERT_RELEASE(
            false, "RemoveLastBasisElementsAsUniqueP is not implemented for EcFftDomain");
    }

    /*
      See the documentation of FftDomainBase::RemoveFirstBasisElementsAsUniquePtr.
    */
    EcFftDomain RemoveFirstBasisElements(size_t n) const
    {
        ASSERT_DEBUG(n <= basis_.size(), "index out of range");
        return EcFftDomain({basis_.begin() + n, basis_.end()}, start_offset_, curve_);
    }

private:
    const std::vector<EcPointT> basis_;
    const EcPointT start_offset_;
    const EcT curve_;
    //const FieldElementT two_tor_;

    using BasisIteratorType = decltype(basis_.begin());
};

template <typename FieldElementT, typename EcPointT = starkware::EcPoint<FieldElementT>>
EcFftDomain<FieldElementT>
MakeEcFftDomain(const EC<FieldElementT> &curve,
                const EcPointT &generator, size_t log_n, const EcPointT &start_offset /*,
                 bool reversed_order = true*/
)
{
    std::vector<EcPointT> basis;
    EcPointT g = generator;
    basis.reserve(log_n);
    for (size_t i = 0; i < log_n - 1; ++i)
    {
        basis.push_back(g);
        g = curve.Double(g);
    }
    basis.push_back(g);
    ASSERT_RELEASE(g.y == FieldElementT::Zero(), "generator order is not Pow2(log_n)");
    // if (reversed_order) {
    //   std::reverse(basis.begin(), basis.end());
    // }
    return EcFftDomain(std::move(basis), start_offset, curve);
}

#endif // NETHERMIND_EC_FFT_DOMAIN_H_
