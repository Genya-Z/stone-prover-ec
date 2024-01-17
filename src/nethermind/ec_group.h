#ifndef NETHERMIND_EC_GROUP_H_
#define NETHERMIND_EC_GROUP_H_

#include "starkware/algebra/elliptic_curve/elliptic_curve.h"

template <typename FieldElementT>
struct EC
{
public:
    using EcPointT = starkware::EcPoint<FieldElementT>;

    bool ContainsPoint(const EcPointT &P) const{
        return P.IsOnCurve(alpha, beta);}

    EcPointT addPoints(const EcPointT &lhs, const EcPointT &rhs) const
    {
        ASSERT_DEBUG(ContainsPoint(lhs),"Point not on curve");
        ASSERT_DEBUG(ContainsPoint(rhs),"Point not on curve");
        if (lhs.x != rhs.x)
            return lhs + rhs;
        else if (lhs.y == -rhs.y)
            ASSERT_RELEASE(false, "EcPoint class can't handle point at infinity.");
        // return 0;
        else
            return lhs.Double(alpha);
    }

    EcPointT Double(const EcPointT &P) const
    {
        ASSERT_DEBUG(ContainsPoint(P),"Point not on curve");
        return P.Double(alpha);
    }

    template <size_t N>
    EcPointT MultiplyByScalar(const EcPointT &P, const starkware::BigInt<N>& scalar) const
    {
        ASSERT_DEBUG(ContainsPoint(P),"Point not on curve");
        return P.MultiplyByScalar(scalar,alpha);
    }

    EcPointT MultiplyByScalar(const EcPointT &P, uint64_t scalar) const
    {
        return MultiplyByScalar(P,starkware::BigInt<1>(scalar));
    }

    /* The image of pt under the two isogeny whose kernel is (a,0)*/
    EcPointT TwoIsogeny(const FieldElementT &a,const EcPointT &pt) const
    {
        ASSERT_DEBUG(ContainsPoint({a,FieldElementT::Zero()}),"Not valid two torsion point");
        ASSERT_DEBUG(ContainsPoint(pt),"Point not on curve");
        FieldElementT idenom = (pt.x - a).Inverse();
        FieldElementT a2 = a*a;
        FieldElementT alpTwoA2 = alpha + a2 + a2;
        return {pt.x + (alpTwoA2 + a2)*idenom,
                pt.y*((pt.x - (a + a))*pt.x - alpTwoA2)*idenom*idenom};
    }

    /* The image of x under the two isogeny whose kernel is (a,0)*/
    FieldElementT TwoIsogeny(const FieldElementT &a,const FieldElementT &x) const
    {
        ASSERT_DEBUG(ContainsPoint({a,FieldElementT::Zero()}),"Not valid two torsion point");
        FieldElementT idenom = (x - a).Inverse();
        FieldElementT a2 = a*a;
        return x + (alpha + a2 + a2 + a2)*idenom;
    }

    /* The codomain of the two isogeny whose kernel is (a,0)*/
    EC<FieldElementT> TwoIsogenyCodomain(const FieldElementT &a) const
    {
        ASSERT_DEBUG(ContainsPoint({a,FieldElementT::Zero()}),"Not valid two torsion point");
        return {-alpha*FieldElementT::FromUint(4) - a*a*FieldElementT::FromUint(15),
                beta*FieldElementT::FromUint(22) + alpha*a*FieldElementT::FromUint(14)};
    }

    EcPointT Random(starkware::Prng* prng) const
    {
        return EcPointT::Random(alpha, beta, prng);
    }

    bool operator==(const EC& rhs) const { return alpha == rhs.alpha && beta == rhs.beta; }
    bool operator!=(const EC& rhs) const { return !(*this == rhs); }

    FieldElementT alpha;
    FieldElementT beta;
};

/*template <typename FieldElementT>
struct TwoIsogeny
{
public:
    using EcPointT = EcPoint<FieldElementT>;

    EcPontT operator()(const EcPointT &pt)
    {
        FieldElementT idenom = (pt.x + c).Inverse();
        return {((pt.x + a) * pt.x + b) * idenom,
                pt.y * ((pt.x + d) * pt.x + e) * idenom * idenom};
    }

    FieldElementT operator()(const FieldElementT &x)
    {
        FieldElementT idenom = (x + c).Inverse();
        return ((x + a) * x + b) * idenom;
    }

    bool CheckValid(EC<FieldElementT> &source, EC<FieldElementT> &target)
    {
        return target.ContainsPoint(this(source.Random()));
    }

    // x' = (x^2 + a*x + b)/(x + c)
    // y' = y*(x^2 + d*x + e)/(x + c)^2

    const FieldElementT a;
    const FieldElementT b;
    const FieldElementT c;
    const FieldElementT d;
    const FieldElementT e;
};*/

#endif // NETHERMIND_EC_FFT_BASES_H_
