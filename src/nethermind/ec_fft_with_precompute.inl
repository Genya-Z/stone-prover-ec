#include "nethermind/ec_fft_with_precompute.h"

//This is quadratic, should ideally be loaded from a file
template <typename FieldElementT>
void compute_omega(const EcFftBases<FieldElementT>& bases,
                   std::vector<FieldElementT> & out)
{
  size_t length = starkware::Pow2(bases[0].BasisSize()-1) - 1;
  std::vector<FieldElementT> tors;
  tors.reserve(length);
  for(size_t j = 0; j < length; j++)
    tors.push_back((bases[0][j+1] - bases[0].StartOffset()).x);
  
  bool icun = bases[0].IsClosedUnderNegation();
  length = starkware::Pow2(bases[0].BasisSize()-(icun ? 1 : 0));
  out.reserve(length);
  for(size_t i = 0; i < length; i++)
  {
    FieldElementT omega = FieldElementT::One();
    FieldElementT val = bases[0][i].x;
    for(FieldElementT tor : tors)
    {
      omega *= (val - tor);
    }
    out.push_back(omega);
  }

}

template <typename FieldElementT>
EcFftWithPrecompute<FieldElementT>::EcFftWithPrecompute(const BasesT bases)
  : bases_(std::move(bases))
{
  compute_omega(bases_,omega_);

  this->icun_ = bases_[0].IsClosedUnderNegation();

  size_t length = starkware::Pow2(bases_[0].BasisSize()-(icun_ ? 1 : 0));
  this->iomega_.reserve(length);
  for(FieldElementT val : omega_)
    this->iomega_.push_back(val.Inverse());

  FieldElementT two_tor = bases_[0].TwoTor();
  this->zeta_.reserve(length);
  for(size_t i = 0; i < length; i++) {
    PointT pt=bases_[0][i];
    this->zeta_.push_back(pt.y/(pt.x - two_tor));
  }

  this->izeta_.reserve(length);
  for(FieldElementT val : zeta_)
    this->izeta_.push_back(val.Inverse());

}

template <typename FieldElementT>
void EcFftWithPrecompute<FieldElementT>::SFft(
    const gsl::span<const FieldElementT> src, const gsl::span<FieldElementT> dst) const
{
  ASSERT_DEBUG(icun_,"Coset not closed under negation.");
  size_t half_size = starkware::Pow2(bases_[0].BasisSize()-1);
  ASSERT_DEBUG(src.size() == 2*half_size,"wrong src size for FFT");
  ASSERT_DEBUG(dst.size() == 2*half_size,"wrong dst size for FFT");
  std::vector<FieldElementT> h0s = FieldElementT::UninitializedVector(half_size);
  std::vector<FieldElementT> h1s = FieldElementT::UninitializedVector(half_size);

  FieldElementT half = FieldElementT::FromUint(2).Inverse(); //absorb into omega?
  for (size_t j = 0; j < half_size; j++) {
    FieldElementT omega=GetOmega()[j];
    FieldElementT izeta=GetIZeta()[j];
    FieldElementT h0=omega*(src[j] + src[2*half_size-j-1])*half;
    FieldElementT h1=izeta*omega*(src[j] - src[2*half_size-j-1])*half;
    h0s[inv_grey(j)]=h0;
    h1s[inv_grey(j)]=h1;
  }
  MFft(h0s,dst.subspan(0,half_size));
  MFft(h1s,dst.subspan(half_size,half_size));
}

template <typename FieldElementT>
void EcFftWithPrecompute<FieldElementT>::SIFft(
    const gsl::span<const FieldElementT> src, const gsl::span<FieldElementT> dst) const
{
  size_t half_size = starkware::Pow2(bases_[0].BasisSize()-1);
  ASSERT_DEBUG(src.size() == 2*half_size,"wrong size for evaluation");
  size_t x_size = icun_ ? half_size : 2*half_size;
  std::vector<FieldElementT> pi0 = FieldElementT::UninitializedVector(x_size);
  std::vector<FieldElementT> pi1 = FieldElementT::UninitializedVector(x_size);
  MIFft(src.subspan(0,half_size),pi0);
  MIFft(src.subspan(half_size,half_size),pi1);
  for(size_t j = 0; j < x_size; j++) {
    FieldElementT iomega=GetIOmega()[j];
    FieldElementT zeta=GetZeta()[j];
    FieldElementT h0=pi0[icun_ ? inv_grey(j) : j];
    FieldElementT zh1=zeta*pi1[icun_ ? inv_grey(j) : j];
    dst[j]=(h0+zh1)*iomega;
    if(icun_ || dst.size() == 2*x_size)
      dst[2*half_size-j-1]=(h0-zh1)*iomega;
  }
}

template <typename FieldElementT>
void EcFftWithPrecompute<FieldElementT>::MFft(
    const gsl::span<const FieldElementT> src, const gsl::span<FieldElementT> dst, size_t level) const
{
  ASSERT_DEBUG(icun_,"Coset not closed under negation.");
  EcFftDomain<FieldElementT> S = bases_[level];
  ASSERT_DEBUG(2*src.size() == S.Size(),"wrong src size");
  ASSERT_DEBUG(dst.size() == src.size(),"wrong dst size");
  if(src.size() == 1) {
    dst[0] = src[0];
    return;
  }
  EcFftDomain<FieldElementT> T = bases_[level + 1];
  //v=self.get_map(i).denom()
  size_t sT = starkware::Pow2(T.BasisSize()-1);
  FieldElementT two_tor = S.TwoTor();
  std::vector<FieldElementT> pi0=FieldElementT::UninitializedVector(sT);
  std::vector<FieldElementT> pi1=FieldElementT::UninitializedVector(sT);
  for(size_t j = 0; j < sT; j++) {
    FieldElementT sj = S.GetFieldElementAt(j,icun_);
    ASSERT_DEBUG(S.ApplyTwoIsogeny(sj)==T.GetFieldElementAt(j,icun_),"Two Isogeny doesn't match");
    FieldElementT vsj=starkware::Pow(sj - two_tor, sT-1);

    FieldElementT sj2 = S.GetFieldElementAt(j + sT,icun_);
    ASSERT_DEBUG(S.ApplyTwoIsogeny(sj2)==T.GetFieldElementAt(j,icun_),"Two Isogeny doesn't match");
    FieldElementT vsj2=starkware::Pow(sj2 - two_tor, sT-1);

    FieldElementT idet=((sj2-sj)*vsj*vsj2).Inverse();
    pi0[j]=(src[j]*sj2*vsj2 - src[j+sT]*sj*vsj)*idet;
    pi1[j]=(-src[j]*vsj2 + src[j+sT]*vsj)*idet;
    //pi0[j],pi1[j]=self.M(pi[j],pi[j+sT])
  }
  MFft(pi0,dst.subspan(0,sT),level + 1);
  MFft(pi1,dst.subspan(sT,sT),level + 1);
}

template <typename FieldElementT>
void EcFftWithPrecompute<FieldElementT>::MIFft(
    const gsl::span<const FieldElementT> src, const gsl::span<FieldElementT> dst, size_t level) const
{
  EcFftDomain<FieldElementT> S = bases_[level];
  ASSERT_DEBUG(2*src.size() == S.Size(),"wrong src size");
  if(src.size() == 1) {
    for(size_t i = 0; i < dst.size(); i++ )
      dst[i]=src[0];
    return;
  }
  EcFftDomain<FieldElementT> T = bases_[level + 1];
  size_t sT = starkware::Pow2(T.BasisSize()-1);
  size_t sTp = icun_ ? sT : 2*sT;
  ASSERT_DEBUG(dst.size() == 2*sTp,"wrong dst size");

  std::vector<FieldElementT> pi0=FieldElementT::UninitializedVector(sTp);
  std::vector<FieldElementT> pi1=FieldElementT::UninitializedVector(sTp);
  MIFft(src.subspan(0,sT),pi0,level + 1);
  MIFft(src.subspan(sT,sT),pi1,level + 1);
  FieldElementT two_tor = S.TwoTor();
  for(size_t j = 0; j < sTp; j++) {
    FieldElementT sj=S.GetFieldElementAt(j,icun_);
    ASSERT_DEBUG(S.ApplyTwoIsogeny(sj)==T.GetFieldElementAt(j,icun_),"Two Isogeny doesn't match");
    FieldElementT vsj=starkware::Pow(sj - two_tor, sT-1);
    dst[j]=(pi0[j] + sj*pi1[j])*vsj;

    FieldElementT sj2 = S.GetFieldElementAt(j + sTp,icun_);
    ASSERT_DEBUG(S.ApplyTwoIsogeny(sj2)==T.GetFieldElementAt(j,icun_),"Two Isogeny doesn't match");
    FieldElementT vsj2=starkware::Pow(sj2 - two_tor, sT-1);
    dst[j+sTp]=(pi0[j] + sj2*pi1[j])*vsj2;
  }
}

#if 0
template <typename BasesT>
void FftWithPrecompute<BasesT>::FftNaturalOrder(
    const gsl::span<const FieldElementT> src, const gsl::span<FieldElementT> dst) const {
  gsl::span<const FieldElementT> curr_src = src;
  size_t precompute_depth = PrecomputeDepth();
  const size_t last_precomputed_layer_size =  // The succeeding layers use FftNoPrecompute.
      Pow2(precompute_depth);

  if (src.size() == 1) {
    dst[0] = src[0];
    return;
  }

  bool full_precompute = src.size() <= last_precomputed_layer_size;
  if (last_precomputed_layer_size > 1) {
    for (size_t i = 0; i < src.size(); i += last_precomputed_layer_size) {
      fft::details::FftUsingPrecomputedTwiddleFactors<FieldElementT>(
          curr_src.subspan(i, last_precomputed_layer_size), twiddle_factors_,
          /*normalize=*/full_precompute, dst.subspan(i, last_precomputed_layer_size));
    }
    curr_src = dst;
  }

  if (!full_precompute) {
    fft::details::FftNoPrecompute<BasesT>(
        curr_src, bases_,
        /*layers_to_skip=*/precompute_depth, dst);
  }
}

template <typename BasesT>
void FftWithPrecompute<BasesT>::FftReversedOrder(
    const gsl::span<const FieldElementT> src, const gsl::span<FieldElementT> dst) const {
  ASSERT_RELEASE(
      twiddle_factors_.size() + 1 == src.size() || src.size() == 1,
      "only full precompute is currently supported");
  fft::details::FftNaturalToReverseWithPrecompute<FieldElementT>(src, twiddle_factors_, dst);
}
#endif //0