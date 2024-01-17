
#ifndef NETHERMIND_LIST_OF_EC_COSETS_H_
#define NETHERMIND_LIST_OF_EC_COSETS_H_

#include <memory>
#include <vector>

#include "starkware/algebra/domains/multiplicative_group.h"
#include "starkware/algebra/domains/list_of_cosets.h"
#include "starkware/error_handling/error_handling.h"
#include "starkware/fft_utils/fft_bases.h"
#include "starkware/math/math.h"

#include "nethermind/ec_fft_bases.h"
#include "nethermind/ec_data.h"

/*
  ListOfEcCosets is a union of cosets of a group.
  Let G be a subgroup of the elliptic curve of size coset_size, then the instance
  represents a set which is a union of the cosets s_0 # G, s_1 # G , ... , s_{n_cosets-1} # G (where
  # is the group opperation).
*/
class ListOfEcCosets : public starkware::ListOfCosetsBase {
 public:
  /*
    Constructs an instance with a group of size coset_size and the number of cosets is n_cosets.
    The cosets offsets are as follows:

    In both cases the offsets are s_0 = c, s_1 = ch, s_2 = c(h^2), ... , s_i = c(h^i).
    Where c is a generator of the field's multiplicative group (in both cases).
    In the multiplicative case h is a generator of a group H such that G is a subgroup of H, and |H|
    is the minimal power of two not smaller than |G|*n_cosets.

    In particular: G is disjoint to all cosets, and in case G is multiplicative and |G|*n_cosets is
    a valid size of a subgroup H, the union of the cosets is a coset of H.

    The order parameter affects the order of Bases(), which represents the domain within each coset.
    The cosest offsets remain unchanged, only the order within each coset.
  */
  using Field = starkware::Field;
  using FieldElement = starkware::FieldElement;
  using FftDomainBase = starkware::FftDomainBase;
  using FftBases = starkware::FftBases;

  static ListOfEcCosets MakeListOfEcCosets(
      size_t coset_size, size_t n_cosets, const EcData& ec_data);

  const FftDomainBase& Group() const { return Bases().At(0); }

  const Field& GetField() const { return field_; }

  //maybe rename FieldElement GroupElement
  const FieldElement& TraceGenerator() const { return trace_generator_; }

  size_t NumCosets() const { return cosets_offsets_.size(); }

  const std::vector<FieldElement>& CosetsOffsets() const { return cosets_offsets_; }

  size_t Size() const { return Group().Size() * cosets_offsets_.size(); }

  const starkware::FftBases& Bases() const { return *fft_bases_; }

  FieldElement ElementByIndex(size_t coset_index, size_t group_index) const;

  // Evaluates the vanishing polynomial of the group.
  FieldElement VanishingPolynomial(const FieldElement& eval_point) const;

 private:
  //template<typename FieldElementT>
  ListOfEcCosets(
      std::unique_ptr<FftBases>&& fft_bases, std::vector<FieldElement> cosets_offsets,
      FieldElement trace_generator, const Field& field);

  //template<typename FieldElementT>
  //ListOfEcCosets(std::vector<FieldElement> cosets_offsets, FieldElement trace_generator, const Field& field);

  std::unique_ptr<FftBases> fft_bases_;
  std::vector<FieldElement> cosets_offsets_;
  FieldElement trace_generator_;
  Field field_;
};

//#include "nethermind/list_of_ec_cosets.cc"

#endif  // NETHERMIND_LIST_OF_EC_COSETS_H_
