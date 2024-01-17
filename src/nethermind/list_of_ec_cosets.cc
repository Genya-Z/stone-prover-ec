#include "nethermind/list_of_ec_cosets.h"
#include "nethermind/ec_data.h"

#include <utility>

#include "starkware/algebra/field_operations.h"
#include "starkware/algebra/utils/invoke_template_version.h"


namespace {

template <typename FieldElementT>
std::vector<starkware::FieldElement> GetCosetsOffsets(
    const size_t n_cosets, const starkware::EcPoint<FieldElementT>& domain_generator,
    const starkware::EcPoint<FieldElementT>& common_offset, const EC<FieldElementT>& ec) {
  // Define result vector.
  std::vector<starkware::FieldElement> result;
  result.reserve(n_cosets);

  // Compute the offsets vector.
  starkware::EcPoint<FieldElementT> offset = common_offset;
  result.emplace_back(offset, ec);
  for (size_t i = 1; i < n_cosets; i++) {
    offset = ec.addPoints(offset, domain_generator);
    result.emplace_back(offset, ec);
  }

  return result;
}

}  // namespace

//template<typename FieldElementT>
inline ListOfEcCosets::ListOfEcCosets(
    std::unique_ptr<FftBases>&& fft_bases, std::vector<FieldElement> cosets_offsets,
    FieldElement trace_generator, const Field& field)
    : fft_bases_(std::move(fft_bases)),
      cosets_offsets_(std::move(cosets_offsets)),
      trace_generator_(std::move(trace_generator)),
      field_(std::move(field)) {}

//template<typename FieldElementT>
/*inline ListOfEcCosets::ListOfEcCosets(std::vector<FieldElement> cosets_offsets,
    FieldElement trace_generator, const Field& field)
    : cosets_offsets_(std::move(cosets_offsets)),
      trace_generator_(std::move(trace_generator)),
      field_(std::move(field)) {}*/

ListOfEcCosets ListOfEcCosets::MakeListOfEcCosets(
    size_t coset_size, size_t n_cosets, const EcData& ec_data) {
  Field field = ec_data.GetField();
  // Compute generator of group containing all cosets.
  ASSERT_RELEASE(n_cosets > 0, "Number of cosets must be positive.");
  ASSERT_RELEASE(coset_size > 1, "Coset size must be > 1.");
  const size_t log_size = starkware::SafeLog2(coset_size);

  auto build_eval_domain = [&](auto field_tag) {
    using FieldElementT = typename decltype(field_tag)::type;

    EC<FieldElementT> ec=ec_data.GetCurve<FieldElementT>();
    //FieldElementT coset_generator = FieldElementT::Uninitialized();
    //EcPoint<FieldElementT> trace_generator = {FieldElementT::Uninitialized(),
    //                                          FieldElementT::Uninitialized()};
    //FieldElementT offset = FieldElementT::Uninitialized();

    // Multiplicative case.
    uint64_t power_of_two_cosets = starkware::Pow2(starkware::Log2Ceil(n_cosets));

    starkware::EcPoint<FieldElementT> coset_generator = ec_data.GetSubGroupGenerator<FieldElementT>(coset_size * power_of_two_cosets);
    ASSERT_RELEASE(power_of_two_cosets%2==0,"power_of_two_coset not multiple of 2");
    starkware::EcPoint<FieldElementT> offset = ec.MultiplyByScalar(coset_generator, power_of_two_cosets/2);
    starkware::EcPoint<FieldElementT> trace_generator = ec.MultiplyByScalar(offset, 2);
    //offset = FieldElementT::One();

    EcFftBases<FieldElementT> bases(trace_generator, log_size, offset, ec);
    return ListOfEcCosets(
        std::make_unique<decltype(bases)>(std::move(bases)),
        GetCosetsOffsets(n_cosets, coset_generator, ec_data.Generator<FieldElementT>(), ec),
        FieldElement(trace_generator, ec), field);
    //return ListOfEcCosets(GetCosetsOffsets(n_cosets, coset_generator, ec_data.Generator<FieldElementT>(), ec),
    //    FieldElement(trace_generator,ec),field);
  };
  return InvokeFieldTemplateVersion(build_eval_domain, field);
}

starkware::FieldElement ListOfEcCosets::ElementByIndex(size_t coset_index, size_t group_index) const {
  ASSERT_RELEASE(coset_index < cosets_offsets_.size(), "Coset index out of range.");

  FieldElement point = Bases().At(0).GetFieldElementAt(group_index);
  //FieldElement offset = cosets_offsets_[coset_index];

  return point /* offset*/;
}
