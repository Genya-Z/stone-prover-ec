#include "nethermind/ec_data.h"

class EcData::WrapperBase {
 public:
  virtual ~WrapperBase() = default;
  //virtual std::unique_ptr<WrapperBase> Clone() const = 0;
  virtual starkware::Field GetField() const = 0;
};

template <typename FieldElementT>
class EcData::Wrapper : public WrapperBase {
 public:
  explicit Wrapper(const EC<FieldElementT>& ec,
                   const starkware::EcPoint<FieldElementT>& gen)
    : ec_(ec),
      gen_(gen) {}

  starkware::Field GetField() const override { return starkware::Field::Create<FieldElementT>(); }

  const EC<FieldElementT> ec_;
  const starkware::EcPoint<FieldElementT> gen_;
};

class EcData::IntWrapperBase {
 public:
  virtual ~IntWrapperBase() = default;
  //virtual std::unique_ptr<IntWrapperBase> Clone() const = 0;
};

template <typename IntType>
class EcData::IntWrapper : public IntWrapperBase {
 public:
  explicit IntWrapper(const IntType& val) : val_(val) {}

  const IntType val_;
};

template <typename FieldElementT, typename IntType /*= decltype(FieldElementT::FieldSize())*/>
EcData::EcData(const EC<FieldElementT>& ec,
               const starkware::EcPoint<FieldElementT>& gen,
               const IntType& order)
        : wrapper_(std::make_unique<Wrapper<FieldElementT>>(ec,gen)),
          int_wrapper_(std::make_unique<IntWrapper<IntType>>(order))
          {}

template <typename FieldElementT>
const starkware::EcPoint<FieldElementT>& EcData::Generator() const {
  auto* ptr = dynamic_cast<const Wrapper<FieldElementT>*>(wrapper_.get());
  ASSERT_RELEASE(ptr != nullptr, "The underlying type of FieldElement is wrong");

  return ptr->gen_;
}

template <typename FieldElementT>
const EC<FieldElementT>& EcData::GetCurve() const {
  auto* ptr = dynamic_cast<const Wrapper<FieldElementT>*>(wrapper_.get());
  ASSERT_RELEASE(ptr != nullptr, "The underlying type of FieldElement is wrong");

  return ptr->ec_;
}

template <typename IntType>
const IntType& EcData::GetOrder() const {
  auto* ptr = dynamic_cast<const IntWrapper<IntType>*>(int_wrapper_.get());
  ASSERT_RELEASE(ptr != nullptr, "The underlying type of Int is wrong");

  return ptr->val_;
}

template <typename FieldElementT>
starkware::EcPoint<FieldElementT> EcData::GetSubGroupGenerator(uint64_t n) const {
  using IntType = decltype(FieldElementT::FieldSize());
  IntType quotient, remainder;
  std::tie(quotient, remainder) =
      IntType::Div(GetOrder<IntType>(), IntType(n));

  ASSERT_RELEASE(remainder == IntType::Zero(), "No subgroup of required size exists");
  return GetCurve<FieldElementT>().MultiplyByScalar(Generator<FieldElementT>(), quotient);
}
