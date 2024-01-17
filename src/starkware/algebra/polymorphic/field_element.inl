// Copyright 2023 StarkWare Industries Ltd.
//
// Licensed under the Apache License, Version 2.0 (the "License").
// You may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// https://www.starkware.co/open-source-license/
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions
// and limitations under the License.

#include "starkware/algebra/field_operations.h"
#include "starkware/algebra/field_element_base.h"
#include "starkware/algebra/polymorphic/field.h"
#include "starkware/error_handling/error_handling.h"
#include "nethermind/ec_group.h"

namespace starkware {

class FieldElement::WrapperBase {
 public:
  virtual ~WrapperBase() = default;
  virtual std::unique_ptr<WrapperBase> Clone() const = 0;
  virtual Field GetField() const = 0;
  //virtual FieldElement operator+(const FieldElement& other) const = 0;
  //virtual FieldElement operator-(const FieldElement& other) const = 0;
  //virtual FieldElement operator-() const = 0;
  virtual FieldElement operator*(const FieldElement& other) const = 0;
  virtual FieldElement operator/(const FieldElement& other) const = 0;
  virtual FieldElement Inverse() const = 0;
  virtual FieldElement Pow(uint64_t exp) const = 0;
  virtual void ToBytes(gsl::span<std::byte> span_out, bool use_big_endian) const = 0;
  virtual bool Equals(const FieldElement& other) const = 0;
  virtual size_t SizeInBytes() const = 0;
  virtual std::string ToString() const = 0;
};

template <typename T>
class FieldElement::Wrapper : public WrapperBase {
 public:
  explicit Wrapper(const T& value) : value_(value) {}

  const T& Value() const { return value_; }

  std::unique_ptr<WrapperBase> Clone() const override {
    return std::make_unique<Wrapper<T>>(value_);
  }

  Field GetField() const override { return Field::Create<T>(); }

  /*FieldElement operator+(const FieldElement& other) const override {
    return FieldElement(value_ + other.As<T>());
  }*/

  /*FieldElement operator-(const FieldElement& other) const override {
    return FieldElement(value_ - other.As<T>());
  }*/

  //FieldElement operator-() const override { return FieldElement(-value_); }

  FieldElement operator*(const FieldElement& other) const override {
    return FieldElement(value_ * other.As<T>());
  }

  FieldElement operator/(const FieldElement& other) const override {
    return FieldElement(value_ / other.As<T>());
  }

  FieldElement Inverse() const override { return FieldElement(value_.Inverse()); }

  FieldElement Pow(uint64_t exp) const override {
    return FieldElement(::starkware::Pow(value_, exp));
  }

  void ToBytes(gsl::span<std::byte> span_out, bool use_big_endian) const override {
    value_.ToBytes(span_out, use_big_endian);
  }

  bool Equals(const FieldElement& other) const override { return value_ == other.As<T>(); }

  size_t SizeInBytes() const override { return value_.SizeInBytes(); }

  std::string ToString() const override { return value_.ToString(); }

 private:
  const T value_;
};

template <typename T>
class FieldElement::EcWrapper : public WrapperBase {
 public:
  explicit EcWrapper(const EcPoint<T>& value, const EC<T>& ec)
           : value_(value), 
             ec_(ec)
  {
    ASSERT_DEBUG(ec.ContainsPoint(value),"Point not on curve");
  }

  const EcPoint<T>& Value() const { return value_; }

  const EC<T>& Curve() const {return ec_; }

  std::unique_ptr<WrapperBase> Clone() const override {
    return std::make_unique<EcWrapper<T>>(value_,ec_);
  }

  Field GetField() const override { return Field::Create<T>(); }

  /*FieldElement operator+(const FieldElement& other) const override {
    return FieldElement(value_ + other.As<T>());
  }*/

  /*FieldElement operator-(const FieldElement& other) const override {
    return FieldElement(value_ - other.As<T>());
  }*/

  //FieldElement operator-() const override { return FieldElement(-value_); }

  FieldElement operator*(const FieldElement& other) const override {
    ASSERT_DEBUG(ec_  == other.GetCurve<T>(),"Points on different curves");
    //return FieldElement(value_ + other.AsEc<T>(), ec_);
    return FieldElement(ec_.addPoints(value_ ,other.AsEc<T>()),ec_);
  }

  FieldElement operator/(const FieldElement& other) const override {
    ASSERT_DEBUG(ec_  == other.GetCurve<T>(),"Points on different curves");
    //return FieldElement(value_ - other.AsEc<T>(),ec_);
    return FieldElement(ec_.addPoints(value_ ,-other.AsEc<T>()),ec_);
  }

  FieldElement Inverse() const override { return FieldElement(-value_,ec_); }

  FieldElement Pow(uint64_t exp) const override {
    return FieldElement(value_.MultiplyByScalar(BigInt<1>(exp), ec_.alpha), ec_);
  }

  void ToBytes(gsl::span<std::byte> span_out, bool use_big_endian) const override {
    value_.x.ToBytes(span_out, use_big_endian);
    value_.y.ToBytes(span_out, use_big_endian);
  }

  bool Equals(const FieldElement& other) const override {
    ASSERT_DEBUG(ec_  == other.GetCurve<T>(),"Points on different curves");
    return value_ == other.AsEc<T>();
  }

  size_t SizeInBytes() const override { return 2*value_.x.SizeInBytes(); }

  std::string ToString() const override
  { return "(" + value_.x.ToString() + ", " + value_.y.ToString() + ")";}

 private:
  const EcPoint<T> value_;
  const EC<T> ec_;
};

template <typename T>
FieldElement::FieldElement(const T& t) : wrapper_(std::make_unique<Wrapper<T>>(t)) {}
/*FieldElement::FieldElement(const T& t) {
  if constexpr(std::is_base_of<FieldElementBase<T>, T>::value) 
    wrapper_=std::make_unique<Wrapper<T>>(t);
  else
    wrapper_=std::make_unique<EcWrapper<T>>(t);
  }*/

template <typename T>
FieldElement::FieldElement(const EcPoint<T>& t,const EC<T>& ec) : wrapper_(std::make_unique<EcWrapper<T>>(t,ec)) {}

template <typename T>
const T& FieldElement::As() const {
  auto* ptr = dynamic_cast<const Wrapper<T>*>(wrapper_.get());
  ASSERT_RELEASE(ptr != nullptr, "The underlying type of FieldElement is wrong");

  return ptr->Value();
}

template <typename T>
const EcPoint<T>& FieldElement::AsEc() const {  
    auto* ptr = dynamic_cast<const EcWrapper<T>*>(wrapper_.get());
    ASSERT_RELEASE(ptr != nullptr, "The underlying type of FieldElement is wrong");

    return ptr->Value();
  
}

template <typename T>
const EC<T>& FieldElement::GetCurve() const {  
    auto* ptr = dynamic_cast<const EcWrapper<T>*>(wrapper_.get());
    ASSERT_RELEASE(ptr != nullptr, "The underlying type of FieldElement is wrong");

    return ptr->Curve();
  
}

}  // namespace starkware
