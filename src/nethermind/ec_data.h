#ifndef NETHERMIND_EC_DATA_H_
#define NETHERMIND_EC_DATA_H_

//#include <cstddef>
//#include <memory>
//#include <ostream>
//#include <string>
//#include <vector>

#include "nethermind/ec_group.h"
#include "starkware/algebra/polymorphic/field.h"
//#include "third_party/gsl/gsl-lite.hpp"

class EcData {
 public:
  template <typename FieldElementT, typename IntType = decltype(FieldElementT::FieldSize())>
  explicit EcData(const EC<FieldElementT>& ec,
                  const starkware::EcPoint<FieldElementT>& gen,
                  const IntType& order);

/*  EcData(const EcData& other);

  EcData& operator=(const EcData& other) &;

  ~EcData() = default;
  EcData(EcData&& other) = default;
  EcData& operator=(EcData&& other) & = default;*/

  starkware::Field GetField() const;

  template <typename FieldElementT>
  const starkware::EcPoint<FieldElementT>& Generator() const;

  template <typename IntType>
  const IntType& GetOrder() const;

  template <typename FieldElementT>
  const EC<FieldElementT>& GetCurve() const;

  template <typename FieldElementT>
  starkware::EcPoint<FieldElementT> GetSubGroupGenerator(uint64_t n) const;

 private:
  EcData(const EcData&);  // Declare the copy constructor as private without providing its implementation

  class WrapperBase;

  template <typename FieldElementT>
  class Wrapper;

  class IntWrapperBase;

  template <typename IntType>
  class IntWrapper;

  std::unique_ptr<WrapperBase> wrapper_;
  std::unique_ptr<IntWrapperBase> int_wrapper_;
};

#include "nethermind/ec_data.inl"

#endif  // NETHERMIND_EC_DATA_H_
