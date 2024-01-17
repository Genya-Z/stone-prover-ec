#include "nethermind/ec_data.h"

starkware::Field EcData::GetField() const { return wrapper_->GetField(); }