#pragma once

#include <expected>
#include <string>

// not convinced of the utility of this header
// over just aliasing string(_view) as the error type

struct Err {
  std::string msg;  // human-readable, for reporting
};

using VoidOrErr = std::expected<void, Err>;
using IntOrErr = std::expected<int, Err>;
using BoolOrErr = std::expected<bool, Err>;

inline Err make_err (const std::string& msg)
{
  return Err{msg};
};
