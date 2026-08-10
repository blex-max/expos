#pragma once

#include <iostream>
#include <string>

// User-facing warning
inline void expos_warn (const std::string& msg)
{
  std::cerr << "warning: " << msg << '\n';
}
