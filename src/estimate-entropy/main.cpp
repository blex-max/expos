#include <iostream>
#include <cstdlib>
#include <string>

#include "lib-stats/string.hpp"

constexpr std::size_t MAX_STDIN_BYTES = 1000;  // max input chars per line

static void check_ascii(const std::string& s) {
  for (const auto c : s) {
    if (c <= 0x20 || c > 0x7E) {
      std::cerr << "Input must be an unbroken ASCII string" << std::endl;
      std::exit(1);
    }
  }
}

static void print_help() {
  std::cout
    << "estimate-entropy\n\n"
    << "Lempel-Ziv (1976) estimated entropy rate per character of input string.\n"
    << "Reads one string per line from stdin and prints one value per line.\n\n"
    << "Usage:\n"
    << "  estimate-entropy < windows.txt\n\n"
    << "Options:\n"
    << "  -h, --help    Print usage\n";
}

// TODO report % string in periodic repeats
// TODO report MFE of string (if DNA)
int main (int argc, char ** argv) {
  bool print_region=false;

  if (argc > 2) {
    std::cerr << "Too many arguments" << std::endl;
    return 1;
  }
  if (argc == 2) {
    std::string arg = argv[1];
    if (arg == "-h" || arg == "--help") {
      print_help();
      return 0;
    } else if (arg == "-l") {
      print_region=true;
    } else {
      std::cerr << "Unrecognized argument: " << arg << std::endl;
      return 1;
    }
  }

  bool input_found = false;
  std::string line;
  while (std::getline(std::cin, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;

    if (line.size() > MAX_STDIN_BYTES) {
      std::cerr << "Input > 1000 chars" << std::endl;
      return 1;
    }
    check_ascii(line);

    if (print_region) {
      std::cout << line;
      std::cout << "\t";
    }
    std::cout << std::to_string(string_stats::entropy_lz76(line)) << std::endl;
    input_found = true;
  }

  if (!input_found) {
    std::cerr << "No input provided" << std::endl;
    return 1;
  }

  return 0;
}

