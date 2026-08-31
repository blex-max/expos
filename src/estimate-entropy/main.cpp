#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

#include "argparse/argparse.hpp"
#include "shared/stats.hpp"

constexpr uint32_t MAX_STDIN_BYTES =
    1000;  // max input chars per line

static bool is_ascii_printable (const std::string& s)
{
  constexpr char MIN_PRINTABLE = 0x20;
  constexpr char MAX_PRINTABLE = 0x7E;
  return std::ranges::all_of (s, [] (const char c) {
    return c > MIN_PRINTABLE && c <= MAX_PRINTABLE;
  });
}

static double entropy_per_char (uint16_t lzCmplx, uint16_t nChar)
{
  return static_cast<double> (lzCmplx * log2 (nChar)) /
         static_cast<double> (nChar);
}


int main (int argc, char** argv)
{
  argparse::ArgumentParser cli ("estimate-entropy");
  cli.add_description (
      "Lempel-Ziv (1976) estimated entropy rate per character "
      "of "
      "input string. Reads one string per line from stdin and "
      "prints one value per line."
  );
  cli.add_argument ("-l", "--show-input")
      .flag()
      .help ("prefix each output line with the input string");
  cli.add_argument ("-n", "--normalise")
      .flag()
      .help (
          "normalise to estimated entropy per character by\n"
          "((complexity * log2(N)) / N )"
      );

  try {
    cli.parse_args (argc, argv);
  }
  catch (const std::exception& ex) {
    std::cerr << ex.what() << "\n" << cli;
    return EXIT_FAILURE;
  }

  const bool showInput = cli.get<bool> ("--show-input");
  const bool normalise = cli.get<bool> ("--normalise");

  bool inputFound = false;
  std::string line;
  while (std::getline (std::cin, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty()) {
      continue;
    }

    if (line.size() > MAX_STDIN_BYTES) {
      std::cerr << "Input > 1000 chars" << std::endl;
      return EXIT_FAILURE;
    }
    if (!is_ascii_printable (line)) {
      std::cerr << "Input must be an unbroken ASCII string"
                << std::endl;
      return EXIT_FAILURE;
    }

    if (showInput) {
      std::cout << line << "\t";
    }
    const auto result = lz76 (line);
    std::cout << std::to_string (
                     (normalise)
                         ? entropy_per_char (result, line.size())
                         : result
                 )
              << std::endl;
    inputFound = true;
  }

  if (!inputFound) {
    std::cerr << "No input provided" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
