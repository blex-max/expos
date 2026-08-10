#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

#include "argparse/argparse.hpp"
#include "shared/stats.hpp"

constexpr std::size_t MAX_STDIN_BYTES =
    1000;  // max input chars per line

static bool is_ascii_printable (const std::string& s)
{
  constexpr char MIN_PRINTABLE = 0x20;
  constexpr char MAX_PRINTABLE = 0x7E;
  return std::ranges::all_of (s, [] (const char c) {
    return c > MIN_PRINTABLE && c <= MAX_PRINTABLE;
  });
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

  try {
    cli.parse_args (argc, argv);
  }
  catch (const std::exception& ex) {
    std::cerr << ex.what() << "\n" << cli;
    return EXIT_FAILURE;
  }

  const bool showInput = cli.get<bool> ("--show-input");

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
    std::cout << std::to_string (entropy_lz76 (line))
              << std::endl;
    inputFound = true;
  }

  if (!inputFound) {
    std::cerr << "No input provided" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
