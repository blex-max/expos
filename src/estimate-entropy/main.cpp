#include <iostream>

#include <cxxopts.hpp>

#include "lib-stats/stats.hpp"

int main (int argc, char ** argv) {
  std::string input;

  cxxopts::Options args (
    "estimate-entropy",
    "\n"
    "Lempel-Ziv (1976) estimated entropy rate per character of input string <STR>"
    "\n"
  );

  args.add_options()
    ("h,help", "Print usage")
    ("str", "input string", cxxopts::value<std::string>());

  args.parse_positional({"str"});
  args.positional_help("<STR>");

  {
    auto parsedargs = args.parse(argc, argv);

    if (parsedargs.count ("help")) {
        std::cout << args.help() << std::endl;
        return 0;
    }

    if (!parsedargs.count ("str")) {
        std::cout << "All positional arguments must be provided" << std::endl;
        return 1;
    } else {
      input = parsedargs["str"].as<std::string>();
    }
  }

  std::cout << std::to_string(entropy_lz76(input)) << std::endl;

  return 0;
}
