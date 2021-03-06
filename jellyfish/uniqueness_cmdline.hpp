
#ifndef __UNIQUENESS_ARGS_HPP__
#define __UNIQUENESS_ARGS_HPP__

#include <jellyfish/yaggo.hpp>

class uniqueness_args {
public:
  bool                           both_strands_flag;
  const char *                   input_arg;
  bool                           input_given;
  const char *                   output_arg;
  bool                           output_given;
  yaggo::string                  db_arg;

  enum {
    USAGE_OPT = 1000
  };

  uniqueness_args(int argc, char *argv[]) :
    both_strands_flag(false),
    input_arg(""), input_given(false),
    output_arg(""), output_given(false)
  {
    static struct option long_options[] = {
      {"both-strands", 0, 0, 'C'},
      {"input", 1, 0, 'i'},
      {"output", 1, 0, 'o'},
      {"help", 0, 0, 'h'},
      {"usage", 0, 0, USAGE_OPT},
      {"version", 0, 0, 'V'},
      {0, 0, 0, 0}
    };
    static const char *short_options = "hVCi:o:";

#define CHECK_ERR(type,val,which) if(!err.empty()) { std::cerr << "Invalid " #type " '" << val << "' for [" which "]: " << err << "\n"; exit(1); }
    while(true) { 
      int index = -1;
      int c = getopt_long(argc, argv, short_options, long_options, &index);
      if(c == -1) break;
      switch(c) {
      case ':': 
        std::cerr << "Missing required argument for "
                  << (index == -1 ? std::string(1, (char)optopt) : std::string(long_options[index].name))
                  << std::endl;
        exit(1);
      case 'h':
        std::cout << usage() << "\n\n" << help() << std::endl;
        exit(0);
      case USAGE_OPT:
        std::cout << usage() << "\nUse --help for more information." << std::endl;
        exit(0);
      case 'V':
        print_version();
        exit(0);
      case '?':
        std::cerr << "Use --usage or --help for some help\n";
        exit(1);
      case 'C':
        both_strands_flag = true;
        break;
      case 'i':
        input_given = true;
        input_arg = optarg;
        break;
      case 'o':
        output_given = true;
        output_arg = optarg;
        break;
      }
    }
    if(argc - optind != 1)
      error("Requires exactly 1 argument.");
    db_arg = yaggo::string(argv[optind]);
    ++optind;
  }
#define uniqueness_args_USAGE "Usage: jellyfish uniqueness [options] db:path"
  const char * usage() const { return uniqueness_args_USAGE; }
  void error(const char *msg) { 
    std::cerr << "Error: " << msg << "\n" << usage()
              << "\nUse --help for more information"
              << std::endl;
    exit(1);
  }
#define uniqueness_args_HELP "Query from a compacted database\n\nQuery a hash. It reads k-mers from the standard input and write the counts on the standard output.\n\n" \
  "Options (default value in (), *required):\n" \
  " -C, --both-strands                       Both strands (false)\n" \
  " -i, --input=file                         Input file\n" \
  " -o, --output=file                        Output file\n" \
  "     --usage                              Usage\n" \
  " -h, --help                               This message\n" \
  " -V, --version                            Version"

  const char * help() const { return uniqueness_args_HELP; }
#define uniqueness_args_HIDDEN "Hidden options:"

  const char * hidden() const { return uniqueness_args_HIDDEN; }
  void print_version(std::ostream &os = std::cout) const {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0"
#endif
    os << PACKAGE_VERSION << "\n";
  }
  void dump(std::ostream &os = std::cout) {
    os << "both_strands_flag:" << both_strands_flag << "\n";
    os << "input_given:" << input_given << " input_arg:" << input_arg << "\n";
    os << "output_given:" << output_given << " output_arg:" << output_arg << "\n";
    os << "db_arg:" << db_arg << "\n";
  }
private:
};

#endif // __UNIQUENESS_ARGS_HPP__"
