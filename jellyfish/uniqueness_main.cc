/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <jellyfish/misc.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/compacted_hash.hpp>
#include <jellyfish/uniqueness_cmdline.hpp>
#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/fstream_default.hpp>
#include <algorithm>
#include <string>
#include <sstream>

class FastaItem{
  public:
    std::string id;
    std::string sequence;
    bool empty;
    bool isNotEmpty(){
      return(!empty);
    }
  private:
};

FastaItem nextFastaFromStream(std::istream &in){
  using namespace std;
  FastaItem res;
  res.empty = false;
  while(in.good()){
    char n = in.peek();
    if(n == '>' && res.sequence.size()){
      return(res);
    }else{
      string line;
      getline(in,line);
      
      if(n=='>'){
        istringstream iss(line);
        string name_with_char;
        iss >> name_with_char;
        res.id = name_with_char.substr(1); //skip the leading '>'
      }else{
        //line.erase(std::remove_if(line.begin(), line.end(), isspace ));
        std::transform(line.begin(), line.end(), line.begin(), ::toupper); //uppercase
        if(line.size()){
          res.sequence.append(line);
        }
      }
    }
  }
  //return last one
  if(res.id.size() && res.sequence.size()){
    return(res);
  }
  res.empty = true;
  return(res);
}

template<typename hash_t>
void print_wig_file(const hash_t &h, std::istream &in, std::ostream &out) {
  using namespace std;
  uint_t mer_len = h.get_mer_len();
  uint_t width = mer_len + 2;
  char mer[width + 1];
  mer[width] = '\0';
  FastaItem fa;
  while((fa = nextFastaFromStream(in)).isNotEmpty()){
    
    out << "fixedStep  chrom=" << fa.id << "  start=0  step=1" << endl;
    
    for(int i = 0; i <= int(fa.sequence.size()) - int(mer_len); i++){
      
      //get mer and check for N
      int lastN = -1;
      for(int j = 0; j < int(width); j++){
        char c = fa.sequence[i+j];
        mer[j] = c;
        if(c == 'N'){
          lastN = j;
        }
      }
      
      //handle output for this mer
      if(lastN >= 0){
        for(int j=i; j< (i+lastN) && j < int(fa.sequence.size()); j++){
          /*uniqueness=1=1 occurence,0.5=2,0.33=3,0.25=4,0=5 or more (or containing sequence ambiguities)*/
          out << 0 << endl;
        }
        i+=lastN;
      }else{
        //no N's
        int count = h[mer];
        /*uniqueness=1=1 occurence,0.5=2,0.33=3,0.25=4,0=5 or more (or containing sequence ambiguities)*/
        switch(count){
        case 1:
          out << 1 << endl;
          break;
        case 2:
          out << 0.5 << endl;
          break;
        case 3:
          out << 0.33 << endl;
          break;
        case 4:
          out << 0.25 << endl;
          break;
        default:
          out << 0 << endl;
          break;
        }
      }
      //handle the last few bases, just give them 0 scores
      if(i >= int(fa.sequence.size()) - int(mer_len)){
        for(int j = i; j < int(fa.sequence.size()); j++){
          out << 0 << endl;
        }
      }
    }
  }
}

int uniqueness_main(int argc, char *argv[])
{
  uniqueness_args args(argc, argv);

  ifstream_default in(args.input_given ? args.input_arg : 0, std::cin);
  if(!in.good())
    die << "Can't open input file '" << args.input_arg << "'" << err::no;
  ofstream_default out(args.input_given ? args.output_arg : 0, std::cout);
  if(!out.good())
    die << "Can't open output file '" << args.output_arg << "'" << err::no;

  mapped_file dbf(args.db_arg.c_str());
  char type[8];
  memcpy(type, dbf.base(), sizeof(type));
  if(!strncmp(type, jellyfish::raw_hash::file_type, sizeof(type))) {
    raw_inv_hash_query_t hash(dbf);
    hash.set_canonical(args.both_strands_flag);
    print_wig_file(hash, in, out);
  } else if(!strncmp(type, jellyfish::compacted_hash::file_type, sizeof(type))) {
    hash_query_t hash(dbf);
    hash.set_canonical(args.both_strands_flag);
    print_wig_file(hash, in, out);
  }

  out.close();

  return 0;
}

