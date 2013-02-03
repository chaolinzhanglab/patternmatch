#ifndef FASTA_FILE_HPP
#define FASTA_FILE_HPP

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <functional>
#include <vector>
#include <sstream>
#include <set>
#include <map>
#include <ctime>
#include <cstdlib>

#include "gzstream.h"
//This is adapted from the implementation 
//of Andrew Smith (http://http://smithlab.usc.edu/plone)


using namespace std;

class FastaFile {

  static const int buffer_size = 1000000;

  string filename;
  vector<string> names;	//sequence name
  vector<string> descs;	//sequence descriptions
  vector<string> sequences;	//sequences

  void ReadFile();

  static char btoc(int i) {
    static char btoc[26] = {
     //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T,  u,  v,  w,  x,  y,  z
      'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A','N','N','N','N','N','N'
    };
    return btoc[i];
  }
  
public:
  
  FastaFile() : filename("") {}
  FastaFile(const FastaFile &);
  FastaFile &operator=(const FastaFile &);
  //~FastaFile();

  FastaFile(string &);
  FastaFile(char *);

  // Accessors
  vector<string> GetNames() const {return names;}
  vector<string> GetDescs () const {return descs;}
  vector<string> GetSequences() const {return sequences;}
	size_t size(){return sequences.size();}

  // Mutators
  void Clean();
  void ToUpper();

  //static methods

  /*
   * GetBaseComposition expects base_composition_vector to be allocated alphabet_size floats.
   * It returns the base composition of sequences in the base_composition_vector array.
   */
  static void GetBaseComposition(vector<string> sequences, float* base_composition_vector, int alphabet_size);
  static string ReverseComplement(string&);
};

FastaFile::FastaFile(const FastaFile &ff) {
  filename = ff.filename;
  names = ff.names;
  sequences = ff.sequences;
}

FastaFile &
FastaFile::operator=(const FastaFile &ff) {
  if (this != &ff) {
    filename = ff.filename;
    names = ff.names;
    sequences = ff.sequences;
  }
  return *this;
}

void FastaFile::ReadFile() {
  //char buffer[buffer_size];

  //ifstream in(filename.c_str());
  igzstream in(filename.c_str());

  if (!in) {
    cerr << "cannot open input file " << filename << endl;
    exit(1);
  }

  string line, s, name = "", desc = "";
  bool first_line = true;
  while (!in.eof()) {
	getline (in, line);

	//in.getline(buffer, buffer_size);
	const char* line_str = line.c_str();
	if (line_str[0] == '>') {
      if (first_line == false && s.length() > 0) {
	names.push_back(name);
	descs.push_back(desc);
	sequences.push_back(s);
      }
      else first_line = false;
		
	  //string str = buffer;
	  string str = line.substr(line.find_first_not_of(">\t "));
	  //cout << str << endl;
	  int pos = str.find_first_of(">\t ");	  
	  if (pos < 0){// no description
		name = str;
	  }
	  else{
	  	name = str.substr(0, pos);
	  	desc = str.substr(pos+1);
	  }
      s = "";
    }
    else s += line;
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
	descs.push_back(desc);
    sequences.push_back(s);
  }
}

FastaFile::FastaFile(string &fn) {
  filename = fn;
  ReadFile();
}

FastaFile::FastaFile(char *fn) {
  filename = string(fn);
  ReadFile();
}

// Not a member of FastaFile !!!
bool Valid(char c) {
  c = ::toupper(c);
  return c=='A' || c=='C' || c=='G'|| c=='T';
}

void
FastaFile::Clean() {
  for (vector<string>::iterator i = sequences.begin(); i != sequences.end(); ++i)
    replace_if(i->begin(), i->end(), not1(ptr_fun(Valid)), 'N');
}

void
FastaFile::ToUpper() {
  for (vector<string>::iterator i = sequences.begin(); i != sequences.end(); ++i)
    transform(i->begin(), i->end(), i->begin(), ::toupper); // not in "namespace std"
}

void
FastaFile::GetBaseComposition(vector<string> sequences, float* f, int alphabet_size) {
  static int btoi[20] = {
  //A, b, C, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, T
    0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3 
  };
  int total = 0;
  fill(f, f + alphabet_size, 0.0);
  for (vector<string>::const_iterator i = sequences.begin(); i != sequences.end(); ++i)
    for (string::const_iterator j = i->begin(); j != i->end(); ++j)
      if (*j != 'N') {f[btoi[*j - 'A']]++; total++;}
  transform(f, f + alphabet_size, f, bind2nd(divides<float>(),total));
}

string
FastaFile::ReverseComplement(string& s)
{
	string r = s;
	for(size_t i = 0; i < r.size(); ++i) r[i]= btoc(r[i] - 'A');
	reverse(r.begin(), r.end());
	return r;
}

#endif
