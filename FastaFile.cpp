#include "FastaFile.hpp"

static char btoc[26] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T,  u,  v,  w,  x,  y,  z
			'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A','N','N','N','N','N','N'};

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
  char buffer[buffer_size];

  ifstream in(filename.c_str());
  if (!in) {
    cerr << "cannot open input file " << filename << endl;
    exit(1);
  }

  string s, name = "";
  bool first_line = true;
  while (!in.eof()) {
    in.getline(buffer, buffer_size);
    if (buffer[0] == '>') {
      if (first_line == false && s.length() > 0) {
	names.push_back(name);
	sequences.push_back(s);
      }
      else first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("> "));
      s = "";
    }
    else s += buffer;
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
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
	for(size_t i = 0; i < r.size(); ++i) r[i]= btoc[r[i] - 'A'];
	reverse(r.begin(), r.end());
	return r;
}
