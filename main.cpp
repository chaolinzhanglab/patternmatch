

// C++ standard headers
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <functional>
#include <cstdlib>
// C headers
#include <popt.h>        // For command line options

// Unix libraries
#include <unistd.h>
#include <sys/stat.h>

/* THIS PROGRAM REQUIRES THE "GNU Scientific Ligrary"  */
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>


#include "FastaFile.hpp"

using namespace std;

#define log2(A) (log(A)/log(2.0))

static const int alphabet_size = 4;

/* INPUT PARAMETERS */
const char *sequences_file_name = static_cast <const char *>(0); // file containing foreground sequences
const char *motif_consensus = static_cast <const char *>(0);
const char *motif_list_file_name = static_cast <const char *>(0); // specify motif consensus or a file with motif concensus
int mismatch_threshold = 0;

const char *output_file_name = static_cast<const char *>(0);        // file in which to print output

int use_both_strands = 0;
int ignore_repeat = 0;
int verbose = 0;

char complement(char c) {
  static char btoc[26] = {
  //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T,  u,  v,  w,  x,  y,  z
   'T','V','G','H','N','N','C','D','N','N','M','N','K','N','N','N','N','Y','S','A','N','B','W','N','R','N' 
  };
  static char btoc_small[26] = {
  //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T,  u,  v,  w,  x,  y,  z
   't','v','g','h','n','n','c','d','n','n','m','n','k','n','n','n','n','y','s','a','n','b','w','n','r','n' 
  };
  if (c <= 'Z' && c >= 'A')
  	return btoc[c - 'A'];
  else
	return btoc_small[c - 'a'];
}
string reverse_complement(string &s) {
  string r;
  transform(s.begin(), s.end(), back_inserter(r), complement);
  reverse(r.begin(), r.end());
  return r;
}

int get_sequence_length(vector<string> &s) {
  int count = 0;
  for (vector<string>::iterator i = s.begin(); i != s.end(); ++i)
    for (string::iterator j = i->begin(); j != i->end(); ++j)
      count += (*j != 'N' && *j != 'n');
  return count;
}

/****************** COMMAND LINE OPTIONS ********************/
poptContext optCon = NULL;
static struct poptOption optionsTable[] = {
  { "motif-consensus", 'c', POPT_ARG_STRING, &motif_consensus, 0, "the consensus sequence of motif" },
  { "number-of-mismatches", 'm', POPT_ARG_INT, &mismatch_threshold, 0, "number of mismatches allowed" },
  { "motif-list-file", 'l', POPT_ARG_STRING, &motif_list_file_name, 0, "motif list file name" },
  { "both-strands", 'b', POPT_ARG_NONE, &use_both_strands, 0, "Search both strands" },
  { "ignore-repeat", 'i', POPT_ARG_NONE, &ignore_repeat, 0, "Ignore repeat (in small letters)"},
  { "verbose", 'v', POPT_ARG_NONE, &verbose, 0, "Verbose mode" },
  { "output-file", 'o', POPT_ARG_STRING, &output_file_name, 0, "output file (default: stdout)" },
   POPT_AUTOHELP { NULL, 0, 0, NULL, 0, }
};

bool  BaseMatch (char q, char c)
{
	//see if the definition of q is included in c
	//the parameter should in in upper case
	
	if (ignore_repeat && q >= 'a' && q<= 'z') //repeat region should be ignored when the switch is on
		return false;
	
	if (c == q || c == 'N')
		return true;
	if ((c == 'R' && (q == 'A' || q == 'G'))
	  ||(c == 'Y' && (q == 'C' || q == 'T')) 
	  ||(c == 'M' && (q == 'A' || q == 'C'))
	  ||(c == 'K' && (q == 'G' || q == 'T'))
	  ||(c == 'S' && (q == 'C' || q == 'G'))
	  ||(c == 'W' && (q == 'A' || q == 'T'))
	  ||(c == 'H' && (q == 'A' || q == 'C' || q == 'T'))
	  ||(c == 'B' && (q == 'C' || q == 'G' || q == 'T'))
	  ||(c == 'V' && (q == 'A' || q == 'C' || q == 'G'))
	  ||(c == 'D' && (q == 'A' || q == 'G' || q == 'T'))
	  )
		return true;
	return false;
}

size_t CountMisMatch (string target, string motif)
{
	size_t w = motif.size ();
	//transform (target.begin (), target.end (), target.begin(),::toupper);
	//transform (motif.begin (), motif.end (), motif.begin(), ::toupper);
	if (target.size() != w)
	{
		cerr << "the length of the motif is wrong" << endl;
		exit (EXIT_FAILURE);
	}

	size_t m = 0;
	string::iterator i;
	string::iterator j;
	for (i = target.begin (), j = motif.begin (); i != target.end () && j != motif.end (); ++i, ++j)
	{
		//if (*i != *j)
		if (!BaseMatch (*i, *j)) // to allow IUB code
			++m;
	}
	return m;
}

size_t CountMisMatchFast (string target, string motif, size_t mismatch_threshold)
{
	size_t w = motif.size ();
	//transform (target.begin (), target.end (), target.begin(),::toupper);
	//transform (motif.begin (), motif.end (), motif.begin(), ::toupper);
	if (target.size() != w)
	{
		cerr << "the length of the motif is wrong:" << target << endl;
		exit (EXIT_FAILURE);
	}

	size_t m = 0;
	string::iterator i;
	string::iterator j;
	for (i = target.begin (), j = motif.begin (); i != target.end () && j != motif.end (); ++i, ++j)
	{
		//if (*i != *j)
		if (!BaseMatch(*i, *j))//to allow IUB code
			++m;
		if (m > mismatch_threshold)
				return w;
	}
	return m;
}



//////////////////////////////////////////
//////////// THE MAIN METHOD /////////////
//////////////////////////////////////////
int main(int argc, char **argv) {
  char c;
  
  /***************** GET COMMAND LINE ARGUMENTS *******************/
  optCon = poptGetContext("storm", argc, (const char **)argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "fasta-file");
  if (argc < 2) {
    poptPrintUsage(optCon, stderr, 0);
    exit(EXIT_SUCCESS);
  }
  if ((c = poptGetNextOpt(optCon)) >= 0) {
    fprintf(stderr, "Error: bad option ?");
    return EXIT_FAILURE;
  }
  if ((c = poptGetNextOpt(optCon)) < -1) {
    fprintf(stderr, "PatternMatch: bad argument %s: %s\n",
            poptBadOption(optCon, POPT_BADOPTION_NOALIAS),
            poptStrerror(c));
    return EXIT_FAILURE;
  }
  if (poptPeekArg(optCon) == NULL) {
    poptPrintUsage(optCon, stderr, 0);
    return EXIT_FAILURE;
  }
  else sequences_file_name = poptGetArg(optCon);
  poptFreeContext(optCon);
  /**********************************************************************/

  // READ IN THE SEQUENCES AND THEIR NAMES
  FastaFile faa = FastaFile(const_cast<char *>(sequences_file_name));
  
  //faa.ToUpper(); 
  faa.Clean();
  vector<string> sequence_names = faa.GetNames();
  vector<string> sequences = faa.GetSequences();
  size_t num_seq = sequence_names.size();

  //read motif list if the file is specified
  vector<pair<string,int> > motif_list;
  
  if (motif_list_file_name != 0)
  {
	if (verbose)
		cerr << "reading motifs from " << motif_list_file_name << "..." << endl;
	ifstream in (motif_list_file_name);
	if (!in)
	{
		cerr << "cannot open motif list file " << motif_list_file_name << endl;
		exit (1);
	}

	string line;
	while (!in.eof())
	{
		getline (in, line);
		if (line.empty())
			continue;

		//still empty lines
		int pos = line.find_first_not_of("\t ");
		if (pos < 0)
			continue;

		//read motif concensus and mismatches
		pos = line.find_first_of("\t ");
		if (pos < 0)
			cerr << "each line of the motif list file must be concensus<space or tab>mismatch" << endl;

		pair<string,int> motif_info;
		motif_info.first = line.substr(0, pos);
		motif_info.second = atoi (line.substr(pos+1).c_str());
		motif_list.push_back (motif_info);

		//cerr << "motif " << motif_iter++ << ":" << motif.first << "\tmismatch=" << motif.second << endl;
	}
  }
  else
  {
	//get motif from the command line
	if (motif_consensus == 0)
	{
		cerr << "motif must be provided either in the command line or in a file" << endl;
		exit (1);
	}
	pair<string,int> motif_info;
	motif_info.first = motif_consensus;
	motif_info.second = mismatch_threshold;
	motif_list.push_back (motif_info);
  }

  //now ready to search motif
  //ostringstream buffer;
  ofstream out;
  if (output_file_name != static_cast<char *> (0))
  {
	//ofstream out;
	out.open (const_cast<char *>(output_file_name));
	//out << buffer.str ();
	//out.close();	
  }
 

  for (size_t m = 0; m < motif_list.size(); m++)
  {
	pair<string,int> motif_info = motif_list[m];
	string motif = motif_info.first;
	int mismatch_threshold = motif_info.second;

    transform (motif.begin (), motif.end (), motif.begin(), ::toupper);	
    //this is necessary because the subroutine reverse_complement can not handle small letters
    size_t motif_width = motif.size ();
  
  	string motif_rc = reverse_complement (motif);
 	if (verbose)
		cerr << "\n\n##search for motif # " << m << ": " << motif << " ..." << endl;
  	for (size_t i = 0; i < num_seq; ++i)
	{
		string name = sequence_names[i];
		string seq = sequences[i];
	
		if (verbose & (i % 500 == 0))
			cerr << "processing sequence # " << i << ":" << name << ", size=" << seq.size () << endl; 
	
		size_t len = seq.size ();
		if (len < motif_width)
			continue;


		for (size_t j = 0; j < len - motif_width + 1; ++j)
		{
			string target = seq.substr (j, motif_width);
			string target2 = target;
			if (!ignore_repeat)
				transform (target2.begin (), target2.end (), target2.begin(), ::toupper);
		
			size_t mismatch = CountMisMatchFast (target2, motif, mismatch_threshold);
			if (mismatch <= (size_t)mismatch_threshold)
			{
				//transform (target.begin (), target.end (), target.begin(), ::toupper);
				if (output_file_name != static_cast<char *> (0))
				{
					out << name << "\t" 
						<< j << "\t"
						<< (j + motif_width) << "\t"
				   		<< motif << ":" << target << "\t"
						<< mismatch << "\t"
						<< "+" << endl;
				}
				else
				{
					cout << name << "\t" 
						<< j << "\t"
						<< (j + motif_width) << "\t"
					   	<< motif << ":" << target << "\t"
						<< mismatch << "\t"
						<< "+" << endl;
				}
			}

			if (use_both_strands == 0)
				continue;
	
			//the negative strand		
			mismatch = CountMisMatchFast (target2, motif_rc, mismatch_threshold);
		
			if (mismatch <= (size_t)mismatch_threshold)
			{
				//transform (target.begin (), target.end (), target.begin(),::toupper);
			
				if (output_file_name != static_cast<char *> (0))
				{
					out << name << "\t" 
						<< j << "\t"
						<< (j + motif_width) << "\t"
				   		<< motif << ":" << reverse_complement (target) << "\t"
						<< mismatch << "\t"
						<< "-" << endl;
				}
				else
				{
					cout << name << "\t" 
						<< j << "\t"
						<< (j + motif_width) << "\t"
					   	<< motif << ":" << reverse_complement (target) << "\t"
						<< mismatch << "\t"
						<< "-" << endl;
				}
			}
		}
	}
  }

  if (output_file_name != static_cast<char *> (0))
  {
	out.close();	
  }
  //else
  //{
	//cout << buffer.str ();
  //}
  return EXIT_SUCCESS;
}

