//g++ -c cgFbat.cpp -pedantic    --> Gives a lot more of the warnings

/* Thomas Hoffmann
 * Fork of condGeneFBATControl.cpp on 07/16/2008.
 *
 * Now we handle vectors of genotypes as well, I don't
 *  know what the hell we were thinking of for doing the
 *  pairwise stuff -- I suppose it was back when
 *  we were thinking about haplotypes...
 *
 * This coding implements the conditional gene tests
 *  that we have been working on.
 * Here we control FBAT to get the haplotype density,
 *  rather than drudging through the FBAT source code.
 * We also keep the data in C++ code, we never export
 *  it back to R. Instead, we write the dataset out from
 *  R, then we run FBAT on the dataset, we read the
 *  output of FBAT into the c++ code, and pass a
 *  reference back to the R code to index the object.
 *  Thus all computations are done in this code.
 */

//#define _cgFbat_DEBUG_
//  g++ cgFbat.cpp -o cgFbat
//  ./cgFbat

#include <R.h>

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
using namespace std;

enum XCODE {
	ADDITIVE, DOMINANT, RECESSIVE, GENOTYPE
};
const int ANALYZE_ALLELE_INDEX = 0; // First marker is the one being analyzed
const int CONDITIONAL_ALLELE_INDEX = 1; // Second Marker is the one that we are conditioning on

/** Misc helper functions **/

// Converts a double to a string
// DEBUGGED
string doubleToString(double d) {
	//string str;
	ostringstream oss;
	oss << d;
	return (oss.str());
}

// Returns in 'ret' all the possible permutations of the 'data'
// DEBUGGED Previously?
void allPerms(vector<int> &data, vector<vector<int> > &ret) {
	if (data.size() == 1) {
		ret.push_back(data);
		return;
	}

	for (unsigned int i = 0; i < data.size(); i++) {
		vector<int> subdata = data;
		subdata.erase(subdata.begin() + i);

		vector<vector<int> > subret;
		allPerms(subdata, subret);
		for (unsigned int s = 0; s < subret.size(); s++) {
			subret[s].push_back(data[i]);
			ret.push_back(subret[s]);
		}
	}
}

/** String tokenizing class **/
// DEBUG: heavily tested
class StrTok {
private:
	vector<string> tokens;
	unsigned int curToken;
public:
	// tokenizes a string based on white spaces
	void tokenize(string &s) {
		// empty the tokens
		tokens.clear();

		// parse out the line
		stringstream ss(s);
		string buf;
		while (ss >> buf)
			tokens.push_back(buf);

		// clear out the current token
		curToken = 0;
	}

	// tokenizes a string based on the delimiters in string 'd'
	void tokenize(string &s, string &d) {
		// empty the tokens
		tokens.clear();

		// parse out the line
		string::size_type fpos = s.find_first_not_of(d, 0); // ignore delimeters at beginning
		string::size_type bpos = s.find_first_of(d, fpos);
		while (fpos != string::npos || bpos != string::npos) {
			tokens.push_back(s.substr(fpos, bpos - fpos));
			fpos = s.find_first_not_of(d, bpos);
			bpos = s.find_first_of(d, fpos);
		}

		// clear out the current token
		curToken = 0;
	}

	// whether there are any more tokens
	bool hasMoreTokens() {
		if (curToken >= tokens.size())
			return (false);
		return (true);
	}

	// get the next token
	string nextToken() {
		// increment the token
		curToken++;

		// if too far past, return an empty string
		if (curToken - 1 >= tokens.size()) {
			Rprintf("StrTok::nextToken() past the end of tokens.\n");
			string empty;
			return (empty);
		}

		// return the current token
		return (tokens[curToken - 1]);
	}

	// get the next token as a _number_
	double nextTokenN() {
		istringstream i(nextToken());
		double x;
		if (!(i >> x))
			return (0.0);
		return (x);
	}

	// expunge the next token
	void nextTokenE() {
		curToken++;
	}

	// get the last token (and update the counter)
	string lastToken() {
		if (curToken >= tokens.size() - 1)
      Rprintf("StrTok::lastToken warning -- already popped off last token.\n");
		curToken = tokens.size() - 1;
		return (nextToken());
	}
	// last token, but as a number
	double lastTokenN() {
		curToken = tokens.size() - 1;
		return (nextTokenN());
	}

	// returns the number of tokens
	unsigned int size() {
		return (tokens.size());
	}

	// converts object to a string, primarily for debugging
	string toString() {
		string string;
		for (unsigned int i = 0; i < tokens.size(); i++)
			string = string + "<" + tokens[i] + ">";
		return (string);
	}
};

/** Reading in lines from output and searching them **/
// DEBUG: Fairly simple class, tested.
class Lines {
public:
	vector<string> l; // all the lines
	string empty; // to prevent compile time warnings

	// loads output from a file
	void load(string &fname) {
		ifstream ifs(fname.c_str(), ios::in);

		string line;
		while (getline(ifs, line))
			l.push_back(line);
	}

	// returns indices to all lines that start with given, starting at start, and ending at end
	void find(string s, vector<int> &locs, int start = 0, int end = -1) {
		if (end == -1)
			end = l.size() - 1;

		// clear out locs
		locs.clear();

		// now check each line
		for (int i = start; i < end; i++) {
			//cout << "'" << s << "'" << ",'" << l[i].substr(0,s.length()) << "'" << endl;
			if (l[i].length() >= s.length() && l[i].substr(0, s.length()) == s)
				locs.push_back(i);
		}
	}

	// convert to a string
	string toString() {
		string s;
		for (unsigned int i = 0; i < l.size(); i++)
			s = s + l[i] + "\n";
		return (s);
	}

	// _safely_ get a line, returns empty string if there is no line
	//  if too slow, add a compile time #define to disable when necessary
	string &operator [](unsigned int index) {
		if (index < 0 || index >= l.size()) {
      Rprintf("Lines index %d is out of bounds [0,%d]\n", index, (l.size() - 1));
			return (empty);
		}

		return (l[index]);
	}

	// returns number of lines
	int size() {
		return (l.size());
	}
};

/** Haplotype phasing **/
//typedef Haplotype char; // 11, 12, 22 (also missing with zeros, but those are useless)
// DEBUGGED: Mostly just a storage class
class Haplotype {
public:
	vector<char> a;
	char empty;

	// convert to a string
	string toString() {
		string s;
		for (unsigned int i = 0; i < a.size(); i++) {
			if (a[i] == 0)
				s += '0';
			if (a[i] == 1)
				s += '1';
			if (a[i] == 2)
				s += '2';
		}
		return (s);
	}

	// number of genotypes in the haplotype
	unsigned int size() {
		return (a.size());
	}

	// returns value of haplotype at the index
	char &operator [](unsigned int index) {
		if (index < 0 || index >= a.size()) {
			Rprintf("Haplotype index %d is out of bounds [0,%d].\n", index, (a.size() - 1));
			return (empty);
		}

		return (a[index]);
	}
};

/** collection of phased haplotypes and their EM weights **/
// DEBUG: Fully debugged, has own debug() routine.
class Genotype {
public:
	// vector over all of the possible phases
	vector<Haplotype> ha;
	vector<Haplotype> hb;
	vector<double> emWeight; // Doesn't generally matter, WARNING, not set properly for parents (since I found that it didn't matter)

	// push a new genotype on
	void push_back(Haplotype &h1, Haplotype &h2, double emWeight) {
		this->ha.push_back(h1);
		this->hb.push_back(h2);
		this->emWeight.push_back(emWeight);
	}

	// convert to a string
	string toString() {
		if (ha.size() != hb.size() || hb.size() != emWeight.size()) {
			Rprintf("Genotype::toString() -- ha, hb, emWeight are not all the same size (%d, %d, %d).\n", ha.size(), hb.size(), emWeight.size());
			string empty;
			return (empty);
		}

		if (ha.size() == 1)
			return (ha[0].toString() + "|" + hb[0].toString());

		string s;
		for (unsigned int i = 0; i < ha.size(); i++)
			s += ha[i].toString() + "|" + hb[i].toString() + " ("
					+ doubleToString(emWeight[i]) + ")" + ", ";
		return (s);
	}

	// the number of phases possible
	int numPhases() {
		return (ha.size());
	}

	// coding of the genotype, assumes bi-allelic
	double xCode(int phase, unsigned int index, char allele, XCODE code) {
		// addition for missing (also changed return type of function to double, rather than int)
		if (ha[phase][index] == 0 || hb[phase][index] == 0)
			return (NAN);

		switch (code) {
		case ADDITIVE:
			return ((int) (ha[phase][index] == allele)
					+ (int) (hb[phase][index] == allele));
		case DOMINANT:
			return ((int) (ha[phase][index] == allele || hb[phase][index]
					== allele));
		case RECESSIVE:
			return ((int) (ha[phase][index] == allele && hb[phase][index]
					== allele));
		case GENOTYPE:
      Rprintf("Genotype::xcode::xCode -- should be using the genotype call, not xCode.\n");
		}
		Rprintf("Genotype::xcode::code misunderstood.\n");
		return (-1);
	}

	// indicator coding of the genotype
	double genotype(int phase, unsigned int index, char allele1, char allele2) {
		// addition for missing (also changed return type of function to double, rather than int)
		if (ha[phase][index] == 0 || hb[phase][index] == 0)
			return (NAN);

		return ((int) ((ha[phase][index] == allele1 && hb[phase][index]
				== allele2) || (ha[phase][index] == allele2 && hb[phase][index]
				== allele1)));
	}

	// is it missing?
	bool missing(int phase, unsigned int index) {
		return( ha[phase][index]==0 || hb[phase][index]==0 );
	}

	// New code for handling strata, taken from our old GxE code! yay! old code!
	static const int gMiss = -1;
	static const int gAA = 0;
	static const int gAB = 1;
	static const int gBB = 2;
	int gCode( int phase, unsigned int index ) {
		if( ha.size()==0 || phase!=0 || index<0 || index>=ha[phase].size() ) // currently only phase=0 is supported, doesn't matter
			return( gMiss );

		int a = ha[phase][index];
		int b = hb[phase][index];

		if( a==0 || b==0 ) return( gMiss );
		if( a==1 && b==1 ) return( gAA );
		if( a==2 && b==2 ) return( gBB );
		return( gAB );
	}

	void debug() {
		Rprintf("BEGIN Genotype::debug\n");

		Genotype g;
		Haplotype h1, h2;
		char DISEASE_ALLELE = 2;
		int numTest = 4; // length of "haplotype", although, really, we'll just be testing each allele in turn
		h1.a.push_back( 0 ); h2.a.push_back( 0 ); // first is missing
		h1.a.push_back( 1 ); h2.a.push_back( 1 ); // second is homozygous minor allele
		h1.a.push_back( 1 ); h2.a.push_back( 2 ); // third is heterozygous
		h1.a.push_back( 2 ); h2.a.push_back( 2 ); // fourth is homozygous disease allele

		// em weights don't really matter, so not going to test here
		ha.push_back( h1 );
		hb.push_back( h2 );
		emWeight.push_back( 1.0 );

		// Now test the functions
		Rprintf("numPhases (should be 1) %d\n", numPhases());

    /*
		const char *genoStrs[] = {"Missing", "Homozygous minor allele", "Heterozygous", "Homozygous major allele"};
		for( int i=0; i<numTest; i++ ) {
			cout << "GENOTYPE: " << genoStrs[i] << endl;
			// first index is zero in all of these, as that is the phase (and there is only one phase for testing, since that's all that is really used in this code anyway.
			cout << "xCode ADDITIVE " << xCode( 0, i, DISEASE_ALLELE, ADDITIVE ) << endl;
			cout << "xCode DOMINANT " << xCode( 0, i, DISEASE_ALLELE, DOMINANT ) << endl;
			cout << "xCode RECESSIVE " << xCode( 0, i, DISEASE_ALLELE, RECESSIVE ) << endl;
			cout << "genotype 1/1 " << genotype( 0, i, 1, 1 ) << endl;
			cout << "genotype 1/2 " << genotype( 0, i, 1, 2 ) << endl;
			cout << "genotype 2/1 " << genotype( 0, i, 2, 1 ) << endl;
			cout << "genotype 2/2 " << genotype( 0, i, 2, 2 ) << endl;
			cout << "missing " << missing( 0, i ) << endl;
			cout << "gCode (gMiss=-1, gAA=0, gAB=1, gBB=2) " << gCode( 0, i ) << endl;
		}

		cout << "END Genotype::debug" << endl;
    */
    Rprintf("((Rest of lines are eliminated.))\n");
	}
};

/** Do haplotypes match? (ignoring phase here) **/
bool unphaseMatch(Haplotype &ha1, Haplotype &hb1, Haplotype &ha2,
		Haplotype &hb2) {
	if (ha1.size() != hb1.size() || ha1.size() != ha2.size() || ha1.size()
			!= hb2.size()) {
		Rprintf("unphaseMatch() -- haplotypes are not the same size!\n");
		return (false);
	}
	for (unsigned int i = 0; i < ha1.size(); i++)
		if (!((ha1[i] == ha2[i] && hb1[i] == hb2[i]) || (ha1[i] == hb2[i]
				&& hb1[i] == ha2[i])))
			return (false);
	return (true);
}

/** Pedigree class, does most of the work **/
// DEBUG: See internally for status
class Pedigree {
public:
	vector<Genotype> g; // compatible offspring genotypes
	vector<double> pg;

	vector<int> observed; // the observed configuration (indices into g)
	vector<double> trait;
	vector<double> traitBackup; // 01/05/2009

	vector<vector<int> > genoDist; // distribution of the genotypes
	vector<double> genoDistP; // probability of each of these

	vector<bool> nonzeroDelX; // nonzero Xc-E(Xc|S); only use these to compute the offset, etc.

	string label; // really just for debugging, contains the pedigree id

	unsigned int pid; // actually it's best to go further and parse it in, to make _sure_ then that things go ahead and link properly later

	Genotype parents[2];

	// parses out FBAT output
	// -- ? should handle longer than 2-length haplotypes?
	// DEBUG: Works on all cases tested, should pretty much completely crash the program (or cryptic error messages) if it goes awry
	void parse(Lines &lines, int start, int end) {
		StrTok tok;

		// New 01/09/09 addition
		// Read in the observed parental haplotypes
		/*  observed genotype configuration
		 *  father         = 0   0   0   0   0   0   0   0   0
		 *                   0   0   0   0   0   0   0   0   0
		 *  mother         = 0   0   0   0   0   0   0   0   0
		 *                   0   0   0   0   0   0   0   0   0
		 *  offspring 1    = 4   1   3   1   4   3   2   4   1
		 *                   4   3   2   2   4   1   1   2   3
		 *  offspring 2    = 4   1   3   1   4   3   2   4   1
		 *                   4   3   2   2   4   1   1   2   3
		 *
		 */
		vector<int> obsLoc;
		lines.find("observed", obsLoc, start, end);
		if( obsLoc.size()==0 ) {
			Rprintf("No parental information ('observed' string not found) available for subject, lines %d to %d", start, end);
		}else{
			// keep going!
			int newStart = obsLoc[0];

			const char* parentStr[] = {"father","mother"};

			for( int p=0; p<2; p++ ) {
				// Parse the father/mother
				obsLoc.clear();
				lines.find(parentStr[p], obsLoc, newStart, end);
				if( obsLoc.size()>0 && obsLoc[0]<(end-1) ) {
					// father/mother was found
					Haplotype ha, hb;

					// the first line
					tok.tokenize(lines[obsLoc[0]]);
					tok.nextTokenE(); // father (or mother)
					tok.nextTokenE(); // =
					while( tok.hasMoreTokens() )
						ha.a.push_back((char)tok.nextTokenN());

					// the second line
					tok.tokenize(lines[obsLoc[0]+1]);
					while( tok.hasMoreTokens() )
						hb.a.push_back((char)tok.nextTokenN());

					parents[p].push_back( ha, hb, 1 ); // EM weight doesn't matter for this.
				}else{
					Rprintf("%s could not be found.\n", parentStr[p]);
				}
			}

			// Parse the mother (almost identical to father)
		}

		// First read in the offspring haplotypes
		// - These will need to be translated into genotypes later
		/*
		 * offspring 1    = 1   1
		 *                  2   2
		 * offspring 2    = 1   2
		 *                  1   2
		 */
		vector<int> offspringLocs;
		lines.find("offspring", offspringLocs, start, end);
		vector<Haplotype> obsA, obsB;
		for (unsigned int i = 0; i < offspringLocs.size(); i++) {
			Haplotype ha, hb;

			// The first line
			tok.tokenize(lines[offspringLocs[i]]);
			tok.nextTokenE(); // offspring
			tok.nextTokenE(); // 1
			tok.nextTokenE(); // =
			while (tok.hasMoreTokens())
				ha.a.push_back((char)tok.nextTokenN());

			// The second line
			tok.tokenize(lines[offspringLocs[i] + 1]);
			while (tok.hasMoreTokens())
				hb.a.push_back((char)tok.nextTokenN());

			if (ha.a.size() != hb.a.size())
				Rprintf("Offspring haplotypes are not the same size!\n");

			obsA.push_back(ha);
			obsB.push_back(hb);
		}

		// now read in the genotypes...
		/*
		 compatible offspring genotype 1 (g1), phase 1; EM P=1.000
		 h1   :   1   2
		 h1   :   1   2

		 compatible offspring genotype 2 (g2), phase 1; EM P=0.333
		 h3   :   2   2
		 h2   :   1   1

		 compatible offspring genotype 2 (g2), phase 2; EM P=0.667
		 h1   :   1   2
		 h4   :   2   1

		 distribution of compatible offspring genotype configurations:

		 #g1     #g2     P[G]
		 1       1       1.000

		 pg=
		 g1      0.500
		 g2      0.500
		 */

		/* NEW: FORMAT 2 (only when solving for a one dimensional nuisance parameter)
		 There are 3 compatible offspring genotypes:

		 compatible offspring genotype 1 (g1) = 2/2

		 compatible offspring genotype 2 (g2) = 2/1

		 compatible offspring genotype 3 (g3) = 1/1

		 distribution of compatible offspring genotype configurations:

		 #g1     #g2     #g3     P[G]
		 0       1       0       0.500
		 1       0       0       0.250
		 0       0       1       0.250

		 pg=
		 g1      0.250
		 g2      0.500
		 g3      0.250
		 */
		vector<int> compatibleHaplotypeLocs;
		lines.find("compatible offspring genotype", compatibleHaplotypeLocs,
				start, end);

		int numGenos = 0;
		if (compatibleHaplotypeLocs.size() > 0) {

			int ii =
					compatibleHaplotypeLocs[compatibleHaplotypeLocs.size() - 1];
			tok.tokenize(lines[ii]);
			tok.nextTokenE(); // compatible
			tok.nextTokenE(); // offspring
			tok.nextTokenE(); // genotype
			numGenos = (int)tok.nextTokenN(); // genotype number
			if (numGenos < 0) {
				Rprintf("Pedigree::numGenos < 0!\n");
				numGenos = 0;
			}
			g.resize(numGenos);

			// 07/23/2008 begin
			//  "compatible offspring genotype configuration not found"
			if (numGenos == 0)
				return;
			// 07/23/2008 end

			for (unsigned int i = 0; i < compatibleHaplotypeLocs.size(); i++) {
				tok.tokenize(lines[compatibleHaplotypeLocs[i]]);
				tok.nextTokenE(); // compatible
				tok.nextTokenE(); // offspring
				tok.nextTokenE(); // genotype
				int curGeno = (int)tok.nextTokenN() - 1;
				if (curGeno < 0 || curGeno >= numGenos) {
					Rprintf("Pedigree::parse -- curGeno isn't proper...\n");
					Rprintf(" Offending line (%d) = '%s'\n", compatibleHaplotypeLocs[i], lines[compatibleHaplotypeLocs[i]].c_str());
					Rprintf(" curGeno=%d numGenos=%d\n", curGeno, numGenos);
					continue; // skip to the next one then...
				}

				Haplotype ha, hb;
				double emWeight = 1.0;

				// addition in case we've got the other method
				tok.nextTokenE(); // (g1)
				if (tok.nextToken() == "=") {
					// then it's the latter of the two methods...
					//cout << "We've got the new method!" << endl;
					string numSlashNum = tok.nextToken();
					//cout << numSlashNum << endl;
					string token = "/";
					tok.tokenize(numSlashNum, token);
					ha.a.push_back((char)tok.nextTokenN());
					hb.a.push_back((char)tok.nextTokenN());
				} else {
					// it's the usual method...

					string p_eq_em = tok.lastToken();
					string delim = "=";
					tok.tokenize(p_eq_em, delim);
					emWeight = tok.lastTokenN();

					tok.tokenize(lines[compatibleHaplotypeLocs[i] + 1]);
					tok.nextTokenE(); // h1
					tok.nextTokenE(); // :
					while (tok.hasMoreTokens())
						ha.a.push_back((char)tok.nextTokenN());

					tok.tokenize(lines[compatibleHaplotypeLocs[i] + 2]);
					tok.nextTokenE(); // h1
					tok.nextTokenE(); // :
					while (tok.hasMoreTokens())
						hb.a.push_back((char)tok.nextTokenN());
				}

				if (curGeno >= 0 && curGeno < numGenos) {
					g[curGeno].ha.push_back(ha);
					g[curGeno].hb.push_back(hb);
					g[curGeno].emWeight.push_back(emWeight);
				} else {
					Rprintf("curGeno past genotype size...\n");
				}
			}

			pg.resize(g.size());
			vector<int> pgLoc;
			lines.find("pg=", pgLoc, start, end);
			if( pgLoc.size() > 0 ){
				for (unsigned int i = 0; i < g.size(); i++) {
					tok.tokenize(lines[pgLoc[0] + i + 1]);
					pg[i] = tok.lastTokenN();
				}
			} else {
			}
		}

		// now need to translate the offspring haplotypes into the genotypes...
		for (unsigned int i = 0; i < obsA.size(); i++) {
			bool found = false;

			// which g matches?
			for (unsigned int j = 0; j < g.size() && !found; j++) {
				// should only need to do the first, since others should be a phase of it
				if (unphaseMatch(obsA[i], obsB[i], g[j].ha[0], g[j].hb[0])) {
					observed.push_back(j);
					found = true;
				}
			}
		}

		// lastly read in the offspring genotype configurations
		vector<int> distLoc;
		lines.find("distribution of compatible", distLoc, start, end);
		if (distLoc.size() > 0) {
			// the header line
			tok.tokenize(lines[distLoc[0] + 2]);

			if (numGenos != (int) tok.size() - 1) {
				Rprintf("numGenos does not match distribution numGenos!\n");
			}

			int nextlineloc = distLoc[0] + 3;
			string nextline = lines[nextlineloc];
			while (nextline != "" && nextlineloc <= end) {
				// make sure everything is formatted correctly...
				tok.tokenize(nextline);
				if (numGenos != (int) tok.size() - 1) {
					nextline = ""; // force exit
					continue;
				}

				// parse in the information
				vector<int> counts;
				for (int t = 0; t < numGenos; t++)
					counts.push_back((int) tok.nextTokenN());
				genoDist.push_back(counts);
				genoDistP.push_back(tok.lastTokenN());

				// and update
				nextlineloc++;
				nextline = lines[nextlineloc];
			}
		}
	}

	// returns all possible enumerations of the phases for the genotypes
	// DEBUG: No longer used, but should be good (the routine just following this one is for debugging this routine)
	void enumPhases(vector<vector<int> > &phases, vector<double> &phaseWeight) {
		if (g.size() == 0)
			return;

		// handle the first one
		vector<int> nPhase;
		nPhase.push_back(0);
		for (int k = 0; k < g[0].numPhases(); k++) {
			nPhase[0] = k;
			phases.push_back(nPhase);
		}

		// now handle the rest
		for (unsigned int i = 1; i < g.size(); i++) {
			// push on the first (and possibly only) phase
			for (unsigned int p = 0; p < phases.size(); p++)
				phases[p].push_back(0);
			int ncopy = phases.size();

			for (int j = 1; j < g[i].numPhases(); j++) {
				for (int c = 0; c < ncopy; c++) {
					phases.push_back(phases[c]);
					phases[phases.size() - 1][i] = j;
				}
			}
		}

		// now get the phaseWeights
		for (unsigned int p = 0; p < phases.size(); p++) {
			double weight = g[0].emWeight[phases[p][0]];
			for (unsigned int j = 0; j < phases[p].size(); j++)
				weight *= g[j].emWeight[phases[p][j]];
			phaseWeight.push_back(weight);
		}
	}

	// really just for debugging enumPhases
	string phasesString() {
		vector<vector<int> > phases;
		vector<double> phasesWeight;
		enumPhases(phases, phasesWeight);
		//cout << "phases.size = " << phases.size() << endl;

		string s = "Phases:\n";
		for (unsigned int p = 0; p < phases.size(); p++) {
			s += " ";
			for (unsigned int j = 0; j < phases[p].size(); j++)
				s += doubleToString(phases[p][j]) + " ";
			s += " (" + doubleToString(phasesWeight[p]) + ")";
			s += "\n";
		}
		return (s);
	}

	// convert to a string
	// DEBUGGED
	string toString() {
		string s;

		s += label + "\n";

		s += " (" + doubleToString(pid) + ")\n";

		s += "Parents\n";
		for( int p=0; p<2; p++ )
			s += " p" + doubleToString(p) + ")" + parents[p].toString() + "\n";

		s += "Observed: ";
		for (unsigned int i = 0; i < observed.size(); i++)
			s += doubleToString(observed[i]) + ",";
		s += "\n";

		if (trait.size() > 0) {
			s += "Observed trait: ";
			for (unsigned int i = 0; i < trait.size(); i++)
				s += doubleToString(trait[i]) + ",";
			s += "\n";
		}

		s += "Genotype probabilities:\n";
		for (unsigned int i = 0; i < g.size(); i++)
			s += " g" + doubleToString(i) + ")" + doubleToString(pg[i]) + ": "
					+ g[i].toString() + "\n";

		s += "Distribution:\n";
		for (unsigned int i = 0; i < genoDistP.size(); i++) {
			for (unsigned int j = 0; j < genoDist[i].size(); j++)
				s += doubleToString(genoDist[i][j]) + " ";
			s += doubleToString(genoDistP[i]) + "\n";
		}

		s += phasesString(); // new, debug

		return (s);
	}

	// Was going to return if a pedigree was informative,
	//  but with the new fbat, we can make it only print the
	//  families which _are_ informative
	bool informative() {
		return (true); // fixed with the newest fbat!
		//return( observed.size() > 0 );
	}

	// DEPRECATED
	bool phaseInferrable() {
		// does any of the genotypes have more than one phase?
		for (unsigned int k = 0; k < g.size(); k++)
			if (g[k].numPhases() > 1)
				return (false);

		// no, they do not
		return (true);
	}

	// for dichotomous traits
	void uimc(double *bm, double *bc0, double *bc1, int *analyze_allele_index,
			int analyze_allele_index_size, int *conditional_allele_index,
			int conditional_allele_index_size, bool onlyComputeConditional,
			double *ret_analyze, double *ret_conditional0,
			double *ret_conditional1) {
		//cout << "ENTERED UIMC" << endl;
		if (observed.size() == 0) {
			for (int a = 0; a < analyze_allele_index_size; a++)
				ret_analyze[a] = 0.0;
			for (int a = 0; a < conditional_allele_index_size; a++) {
				ret_conditional0[a] = 0.0;
				ret_conditional1[a] = 0.0;
			} // a
			return;
		}

		// compute the observed
		//double sum_xijm[analyze_allele_index_size];
		vector<double> sum_xijm; sum_xijm.resize(analyze_allele_index_size);
		if (!onlyComputeConditional) {
			// new, only compute the conditional if asked to save computation time
			for (int a = 0; a < analyze_allele_index_size; a++) {
				sum_xijm[a] = 0.0;
				for (unsigned int j = 0; j < observed.size(); j++)
					if (trait[j] == 1)
						sum_xijm[a] += g[observed[j]].xCode(0,
								analyze_allele_index[a], 2, ADDITIVE);
			}
		}

		//double sum_xijc0[conditional_allele_index_size];
		//double sum_xijc1[conditional_allele_index_size];
		vector<double> sum_xijc0; sum_xijc0.resize(conditional_allele_index_size);
    vector<double> sum_xijc1; sum_xijc1.resize(conditional_allele_index_size);
		for (int a = 0; a < conditional_allele_index_size; a++) {
			sum_xijc0[a] = sum_xijc1[a] = 0.0;
			for (unsigned int j = 0; j < observed.size(); j++) {
				if (trait[j] == 1) {
					sum_xijc0[a] += g[observed[j]].genotype(0,
							conditional_allele_index[a], 2, 2);
					sum_xijc1[a] += g[observed[j]].genotype(0,
							conditional_allele_index[a], 1, 2);
				}
			}
		}

		// Now calculate the other pieces
		//double numM[analyze_allele_index_size];
		vector<double> numM; numM.resize(analyze_allele_index_size);
		for (int a = 0; a < analyze_allele_index_size; a++)
			numM[a] = 0.0;

		//double numC0[conditional_allele_index_size];
		//double numC1[conditional_allele_index_size];
    vector<double> numC0; numC0.resize(conditional_allele_index_size);
    vector<double> numC1; numC1.resize(conditional_allele_index_size);
		for (int a = 0; a < conditional_allele_index_size; a++)
			numC0[a] = numC1[a] = 0.0;

		double den = 0.0;

		for (unsigned int x = 0; x < genoDist.size(); x++) {
			//cout << "BEFORE CREATING GENOTYPE PERMUTATIONS" << endl;
			// create all genotype permutations
			vector<int> genoToPerm;
			for (unsigned int j = 0; j < genoDist[x].size(); j++)
				for (int k = 0; k < genoDist[x][j]; k++)
					genoToPerm.push_back(j);
			vector<vector<int> > genoPerm;
			allPerms(genoToPerm, genoPerm);
			//cout << "AFTER CREATING GENOTYPE PERMUTATIONS" << endl;

			//cout << "BEFORE COMPUTING TEST STATISTIC" << endl;
			// calculate the test statistic
			double permWeight = (double) 1.0 / (double) genoPerm.size();
			for (unsigned int p = 0; p < genoPerm.size(); p++) {
				//double sum_xstarjm[analyze_allele_index_size];
				vector<double> sum_xstarjm; sum_xstarjm.resize(analyze_allele_index_size);
				if (!onlyComputeConditional) {
					for (int a = 0; a < analyze_allele_index_size; a++) {
						sum_xstarjm[a] = 0.0;
						for (unsigned int j = 0; j < genoPerm[p].size(); j++) {
							if (trait[j] == 1)
								sum_xstarjm[a] += g[genoPerm[p][j]].xCode(0,
										analyze_allele_index[a], 2, ADDITIVE);
						} // j
					} // a
				}

				//double sum_xstarjc0[conditional_allele_index_size];
				//double sum_xstarjc1[conditional_allele_index_size];
				vector<double> sum_xstarjc0; sum_xstarjc0.resize(conditional_allele_index_size);
        vector<double> sum_xstarjc1; sum_xstarjc1.resize(conditional_allele_index_size);
				for (int a = 0; a < conditional_allele_index_size; a++) {
					sum_xstarjc0[a] = 0.0;
					sum_xstarjc1[a] = 0.0;
					for (unsigned int j = 0; j < genoPerm[p].size(); j++) {
						if (trait[j] == 1) {
							sum_xstarjc0[a] += g[genoPerm[p][j]].genotype(0,
									conditional_allele_index[a], 2, 2);
							sum_xstarjc1[a] += g[genoPerm[p][j]].genotype(0,
									conditional_allele_index[a], 1, 2);
						}
					} // j
				} // a

				double temp = 0.0; // exp( b0*x0 + ... )
				if (!onlyComputeConditional)
					for (int a = 0; a < analyze_allele_index_size; a++)
						temp += bm[a] * sum_xstarjm[a];
				for (int a = 0; a < conditional_allele_index_size; a++)
					temp += (bc0[a] * sum_xstarjc0[a] + bc1[a]
							* sum_xstarjc1[a]);
				temp = exp(temp) * genoDistP[x] * permWeight;

				if (!onlyComputeConditional)
					for (int a = 0; a < analyze_allele_index_size; a++)
						numM[a] += sum_xstarjm[a] * temp;

				for (int a = 0; a < conditional_allele_index_size; a++) {
					numC0[a] += sum_xstarjc0[a] * temp;
					numC1[a] += sum_xstarjc1[a] * temp;
				}

				den += temp;
			} // p
			//cout << "AFTER COMPUTING TEST STATISTIC" << endl;
		} // x

		if (!onlyComputeConditional)
			for (int a = 0; a < analyze_allele_index_size; a++)
				ret_analyze[a] = sum_xijm[a] - (numM[a] / den);

		for (int a = 0; a < conditional_allele_index_size; a++) {
			ret_conditional0[a] = sum_xijc0[a] - (numC0[a] / den);
			ret_conditional1[a] = sum_xijc1[a] - (numC1[a] / den);
		} // a

		//cout << "ret_analyze[0] " << ret_analyze[0] << " " << ret_conditional0[0] << " " << ret_conditional1[0] << " " << sum_xijm[0] << " " << sum_xijc0[0] << " " << sum_xijc1[0] << endl;

		// COMPLETED WRITING, BUT _NOT_ DEBUGGED (OR EVEN COMPILED FOR THAT MATTER!)
	} // uimc

	// stored by _row_
	int imcIndex(int r, int c, int R) {
		return (r + c * R);
	}
	// dichotomous trait, for computing the model-based variance
	void imc(double *bm, double *bc0, double *bc1, int *analyze_allele_index,
			int analyze_allele_index_size, int *conditional_allele_index,
			int conditional_allele_index_size, double *I) {
		int R = analyze_allele_index_size + 2 * conditional_allele_index_size;

		if (observed.size() == 0) {
			int Rsq = R * R;
			for (int r = 0; r < Rsq; r++)
				I[r] = 0.0;
			return;
		}

		double ai0 = 0.0;
		//double ai1[R];
    vector<double> ai1; ai1.resize(R);
		//double ai2[R][R];
    vector<vector<double> > ai2; ai2.resize(R); for(int r=0; r<R; r++) ai2[r].resize(R);
		for (int k1 = 0; k1 < R; k1++) {
			ai1[k1] = 0.0;
			for (int k2 = 0; k2 < R; k2++)
				ai2[k1][k2] = 0.0;
		}

		for (unsigned int x = 0; x < genoDist.size(); x++) {
			// create all genotype permutations
			vector<int> genoToPerm;
			for (unsigned int j = 0; j < genoDist[x].size(); j++)
				for (int k = 0; k < genoDist[x][j]; k++)
					genoToPerm.push_back(j);
			vector<vector<int> > genoPerm;
			allPerms(genoToPerm, genoPerm);

			// calculate the test statistic pieces
			double permWeight = (double) 1.0 / (double) genoPerm.size();
			for (unsigned int p = 0; p < genoPerm.size(); p++) {
				//double sum_xstarjm[analyze_allele_index_size];
				vector<double> sum_xstarjm; sum_xstarjm.resize(analyze_allele_index_size);
				if (true) { // !onlyComputeConditional ) {
					for (int a = 0; a < analyze_allele_index_size; a++) {
						sum_xstarjm[a] = 0.0;
						for (unsigned int j = 0; j < genoPerm[p].size(); j++) {
							if (trait[j] == 1)
								sum_xstarjm[a] += g[genoPerm[p][j]].xCode(0,
										analyze_allele_index[a], 2, ADDITIVE);
						} // j
					} // a
				}

				//double sum_xstarjc0[conditional_allele_index_size];
				//double sum_xstarjc1[conditional_allele_index_size];
				vector<double> sum_xstarjc0; sum_xstarjc0.resize(conditional_allele_index_size);
        vector<double> sum_xstarjc1; sum_xstarjc1.resize(conditional_allele_index_size);
				for (int a = 0; a < conditional_allele_index_size; a++) {
					sum_xstarjc0[a] = 0.0;
					sum_xstarjc1[a] = 0.0;
					for (unsigned int j = 0; j < genoPerm[p].size(); j++) {
						if (trait[j] == 1) {
							sum_xstarjc0[a] += g[genoPerm[p][j]].genotype(0,
									conditional_allele_index[a], 2, 2);
							sum_xstarjc1[a] += g[genoPerm[p][j]].genotype(0,
									conditional_allele_index[a], 1, 2);
						}
					} // j
				} // a

				double temp = 0.0; // exp( b0*x0 + ... )
				if (true) // !onlyComputeConditional )
					for (int a = 0; a < analyze_allele_index_size; a++)
						temp += bm[a] * sum_xstarjm[a];
				for (int a = 0; a < conditional_allele_index_size; a++)
					temp += (bc0[a] * sum_xstarjc0[a] + bc1[a]
							* sum_xstarjc1[a]);
				temp = exp(temp) * genoDistP[x] * permWeight;

				// fill in xstar from xstarjm, xstarjc0, xstarjc1
				//  for ease of computation (or code rather)
				//double sum_xstar[analyze_allele_index_size + 2 * conditional_allele_index_size];
        vector<double> sum_xstar; sum_xstar.resize(analyze_allele_index_size + 2 * conditional_allele_index_size);
				int ca = 0;

				for (int a = 0; a < analyze_allele_index_size; a++) {
					sum_xstar[ca] = sum_xstarjm[a];
					ca++;
				}
				for (int a = 0; a < conditional_allele_index_size; a++) {
					sum_xstar[ca] = sum_xstarjc0[a];
					ca++;
				}
				for (int a = 0; a < conditional_allele_index_size; a++) {
					sum_xstar[ca] = sum_xstarjc1[a];
					ca++;
				}

				ai0 += temp;
				for (int k1 = 0; k1 < R; k1++) {
					ai1[k1] += sum_xstar[k1] * temp;
					for (int k2 = 0; k2 < R; k2++) { // somewhat inefficient
						ai2[k1][k2] += sum_xstar[k1] * sum_xstar[k2] * temp;
					} // k2
				} // k1
			} // p
		} // x

		// finally compute I!
		double ai0sq = ai0 * ai0;
		for (int k1 = 0; k1 < R; k1++)
			for (int k2 = 0; k2 < R; k2++)
				I[imcIndex(k1, k2, R)]
						= (ai2[k1][k2] * ai0 - ai1[k1] * ai1[k2]) / ai0sq;
	} // imc

	// migrating only to empirical variance
	// dichotomous / continuous trait
	void robustStat(int *analyze_allele_index, int analyze_allele_index_size,
			int *conditional_allele_index, int conditional_allele_index_size,
			double *ret_analyze) {
		for (int aa = 0; aa < analyze_allele_index_size; aa++)
			ret_analyze[aa] = 0.0;

		for (unsigned int j = 0; j < observed.size(); j++) {
			if (!isnan(trait[j])) {
				double weight = 0.0;
				//double ex[analyze_allele_index_size];
        vector<double> ex; ex.resize(analyze_allele_index_size);
				for (int aa = 0; aa < analyze_allele_index_size; aa++)
					ex[aa] = 0.0;

				for (unsigned int k = 0; k < g.size(); k++) {
					// make sure all alleles match
					bool good = true;
					for (int ca = 0; good && ca < conditional_allele_index_size; ca++)
						if (g[k].genotype(
								0,
								conditional_allele_index[ca],
								g[observed[j]].ha[0][conditional_allele_index[ca]],
								g[observed[j]].hb[0][conditional_allele_index[ca]])
								== 0)
							good = false;

					if (good) {
						// all alleles match
						weight += pg[k];
						for (int aa = 0; aa
								< analyze_allele_index_size; aa++)
							ex[aa] += pg[k] * g[k].xCode(0,
									analyze_allele_index[aa], 2, ADDITIVE);
					}
				} // k

				for (int aa = 0; aa < analyze_allele_index_size; aa++) {
					ex[aa] /= weight;
					ret_analyze[aa] += trait[j] * (g[observed[j]].xCode(0,
							analyze_allele_index[aa], 2, ADDITIVE) - ex[aa]);
				} // aa
			} // !isnan(...)
		} // j
	} // robustStat

	void contsX(int *analyze_allele_index, int analyze_allele_index_size,
			int *conditional_allele_index, int conditional_allele_index_size,
      int gIndex, vector<double> &x) {
    //int gIndex, double *x) {
      for (int a = 0; a < analyze_allele_index_size; a++) {
			x[a] = g[gIndex].xCode(0, analyze_allele_index[a], 2, ADDITIVE);
		}
		for (int a = 0; a < conditional_allele_index_size; a++) {
			x[a + analyze_allele_index_size] = g[gIndex].genotype(0,
					conditional_allele_index[a], 2, 2);
			x[a + analyze_allele_index_size + conditional_allele_index_size]
					= g[gIndex].genotype(0, conditional_allele_index[a], 1, 2);
		}
	}

	double contsBetaX(double *b, double *x, int R) {
		double bx = 0.0;
		for (int r = 0; r < R; r++)
			bx += (b[r] * x[r]);
		return (bx);
	}

	double contsA(double alpha, double sigmaSq, double *b, double *x,
			double bx, int j, int *analyze_allele_index,
			int analyze_allele_index_size, int *conditional_allele_index,
			int conditional_allele_index_size, bool ignoreBtX) {
		//
		if (isnan(trait[j]))
			return (0.0);

		if (ignoreBtX)
			return ((trait[j] - alpha) * bx / sigmaSq);
		return ((-0.5 * bx * bx + (trait[j] - alpha) * bx) / sigmaSq);
	} // contsA

	void contsB(double alpha, double sigmaSq, double *b, double *x, double bx,
			int j, int *analyze_allele_index, int analyze_allele_index_size,
			int *conditional_allele_index, int conditional_allele_index_size,
      bool ignoreBtX, vector<double> &res) {
      //bool ignoreBtX, double *res) {
    int R = analyze_allele_index_size + conditional_allele_index_size * 2;

		if (isnan(trait[j])) {
			for (int r = 0; r < R; r++)
				res[r] = 0.0;
			return;
		}

		for (int r = 0; r < R; r++) {
			if (ignoreBtX) {
				res[r] = ((trait[j] - alpha) * x[r]) / sigmaSq;
			} else {
				res[r] = (-bx * x[r] + (trait[j] - alpha) * x[r]) / sigmaSq;
			}
		}
	} // contsB

	void contsC(double alpha, double sigmaSq, double *b, double *x, double bx,
			int j, int *analyze_allele_index, int analyze_allele_index_size,
			int *conditional_allele_index, int conditional_allele_index_size,
			bool ignoreBtX, double *res) {
		//
		int R = analyze_allele_index_size + conditional_allele_index_size * 2;

		if (isnan(trait[j])) {
			int Rsq = R * R;
			for (int r = 0; r < Rsq; r++)
				res[r] = 0.0;
			return;
		} // isnan(...)

		for (int k1 = 0; k1 < R; k1++) {
			for (int k2 = 0; k2 < R; k2++) {
				if (ignoreBtX) {
					res[imcIndex(k1, k2, R)] = 0.0;
				} else {
					res[imcIndex(k1, k2, R)] = -x[k1] * x[k2] / sigmaSq;
				}
			} // k1
		} // k2
	} // contsC

	// Deprecated - continuous, normal model
	void contsUimc(double alpha,
			double sigmaSq,
			double *b, // bm, bc0, bc1
			int *analyze_allele_index, int analyze_allele_index_size,
			int *conditional_allele_index, int conditional_allele_index_size,
			bool onlyComputeConditional, bool ignoreBtX, double *ret_b) {
		//
		int R = analyze_allele_index_size + conditional_allele_index_size * 2;

		if (observed.size() == 0) {
			for (int r = 0; r < R; r++)
				ret_b[r] = 0.0;
			return;
		}

		// calculate the observed part of the vector
		//double sumB[R];
		vector<double> sumB; sumB.resize(R);
		for (int r = 0; r < R; r++)
			sumB[r] = 0.0;
		//double X[R];
		//double sumBAddi[R];
    vector<double> X; X.resize(R);
    vector<double> sumBAddi; sumBAddi.resize(R);
		for (unsigned int j = 0; j < observed.size(); j++) {
			contsX(analyze_allele_index, analyze_allele_index_size,
					conditional_allele_index, conditional_allele_index_size,
					observed[j], X);
			contsB(alpha, sigmaSq, b, (double *)X.data(), contsBetaX(b, (double *)X.data(), R), j,
					analyze_allele_index, analyze_allele_index_size,
					conditional_allele_index, conditional_allele_index_size,
					ignoreBtX, sumBAddi);

			//for( unsigned int a=0; a<analyze_allele_index_size; a++ ) {
			//  sumB[a] += sumBAddi[a];
			//} // a

			for (int r = 0; r < R; r++)
				sumB[r] += sumBAddi[r];
		} // j

		// now calculate the expected piece
		//double num[R];
		vector<double> num; num.resize(R);
		for (int r = 0; r < R; r++)
			num[r] = 0.0;
		double den = 0.0;
		for (unsigned int x = 0; x < genoDist.size(); x++) {
			// create all genotype permutations
			vector<int> genoToPerm;
			for (unsigned int j = 0; j < genoDist[x].size(); j++)
				for (int k = 0; k < genoDist[x][j]; k++)
					genoToPerm.push_back(j);
			vector<vector<int> > genoPerm;
			allPerms(genoToPerm, genoPerm);

			// calculate the test statistic
			double permWeight = (double) 1.0 / (double) genoPerm.size();
			for (unsigned int p = 0; p < genoPerm.size(); p++) {
				double sumAStar = 0.0;
				//double sumBStar[R];
        vector<double> sumBStar; sumBStar.resize(R);
				for (int r = 0; r < R; r++)
					sumBStar[r] = 0.0;

				for (unsigned int j = 0; j < genoPerm[p].size(); j++) {
					//double XStar[R];
					vector<double> XStar; XStar.resize(R);
					contsX(analyze_allele_index, analyze_allele_index_size,
							conditional_allele_index,
							conditional_allele_index_size, genoPerm[p][j],
							XStar);
					double bXStar = contsBetaX(b, (double *)XStar.data(), R);

					sumAStar += contsA(alpha, sigmaSq, b, (double *)XStar.data(), bXStar, j,
							analyze_allele_index, analyze_allele_index_size,
							conditional_allele_index,
							conditional_allele_index_size, ignoreBtX);

					contsB(alpha, sigmaSq, b, (double *)XStar.data(), bXStar, j,
							analyze_allele_index, analyze_allele_index_size,
							conditional_allele_index,
							conditional_allele_index_size, ignoreBtX, sumBAddi);
					for (int r = 0; r < R; r++)
						sumBStar[r] += sumBAddi[r];
				} // j

				double temp = exp(sumAStar) * genoDistP[x] * permWeight;
				for (int r = 0; r < R; r++)
					num[r] += sumBStar[r] * temp;
				den += temp;
			} // p
		} // x

		// finally calculate the return pieces!
		for (int r = 0; r < R; r++)
			ret_b[r] = sumB[r] - (num[r] / den);
	} // contsUimc

	// Deprecated, continuous normal model
	void contsImc(double alpha,
			double sigmaSq,
			double *b, // bm, bc0, bc1
			int *analyze_allele_index, int analyze_allele_index_size,
			int *conditional_allele_index, int conditional_allele_index_size,
			bool ignoreBtX, double *ret_I) {
		//
		int R = analyze_allele_index_size + conditional_allele_index_size * 2;
		int Rsq = R * R;

		if (observed.size() == 0) {
			for (int rr = 0; rr < Rsq; rr++)
				ret_I[rr] = 0.0;
			return;
		}

		// calculate the observed part of the vector
		//double X[R];
		//double sumC[Rsq];
		//double sumCAddi[Rsq];
		//double sumBAddi[R];
		vector<double> X; X.resize(R);
    vector<double> sumC; sumC.resize(R);
    vector<double> sumCAddi; sumCAddi.resize(Rsq);
    vector<double> sumBAddi; sumBAddi.resize(R);
		for (int rr = 0; rr < Rsq; rr++)
			sumC[rr] = 0.0;
		for (unsigned int j = 0; j < observed.size(); j++) {
			contsX(analyze_allele_index, analyze_allele_index_size,
					conditional_allele_index, conditional_allele_index_size,
					observed[j], X);
			contsC(alpha, sigmaSq, b, (double *)X.data(), contsBetaX(b, (double *)X.data(), R), j,
					analyze_allele_index, analyze_allele_index_size,
					conditional_allele_index, conditional_allele_index_size,
					ignoreBtX, (double *)sumCAddi.data());

			for (int rr = 0; rr < Rsq; rr++)
				sumC[rr] += sumCAddi[rr];
		} // j

		// now calculate the expected piece
		double Ai0 = 0.0;
		//double Ai1[R];
    vector<double> Ai1; Ai1.resize(R);
		//double Ai2[R][R];
    vector<vector<double> > Ai2; Ai2.resize(R); for(int r=0; r<R; r++) Ai2[r].resize(R);
		for (int k1 = 0; k1 < R; k1++) {
			Ai1[k1] = 0.0;
			for (int k2 = 0; k2 < R; k2++)
				Ai2[k1][k2] = 0.0;
		}

		for (unsigned int x = 0; x < genoDist.size(); x++) {
			// create all genotype permutations
			vector<int> genoToPerm;
			for (unsigned int j = 0; j < genoDist[x].size(); j++)
				for (int k = 0; k < genoDist[x][j]; k++)
					genoToPerm.push_back(j);
			vector<vector<int> > genoPerm;
			allPerms(genoToPerm, genoPerm);

			// calculate the test statistic
			double permWeight = (double) 1.0 / (double) genoPerm.size();
			for (unsigned int p = 0; p < genoPerm.size(); p++) {
				double sumAStar = 0.0;
				//double sumBStar[R];
				//double sumCStar[Rsq];
        vector<double> sumBStar; sumBStar.resize(R);
        vector<double> sumCStar; sumCStar.resize(Rsq);
				for (int r = 0; r < R; r++)
					sumBStar[r] = 0.0;
				for (int rr = 0; rr < Rsq; rr++)
					sumCStar[rr] = 0.0;

				for (unsigned int j = 0; j < genoPerm[p].size(); j++) {
					//double XStar[R];
					vector<double> XStar; XStar.resize(R);
					contsX(analyze_allele_index, analyze_allele_index_size,
							conditional_allele_index,
							conditional_allele_index_size, genoPerm[p][j],
							XStar);
					double bXStar = contsBetaX(b, (double *)XStar.data(), R);

					sumAStar += contsA(alpha, sigmaSq, b, (double *)XStar.data(), bXStar, j,
							analyze_allele_index, analyze_allele_index_size,
							conditional_allele_index,
							conditional_allele_index_size, ignoreBtX);

					contsB(alpha, sigmaSq, b, (double *)XStar.data(), bXStar, j,
							analyze_allele_index, analyze_allele_index_size,
							conditional_allele_index,
							conditional_allele_index_size, ignoreBtX, sumBAddi);
					for (int r = 0; r < R; r++)
						sumBStar[r] += sumBAddi[r];

					contsC(alpha, sigmaSq, b, (double *)XStar.data(), bXStar, j,
							analyze_allele_index, analyze_allele_index_size,
							conditional_allele_index,
							conditional_allele_index_size, ignoreBtX, (double *)sumCAddi.data());
					for (int rr = 0; rr < Rsq; rr++)
						sumCStar[rr] = +sumCAddi[rr];
				} // j

				double temp = exp(sumAStar) * genoDistP[x] * permWeight;

				Ai0 += temp;

				for (int r = 0; r < R; r++)
					Ai1[r] += sumBStar[r] * temp;

				for (int k1 = 0; k1 < R; k1++)
					for (int k2 = 0; k2 < R; k2++)
						Ai2[k1][k2] += (sumCStar[imcIndex(k1, k2, R)]
								+ sumBStar[k1] * sumBStar[k2]) * temp;
			} // p
		} // x

		// finally calculate the return piece
		double Ai0sq = Ai0 * Ai0;
		for (int k1 = 0; k1 < R; k1++)
			for (int k2 = 0; k2 < R; k2++)
				ret_I[imcIndex(k1, k2, R)] = (Ai2[k1][k2] * Ai0 - Ai1[k1]
						* Ai1[k2]) / Ai0sq;
	} // contsImc

	// a hash for the sufficient statistic group (index to allele, can only be done on an allele-wide basis)
	// DEBUG: Taken from previous code that was debugged, should be proper.
	int ssGroup( int index ) {
		//return( 0 );
		// Hash based on code from GxE work!
		int nG[3] = {0,0,0};
		for( unsigned int j=0; j<observed.size(); j++ ) {
			int gc = g[ observed[j] ].gCode( 0, index );
			if( gc != Genotype::gMiss )
				nG[ gc ]++;
		}
		//cout << "nG (" << nG[0] << "," << nG[1] << "," << nG[2] << ") " << toString() << endl; // DEBUG

		int gP1 = parents[0].gCode( 0, index );
		int gP2 = parents[1].gCode( 0, index );

		// gMiss=-1, so need to add 1 to parents code to hash correctly
		int pCode = (int)( (gP2+1)*1e6 + (gP1+1)*1e7 );
		if( gP1>gP2 )
			pCode = (int)( (gP1+1)*1e6 + (gP2+1)*1e7 ); // should have been the other way around
		if( gP2 != Genotype::gMiss && gP1 != Genotype::gMiss )
			return( pCode ); // coding only based on parents if both present
		return( (int)( nG[0]*1e0 + nG[1]*1e2 + nG[2]*1e4 ) + pCode ); // otherwise need coding of offspring and parent in here
	}
};

// DEBUG: Mostly should be good, Data:: LinkTrait is a little iffy in how it works -- what happens when people have missing traits?
class Data {
public:
	vector<Pedigree> ped;

	// create from FBAT output
	void create(string &filename) {
		// Read in the output
		Lines lines;
		lines.load(filename);

		// Find each pedigree
		vector<int> pedigreeLocs;
		lines.find("ped", pedigreeLocs);

		// size up the pedigrees
		ped.resize(pedigreeLocs.size());

		StrTok tok;
		for (unsigned int i = 0; i < pedigreeLocs.size(); i++) {
			ped[i].label = lines[pedigreeLocs[i]];
			tok.tokenize(ped[i].label);
			tok.nextTokenE();
			ped[i].pid = (unsigned int) tok.nextTokenN();

			//cout << lines[pedigreeLocs[i]] << endl;
			if (i == pedigreeLocs.size() - 1) {
				ped[i].parse(lines, pedigreeLocs[i], lines.size() - 1);
			} else {
				ped[i].parse(lines, pedigreeLocs[i], pedigreeLocs[i + 1] - 1);
			}
		} // i
	} // create

	// clear the object
	void clear() {
		ped.resize(0);
	} // clear

	// convert to a string
	string toString() {
		string s;
		for (unsigned int i = 0; i < ped.size(); i++) {
			if (ped[i].informative())
				s += ped[i].toString() + "\n";
		}
		return (s);
	} // toString

	// link the trait to the pedigree individual
	// DEBUG: MOST CONCERNING ROUTINE -- what happens when a trait is missing?
	void linkTrait(int *pid, double *trait, int n) {
		// somewhat inefficient code to do this... but it isn't what is really taking the most time (i.e. other optimizations would probably be better)
		for( unsigned int iped = 0; iped < ped.size(); iped++) {
			for( unsigned int ipid = 0; (int)ipid < n; ipid++)
				if ( (int)ped[iped].pid == pid[ipid]) // && !isnan(trait[ipid]) ) // added second piece 7/28/2008
					ped[iped].trait.push_back(trait[ipid]);
		} // iped

		// now make sure that the offspring and traits match
		// - actually, it's really only a problem if the observed were bigger than the trait,
		//   otherwise it doesn't really matter.
		for (unsigned int p = 0; p < ped.size(); p++) {
			if (ped[p].observed.size() != ped[p].trait.size()
					&& ped[p].observed.size() != 0) {
				//if( ped[p].observed.size() > ped[p].trait.size() ) {
        Rprintf("data::linkTrait::observed.size()(%d) != trait.size()(%d) for pedigree %d\n", ped[p].observed.size(), ped[p].trait.size(), ped[p].pid);
				Rprintf("%s\n", ped[p].toString().c_str());
			}
		} // p
	} // linkTrait

	// centering the trait
	void centerTrait(double* center, bool mean) {
		if (mean) {
			// _replaces_ center with mean
			*center = 0.0;
			int count = 0;

			for (unsigned int p = 0; p < ped.size(); p++) {
				for (unsigned int j = 0; j < ped[p].trait.size(); j++) {
                                        if (!isnan(ped[p].trait[j])) {
//                                        if ( !isnan(ped[p].trait[j]) && ped[p].g[j].genotype(0,0,1,1)==1 ) {
						*center += ped[p].trait[j];
						count++;
					}
				} // j
			} // p
			*center /= count;
		} // if(mean)

		// now suptract off the correction
		for (unsigned int p = 0; p < ped.size(); p++)
			for (unsigned int j = 0; j < ped[p].trait.size(); j++)
				ped[p].trait[j] -= *center; // nan's are still nan...
	} // centerTrait

	double proportionInferrable() {
		double inferrable = 0.0;
		for (unsigned int p = 0; p < ped.size(); p++)
			inferrable += (int) ped[p].phaseInferrable();
		inferrable /= (double) ped.size();
		return (inferrable);
	} // proportionInferrable

	void removeUnphased() {
		// removes pedigrees where phase cannot be determined
		for (unsigned int p = 0; p < ped.size(); p++) {
			if (!ped[p].phaseInferrable()) {
				// remove the pedigree
				ped.erase(ped.begin() + p);

				// and decrement
				p--;
			}
		}
	} // removeUnphased

	// each allele goes in a seperate "column" after being dealt with
	//  using the matrix command in R
	void uimc(double *bm, double *bc0, double *bc1, int *analyze_allele_index,
			int analyze_allele_index_size, int *conditional_allele_index,
			int conditional_allele_index_size, bool onlyComputeConditional,
			double *ret_analyze, double *ret_conditional0,
			double *ret_conditional1) {
		//double temp_ret_analyze[analyze_allele_index_size];
		//double temp_ret_conditional0[conditional_allele_index_size];
		//double temp_ret_conditional1[conditional_allele_index_size];
		vector<double> temp_ret_analyze; temp_ret_analyze.resize(analyze_allele_index_size);
    vector<double> temp_ret_conditional0; temp_ret_conditional0.resize(conditional_allele_index_size);
    vector<double> temp_ret_conditional1; temp_ret_conditional1.resize(conditional_allele_index_size);

		unsigned int P = ped.size();

		for (unsigned int p = 0; p < ped.size(); p++) {
			ped[p].uimc(bm, bc0, bc1, analyze_allele_index,
					analyze_allele_index_size, conditional_allele_index,
					conditional_allele_index_size, onlyComputeConditional,
					(double *)temp_ret_analyze.data(), (double *)temp_ret_conditional0.data(),
					(double *)temp_ret_conditional1.data());
			for (int a = 0; a < analyze_allele_index_size; a++)
				ret_analyze[p + a * P] = temp_ret_analyze[a];

			for (int a = 0; a < conditional_allele_index_size; a++) {
				ret_conditional0[p + a * P] = temp_ret_conditional0[a];
				ret_conditional1[p + a * P] = temp_ret_conditional1[a];
			} // a
		} // p
	} // uimc

	void imc(double *bm, double *bc0, double *bc1, int *analyze_allele_index,
			int analyze_allele_index_size, int *conditional_allele_index,
			int conditional_allele_index_size, double *I) {
		int R = analyze_allele_index_size + 2 * conditional_allele_index_size;
		int Rsq = R * R;
		//double Iplus[Rsq];
    vector<double> Iplus; Iplus.resize(Rsq);
		for (unsigned int p = 0; p < ped.size(); p++) {
			ped[p].imc(bm, bc0, bc1, analyze_allele_index,
					analyze_allele_index_size, conditional_allele_index,
					conditional_allele_index_size, (double *)Iplus.data());

			for (int r = 0; r < Rsq; r++)
				if (!isnan(Iplus[r]))
					I[r] += Iplus[r];
		}
	} // imc

	void robustStat(int *analyze_allele_index, int analyze_allele_index_size,
			int *conditional_allele_index, int conditional_allele_index_size,
			double *ret_analyze) {
		//double temp_ret_analyze[analyze_allele_index_size];
		vector<double> temp_ret_analyze; temp_ret_analyze.resize(analyze_allele_index_size);
		unsigned int P = ped.size();
		for (unsigned int p = 0; p < ped.size(); p++) {
			ped[p].robustStat(analyze_allele_index, analyze_allele_index_size,
					conditional_allele_index, conditional_allele_index_size,
					(double *)temp_ret_analyze.data());

			for (int a = 0; a < analyze_allele_index_size; a++)
				ret_analyze[p + a * P] = temp_ret_analyze[a];
		}
	} // robustStat

	void contsUimc(double alpha,
			double sigmaSq,
			double *b, // bm, bc0, bc1
			int *analyze_allele_index, int analyze_allele_index_size,
			int *conditional_allele_index, int conditional_allele_index_size,
			bool onlyComputeConditional, bool ignoreBtX, double *ret_b) {
		int R = analyze_allele_index_size + conditional_allele_index_size * 2;
		//double temp_ret_b[R];
    vector<double> temp_ret_b; temp_ret_b.resize(R);

		unsigned int P = ped.size();

		for (unsigned int p = 0; p < ped.size(); p++) {
			ped[p].contsUimc(alpha, sigmaSq, b, analyze_allele_index,
					analyze_allele_index_size, conditional_allele_index,
					conditional_allele_index_size, onlyComputeConditional,
					ignoreBtX, (double *)temp_ret_b.data());
			for (int r = 0; r < R; r++)
				ret_b[p + r * P] = temp_ret_b[r];
		} // p
	} // contsUimc

	void contsImc(double alpha,
			double sigmaSq,
			double *b, // bm, bc0, bc1
			int *analyze_allele_index, int analyze_allele_index_size,
			int *conditional_allele_index, int conditional_allele_index_size,
			bool ignoreBtX, double *ret_I) {
		//
		int R = analyze_allele_index_size + conditional_allele_index_size * 2;
		int Rsq = R * R;
		/////unsigned int P = ped.size();
		//double temp_ret_I[Rsq];
    vector<double> temp_ret_I; temp_ret_I.resize(Rsq);

		for (int rr = 0; rr < Rsq; rr++)
			ret_I[rr] = 0.0;

		for (unsigned int p = 0; p < ped.size(); p++) {
			ped[p].contsImc(alpha, sigmaSq, b, analyze_allele_index,
					analyze_allele_index_size, conditional_allele_index,
					conditional_allele_index_size, ignoreBtX, (double *)temp_ret_I.data());

			for (int rr = 0; rr < Rsq; rr++)
				if (!isnan(temp_ret_I[rr]))
					ret_I[rr] += temp_ret_I[rr];
		}
	} // contsImc

	int numInfFam() {
		int num = 0;
		for (unsigned int p = 0; p < ped.size(); p++) {
			if (ped[p].observed.size() > 0)
				num++; // then they are informative
		} // p
		return (num);
	} // numInf

	void pids(int *pid) {
		for (unsigned int p = 0; p < ped.size(); p++)
			pid[p] = ped[p].pid;
	}
};

/** Code related to talking / exporting to R **/
vector<Data> data; // vector of data stored in the shared library, not in R
vector<bool> dataUsed; // whether the current vector slot is in use
// allocate by extending data, or using an empty slot
unsigned int dataAllocate() {
	for (unsigned int i = 0; i < dataUsed.size(); i++) {
		if (dataUsed[i] == false) {
			dataUsed[i] = true;
			return (i);
		}
	}
	data.resize(data.size() + 1);
	dataUsed.push_back(true);
	return (data.size() - 1);
}
// free up the data
void dataFree(int dataIndex) {
	data[dataIndex].clear();
	dataUsed[dataIndex] = false;
}

// EXPORTED FUNCTION TO R
extern "C" {
// loads and stores the output in memory here
// -- reference is a returned value to index in for the rest of the routines
void condGeneFBATControl_load(char** filename, int *reference) {
	int dataIndex = dataAllocate();
	string strFilename = *filename;
	data[dataIndex].create(strFilename);

	*reference = dataIndex;
}

// frees memory allocated here
void condGeneFBATControl_free(int *reference) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_free::Reference %d no longer exists.\n", (*reference));
    return;
	}

	dataFree(*reference);
}

// printing (primarily for debugging)
void condGeneFBATControl_print(int *reference) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_print::Reference %d no longer exists.\n", (*reference));
    return;
	}
	Rprintf("%s\n", data[*reference].toString().c_str());
}

// linking in the trait
void condGeneFBATControl_linkTrait(int *reference, int *pid, double *trait,
		int *n) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_free::linkTrait %d no longer exists.\n", (*reference));
    return;
	}
	data[*reference].linkTrait(pid, trait, *n);
}

// informative
void condGeneFBATControl_proportionInformative(int *reference, double *ret) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_proportionInformative::Reference %d no longer exists.\n", (*reference));
    return;
	}

	*ret = data[*reference].proportionInferrable();
}

// removing unphased individuals
void condGeneFBATControl_removeUnphased(int *reference) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_removeUnphased::Reference %d no longer exists.\n", (*reference));
    return;
	}

	data[*reference].removeUnphased();
}

// how many families are there?
void condGeneFBATControl_numFam(int *reference, int *numFam) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_numFam::Reference %d no longer exists.\n", (*reference));
    return;
	}

	*numFam = data[*reference].ped.size();
}

// centering the trait
void condGeneFBATControl_centerTrait(int *reference, double *center, int *mean) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_centerTrait::Reference %d no longer exists.\n", (*reference));
    return;
	}

	data[*reference].centerTrait(center, (bool) ((*mean) == 1));
}

void condGeneFBATControl_uimc(int *reference, double *bm, double *bc0,
		double *bc1, int *analyze_allele_index, int *analyze_allele_index_size,
		int *conditional_allele_index, int *conditional_allele_index_size,
		int *onlyComputeConditional, double *ret_analyze,
		double *ret_conditional0, double *ret_conditional1) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_uimc %d no longer exists.\n", (*reference));
    return;
	}

	data[*reference].uimc(bm, bc0, bc1, analyze_allele_index,
			*analyze_allele_index_size, conditional_allele_index,
			*conditional_allele_index_size, ((*onlyComputeConditional) != 0),
			ret_analyze, ret_conditional0, ret_conditional1);
}

void condGeneFBATControl_imc(int *reference, double *bm, double *bc0,
		double *bc1, int *analyze_allele_index, int *analyze_allele_index_size,
		int *conditional_allele_index, int *conditional_allele_index_size,
		double *ret_I) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_imc %d no longer exists.\n", (*reference));
    return;
	}

	data[*reference].imc(bm, bc0, bc1, analyze_allele_index,
			*analyze_allele_index_size, conditional_allele_index,
			*conditional_allele_index_size, ret_I);
}

void condGeneFBATControl_robustStat(int *reference, int *analyze_allele_index,
		int *analyze_allele_index_size, int *conditional_allele_index,
		int *conditional_allele_index_size, double *ret_analyze) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_robustStat::Reference %d no longer exists.\n", (*reference));
    return;
	}

	data[*reference].robustStat(analyze_allele_index,
			*analyze_allele_index_size, conditional_allele_index,
			*conditional_allele_index_size, ret_analyze);
}

void condGeneFBATControl_contsUimc(int *reference, double *alpha,
		double *sigma, double *b, int *analyze_allele_index,
		int *analyze_allele_index_size, int *conditional_allele_index,
		int *conditional_allele_index_size, int *onlyComputeConditional,
		int *ignoreBtX, double *ret_b) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_contsUimc::Reference %d no longer exists.\n", (*reference));
    return;
	}

	/*
	 int R = (*analyze_allele_index_size) + (*conditional_allele_index_size) * 2;
	 cout << "beta ";
	 for( int r=0; r<R; r++ )
	 cout << b[r] << " ";
	 cout << endl;
	 */

	data[*reference].contsUimc(*alpha, (*sigma) * (*sigma), b,
			analyze_allele_index, *analyze_allele_index_size,
			conditional_allele_index, *conditional_allele_index_size,
			((*onlyComputeConditional) != 0), ((*ignoreBtX) != 0), ret_b);
} // condGeneFBATControl_contsUimc

void condGeneFBATControl_contsImc(int *reference, double *alpha, double *sigma,
		double *b, int *analyze_allele_index, int *analyze_allele_index_size,
		int *conditional_allele_index, int *conditional_allele_index_size,
		int *ignoreBtX, double *ret_I) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_contsImc::Reference %d no longer exists.\n", (*reference));
    return;
	}

	data[*reference].contsImc(*alpha, (*sigma) * (*sigma), b,
			analyze_allele_index, *analyze_allele_index_size,
			conditional_allele_index, *conditional_allele_index_size,
			((*ignoreBtX) != 0), ret_I);
} // condGeneFBATControl_contsImc


void condGeneFBATControl_estEqNuis(
		int *referenceCondition, int *referenceConditionSize,
		double *offset,     // if nonzero, indicates dichotomous
		double *ret_lhs, // matrix of size (referenceCondition*2)^2
		double *ret_rhs ) {
	// Make sure the references are good
	for( int c=0; c<*referenceConditionSize; c++ ) {
		if( referenceCondition[c] < 0 || referenceCondition[c] >= (int) data.size()) {
      Rprintf("condGeneFBATControl_estEqNuis::Reference %d no longer exists.\n", referenceCondition[c]);
      return;
		}
	}

	// is it a quantitative trait?
	bool qtl = (*offset == 0);

	// Notes:
	// - the allele index will always be zero

	// make some constants
	int nc = *referenceConditionSize;
	int nc2 = 2 * nc;
	int np = data[ referenceCondition[0] ].ped.size();

	// zero out the return piece (actually this is done in the R code, so not really necessary)
	for( int c=0; c<nc2*nc2; c++ )
		ret_lhs[c] = 0.0;
	for( int c=0; c<nc2; c++ )
		ret_rhs[c] = 0.0;

	// across each pedigree
	for( int i=0; i<np; i++ ) {
		//double xijc[nc2], exijc[nc2], eij=0;
		vector<double> xijc; xijc.resize(nc2);
    vector<double> exijc; exijc.resize(nc2);
    double eij = 0;

		// need to know the number of offspring
		unsigned int noff = data[referenceCondition[0]].ped[i].observed.size();
		for( int c=1; c<nc; c++)
			if( data[referenceCondition[c]].ped[i].observed.size() < noff )
				noff = data[referenceCondition[c]].ped[i].observed.size();

		// go across each of the offspring
		for( unsigned int j=0; j<noff; j++ ) {
			// go across all of the conditioning alleles
			bool traitSet = false;
			eij = 0.0;

			// OK, this was a horrible choice of words... I don't really mean "informative" as in X-E(X|S)!=0
			//  here. Instead, what we really mean is that this marker doesn't have any missing data
			//  associated with it...
			// In fact, we'd like to have a separate way of indicating informative...
			// It might be best to find all the "informative" and replace them with "nonmissing" as the varname
			bool informative = true; // this is just all bad -- really just X-E(X|S)=0 for the current piece

			for( int c=0; c<nc; c++ ) {
				// Some pedigrees might be noninformative, so no contribution for a particular allele
				Pedigree *ped = &data[referenceCondition[c]].ped[i];

				// 01.20.2009
				if( ped->observed.size() != ped->nonzeroDelX.size() ) {
					ped->nonzeroDelX.resize( ped->observed.size() );
					for( unsigned int jj=0; jj<ped->nonzeroDelX.size(); jj++ )
						ped->nonzeroDelX[jj] = false;
				}

				////if( ped[i].observed.size() < 0 ) {
				//if( ped->observed.size() <= 0 ) {
				if( j >= ped->observed.size() ) {
					// noninformative, set everything to be zero
					xijc[c*2+0] = xijc[c*2+1] = 0.0;
					exijc[c*2+0] = exijc[c*2+1] = 0.0;
					informative = false;
				}else{
					// compute the observed
					xijc[c*2+0] = ped->g[ ped->observed[j] ].genotype( 0, 0, 2, 2 );
					xijc[c*2+1] = ped->g[ ped->observed[j] ].genotype( 0, 0, 1, 2 );

					// compute the expected [[yes, this really shouldn't be necessary for all offspring, inefficient]]
					exijc[c*2+0] = exijc[c*2+1] = 0.0;
					for( unsigned int gg=0; gg<ped->g.size(); gg++ ) {
						exijc[c*2+0] += ped->g[gg].genotype( 0, 0, 2, 2 ) * ped->pg[gg];
						exijc[c*2+1] += ped->g[gg].genotype( 0, 0, 1, 2 ) * ped->pg[gg];
					}// gg

					// is the trait set?
					if( !traitSet && !isnan(ped->trait[j]) ) {
						eij = ped->trait[ j ];
						//if( !qtl ) eij *= exp( - (*offset) );
						traitSet = true;
					}//fi( !traitSet )
				}//fi( ped[i].observed.size() < 0 )
			}// c

			if( traitSet ) { //& informative ) { // make sure someone had a contribution -- not really necessary
				bool nonzeroDelX = false; // 01.20.2009
				// go across all of the conditioning alleles again
				for( int c1=0; c1<nc2; c1++ ) {
					double delxijc1 = xijc[c1] - exijc[c1];
					nonzeroDelX = nonzeroDelX || delxijc1!=0; // 01.20.2009
					ret_rhs[c1] += delxijc1 * eij;
					for( int c2=0; c2<nc2; c2++ ) {
						//double delxijc2 = xijc[c2] - exijc[c2];
						//ret_lhs[c1 + c2*nc2] += delxijc1*xijc[c1] * delxijc2*xijc[c2];

						// m11 m21 m31 ... m21 m22 m32 ... (order it should be strung out so R puts it in a matrix properly)
						ret_lhs[c1 + c2*nc2] += delxijc1*xijc[c2];
						//ret_lhs[c1*nc2 + c2] += delxijc1*xijc[c2];
					}// c2
				}// c1

				// 01.20.2009 update
				for( int c=0; c<nc; c++ ) {
					Pedigree *ped = &data[referenceCondition[c]].ped[i];
					ped->nonzeroDelX[j] = nonzeroDelX;
				}

				//if( !nonzeroDelX )
				//	cout << "NONINFORMATIVE (" << i << "," << j << ")" << endl;
			}//fi( traitSet )
		}// j
	}// i
}

// for updating Y until convergence is reached, no return value
// Does it based on the mean of the residuals (does NOT adjust for strata)
void condGeneFBATControl_estEqNuisUpdate(
		int *referenceCondition, int *referenceConditionSize,
		double *bc ) {
	// Make sure the references are good
	for( int c=0; c<*referenceConditionSize; c++ ) {
		if( referenceCondition[c] < 0 || referenceCondition[c] >= (int) data.size()) {
      Rprintf("condGeneFBATControl_estEqNuis %d no longer exists.\n", referenceCondition[c]);
      return;
		}
	}

	// is it a quantitative trait?
	// really has to be if in this routine!

	// Notes:
	// - the allele index will always be zero

	// make some constants
	int nc = *referenceConditionSize;
	int nc2 = 2 * nc;
	int np = data[ referenceCondition[0] ].ped.size();

	// across each pedigree
	double eijSum = 0.0; int eijCount = 0;
	for( int i=0; i<np; i++ ) {
		//double xijc[nc2], eij=0;
		vector<double> xijc; xijc.resize(nc2);
    double eij = 0;

		// need to know the number of offspring
		unsigned int noff = data[referenceCondition[0]].ped[i].observed.size();
		for( int c=1; c<nc; c++)
			if( data[referenceCondition[c]].ped[i].observed.size() < noff )
				noff = data[referenceCondition[c]].ped[i].observed.size();

		// go across each of the offspring
		for( unsigned int j=0; j<noff; j++ ) {
			// go across all of the conditioning alleles
			bool traitSet = false;
			eij = 0.0;
			bool informative = true;
			for( int c=0; c<nc; c++ ) {
				// Some pedigrees might be noninformative, so no contribution for a particular allele
				Pedigree *ped = &data[referenceCondition[c]].ped[i];
				////if( ped[i].observed.size() < 0 ) {
				//if( ped->observed.size() <= 0 ) {
				if( j>=ped->observed.size() ) {
					// noninformative, set everything to be zero
					informative = false;
				//}else if( !ped->nonzeroDelX[j] ) { // 01.20.2009
				//	informative = false;
				}else{
					// compute the observed
					xijc[c*2+0] = ped->g[ ped->observed[j] ].genotype( 0, 0, 2, 2 );
					xijc[c*2+1] = ped->g[ ped->observed[j] ].genotype( 0, 0, 1, 2 );

					// is the trait set?
					if( !traitSet && !isnan(ped->trait[j]) ) {
						eij = ped->trait[ j ];
						traitSet = true;
					}//fi( !traitSet )
				}//fi( ped[i].observed.size() < 0 )
			}// c

			if( traitSet ) {//& informative ) {
				double bx = 0.0;
				for( int c1=0; c1<nc2; c1++ )
					bx += bc[c1] * xijc[c1];
				eij -= bx;
				eijSum += eij;
				eijCount++;
			}//fi( traitSet )
		}// j
	}// i

	// compute mean of eij
	eijSum /= eijCount;

	// subtract mean of eij (eijSum) from every trait
	// NOTE: this is only doing the traits on these alleles
	// WARNING: TODO: OH HELL: Make SURE that we search the conditioning alleles for the trait first!!!
	for( int i=0; i<np; i++ ) {
		for( int c=0; c<nc; c++ )
			for( unsigned int j=0; j<data[referenceCondition[c]].ped[i].trait.size(); j++ )
				data[referenceCondition[c]].ped[i].trait[j] -= eijSum;
	}
}

class SSBucketMember{
public:
  int hash;
  vector<unsigned int> memberPedIndexI;
  vector<unsigned int> memberPedIndexJ;

  unsigned int size() {
    // if( memberPedIndexI.size() != memberPedIndexJ.size() ) cout << "Help me..." << endl;
    return( memberPedIndexI.size() );
  }//size

  void add( unsigned int i, unsigned int j ) {
    memberPedIndexI.push_back( i );
    memberPedIndexJ.push_back( j );
  }//add

  void clear() {
    memberPedIndexI.clear();
    memberPedIndexJ.clear();
  }//clear

  // Helper Function - Converts a double to a string
  string doubleToString(double d) {
    //string str;
    ostringstream oss;
    oss << d;
    return (oss.str());
  }

  string toString() {
    string str;
    str = doubleToString(hash);
    str += ": ";
    for( unsigned int k=0; k<memberPedIndexI.size(); k++ ) {
      str += "(";
      str += doubleToString(memberPedIndexI[k]);
      str += ",";
      str += doubleToString(memberPedIndexJ[k]);
      str += ") ";
    }
    return( str );
  }//toString
};//SSBucketMember

class SSBucket{
public:
  vector<SSBucketMember> bucket;
  SSBucketMember empty;

  void clear() {
    bucket.clear();
  }// clear

  unsigned int size() {
    return( bucket.size() );
  }// size

  string toString( const char* bucketName=NULL ) {
    string s = "SSBucket";
    if( bucketName!=NULL ) {
      s += "("; s+= bucketName; s+=")";
    }
    s += "\n";

    for( unsigned int b=0; b<bucket.size(); b++ ){
      s += " ";
      s += bucket[b].toString();
      s += "\n";
    }

    return( s );
  }

  SSBucketMember &operator []( unsigned int index ){
    if( index<0 || index>=size() ) {
      Rprintf("Bucket member %d is out of bounds [0,%d]\n", index, (size()-1));
      return( empty );
    }

    return( bucket[index] );
  }// []

  // returns index to hash, creating it if necessary
  unsigned int get( int hash ) {
    for( unsigned int b=0; b<size(); b++ )
      if( bucket[b].hash == hash )
        return(b);
    // couldn't be found, create it
    SSBucketMember newBucket;
    newBucket.hash = hash;
    bucket.push_back( newBucket );
    return( bucket.size() - 1 );
  }// get

  // for adding an entry
  void add( int hash, unsigned int pedIndexI, unsigned int pedIndexJ ) {
    unsigned int index = get( hash );
    bucket[index].add( pedIndexI, pedIndexJ );
  }// add

  // Reduces the number of hashes, only call _after finalized_
  void reduce() {
    for( unsigned int b=0; b<bucket.size(); b++ )
      bucket[b].hash = b;
  }// reduceHash

  // Will call reduce(), so needs to be _finalized_
  SSBucket merge( SSBucket &ssb ) {
    // The output
    SSBucket mergeBucket;

    // reduce hashes
    reduce();
    ssb.reduce();

    // multiplier
    int multiplier = size() + 2; // 2 is for safety, 1 should be sufficient, maybe 0 works
    if( ssb.size() > size() )
      multiplier = ssb.size() + 2;

    for( unsigned int b=0; b<bucket.size(); b++ ) {
      for( unsigned int k=0; k<bucket[b].size(); k++ ) {
        // the entry we are on
        unsigned int i = bucket[b].memberPedIndexI[k];
        unsigned int j = bucket[b].memberPedIndexJ[k];
        // find it in ssb
        int ssb_hash = -1;
        for( unsigned int ssb_b=0; ssb_b<ssb.bucket.size() && ssb_hash==-1; ssb_b++ ) {
          for( unsigned int ssb_k=0; ssb_k<ssb.bucket[ssb_b].size() && ssb_hash==-1; ssb_k++ ) {
            if( i == ssb.bucket[ssb_b].memberPedIndexI[ssb_k] && j == ssb.bucket[ssb_b].memberPedIndexJ[ssb_k] )
              ssb_hash = ssb.bucket[ssb_b].hash;
          }
        }
        // Err if there is someone in one bucket and not in the other
        if( ssb_hash == -1 ) {
          Rprintf("SSBucket::mergeHash:: The individual (%d,%d) is in one SSBucket, but not in the other!\n", i, j);
        }
        // and add it to the new bucket
        int newHash = ssb_hash + multiplier * bucket[b].hash;
        mergeBucket.add( newHash, i, j );
      }// k
    }// b

    // now, reduce again
    reduce();

    // and return the new SSBucket
    return( mergeBucket );
  }// mergeHash
};// SSBucket

// for updating Y until convergence is reached, no return value
// This routine adjusts for the strata of sufficient statistics
void condGeneFBATControl_estEqNuisUpdate2(
		int *referenceCondition, int *referenceConditionSize,
		double *bc ) {
	// Make sure the references are good
	for( int c=0; c<*referenceConditionSize; c++ ) {
		if( referenceCondition[c] < 0 || referenceCondition[c] >= (int) data.size()) {
      Rprintf("condGeneFBATControl_free::Reference %d no longer exists.\n", referenceCondition[c]);
      return;
		}
	}

	// a few constants
	int nc = *referenceConditionSize;
	int nc2 = nc*2;
	int np = data[ referenceCondition[0] ].ped.size();

	//cout << "about to create bucket" << endl;
	// Create an SSBucket of all of the pedigree members for _each_ allele (then merge them)
	//SSBucket bucket[nc];
  vector<SSBucket> bucket; bucket.resize(nc);
	for( int b=0; b<nc; b++ ) {
		//cout << "b=" << b << ";";
		for( int i=0; i<np; i++ ) {
			//cout << " " << i;
			//////int c = b;

			// get the number of offspring
			unsigned int noff = data[referenceCondition[0]].ped[i].observed.size();
			for( int ck=1; ck<nc; ck++ )
				if( data[referenceCondition[ck]].ped[i].observed.size() < noff )
					noff = data[referenceCondition[ck]].ped[i].observed.size();
			//cout << "noff = " << noff << endl;

			if( noff > 0 ) {
				// compute the group
				//cout << "b=" << b << " i=" << i << endl;
				//cout << "data.size() = " << data.size() << endl;
				//cout << "data[referenceCondition[b]].ped.size()" << data[referenceCondition[b]].ped.size() << endl;

				int ss = data[referenceCondition[b]].ped[i].ssGroup(0); // always will be index 0
				//cout << "ss = " << ss << endl;

				// go across the offspring, and add them to the hash table
				for( unsigned int j=0; j<noff; j++ )
					bucket[b].add( ss, i, j );
			}
		}
		//cout << endl;
	}
	//cout << "bucket created" << endl;

	// Merge the SSBuckets
	SSBucket mbucket = bucket[0];
	for( int b=1; b<nc; b++ )
		mbucket = mbucket.merge( bucket[b] );
	//cout << "bucket merged" << endl;

	//cout << mbucket.toString() << endl;

	//cout << "mbucket.size() " << mbucket.size(); // the number of strata...

	// Now for the merged bucket, compute the mean of the eij, and then subtract it off of the traits
	for( unsigned int h=0; h<mbucket.size(); h++ ) {
		double eijSum = 0.0;
		double eijCount = 0;

		for( unsigned int k=0; k<mbucket[h].size(); k++ ) {
			int i = mbucket[h].memberPedIndexI[k];
			int j = mbucket[h].memberPedIndexJ[k];

			// go across all of the conditioning alleles
			bool traitSet = false;
			bool informative = true;
			//double xijc[nc2], eij = 0.0;
      vector<double> xijc; xijc.resize(nc2);
      double eij = 0.0;
			for( int c=0; c<nc; c++ ) {
				// Some pedigrees might be noninformative, so no contribution for a particular allele
				Pedigree *ped = &data[referenceCondition[c]].ped[i];
				////if( ped[i].observed.size() < 0 ) {
				//if( ped->observed.size() <= 0 ) {
				if( j >= (int)ped->observed.size() ) {
					// noninformative, set everything to be zero
					informative = false;
				//}else if( !ped->nonzeroDelX[j] ) { // 01.20.2009
				//	informative = false;
				}else{
					// compute the observed
					xijc[c*2+0] = ped->g[ ped->observed[j] ].genotype( 0, 0, 2, 2 );
					xijc[c*2+1] = ped->g[ ped->observed[j] ].genotype( 0, 0, 1, 2 );

					// is the trait set?
					if( !traitSet && !isnan(ped->trait[j]) ) {
						eij = ped->trait[ j ];
						traitSet = true;
					}//fi( !traitSet )
				}//fi( ped[i].observed.size() < 0 )
			}

			if( traitSet ){//& informative ) {
				double bx = 0.0;
				for( int c1=0; c1<nc2; c1++ )
					bx += bc[c1] * xijc[c1];
				eij -= bx;
				eijSum += eij;
				eijCount++;
				//cout << "eij " << eij << endl;
			}//fi( traitSet )
		}

		eijSum /= eijCount;
		//cout << "eijSum (mean) = " << eijSum << endl;

		// and subtract the mean from every trait
		for( unsigned int k=0; k<mbucket[h].size(); k++ ) {
			int i = mbucket[h].memberPedIndexI[k];
			int j = mbucket[h].memberPedIndexJ[k];

			for( int c=0; c<nc; c++ ) {
				if( j < (int)data[referenceCondition[c]].ped[i].trait.size() ) {
					/*
					cout << "Trait (" << i << "," << j << "," << c << ") before=" <<
						data[referenceCondition[c]].ped[i].trait[j] <<
						" after=" <<
						(data[referenceCondition[c]].ped[i].trait[j] - eijSum) <<
						endl;
						*/

					data[referenceCondition[c]].ped[i].trait[j] -= eijSum;
				}
			}
		}// k
	}// h
}

void condGeneFBATControl_estEq(
		int *referenceAnalyze, int *referenceAnalyzeSize,
		int *referenceCondition, int *referenceConditionSize,
		double *bc,
		double *offset, // if nonzero, indicates dichotomous
		double *ret_uij, // ret_uijm, ret_uijc
		double *ret_xmxc, double *ret_xcxc ){
	// is it a quantitative trait
	bool qtl = (*offset==0);

	// make some constants
	int na = *referenceAnalyzeSize;
	int nc = *referenceConditionSize;
	int nc2 = nc * 2;
	int np = data[ referenceCondition[0] ].ped.size();

	//cout << "na = " << na << " nc=" << nc << " nc2=" << nc2 << " np=" << np << endl;
	//for( int a=0; a<na; a++ )
	//	cout << " referenceAnalyze[" << a << "] = " << referenceAnalyze[a] << endl;

	// make sure the references are good
	for( int a=0; a<na; a++ ) {
		if( referenceAnalyze[a] < 0 || referenceAnalyze[a] >= (int)data.size() ) {
			Rprintf("condGeneFBATControl_estEq %d no longer exists.\n", referenceAnalyze[a]);
			return;
		}
	}
	for( int c=0; c<nc; c++ ) {
		if( referenceCondition[c] < 0 || referenceCondition[c] >= (int) data.size()) {
			Rprintf("condGeneFBATControl_estEq %d no longer exists.\n", referenceCondition[c]);
			return;
		}
	}

	// zero out the return
	for( int ij=0; ij<np*(na+nc2); ij++ )
		ret_uij[ij] = 0.0;
	for( int ij=0; ij<na*nc2; ij++ ) // new additions...
		ret_xmxc[ij] = 0.0;
	for( int ij=0; ij<nc2*nc2; ij++ )
		ret_xcxc[ij] = 0.0;

	// across each pedigree
	for( int i=0; i<np; i++ ) {
		// pieces we need to compute (for each offspring)
		//double xija[na], exija[na],  xijc[nc2], exijc[nc2],  eij=0.0;
		vector<double> xija; xija.resize(na);
    vector<double> exija; exija.resize(na);
    vector<double> xijc; xijc.resize(nc2);
    vector<double> exijc; exijc.resize(nc2);
    double eij = 0.0;

		// need to know the number of offspring? why?
		unsigned int noff = data[referenceCondition[0]].ped[i].observed.size();
		for( int c=1; c<nc; c++ )
			if( data[referenceCondition[c]].ped[i].observed.size() > noff )
				noff = data[referenceCondition[c]].ped[i].observed.size();
		for( int a=0; a<na; a++ )
			if( data[referenceAnalyze[a]].ped[i].observed.size() > noff )
				noff = data[referenceAnalyze[a]].ped[i].observed.size();

		// go across each of the offspring
		for( unsigned int j=0; j<noff; j++ ) {
			bool traitSet = false;
			eij = 0.0;

			bool informative = true;
			// go across the analyze alleles
			for( int a=0; a<na; a++ ) {
				Pedigree *ped = &data[referenceAnalyze[a]].ped[i];
				//if( ped->observed.size() <= 0 ) {
				if( j >= ped->observed.size() ) {
					// noninformative
					xija[a] = exija[a] = 0.0;
					informative = false;
				}else{
					// compute the observed
					xija[a] = ped->g[ ped->observed[j] ].xCode( 0, 0, 2, ADDITIVE );

					// compute the expected
					exija[a] = 0.0;
					for( unsigned int gg=0; gg<ped->g.size(); gg++ )
						exija[a] += ped->g[gg].xCode( 0, 0, 2, ADDITIVE ) * ped->pg[gg];

          // NO -- this should _NOT_ be done, as it's only adjusted for the conditioning alleles
					// set the trait?
					if( !traitSet && !isnan(ped->trait[j]) ) {
						eij = ped->trait[j];
						traitSet = true;
						if(!qtl)
              Rprintf("DEATH KNELL -- DICHOTOMOUS IS NOT YET SUPPORTED!\n");
					}
				}//fi( ped->observed.size() <= 0 )
			}// a

			// go across the condition alleles
			for( int c=0; c<nc; c++ ) {
				// Some pedigrees might be noninformative, so no contribution for a particular allele
				Pedigree *ped = &data[referenceCondition[c]].ped[i];
				////if( ped[i].observed.size() < 0 ) {
				//if( ped->observed.size() <= 0 ) {
				if( j >= ped->observed.size() ) {
					// noninformative, set everything to be zero
					xijc[c*2+0] = xijc[c*2+1] = 0.0;
					exijc[c*2+0] = exijc[c*2+1] = 0.0;
					informative = false;
				//}else if( !ped->nonzeroDelX[j] ) { // 01.20.2009
				//	informative = false; // VERY VERY VERY IFFY, DO WE INCLUDE OR NOT INCLUDE THESE?
				}else{
					// compute the observed
					xijc[c*2+0] = ped->g[ ped->observed[j] ].genotype( 0, 0, 2, 2 );
					xijc[c*2+1] = ped->g[ ped->observed[j] ].genotype( 0, 0, 1, 2 );

					// compute the expected [[yes, this really shouldn't be necessary for all offspring, inefficient]]
					exijc[c*2+0] = exijc[c*2+1] = 0.0;
					for( unsigned int gg=0; gg<ped->g.size(); gg++ ) {
						exijc[c*2+0] += ped->g[gg].genotype( 0, 0, 2, 2 ) * ped->pg[gg];
						exijc[c*2+1] += ped->g[gg].genotype( 0, 0, 1, 2 ) * ped->pg[gg];
					}// gg

					// is the trait set?
          //if( !traitSet && !isnan(ped->trait[j]) ) {
          if( !isnan(ped->trait[j]) ) {
						eij = ped->trait[ j ];
						if( !qtl ) eij *= exp( - (*offset) );
						traitSet = true;
					}//fi( !traitSet )
				}//fi( ped[i].observed.size() < 0 )
			}// c

			// now compute the uij
			if( traitSet ){//& informative ) { // make sure someone had a contribution -- not really necessary (wasteful otherwise)
				// adjust the trait!
				for( int cc=0; cc<nc2; cc++ )
					eij -= bc[cc]*xijc[cc];

				// across the analyze alleles, compute uij
				if( !isnan(eij) ) {
					for( int a=0; a<na; a++ )
						ret_uij[ i + a*np ] += ( xija[a] - exija[a] ) * eij;
					// across the condition alleles
					for( int cc=0; cc<nc2; cc++ )
						ret_uij[ i + (na + cc)*np ] += ( xijc[cc] - exijc[cc] ) * eij;

					// FOR TWO STEP APPROACH ADDITION 10/22/2008 -- simply BAD, exported elsewhere!!!
					// -- a little overkill? Really should be able to just get away with the first one...
					//for( int a=0; a<na; a++ ) {
					//	Pedigree *ped = &data[referenceAnalyze[a]].ped[i];
					//	if( !isnan(ped->trait[j]) )
					//		ped->trait[j] = eij;
					//}
					//for( int c=0; c<nc; c++ ) {
					//	Pedigree *ped = &data[referenceCondition[c]].ped[i];
					//	if( !isnan(ped->trait[j]) )
					//		ped->trait[j] = eij;
					//}
					// NOITIDDA HCAORPPA PETS OWT ROF

					// compute delta x_ijm * xij_c^T
					for( int a=0; a<na; a++ )
						for( int cc=0; cc<nc2; cc++ )
							ret_xmxc[ a + cc*na ] += - ( xija[a] - exija[a] ) * xijc[cc]; // may be backwards
							//ret_xmxc[ cc + a*nc2 ] += - ( xija[a] - exija[a] ) * xijc[cc]; // may be backwards
					// compute delta x_ijc * xij_c^T
					for( int cc1=0; cc1<nc2; cc1++ )
						for( int cc2=0; cc2<nc2; cc2++ )
							ret_xcxc[ cc1 + cc2*nc2 ] += - ( xijc[cc1] - exijc[cc1] ) * xijc[cc2]; // symmetric? -- NO!
							//ret_xcxc[ cc2 + cc1*nc2 ] += - ( xijc[cc1] - exijc[cc1] ) * xijc[cc2]; // symmetric? -- NO!
				}// fi(!isnan(eij))
			}//fi( traitSet )
		}// j
	}// i
}

void condGeneFBATControl_numInfFam(int *reference, int *numInf) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_numInfFam %d no longer exists.\n", (*reference));
    return;
	}

	*numInf = data[*reference].numInfFam();
} // condGeneFBATControl_numInfFam

void condGeneFBATControl_pids(int *reference, int *pid) {
	if (*reference < 0 || *reference >= (int) data.size()) {
    Rprintf("condGeneFBATControl_pids %d no longer exists.\n", (*reference));
    return;
	}

	data[*reference].pids(pid);
}

void condGeneFBATControl_dUdBc(
	int *referenceAnalyze, int *referenceAnalyzeSize,
	int *referenceCondition, int *referenceConditionSize,
	int *analyze_allele_index, int *analyzeAlleleIndexSize,
	int *conditional_allele_index, int *conditionAlleleIndexSize, // for referenceAnalyze
	int *conditionAlleleIndex2, int *conditionAlleleIndexSize2, // for referenceCondition
	double *bc,
	double *ret_m, double *ret_c ) {

	int analyze_allele_index_size = *analyzeAlleleIndexSize;
	int conditional_allele_index_size = *conditionAlleleIndexSize;
	int conditional_allele_index_size2 = *conditionAlleleIndexSize2;

	int na = *referenceAnalyzeSize;
	int nc = *referenceConditionSize;
	int nc2 = nc*2;
	// make sure the references are good
	for( int a=0; a<na; a++ ) {
		if( referenceAnalyze[a] < 0 || referenceAnalyze[a] >= (int)data.size() ) {
			Rprintf("condGeneFBATControl_dUmdBc %d no longer exists. Terminating.\n", referenceAnalyze[a]);
			return;
		}
	}
	for( int c=0; c<nc; c++ ) {
		if( referenceCondition[c] < 0 || referenceCondition[c] >= (int)data.size() ) {
			Rprintf("condGeneFBATControl_dUmdBc %d no longer exists. Terminating.\n", referenceCondition[c]);
			return;
		}
	}
	int np = data[ referenceAnalyze[0] ].ped.size();

	// split up the bc
	//double bc0[nc], bc1[nc];
  vector<double> bc0; bc0.resize(nc);
  vector<double> bc1; bc1.resize(nc);
	for( int kc=0; kc<nc; kc++ ) {
    // CHECK THIS AGAIN -- THIS DOESN'T LOOK QUITE RIGHT, BUT SHOULD WORK FOR 2 MARKERS
		bc0[kc] = bc[kc];
		bc1[kc] = bc[kc+nc];
	}

	// zero out the return
	for( int ij=0; ij<na*nc2;  ij++ ) ret_m[ij] = 0.0;
	for( int ij=0; ij<nc2*nc2; ij++ ) ret_c[ij] = 0.0;

	// go across each pedigree
	for( int i=0; i<np; i++ ) {
		// across each analysis allele
		for( int a=0; a<na; a++ ) {
			//double Ai11[na][nc2];
			//double Ai10[na];
			//double Ai01[nc2];
			vector<vector<double> > Ai11; Ai11.resize(na); for(int aa=0; aa<na; aa++) Ai11[aa].resize(nc2);
      vector<double> Ai10; Ai10.resize(na);
      vector<double> Ai01; Ai01.resize(nc2);
			double Ai00 = 0.0;

			// Zero them out...
			for( int ka=0; ka<na; ka++ )
				Ai10[ka] = 0.0;
			for( int kc=0; kc<nc2; kc++ ) {
				Ai01[kc] = 0.0;
				for( int ka=0; ka<na; ka++ )
					Ai11[ka][kc] = 0.0;
			} // kc

			Pedigree *pp = &data[referenceAnalyze[a]].ped[i];
			for( unsigned int x=0; x<pp->genoDist.size(); x++ ) {
				// create all genotype permutations
				vector<int> genoToPerm;
				for (unsigned int j = 0; j < pp->genoDist[x].size(); j++)
					for (int k = 0; k < pp->genoDist[x][j]; k++)
						genoToPerm.push_back(j);
				vector<vector<int> > genoPerm;
				allPerms(genoToPerm, genoPerm);

				// calculate the test statistic pieces
				double permWeight = (double) 1.0 / (double) genoPerm.size();
				for( unsigned int p = 0; p < genoPerm.size(); p++ ) {
					//double sum_xstarjm[analyze_allele_index_size];
					vector<double> sum_xstarjm; sum_xstarjm.resize(analyze_allele_index_size);
					for( int a = 0; a < analyze_allele_index_size; a++ ) {
						sum_xstarjm[a] = 0.0;
						for (unsigned int j = 0; j < genoPerm[p].size(); j++) {
							if( pp->trait[j] == 1 )
								sum_xstarjm[a] += pp->g[genoPerm[p][j]].xCode( 0,analyze_allele_index[a], 2, ADDITIVE );
						} // j
					} // a

					//double sum_xstarjc0[conditional_allele_index_size];
					//double sum_xstarjc1[conditional_allele_index_size];
					vector<double> sum_xstarjc0; sum_xstarjc0.resize(conditional_allele_index_size);
          vector<double> sum_xstarjc1; sum_xstarjc1.resize(conditional_allele_index_size);
					for (int a = 0; a < conditional_allele_index_size; a++) {
						sum_xstarjc0[a] = 0.0;
						sum_xstarjc1[a] = 0.0;
						for (unsigned int j = 0; j < genoPerm[p].size(); j++) {
							if( pp->trait[j] == 1 ) {
								sum_xstarjc0[a] += pp->g[genoPerm[p][j]].genotype(0, conditional_allele_index[a], 2, 2);
								sum_xstarjc1[a] += pp->g[genoPerm[p][j]].genotype(0, conditional_allele_index[a], 1, 2);
							}
						} // j
					} // a

					double temp = 0.0; // exp( b0*x0 + ... )
					for (int a = 0; a < conditional_allele_index_size; a++)
						temp += (bc0[a] * sum_xstarjc0[a] + bc1[a] * sum_xstarjc1[a]);
					temp = exp(temp) * pp->genoDistP[x] * permWeight;

                    // finally compute the pieces
					Ai00 += temp;
					for( int ka=0; ka<na; ka++ )
						Ai10[ka] += sum_xstarjm[ka] * temp;
					for( int kc=0; kc<nc; kc++ ) {
						// WARNING -- MAY BE WRONG ORDER!!!
						Ai01[kc]    += sum_xstarjc0[kc] * temp;
						Ai01[nc+kc] += sum_xstarjc1[kc] * temp;
						for( int ka=0; ka<na; ka++ ) {
							Ai11[ka][kc]    += sum_xstarjm[ka] * sum_xstarjc0[kc] * temp;
							Ai11[ka][nc+kc] += sum_xstarjm[ka] * sum_xstarjc1[kc] * temp;
						}
					}

					// Fill in the return!
					double Ai00sq = Ai00 * Ai00;
					for( int ka=0; ka<na; ka++ ) {
						for( int kc=0; kc<nc2; kc++ ) {
							ret_m[ ka*nc2 + kc ] -= ( Ai11[ka][kc]*Ai00 - Ai10[ka]*Ai01[kc] ) / Ai00sq;
						} // kc
					} // ka
				} // p
			} // x
		} // a

		// Now we've got to a similar thing across the distribution for the conditional allele
		Pedigree *pp = &data[referenceCondition[0]].ped[i];
		//double Bi2[nc2][nc2];
    vector<vector<double> > Bi2; Bi2.resize(nc2); for(int bb=0; bb<nc2; bb++) Bi2[bb].resize(nc2);
		//double Bi1[nc2];
    vector<double> Bi1; Bi1.resize(nc2);
		double Bi0 = 0.0;
		for( int kc1=0; kc1<nc2; kc1++ ) {
			Bi1[kc1] = 0.0;
			for( int kc2=0; kc2<nc2; kc2++ )
				Bi2[kc1][kc2] = 0.0;
		} // kc1
		for( unsigned int x=0; x<pp->genoDist.size(); x++ ) {
			// create all genotype permutations
			vector<int> genoToPerm;
			for (unsigned int j = 0; j < pp->genoDist[x].size(); j++)
				for (int k = 0; k < pp->genoDist[x][j]; k++)
					genoToPerm.push_back(j);
			vector<vector<int> > genoPerm;
			allPerms(genoToPerm, genoPerm);

			// calculate the test statistic pieces
			double permWeight = (double) 1.0 / (double) genoPerm.size();
			for (unsigned int p = 0; p < genoPerm.size(); p++) {
				//double sum_xstarjc0[conditional_allele_index_size2];
				//double sum_xstarjc1[conditional_allele_index_size2];
				vector<double> sum_xstarjc0; sum_xstarjc0.resize(conditional_allele_index_size2);
        vector<double> sum_xstarjc1; sum_xstarjc1.resize(conditional_allele_index_size2);
				for (int a = 0; a < conditional_allele_index_size2; a++) {
					sum_xstarjc0[a] = 0.0;
					sum_xstarjc1[a] = 0.0;
					for (unsigned int j = 0; j < genoPerm[p].size(); j++) {
						if( pp->trait[j] == 1 ) {
							sum_xstarjc0[a] += pp->g[genoPerm[p][j]].genotype(0,conditionAlleleIndex2[a], 2, 2);
							sum_xstarjc1[a] += pp->g[genoPerm[p][j]].genotype(0,conditionAlleleIndex2[a], 1, 2);
						}
					} // j
				} // a

				double temp = 0.0; // exp( b0*x0 + ... )
				for (int a = 0; a < conditional_allele_index_size2; a++)
					temp += (bc0[a] * sum_xstarjc0[a] + bc1[a]*sum_xstarjc1[a]);
				temp = exp(temp) * pp->genoDistP[x] * permWeight;

                // finally compute the pieces
				Bi0 += temp;
				for( int kc1=0; kc1<nc; kc1++ ) { // different kc1 than below
					Bi1[kc1]    += sum_xstarjc0[ kc1 ] * temp;
					Bi1[nc+kc1] += sum_xstarjc1[ kc1 ] * temp;
					for( int kc2=0; kc2<nc; kc2++ ) { // different kc1 than below
						//Bi2[kc1   ][kc2   ] += sum_xstarjc0[ kc1    ] * sum_xstarjc1[ kc2    ] * temp;
						//Bi2[nc+kc1][kc2   ] += sum_xstarjc0[ nc+kc1 ] * sum_xstarjc1[ kc2    ] * temp;
						//Bi2[kc1   ][nc+kc2] += sum_xstarjc0[ kc1    ] * sum_xstarjc1[ nc+kc2 ] * temp;
						//Bi2[nc+kc1][nc+kc2] += sum_xstarjc0[ nc+kc1 ] * sum_xstarjc1[ nc+kc2 ] * temp;
            ///////////////////////////////////////////////////
            // HERE IS A MAJOR CHANGE FROM THE OTHER CODE!!! //
            ///////////////////////////////////////////////////
            Bi2[kc1   ][kc2   ] += sum_xstarjc0[kc1] * sum_xstarjc0[kc2] * temp;
            Bi2[nc+kc1][kc2   ] += sum_xstarjc1[kc1] * sum_xstarjc0[kc2] * temp;
            Bi2[kc1   ][nc+kc2] += sum_xstarjc0[kc1] * sum_xstarjc1[kc2] * temp;
            Bi2[nc+kc1][nc+kc2] += sum_xstarjc1[kc1] * sum_xstarjc1[kc2] * temp;
					} // kc2
				} // kc1

				// And finally fill in the return!
				double Bi0sq = Bi0*Bi0;
				for( int kc1=0; kc1<nc2; kc1++ ) {
					for( int kc2=0; kc2<nc2; kc2++ ) {
						ret_c[ kc1*nc2 + kc2 ] += ( Bi2[kc1][kc2]*Bi0 - Bi1[kc1]*Bi1[kc2] ) / Bi0sq;
					} // kc2
				} // kc1
			} // p
		} // x
	} // i
}

/*
void condGeneFBATControl_dUdBc(
  int *referenceAnalyze, int *referenceAnalyzeSize,
  int *referenceCondition, int *referenceConditionSize,
  int *analyze_allele_index, int *analyzeAlleleIndexSize,
  int *conditional_allele_index, int *conditionAlleleIndexSize, // for referenceAnalyze
  int *conditionAlleleIndex2, int *conditionAlleleIndexSize2, // for referenceCondition
  double *bc,
  double *ret_m, double *ret_c ) {

  int analyze_allele_index_size = *analyzeAlleleIndexSize;
  int conditional_allele_index_size = *conditionAlleleIndexSize;
  int conditional_allele_index_size2 = *conditionAlleleIndexSize2;

  int na = *referenceAnalyzeSize;
  int nc = *referenceConditionSize;
  int nc2 = nc*2;
  // make sure the references are good
  for( int a=0; a<na; a++ ) {
    if( referenceAnalyze[a] < 0 || referenceAnalyze[a] >= (int)data.size() ) {
      cout << "condGeneFBATControl_dUmdBc " << referenceAnalyze[a] << " no longer exists. Terminating." << endl;
      return;
    }
  }
  for( int c=0; c<nc; c++ ) {
    if( referenceCondition[c] < 0 || referenceCondition[c] >= (int)data.size() ) {
      cout << "condGeneFBATControl_dUmdBc " << referenceCondition[c] << " no longer exists. Terminating." << endl;
      return;
    }
  }
  int np = data[ referenceAnalyze[0] ].ped.size();

  // split up the bc
  double bc0[nc], bc1[nc];
  for( int kc=0; kc<nc; kc++ ) {
    bc0[kc] = bc[kc];
    bc1[kc] = bc[kc+nc];
  }

  // zero out the return
  for( int ij=0; ij<na*nc2;  ij++ ) ret_m[ij] = 0.0;
  for( int ij=0; ij<nc2*nc2; ij++ ) ret_c[ij] = 0.0;

  // go across each pedigree
  for( int i=0; i<np; i++ ) {
    // Pointer to the current pedigree
    Pedigree *pp = &data[referenceAnalyze[a]].ped[i];

    // Storage for Ai__
    double Ai11[na][nc2];
    double Ai10[na];
    double Ai01[nc2];
    double Ai00 = 0.0;
    for( int ka=0; ka<na; ka++ )
      Ai10[ka] = 0.0;
    for( int kc=0; kc<nc2; kc++ ) {
      Ai01[ka] = 0.0;
      for( int ka=0; ka<na; ka++ )
        Ai11[ka][kc] = 0.0;
    }// kc

    // Storage for Bi_
    double Bi2[nc2][nc2];
    double Bi1[nc2];
    double Bi0 = 0.0;
    for( int kc1=0; kc1<nc2; kc1++ ) {
      Bi1[kc1] = 0.0;
      for( int kc2=0; kc2<nc2; kc2++ )
        Bi2[kc1][kc2] = 0.0;
    }// kc1

    // now, for each genotype offspring possibility
    for( unsigned int x=0; x<pp->genoDist.size(); x++ ) {
      // create all genotype permutations
      vector<int> genoToPerm;
      for( unsigned int j=0; j<pp->genoDist[x].size(); j++ )
        for( int k=0; k<genoDist[x][j]; j++ )
          genoToPerm.push_back(j);
      vector< vector<int> > genoPerm;
      allPerms( genoToPerm, genoPerm );

      // now, go across the permutations
      double permWeight = 1.0 / (double)genoPerm.size();
      for( unsigned int p=0; p<genoPerm.size(); p++ ) {
        // ----------------------
        // computations for Ai_
        double sum_xstarjm[analyze_allele_index_size];
        for( int a=0; a<analyze_allele_index_size; a++ ) {
          sum_xstarjm[a] = 0.0;
          for( unsigned int j=0; j<genoPerm[p].size(); j++ ) {
            if( pp->trait[j] == 1 ) {
              sum_xstarjm[a] += pp->g[genoPerm[p][j]].xCode( 0, analyze_allele_index[a], 2, ADDITIVE );
            }//fi
          }//j
        }//a
        double sum_xstarjc0[conditional_allele_index_size];
        double sum_xstarjc1[conditional_allele_index_size];
        for( int a=0; a<conditional_allele_index_size; a++ ){
          sum_xstarjc0[a] = sum_xstarjc1[a] = 0.0;
          for( unsigned int j=0; j<genoPerm[p].size(); j++ ) {
            if( pp->trait[j] == 1 ) {
              sum_xstarjc0[a] += pp->g[genoPerm[p][j]].genotype( 0, conditional_allele_index[a], 2, 2 );
              sum_xstarjc1[a] += pp->g[genoPerm[p][j]].genotype( 0, conditional_allele_index[a], 1, 2 );
            }//fi
          }//j
        }//a

        // exp( bc1*xc1 + ... )
        double temp = 0.0;
        for( int a=0; a<conditional_allele_index_size; a++ )
          temp += (bc0[a]*sum_xstarjc0[a] + bc1[a]*sum_xstarjc1[a]);
        temp = exp(temp)*pp->genoDistP[x] * permWeight;

        // compute the pieces
        Ai00 += temp;
        for( int ka=0; ka<na; ka++ )
          Ai10[ka] += sum_xstarjm[ka] * temp;
        for( int kc=0; kc<nc; kc++ ) {
          Ai01[kc] += sum_xstarjc0[kc] * temp;
          Ai01[nc+kc] += sum_xstarjc[kc] * temp;
          for( int ka=0; ka<na; ka++ ) {
            Ai11[ka][kc]   += sum_xstarjm[ka] * sum_xstarjc0[kc] * temp;
            Ai11[ka][nc+kc]+= sum_xstarjm[ka] * sum_xstarjc1[kc] * temp;
          }//ka
        }//kc

        // fill in the return for Ai_
        double Ai00sq = Ai00 * Ai00;
        for( int ka=0; ka<na; ka++ )
          for( int kc=0; kc<nc2; kc++ )
            ret_m[ ka*nc2 + kc ] -= ( Ai11[ka][kc]*Ai00 - Ai10[ka]*Ai01[kc] ) / Ai00sq;

        //-----------------------
        // Computations for Bi_
        // TRY COMMENT OUT BEGIN --->
        for( int a=0; a<conditional_allele_index_size2; a++ ) {
          sum_xstarjc0[a] = sum_xstarjc1[a] = 0.0;
          for( unsigned int j=0; j<genoPerm[p].size(); j++ ) {
            if( pp->trait[j] == 1 ) {
              sum_xstarjc0[a] += pp->g[genoPerm[p][j]].genotype( 0, conditionAlleleIndex2[a], 2, 2 );
              sum_xstarjc1[a] += pp->g[genoPerm[p][j]].genotype( 0, conditionAlleleIndex2[a], 1, 2 );
            }//fi
          }//j
        }//a

        // exp( bc1*xc1 + ... ) --> This is really, really confusing...
        temp = 0.0;
        for( int a=0; a<conditional_allele_index_size2; a++ )
          temp += (bc0[a]*sum_xstarjc0[a] + bc1[a]*sum_xstarjc1[a]);
        temp = exp(temp) * pp->genoDistP[x] * permWeight;
        // TRY COMMENT OUT END <--

        // compute the pieces
        Bi0 += temp;
        for( int kc1=0; kc1<nc; kc1++ ) {
          Bi1[kc1]   += sum_xstarjc0[kc1] * temp;
          Bi1[nc+kc1]+= sum_xstarjc1[kc1] * temp;
          for( int kc2=0; kc2<nc; kc2++ ) {
            ///////////////////////////////////////////////////
            // HERE IS A MAJOR CHANGE FROM THE OTHER CODE!!! //
            ///////////////////////////////////////////////////
            Bi2[kc1][kc2]      += sum_xstarjc0[kc1] * sum_xstarjc0[kc2];
            Bi2[nc+kc1][kc2]   += sum_xstarjc1[kc1] * sum_xstarjc0[kc2];
            Bi2[kc1][nc+kc2]   += sum_xstarjc0[kc1] * sum_xstarjc1[kc2];
            Bi2[nc+kc1][nc+kc2]+= sum_xstarjc1[kc1] + sum_xstarjc1[kc2];
          }//kc2
        }//kc1

        // fill in the return for Bi_
        double Bi0sq = Bi0 * Bi0;
        for( int kc1=0; kc1<nc2; kc1++ )
          for( int kc2=0; kc2<nc2; kc2++ )
            ret_c[ kc1*nc2 + kc2 ] += ( Bi2[kc1][kc2]*Bi0 - Bi1[kc1]*Bi1[kc2] ) / Bi0sq;
      }//p
    }//x
  }//i
}
*/

// DEPRECATED
// Alleles are treated as conditioning alleles, so both test and condition are two parameter
// Proportion of variability explained (like R^2) = 1 - SS_err / SS_tot = 1 - ( sum[yi-yhati] ) / ( sum[ yi - ybar ] )
void condGeneFBATControl_varExplConts(
		int *reference, int *referenceSize,
		double *betaEst,
		double *ret_varExpl ) {
	// ASSUMES DATA IS ALREADY MEAN CENTERED...

	Rprintf("condGeneFBATControl_varExplConts is deprecated.\n");

	int nc = *referenceSize;
	//////int nc2 = nc * 2;
	int np = data[ reference[0] ].ped.size();

	//cout << "nc: " << nc << " nc2: " << nc2 << " np: " << np << endl;

	// make sure the references are good
	for( int c=0; c<nc; c++ ) {
		if( reference[c] < 0 || reference[c] >= (int)data.size() ) {
			Rprintf("condGeneFbatControl_varExplConts %d no longer exists\n", reference[c]);
		}
	}

	double ymean = 0.0; // in case data isn't really centered...
	double varModel = 0.0;
	double varMean = 0.0;

	// first compute ymean
	// - go across each pedigree
	int ymeanCount = 0;
	for( int i=0; i<np; i++ ) {
		// determine the number of offspring
		unsigned int noff = data[ reference[0] ].ped[i].observed.size();
		for( int c=1; c<nc; c++ )
			if( data[reference[c]].ped[i].observed.size() < noff )
				noff = data[reference[c]].ped[i].observed.size();

		// go across each offspring
		for( unsigned int j=0; j<noff; j++ ) {
			double y = data[ reference[0] ].ped[i].trait[j];
			if( !isnan(y) ) {
				ymean += y;
				ymeanCount++;
			}
    }//j
	}//i
	if( ymeanCount < 1 ) {
		*ret_varExpl = 0.0;
		Rprintf("No variation in trait!\n");
		return;
	}
	ymean /= (double)ymeanCount;
	//ymean = 0.0;

  /*
	// DEBUG ONLY -- MEAN OF EACH GENOTYPE
	int geno_a[] = {2,2,1}, geno_b[] = {2,1,1};
	for( int g=0; g<3; g++ ) {
		int c = 0;
		double ym = 0.0; int ymCount = 0;
		for( int i=0; i<np; i++ ) {
			// determine the number of offspring
			int noff = data[ reference[0] ].ped[i].observed.size();
			for( int c=1; c<nc; c++ )
				if( data[reference[c]].ped[i].observed.size() < noff )
					noff = data[reference[c]].ped[i].observed.size();

			Pedigree *ped = &data[reference[c]].ped[i];
			// go across each offspring
			for( int j=0; j<noff; j++ ) {
				double y = data[ reference[0] ].ped[i].trait[j];
				if( !isnan(y) && ped->g[ped->observed[j]].genotype(0,0,geno_a[g],geno_b[g])==1 ) {
					ym += y;
					ymCount++;
				}
			}
		}
		ym /= (double)ymCount;
		cout << "MEAN TRAIT (" << geno_a[g] << "/" << geno_b[g] << "): " << ym << endl;
	}
	// DEBUG END
  */

	// Now compute the variances
	// - go across each pedigree to compute the variances
	for( int i=0; i<np; i++ ) {
		// determine the number of offspring
		unsigned int noff = data[ reference[0] ].ped[i].observed.size();
		for( int c=1; c<nc; c++ )
			if( data[reference[c]].ped[i].observed.size() < noff )
				noff = data[reference[c]].ped[i].observed.size();

    // ????
		if( data[reference[0] ].ped[i].trait.size() < noff )
      noff = data[reference[0]].ped[i].trait.size();

		// go across each offspring
		for( unsigned int j=0; j<noff; j++ ) {
			bool informative = true;

			// go across the "conditioning" alleles
			double bx = 0.0;
			for( int c=0; c<nc; c++ ) {
				Pedigree *ped = &data[reference[c]].ped[i];
				if( j >= ped->observed.size() ) {
					informative = false;
				}else{
					// compute the observed
					double xijc0 = ped->g[ ped->observed[j] ].genotype( 0, 0, 2, 2 );
					double xijc1 = ped->g[ ped->observed[j] ].genotype( 0, 0, 1, 2 );

					bx += betaEst[c*2]*xijc0 + betaEst[c*2+1]*xijc1;
					//bx += betaEst[c]*xijc0 + betaEst[c+nc]*xijc1;
				}
			}

			if( informative ) {
				// finally compute the variance pieces
				double y = data[ reference[0] ].ped[i].trait[j];
				if( !isnan(y) ) {
					//ymean += y;
					varModel += ( y - bx ) * ( y - bx );
					varMean  += ( y - ymean ) * ( y - ymean );
					//varModel += fabs( y - bx );
					//varMean  += fabs( y - ymean );

					/*
					cout << "y: " << y << " ymean: " << ymean << " yfit: " << bx << " |y-ymean|: " << fabs(y-ymean) << " |y-yfit|: " << fabs(y-bx) << endl;
					//cout << "beta (1 parm): " << betaEst[0] << "," << betaEst[1] << endl;
					cout << " genotype (2parm): "
						<< data[reference[0]].ped[i].g[ data[reference[0]].ped[i].observed[j] ].genotype( 0, 0, 2, 2 ) << ","
						<< data[reference[0]].ped[i].g[ data[reference[0]].ped[i].observed[j] ].genotype( 0, 0, 1, 2 ) << ","
						<< data[reference[0]].ped[i].g[ data[reference[0]].ped[i].observed[j] ].genotype( 0, 0, 1, 1 ) << "  " << endl;
						*/

					/*
					cout << "y: " << y << " ymean: " << ymean << " yfit: " << bx << " |y-ymean|: " << fabs(y-ymean) << " |y-yfit|: " << fabs(y-bx) << endl;
					cout << "beta (2 parm): " << betaEst[0] << "," << betaEst[1] << "," << betaEst[2] << "," << betaEst[3] << endl;

					cout << " genotype (2parm): "
						<< data[reference[0]].ped[i].g[ data[reference[0]].ped[i].observed[j] ].genotype( 0, 0, 2, 2 ) << ","
						<< data[reference[0]].ped[i].g[ data[reference[0]].ped[i].observed[j] ].genotype( 0, 0, 1, 2 ) << ","
						<< data[reference[0]].ped[i].g[ data[reference[0]].ped[i].observed[j] ].genotype( 0, 0, 1, 1 ) << "  "
						<< data[reference[1]].ped[i].g[ data[reference[1]].ped[i].observed[j] ].genotype( 0, 0, 2, 2 ) << ","
						<< data[reference[1]].ped[i].g[ data[reference[1]].ped[i].observed[j] ].genotype( 0, 0, 1, 2 ) << ","
						<< data[reference[1]].ped[i].g[ data[reference[1]].ped[i].observed[j] ].genotype( 0, 0, 1, 1 ) << endl; */
				}
			}
		}
	}

	//cout << "varModel: " << varModel << ", varMean: " << varMean << endl;
	//cout << "varModel / varMean : " << varModel / varMean << endl;

	*ret_varExpl = 1.0 - ( varModel / varMean );
}

// references are needed so that we know which traits won't be used because they have missing genotype data
void condGeneFBATControl_varContsMean(
		int *reference, int *referenceSize,
		double *betaEst,
		double *ret_var ) {
	// make sure the references are good
	for( int c=0; c<*referenceSize; c++ ) {
		if( reference[c] < 0 || reference[c] >= (int)data.size() ) {
			Rprintf("condGeneFbatControl_varExplConts %d no longer exists\n", reference[c]);
			return;
		}
	}

	// fill in some constants
	int nc = *referenceSize;
	//////int nc2 = nc*2;
	int np = data[ reference[0] ].ped.size();

	// fill in the mean, but select those that are informative to the current markers...
	double ymean = 0.0;
	vector<int> good_i, good_j;
	for( int i=0; i<np; i++ ) {
		// determine the number of offspring
		unsigned int noff = data[reference[0]].ped[i].observed.size();
		for( int c=1; c<nc; c++ )
			if( data[reference[c]].ped[i].observed.size() < noff )
				noff = data[reference[c]].ped[i].observed.size();

		// go across each offspring
		for( unsigned int j=0; j<noff; j++ ) {
			bool informative = true;

			// go across the conditioning alleles
			for( int c=0; c<nc; c++ ) {
				Pedigree *ped = &data[reference[c]].ped[i];
				if( j >= ped->observed.size() ) {
					informative = false;
				}else if( !ped->nonzeroDelX[j] ) { // 01.20.2009
					informative = false;
				}else{ // doesn't matter here
				}
			}

			if( informative ) {
				double y = data[reference[0]].ped[i].trait[j];
				if( !isnan(y) ) {
					//cout << "y = " << y << endl;
					ymean += y;
					good_i.push_back(i);
					good_j.push_back(j);
				}
			}// fi(informative)
		}// j
	}// i}

	ymean /= good_i.size();

	//cout << "good_i.size() " << good_i.size() << endl;

	//cout << "YMEAN = " << ymean << endl;

	// finally compute the variance
	double var = 0.0;
	for( unsigned int k=0; k<good_i.size(); k++ ) {
		double y = data[reference[0]].ped[good_i[k]].trait[good_j[k]];
		var += ( y - ymean ) * ( y - ymean );
	}

	//cout << "VAR = " << var << endl;

	*ret_var = var;
}// condGeneFBATControl_varContsMean

void condGeneFBATControl_varContsModel(
		int *reference, int *referenceSize,
		double *betaEst, // actually betaEst isn't used...
		double *ret_var ) {
	// make sure the references are good
	for( int c=0; c<*referenceSize; c++ ) {
		if( reference[c] < 0 || reference[c] >= (int)data.size() ) {
			Rprintf("condGeneFbatControl_varExplConts %d no longer exists.\n", reference[c]);
			return;
		}
	}

	// fill in some constants
	int nc = *referenceSize;
	//////int nc2 = nc*2;
	int np = data[ reference[0] ].ped.size();

	// fill in the variance
	double varModel = 0.0;
	for( int i=0; i<np; i++ ) {
		// determine the number of offspring
		unsigned int noff = data[reference[0]].ped[i].observed.size();
		for( int c=1; c<nc; c++ )
			if( data[reference[c]].ped[i].observed.size() < noff )
				noff = data[reference[c]].ped[i].observed.size();

		// go across each offspring
		for( unsigned int j=0; j<noff; j++ ) {
			bool informative = true;

			// go across the conditioning alleles
			double bx = 0.0;
			for( int c=0; c<nc; c++ ) {
				Pedigree *ped = &data[reference[c]].ped[i];
				if( j >= ped->observed.size() ) {
					informative = false;
//				}else if( !ped->nonzeroDelX[j] ) { // 01.20.2009
//					informative = false;
				}else{
					// compute the observed
					double xijc0 = ped->g[ ped->observed[j] ].genotype( 0, 0, 2, 2 );
					double xijc1 = ped->g[ ped->observed[j] ].genotype( 0, 0, 1, 2 );

					bx += betaEst[c*2]*xijc0 + betaEst[c*2+1]*xijc1;
				}
			}

			if( informative ) {
				// finally compute the variance pieces
				double y = data[reference[0]].ped[i].trait[j];
				if( !isnan(y) )
					varModel += ( y - bx ) * ( y - bx );
			}// fi(informative)
		}// j
	}// i

	*ret_var = varModel;
}

void condGeneFBATControl_backupTrait( int *reference, int *referenceSize ) {
	int R = *referenceSize;
	for( int c=0; c<R; c++ ) {
		// first make sure it's a good reference
		if( reference[c] < 0 || reference[c] >= (int)data.size() ) {
      Rprintf("condGeneFBATControl_backupTrait::Reference %d no longer exists.\n", reference[c]);
      return;
		}

		// now go ahead and backup
		int P = data[reference[c]].ped.size();
		for( int p=0; p<P; p++ )
			data[reference[c]].ped[p].traitBackup = data[reference[c]].ped[p].trait;
	}
}
void condGeneFBATControl_restoreTrait( int *reference, int *referenceSize ) {
	int R = *referenceSize;
	for( int c=0; c<R; c++ ) {
		// first make sure it's a good reference
		if( reference[c] < 0 || reference[c] >= (int)data.size() ) {
			//cout << "condGeneFbatControl_restoreTrait " << reference[c] << " no longer exists." << endl;
      Rprintf("condGeneFBATControl_restoreTrait::Reference %d no longer exists.\n", reference[c]);
      return;
		}

		// now go ahead and restore!
		int P = data[reference[c]].ped.size();
		for( int p=0; p<P; p++ )
			data[reference[c]].ped[p].trait = data[reference[c]].ped[p].traitBackup;
	}
}


} // extern "C"



#ifdef _cgFbat_DEBUG_

int main() {
	Genotype g;
	g.debug();

	return(0);
}

#endif
