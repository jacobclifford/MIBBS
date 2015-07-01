#ifndef EXPR_PREDICTOR_H
#define EXPR_PREDICTOR_H



#include <cctype>
#include <cstring>
#include <algorithm>
#include "Tools.h"
#include <stdio.h>


// alphabet
const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' };
const char ALPHABET2[] = { 'A', 'C', 'G', 'T', 'N', '-' };
const int NBASES = 4;	// number of bases
const int ALPHABET_SIZE = 6;
const int MISSING = 4;		// missing nt.
const int GAP = 5;		// '-'

// test if a is a nucleotide {A,G,C,T}
bool isNt( int a );
// the complement of a nt.
int complement( int a );
// conversion of a symbol (nt. or gap or N) to its index in ALPHABET
int symbolToInt( char c );
// character of the strand (+ or -)
char strand2char( bool strand );
// strand (+ or -) from the character
bool char2strand( char c );
// nucleotide distribution from GC content
vector< double > createNtDistr( double gcContent ); 

// file formats for Sequence
enum FileFormats { FASTA, PHYLIP };

class Sequence {
public:
    // constructors
    Sequence() : nts() {}     // nts is a private data member vector< int > nts, nts(i) is the ith nucleotide
    Sequence( const vector< int >& _nts ) : nts( _nts ) {}
    Sequence( const string& str );
    Sequence( const Sequence& other, int start, int length, bool strand = true ); 
    void copy( const Sequence& other ) { nts = other.nts; }
    Sequence( const Sequence& other ) { copy( other ); }
    
    // assignment
    Sequence& operator=( const Sequence& other ) { copy( other ); return *this; }

    // sequence comparison
    bool operator==( const Sequence& other ) {
        if ( nts == other.nts ) return true;
        else return false;
    }
    
    // access methods
    int size() const { return nts.size(); }
    const int& operator[]( int pos ) const {
        assert( pos >= 0 && pos < size() );
        return nts[ pos ];	
    }
    int& operator[]( int pos ) {
        assert( pos >= 0 && pos < size() );
        return nts[ pos ];	
    }
    
    // insert a nt. or a sequence segment at the end of sequence
    int push_back( int nt );
    int push_back( const Sequence& elem );	
    
    // the reverse complement of a sequence
    Sequence compRevCompl() const;
    
    // information of the sequence
    void getNtCounts( vector< int >& counts ) const;
    bool containsMissing() const;
            
    // clear the sequence
    void clear() { nts.clear(); }
    
    // load Sequence from a file
    int load( const string& file, string& name, int format = FASTA );
    int load( const string& file, int format = FASTA );	
    
    // output
    friend ostream& operator<<( ostream& os, const Sequence& seq );
private:
    vector< int > nts;		// nts[ i ]: the i-th nt.
};

// read sequences from a file
int readSequences( const string& file, vector< Sequence >& seqs, vector< string >& names, int format = FASTA );
int readSequences( const string& file, vector< Sequence >& seqs, int format = FASTA );

// write sequences to a file
int writeSequences( const string& file, const vector< Sequence >& seqs, const vector< string >& names, int format = FASTA );
int writeSequences( const string& file, const vector< Sequence >& seqs, int format = FASTA );

// output max string length in a vector of strings
int maxstring( string& file );

const double PSEUDO_COUNT = 0.25;

// construct the position weight matrix from the count matrix
Matrix compWtmx( const Matrix& countMatrix, double pseudoCount );
class Motif {
public:	
    // constructors
    Motif() : pwm(), background() {}	
    //Motif( const string& file ) ;
    Motif( const Matrix& _pwm, const vector< double >& _background ); 
    Motif( const Matrix& countMatrix, double pseudoCount, const vector< double >& _background );		// countMatrix in Transfac format
    void copy( const Motif& other ) { pwm = other.pwm; background = other.background; LLRMat = other.LLRMat; maxSite = other.maxSite; maxLLR = other.maxLLR; infomat=other.infomat;}
    Motif( const Motif& other ) { copy( other ); }


    // assignment
    Motif& operator=( const Motif& other ) { copy( other ); return *this; }
    
    // access methods
    int length() const { return pwm.nRows(); }
    const Matrix& getPwm() const { return pwm; } 
	  Matrix& getPwm2() { return pwm; } 
    const vector< double >& getBackground() const { return background; }
    const Matrix& getLLRMat() const { return LLRMat; }
    const Sequence& getMaxSite() const { return maxSite; }
    double getMaxLLR() const { return maxLLR; }
    double getgenomicMean();
	Matrix getInfoMat() ;   // returns the information matrix, each element is f*log(f)
    // compute the log-likelihood ratio of a sequence element
	Matrix getinfomat2() { return infomat; }
    double LLR( const Sequence& elem ) const; 
    double LLR2(  Sequence& elem ) ;
    // compute the energy of a sequence element, relative to the strongest site (thus always >= 0)
    double energy( const Sequence& elem ) const;
     double INFO( const Sequence& elem ) const; 
     double information( const Sequence& elem ) const; 
    double getMean();
    double getRMSE_X_djor();
    double getInformation();
	double getInformation( int nseqs);
    // sample a site from PWM
    void sample( const gsl_rng* rng, Sequence& elem, bool strand = true ) const;
    double Pvalue( double alpha ) ;   
    double PvalueSum( double alpha );              
    // load a Motif
    int load( const string& file, const vector< double >& background, string& name );
    int load( const string& file, const vector< double >& background );
    void setEth(double e){ eth = e;  }
    double getEth(){ return eth;  }
    // output
    friend ostream& operator<<( ostream& os, const Motif& motif );	
private:
	//Matrix counts;
    Matrix pwm;	// the position weight matrix: f_i(b), the frequency of b at position i
	double eth;
   // Matrix pwm2;
    vector< double > background;	// background distribution
    Matrix LLRMat;	// LLR matrix, M(i,b) = log( f_i(b) / p(b) ), where f_i(b) is the frequency of b at position i of PWM, and p(b) the frequency of b at the background
    Sequence maxSite;	// the sequence of the strongest site
    double maxLLR;	// LLR of the strongest site
     Matrix infomat;   
	  Matrix countmat;   
    // initialization: compute the LLR matrix and maxLLR
    void init();	
};

// read the motifs (PWMs) in a FASTA-like format: pseudocounts within the file
int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs, vector< string >& names );
int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs );

void alignNEW(string& file,int nrandStarts, int nUpdates , int wid, int l   );
Matrix countmatrix( const string& file );
Matrix countmatrixArich( const string& file );
Matrix countmatrixflip( const string& file );
Matrix countmatrixS( const string& file );
Matrix countmatrix2( string& file, string& file2  );
void scorematrix( Matrix e , string& file, double maxLLR, string& file2 ); // e is LLRmat, file contains sequences LLRmat acts on, maxllr is maxLLR, file2 is output energies.
Matrix countmatrixerr( const string& file ) ;
// Sensitivity compares
void Sensitivity( const string& dc , const string& dc2 );  
// SensitiivityFN compares the data with the iferred sites.  given an inferred site, it searches all the inferred sies for a hit (TP ), otherwise it's labled as (TN )
void SensitivityFN( const string& core , const string& dc2 ); // takes core sites and checkes against inferred sites (dc2), only works for dorsalcondit
void SensitivityFN2( const string& core , const string& dc2, const string& tpout , const string& tnout  );// same as SensitivityFN but added the files for output.

class Site {
public:
    // constructors
    Site() : start( 0 ), strand( true ), factorIdx( 0 ), energy( 0 ), wtRatio( 1 ) , dist( 0 ) , information(0){}
    Site( int _start, bool _strand, int _factorIdx ) : start( _start ), strand( _strand ), factorIdx( _factorIdx ), energy( 0 ), wtRatio( 1 ), dist( 0 ), information(0) {}
    Site( int _start, bool _strand, int _factorIdx, double _energy, double _information ) : start( _start ), strand( _strand ), factorIdx( _factorIdx ), energy( _energy ), information( _information) { wtRatio = exp(- energy ); dist = 0; }	
 Site( int _start, bool _strand, int _factorIdx, double _energy ) : start( _start ), strand( _strand ), factorIdx( _factorIdx ), energy( _energy ) { wtRatio = exp(- energy ); dist = 0; }	
	 Site( int _start, bool _strand, double _information ) : start( _start ), strand( _strand ), information( _information){ wtRatio = 1; dist = 0; energy=0; }
    void copy( const Site& other ) { start = other.start; strand = other.strand; factorIdx = other.factorIdx; energy = other.energy; wtRatio = other.wtRatio; dist = other.dist; information=other.information;}
    Site( const Site& other ) { copy( other ); }
    
    // assignment
    Site& operator=( const Site& other ) { copy( other ); return *this; }	
            
    friend ostream& operator<<( ostream& os, const Site& site );	
    
    int start;		// start position: 0-based  // really ? it seems 1 based, from couting sequences of sites.
    bool strand;	// 1: positive; 0: negative  //negative means reverse so 0
    int factorIdx;	// the index of the associated TF, starting from 0
    double energy;	// the energy relative to the strongest site (nonnegative)
    double wtRatio;	// the binding weight ratio (<= 1) of site vs the strongest site: K(S) / K(S_max) = q(S) / q(S_max)
    int dist;           // determines how far in basepairs a dorsal site is from nearest twist site
    double information;
};

// test if two sites overlap
bool siteOverlap( const Site& a, const Site& b, const vector< Motif >& motifs );

// find the strongest site among sitevector
//Site SeqAnnotator::siteMax( vector< Site >& sites);

// delete all sites within sitest that overlap with sites from sitesp
void sitestoverlap( vector< Site >& sitest, vector< Site >& sitesp );
// representation of Sequence as a Site vector
typedef vector< Site > SiteVec;		

// read sites (of potentially multiple sequences)
int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, vector< string >& names, bool readEnergy = false );
int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, bool readEnergy = false );

int readSites2Matrix( const string& file );
bool siteSortPredicate(const Site& d1, const Site& d2);
bool siteSortPredicate2(double d1, double d2);
int Random(int n);




enum ModelType {
    LOGISTIC,   // logistic regression
    DIRECT,     // direct interaction between TF and BTM, repressor works through BTM
    QUENCHING,      // repressor stops activator from interacting with BTM
    CHRMOD_UNLIMITED,     // repressor works by chromatin modification (making it unaccessible), unlimited activation
    CHRMOD_LIMITED,        // repressor works by chromatin modification (making it unaccessible), limited activation
    BINS
};

ModelType getModelOption( const string& modelOptionStr );
string getModelOptionStr( ModelType modelOption ); 

enum FactorIntType {
    BINARY,     // Binary model of interaction
    GAUSSIAN,    // Gaussian model of interaction
    BINSF
}; 

string getIntOptionStr( FactorIntType intOption );

enum ObjType {
    SSE,    // sum of squared error
    CORR,   // Pearson correlation
    CROSS_CORR  // cross correlation (maximum in a range of shifts)
};

ObjType getObjOption( const string& objOptionStr );
string getObjOptionStr( ObjType objOption );

enum SearchType {
    UNCONSTRAINED,  // unconstrained search
    CONSTRAINED     // constrained search
};

string getSearchOptionStr( SearchType searchOption );

class FactorIntFunc {
public:
    // compute the factor interaction, given the normal interaction (when they are close enough)
    virtual double compFactorInt( double normalInt, double dist, bool orientation ) const = 0;

    // the maximum distance beyond which there is no interaction
    virtual double getMaxDist() const = 0;    
};

class FactorIntFuncBinsf : public FactorIntFunc {
public:
    // constructors // if the 2 and 3 parameters are switched, a compile error: default argument missing for parameter 3 of ‘FactorIntFuncBinsf::FactorIntFuncBinsf(double, double, int)
FactorIntFuncBinsf(double _distThr,  int _nbins, double _orientationEffect = 1.0 ) : distThr( _distThr ),  nbins(_nbins),orientationEffect( _orientationEffect ) {}

    // compute the factor interaction
    double compFactorInt( double normalInt, double dist, bool orientation ) const;

    // the maximum distance beyond which there is no interaction
    double getMaxDist() const {
        return distThr;
    } 
private:
    int nbins;
    double distThr;		// if distance < thr, the "normal" value; otherwise 1 (no interaction)
    double orientationEffect;	// the effect of orientation: if at different strands, the effect should be multiplied this value	
};


class FactorIntFuncBinary : public FactorIntFunc {
public:
    // constructors
    FactorIntFuncBinary( double _distThr, double _orientationEffect = 1.0 ) : distThr( _distThr ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

    // compute the factor interaction
    double compFactorInt( double normalInt, double dist, bool orientation ) const;

    // the maximum distance beyond which there is no interaction
    double getMaxDist() const {
        return distThr;
    } 
    
private:
   
    double distThr;		// if distance < thr, the "normal" value; otherwise 1 (no interaction)
    double orientationEffect;	// the effect of orientation: if at different strands, the effect should be multiplied this value	
};

class FactorIntFuncGaussian : public FactorIntFunc {
public: 
    // constructors
    FactorIntFuncGaussian( double _distThr, double _sigma ) : distThr( _distThr ), sigma( _sigma ) {
        assert( distThr > 0 && sigma > 0 );
    }

    // compute the factor interaction
    double compFactorInt( double normalInt, double dist, bool orientation ) const; 

    // the maximum distance beyone which there is no interaction
    double getMaxDist() const {
        return distThr;
    } 
private: 
    double distThr;     // no interaction if distance is greater than thr. 
    double sigma;       // standard deviation of 
};

class FactorIntFuncGeometric : public FactorIntFunc {
public:
    // constructors
    FactorIntFuncGeometric( double _distThr, double _spacingEffect, double _orientationEffect ) : distThr( _distThr ), spacingEffect( _spacingEffect ), orientationEffect( _orientationEffect ) { assert( distThr > 0 ); }

    // compute the factor interaction
    double compFactorInt( double normalInt, double dist, bool orientation ) const;

    // the maximum distance beyond which there is no interaction
    double getMaxDist() const {
        return distThr;
    } 
private:
    double distThr;		// if distance < thr, the "normal" value; otherwise decay with distance (by parameter spacingEffect)
    double spacingEffect;		// the effect of spacing
    double orientationEffect;	// the effect of orientation: if at different strands, the effect should be multiplied this value
};


class ExprPar {
public:
    // constructors 
    ExprPar() : factorIntMat(), theV(), theVr() {}  //  is this the constructor called from the initialization list of ExprFunc (..,  par, .) : par(_par)
    ExprPar( int _nFactors );		// default values of parameters
  //  ExprPar( int _nFactors);
ExprPar( int _nFactors, vector< vector< vector<double> > > _theV , vector< vector< vector<double> > > _theVr);
    ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, double _basalTxp );
    ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators );	// construct from a "flat" vector of free parameters (assuming they are in the correct/uniform scale)
    void copy( const ExprPar& other ) { maxBindingWts = other.maxBindingWts; factorIntMat = other.factorIntMat; txpEffects = other.txpEffects; repEffects = other.repEffects; theV = other.theV; theVr = other.theVr;  }
    ExprPar( const ExprPar& other ) { copy( other ); }
	
    // assignment
    ExprPar& operator=( const ExprPar& other ) { copy( other ); return *this; }	
	
    // access methods
    int nFactors() const { return maxBindingWts.size(); }
	
    // get the free parameters (in the correct/uniform scale)
    void getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const; 
	   void getFreePars2( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) ;
    // print the parameters
    void print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const;
    //void print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat) const //, const vector< vector< vector<double > > > _theV ) const;
 		
    // load the parameter values from a file, assuming the parameter has the correct dimensions (and initialized)
    int load( const string& file, IntMatrix& _coopMat , const vector< bool >& _repIndicators ); 

    // adjust the values of parameters: if the value is close to min or max allowed value, slightly change it s.t. it is away from the boundary
    void adjust(); 
    
    // parameters
    vector< double > maxBindingWts;			// binding weight of the strongest site for each TF: K(S_max) [TF_max]
    Matrix factorIntMat; 		// (maximum) interactions between pairs of factors: omega(f,f')
    vector< double > txpEffects;    // transcriptional effects: alpha for Direct and Quenching model, exp(alpha) for Logistic model (so that the same default values can be used). Equal to 1 if a TF is not an activator under the Quenching model
    vector< double > repEffects;    // repression effects: beta under ChrMod models (the equlibrium constant of nucleosome association with chromatin). Equal to 0 if a TF is not a repressor. 
    double basalTxp;        // basal transcription: q_p for Direct and Quenching model, exp(alpha_0) for Logistic model (so that the same default value can be used)
    vector< vector< vector<double > > > theV ;
 vector< vector< vector<double > > > theVr ;
vector< vector< vector<double > > > getV()
{ return theV; }
//     double expRatio; 		// constant factor of measurement to prediction 

    static ModelType modelOption;     // model option
    static SearchType searchOption;    // search option: 0 - unconstrained search; 1 - constrained search
    static int estBindingOption;    // whether to estimate binding parameters
    static int nbins;
    static double default_weight;	// default binding weight
    static double default_interaction;		// default factor interaction
    static double default_effect_Logistic;   // default transcriptional effect under Logistic model
    static double default_effect_Thermo;     // default transcriptional effect under thermo. models
    static double default_repression;   // default repression
    static double default_basal_Logistic;       // default basal transcription under Logistic model
    static double default_basal_Thermo;         // default basal transcriptional under thermo. models
    static double min_weight;		// min. binding weight
    static double max_weight;		// max. binding weight
    static double min_interaction;	    // min. interaction
    static double max_interaction;	    // max. interaction
    static double min_interactionr;	    // min. interaction
    static double max_interactionr;	    // max. interaction
    static double min_effect_Logistic;   // min. transcriptional effect under Logistic model
    static double max_effect_Logistic;   // max. transcriptional effect under Logistic model
//     static double min_effect_Direct;   // min. transcriptional effect under Direct model
    static double min_effect_Thermo;    // min. transcriptional effect under thermo. models
    static double max_effect_Thermo;   // max. transcriptional effect under thermo. models   
    static double min_repression;   // min. repression
    static double max_repression;   // max. repression
    static double min_basal_Logistic;    // min. basal transcription under Logistic model
    static double max_basal_Logistic;    // max. basal transcription under Logistic model
    static double min_basal_Thermo;   // min. basal transcription under thermo. models
    static double max_basal_Thermo;   // max. basal transcription under thermo. models
    static double delta;        // a small number for testing the range of parameter values
// 	static double wt_step;		// step of maxExprWt (log10)
// 	static double int_step;		// step of interaction (log10)
// 	static double ratio_step;	// step of expRatio
};


//int ExprPar::bn =5;
class ExprFunc {
public:
    // constructors
 // ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par  ,const vector< int >& _B );

  //  ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par );
ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par,  const vector< int >& _B ,  const vector< int >& _Br);

//void predictOcc( const SiteVec& _sites, int length, const vector< double >& factorConcs, vector< double >& factorOcc );  // k is the current sequence, followed the predictExpr template

void predictOcc( const SiteVec& _sites, int length, const vector< double >& factorConcs, vector< double >& fOcc );

void predictTwist( vector< vector< vector< Site > > >& _d );
// sort must be called before this is invoked.
double predictZ( const SiteVec& _sites, const vector< double >& factorConcs );
    // access methods
    const vector< Motif >& getMotifs() const {
        return motifs;
    }
  //  void printpar( ostream& os);
    // predict the expression value of a given sequence (its site representation, sorted by the start positions) under given TF concentrations
    double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs );
   static vector<int > Bc;  // bin counter  Bc = B10
static vector<int > Bct;  // bin counter  Bct = B02
static vector<int > B12;  // bin counter  1= DU, 2 = Tw
static vector<int > B03;  // bin counter  0 = DC, 3 = CB
static vector<int > B13;  // bin counter  1 = DU, 3 = CB
static vector<int > B23;  // bin counter  2 = tw, 3 = CB
    static ModelType modelOption;     // model option   
void setsites( SiteVec& s) { sites = s ;}
private:
    // TF binding motifs
    const vector< Motif >& motifs; 
	 
  //  double coopDistThr;
   // int nbins; 
    // control parameters
    const FactorIntFunc* intFunc;   // function to compute distance-dependent TF-TF interactions     
    const vector< bool >& actIndicators;    // 1 if the TF is in the activator set
    int maxContact;     // the maximum contact     
    const vector< bool >& repIndicators;    // 1 if the TF is in the repressor set
    const IntMatrix& repressionMat;    // repression matrix: R(f,f') = 1 if f can repress f'
    double repressionDistThr;   // distance threshold for repression: d_R
   const vector<int >& B;     // bin borders	
const vector<int >& Br;
	///////////////////////
	vector< double > Z;
    vector< double > Zt;
				
	// compute the partition function 
	double compPartFunc();
vector< double > factorOcc;
///////////////////////
    // model parameters
    const ExprPar& par;
		    
    // the sequence whose expression is to be predicted
    SiteVec sites;
    vector< int > boundaries;   // left boundary of each site beyond which there is no interaction    

    // intermediate computational results
    vector< double > bindingWts; 
		
   
   
    
    // compute the TF-TF interaction between two occupied sites
    double compFactorInt( const Site& a, const Site& b ) const;

    // test if one site represses another site
    bool testRepression( const Site& a, const Site& b ) const;
};
class SeqAnnotator   {
public:
    // constructors
 SeqAnnotator( const vector< Motif >& _motifs, const vector< double >& _energyThrs , Motif& _d2) : motifs( _motifs ), energyThrs( _energyThrs ), d2( _d2), B(),  seqSitesm1delete1(  ) { assert( motifs.size() == energyThrs.size() ); }

    SeqAnnotator( const vector< Motif >& _motifs, const vector< double >& _energyThrs ) : motifs( _motifs ), energyThrs( _energyThrs ), B(),  seqSitesm1delete1(  ) { assert( motifs.size() == energyThrs.size() ); }
         SeqAnnotator( const vector< Motif >& _motifs, const vector< double >& _energyThrs ,const  vector<int> _B ) : motifs( _motifs ), energyThrs( _energyThrs ), B(_B),seqSitesm1delete1( ) { assert( motifs.size() == energyThrs.size() ); }   
    // annotate a sequence
    int annot( const Sequence& seq, SiteVec& sites ) const;	
	int annotI( const Sequence& seq, SiteVec& sites ) const;

    int annoty( const Sequence& seq, SiteVec& sites, Matrix f, Matrix e, vector< double > p) ;

int annoty2( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ; //classifies meso and neuro, in meso uses (n,2)dt and dd and compares which Z is higher..
int annoty3( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites) ;//takes annotated sitevec from either annoty2 or annoty4 and uses those sites to add snail sites for quenching.
int annotyd( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites, vector< SiteVec >& dd) ;//contains deletion structure for robustness analysis.
int annoty4( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ; // checks for multiple copy of maxsites.
////
int annotydorsalZ(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ;

////
//int annotydorsalZ2(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ;
int annotydorsalZ2( vector< Sequence > seq,vector< vector< vector< Site > > >& sitespp , Matrix f, Matrix e, ExprFunc& p, int column, vector< string >  names) ;
int annotydorsalZarend( vector< Sequence > seq, vector< vector< vector< Site > > >& sitespp , Matrix f, Matrix e, ExprFunc& p, int column, vector< string >  names) ;

int annotydorsal(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ; // combines annoty4 and annoty2
// so now annotydorsal passes sitevec to annoty3
int annotydorsalold(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) ;
ExprPar SeqPar;
int maxZ( vector< double >& Z );
double maxZs( vector< double >& Z );
Site siteMax( vector< Site >& sites);
Site siteMaxI( vector< Site >& sitess);
bool DeltaZ( vector< vector< Site > >& sitespa , vector< double >& cons , double Zbefore, ExprFunc& p, double Obefore);
// delete all sites within sitest that overlap with sites from sitesp
void sitestoverlap( vector< Site >& sitest, vector< Site >& sitesp );
//void annotprint();
    // compute the energy of sites (update the input sites)
    int compEnergy( const Sequence& seq, SiteVec& sites ) const;
int Delete1( vector< Site > &t ,vector< vector< Site > >& tt);
Motif d2;
void killmo(){ motifs.clear();}
void addmo( Motif a ) { motifs.push_back( a );  cout << motifs << endl;}
private:
    vector< Motif > motifs;	// all the TF binding motifs
    vector< double > energyThrs;	// energy thresholds for all the motifs
 const  vector<int> B;
vector< vector< Site > > seqSitesm1delete1;
};

class ExprPredictor {
public:
    // constructors
    ExprPredictor( const vector< vector< SiteVec > >& _seqSitesb, const vector< vector< double > >& _bindingData,  vector< SiteVec > & _seqSitesa, vector< SiteVec > & _seqSites, const vector< int >& _seqLengths,  Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const FactorIntFunc* _intFunc, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const vector<int>& _binBord, const vector<int>& _binBordr , const vector < bool >& _indicator_bool, SeqAnnotator _anny, const string& _file ,  vector< SiteVec > & _seqSitesbot ,  vector< SiteVec >& _seqSitesm1,  vector< SiteVec >& _seqSitesm2, vector< SiteVec >& _seqSitesf2 ,vector< SiteVec >& _seqSitesbotf2, vector< SiteVec >& _seqSitesm1f2, vector< SiteVec >& _seqSitesm2f2, vector< SiteVec >& _seqSitesf3 ,  vector< SiteVec >& _seqSitesbotf3,vector< SiteVec >& _seqSitesm1f3 , vector< SiteVec >& _seqSitesm2f3,  vector< vector< vector< Site > > >& _d);
///////////////////////////////////////729
//Motif() : pwm(), background() {}
//ExprPredictor()  {}
    // access methods 
    int nSeqs() const {
        return seqSites.size();
    }
	int nSeqsa() const {
        return seqSitesa.size();
    }
////////////////////////////////////////
	int nSeqsb() const {
        return seqSitesb.size();
    }
    int nFactors() const { 
        return motifs.size(); 
    }
    int nConds() const {
         return exprData.nCols();
    }
    const IntMatrix& getCoopMat() const {
        return coopMat;
    }
    const vector< bool >& getActIndicators() const {
        return actIndicators;
    }
    const vector< bool >& getRepIndicators() const {
        return repIndicators;
    }
    const IntMatrix& getRepressionMat() const {  // the  first "const" is the return type's quantity, the second const refers to the function
        return repressionMat;
    }
    const ExprPar& getPar() const { return par_model; }
    double getObj() const { return obj_model; }
    void printPar2(  );
const string& file;

    // the objective function to be minimized
    double objFunc( const ExprPar& par ) const;
    double objFunc2(  ExprPar& par ) ;
	double objFuncborder(  ExprPar& par ) ;
double objFuncborder2(  ExprPar& par );
    // training the model
    int train( const ExprPar& par_init ); 	// training with the initial values given
    int train( const ExprPar& par_init, const gsl_rng* rng );   // training with the initial values and allowing random starts
    int train();	// automatic training: first estimate the initial values, then train

 int train3( const ExprPar& par_init, const gsl_rng* rng );   // training with the initial values and allowing random starts
/////////////////////////////////
      SeqAnnotator anny;      
/////////////////////////////////
    // predict expression values of a sequence (across the same conditions)
    int predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs ) ; 
// 11-7-2011 removed the const from predict
 void predict2( vector< vector< double > >& occs ) ; 
    // test the model, perfOption = 0: RMSE
// 	double test( const vector< Sequence  >& testSeqs, const Matrix& testExprData, Matrix& predictions ) const;    
void compOccMat(  const gsl_vector* x, void* params ) ;//was par
void printFile( ostream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
void printFile2( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
void printFiled( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
void printFilePar_KfoldCV( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;

int createccdata();
int createccdata2();
int load( const string& fila );  
 // void printFile( ostream& os, ExprPredictor jRMSE ) const;
    //void printFile( ostream& os, ExprPredictor jRMSE,  const vector< string >& motifNames, const IntMatrix& coopMat ) const
    static ModelType modelOption;     // model option
    static int estBindingOption;    // whether estimate binding parameters
    static ObjType objOption;       // option of the objective function
///////////////////////////////////////////
static vector< Sequence > seqsy;
static vector< Sequence > seqsya;
static vector< string > seqNmes;
static vector< string > seqNmesa;
//////////////////////////////////////////
 vector< int > cell ;
    // the similarity between two expression patterns, using cross-correlation
    static double exprSimCrossCorr( const vector< double >& x, const vector< double >& y ); 
    static int maxShift;    // maximum shift when computing cross correlation
    static double shiftPenalty;     // the penalty for shift (when weighting different positions)

    // the parameters for the optimizer
    static int nAlternations;   // number of alternations (between two optimization methods)
    static int nRandStarts;     // number of random starts
    static double min_delta_f_SSE;      // the minimum change of the objective function under SSE
    static double min_delta_f_Corr;     // the minimum change of the objective function under correlation
    static double min_delta_f_CrossCorr;    // the minimum change of the objective function under cross correlation
    static int nSimplexIters;       // maximum number of iterations for Simplex optimizer
    static int nGradientIters;      // maximum number of iterations for Gradient optimizer
   // static int nbins; 
void printBc( ) const;    
const vector< bool >& repIndicators;    
    // print the parameter values (the ones that are estimated) in a single line
    void printPar( const ExprPar& par ) const;
 vector < bool > indicator_bool;
   vector < double > fix_pars;
    vector < double > free_pars;
////////////////////////////////////////////7/7/2012 holds the sites for all enhancer for Adami
vector< SiteVec >& seqSitesa;
////////////////////////////////////////////
 vector< SiteVec >& seqSites;		// the extracted sites for all sequences
 vector< SiteVec >& seqSitesbot;
vector< SiteVec >& seqSitesm1;
vector< SiteVec >& seqSitesm2;

//////////////dorsal + noise
 vector< SiteVec >& seqSitesf2;		// the extracted sites for all sequences
 vector< SiteVec >& seqSitesbotf2;
vector< SiteVec >& seqSitesm1f2;
vector< SiteVec >& seqSitesm2f2;
///////////////dorsal - noise
 vector< SiteVec >& seqSitesf3;		// the extracted sites for all sequences
 vector< SiteVec >& seqSitesbotf3;
vector< SiteVec >& seqSitesm1f3;
vector< SiteVec >& seqSitesm2f3;

int clasvar;
static Matrix exprData2;		// expressions of the corresponding sequences across multiple conditions         
    ExprFunc* createExprFunc2(  ExprPar& par ) ;
    static Matrix factorExprData2;
void compvar(vector< double >& vars);
vector< vector< Site > > seqSitesm1d1;
 vector< vector< vector< Site > > >& d;
 double compAvgCorr2( ExprPar& par  , int thres) ;
double compAvgCorr3( ExprPar& par  , int thres) ;
void printFile22( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
void printFile23( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
void printFile24( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
void printFile25( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const;
private:
 //const vector< SiteVec >& seqSites;		// the extracted sites for alntl sequences
    // training data
    const vector< int >& seqLengths;           // lengths of all sequences
  Matrix& exprData;		// expressions of the corresponding sequences across multiple conditions         
    const vector< Motif >& motifs;		// TF binding motifs
   const Matrix& factorExprData;		// [TF] of all factors over multiple conditions	    
//vector< vector< Sequence > > seqs2;
    const vector< vector< SiteVec > >& seqSitesb;		// the extracted sites for all sequences
   // const vector< int >& seqLengths;           // lengths of all sequences
  //  const Matrix& bData;		// expressions of the corresponding sequences across multiple conditions   
const vector< vector< double > >& bindingData;		// binding intensities of sequences: one row per experiment (in-factor)
    // control parameters 
    const FactorIntFunc* intFunc;   // function to compute distance-dependent TF-TF interactions   
    const IntMatrix& coopMat;       // cooperativity matrix: C(f,f') = 1 if f and f' bind cooperatively    
    const vector< bool >& actIndicators;   // 1 if the TF is in the activator set
    int maxContact;     // the maximum contact     
    
    // 1 if the TF is in the repressor set
    const IntMatrix& repressionMat;    // repression matrix: R(f,f') = 1 if f can repress f'
    double repressionDistThr;   // distance threshold for repression: d_R
    
    // model parameters and the value of the objective function
    ExprPar par_model;
    double obj_model;	
    double obj2;
  const  vector<int>& binBord;
const  vector<int>& binBordr;
   
    // randomly sample parameter values (only those free parameters), the parameters should be initialized
    int randSamplePar( const gsl_rng* rng, ExprPar& par ) const; 

    // check if some parameter combination is valid
    bool testPar( const ExprPar& par ) const; 

  // create the expression function
    ExprFunc* createExprFunc( const ExprPar& par ) const;
 //  double coopDistThr;
    
    // objective functions
    double compRMSE( const ExprPar& par ) const;		// root mean square error between predicted and observed expressions
    double compAvgCorr( const ExprPar& par ) const;     // the average Pearson correlation
// double compAvgCorr2( ExprPar& par ) ;
double compAvgCorrborder( ExprPar& par ) ;
    double compAvgCrossCorr( const ExprPar& par ) const;    // the average cross correlation -based similarity
double compAvgCorrborder2(  ExprPar& par ) ;
    // minimize the objective function, using the current model parameters as initial values
    int simplex_minimize( ExprPar& par_result, double& obj_result ) ;//const;	// simplex	
    int gradient_minimize( ExprPar& par_result, double& obj_result ) ;// const;	// gradient: BFGS or conjugate gradient
//  	int SA_minimize( ExprPar& par_result, double& obj_result ) const;	// simulated annealing 		
};

// the objective function and its gradient of ExprPredictor::simplex_minimize or gradient_minimize
double gsl_obj_f( const gsl_vector* v, void* params );
void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad ); 
void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad ); 

#endif
