#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include "ExprPredictor.h"

#include <sstream>

bool isNt( int a )
{
    if ( a < 0 || a > 3 ) return false;
    else return true;	
}

int complement( int a )
{
    assert( a >= 0 && a < ALPHABET_SIZE );
            
    if ( a == 0 ) return 3;
    if ( a == 1 ) return 2;
    if ( a == 2 ) return 1;
    if ( a == 3 ) return 0;	
    if ( a == MISSING ) return MISSING;
    if ( a == GAP ) return GAP;	
}

int symbolToInt( char c )  
{
    char lower = toupper( c );
    char upper = tolower( c );
    for ( int i = 0; i < ALPHABET_SIZE; i++ ) {
        if ( ALPHABET[ i ] == upper || ALPHABET[i] == lower ) return i;	
    }
    
    return -1;
}

char strand2char( bool strand )
{
    if ( strand ) return '+';
    else return '-';	
}

bool char2strand( char c )
{
    assert( c == '+' || c == '-' );
    
    if ( c == '+' ) return true;
    else return false;
}

vector< double > createNtDistr( double gcContent )
{
    assert( gcContent >= 0 && gcContent <= 1.0 );
    
    vector< double > freqs( 4 );
    freqs[0] = ( 1.0 - gcContent ) / 2.0;
    freqs[1] = gcContent / 2.0;
    freqs[2] = freqs[1];
    freqs[3] = freqs[0];

    return freqs;
}

Sequence::Sequence( const string& str )
{
    for ( int i = 0; i < str.size(); i++ ) {
        int nt = symbolToInt( str[ i ] );	
        if ( nt >= 0 && nt < ALPHABET_SIZE ) {
            nts.push_back( nt );
        } else {
            cerr << "Illegal symbol: " << nt << " in " << str << endl;
            exit( 0 );	
        }       
    }
}

Sequence::Sequence( const Sequence& other, int start, int length, bool strand )
{
    assert( start >= 0 && length >= 0 && ( start + length ) <= other.size() );	

    for ( int i = 0; i < length; i++ ) {
        if ( strand ) {	nts.push_back( other[ start + i ] ); }
        else { nts.push_back( complement( other[ start + length - 1 - i ] ) ); }
    }	
}

int Sequence::push_back( int nt )
{	
    assert( nt >= 0 && nt < ALPHABET_SIZE );
    nts.push_back( nt );
    
    return 0;
}

int Sequence::push_back( const Sequence& elem )
{
    for ( int i = 0; i < elem.size(); i++ ) push_back( elem[ i ] );	
    return 0;
}

Sequence Sequence::compRevCompl() const
{
    return Sequence( *this, 0, size(), false );	
}

void Sequence::getNtCounts( vector< int >& counts ) const
{
    counts.clear();
    for ( int i = 0; i < NBASES; i++ ) {
        counts.push_back( 0 );
    }
    
    for ( int i = 0; i < nts.size(); i++ ) {
        if ( nts[ i ] != GAP ) counts[ nts[ i ] ]++;	
    }
}

bool Sequence::containsMissing() const
{
    for ( int i = 0; i < nts.size(); i++ ) {
        if ( nts[ i ] == MISSING ) return true;	
    }	
    
    return false;
}

int Sequence::load( const string& file, string& name, int format )
{
    vector< Sequence > seqs;
    vector< string > names;
    int rval = readSequences( file, seqs, names, format );
    if ( rval == RET_ERROR ) return RET_ERROR;
    
    copy( seqs[ 0 ] );
    name = names[ 0 ];
    return rval;
}

int Sequence::load( const string& file, int format )
{
    string name;
    int rval = load( file, name, format );
    
    return rval;	
}

ostream& operator<<( ostream& os, const Sequence& seq )
{
    
    for ( int i = 0; i < seq.size(); i++ ) {
        os << ALPHABET[ seq[ i ] ];	
    }	
                    
    return os;
}

int readSequences( const string& file, vector< Sequence >& seqs, vector< string >& names, int format )
{
    
    if ( format != FASTA ) { return RET_ERROR; }
    seqs.clear();
    names.clear();
     
    
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }

    string line;
    Sequence seq;
    
    
    if ( format == FASTA ) {
        while ( getline( fin, line ) ) {
            
            
            if ( line[ 0 ] == '>' ) { 	
                if ( seq.size() ) {
                    seqs.push_back( seq );
                    seq.clear();	
                }
                        
                stringstream ss( line.substr( 1 ) );
                string name; 
                ss >> name;
                names.push_back( name );
            } else { 
                
                int start = line.find_first_not_of( " \t\r" );
                int last = line.find_last_not_of( " \t\r" );
                if ( start == string::npos || last == string::npos ) continue;
                        
                
                for ( int i = start; i <= last; i++ ) {
                    int nt = symbolToInt( line[ i ] );	
		    if( i == start){ for(int j=0; j < 10; j++){seq.push_back(4);} }
		    
                    if ( nt >= 0 && nt < ALPHABET_SIZE ) {
                        seq.push_back( nt );
                    } 
			else {
			
                        cerr << "Illegal symbol: " << nt << " in " << file << endl;
                        return RET_ERROR;	
                    } 
		    if( i == last){ for(int j=0; j < 10; j++){seq.push_back(4);} }
		  
                }
            }			
        }
            
        
        if( seq.size() ) seqs.push_back( seq );
                        
        return 0;
    }	
}
int readSequences( const string& file, vector< Sequence >& seqs, int format )
{
    vector< string > names;
    int rval = readSequences( file, seqs, names, format );	
    return rval;
}

int writeSequences( const string& file, const vector< Sequence >& seqs, const vector< string >& names, int format )
{
    assert( seqs.size() == names.size() );
    
    
    if ( format != FASTA ) { return RET_ERROR; }
            
    ofstream fout( file.c_str() );
    
    if ( format == FASTA ) {
        for ( int i = 0; i < seqs.size(); i++ ) {
            fout << ">" << names[ i ] << endl;
            fout << seqs[ i ] << endl;
        }
    }
    
    return 0;
}

int writeSequences( const string& file, const vector< Sequence >& seqs, int format )
{
    
    vector< string > names;
    for ( int i = 0; i < seqs.size(); i++ ) {
        char buffer[ 10 ];
        sprintf( buffer, "%i", i );
        names.push_back( string( buffer ) );	
    }	
    
    
    return writeSequences( file, seqs, names, format );
}

Matrix compWtmx( const Matrix& countMatrix, double pseudoCount )
{
    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
    
    int l = countMatrix.nRows();		
    Matrix pwm( l, 4 );

    
    
    for ( int i = 0; i < l; i++ ) {
        double n = 0;       
        for ( int j = 0; j < 4; j++ ) {
            n += countMatrix( i, j );
        }
        for ( int j = 0; j < 4; j++ ) {
            pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
        }	
    }

    return pwm;		
}


Motif::Motif( const Matrix& _pwm, const vector< double >& _background ) : pwm( _pwm ), background( _background ), LLRMat( pwm.nRows(), 4 ), eth( 0 ), infomat(pwm.nRows(),4), countmat( pwm.nRows(), 4  )
{
    assert( background.size() == 4 );	
    
    init();
}

Motif::Motif( const Matrix& _pwm, Motif& _motif ): pwm( _pwm), LLRMat( _pwm.nRows(), 4 ), eth( 0 ), infomat(_pwm.nRows(),4), countmat( _pwm.nRows(), 4  )
{
    
    
    
 int l = _motif.pwm.nRows();
    Matrix LLRcb=_motif.getLLRMatpos();
	Matrix pwmCB =_motif.getPwm();
	minLLR =_motif.getMinLLR();
	maxLLR =_motif.getMaxLLR();
    
				
    for ( int i = 0; i < l; i++ ) {
        for ( int j = 0; j < 4; j++ ) {			
            LLRMat( i, j ) = -(-maxLLR/double(l)  +LLRcb(i,j)+  log( pwmCB( i, j )/_pwm(i,j) )  );  
        }
    }

  


    
    for ( int i = 0; i < l; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
            infomat(i,j)=  log( 4*_pwm(i,j));	
        }	
    }

}

Motif::Motif( const Matrix& countMatrix, double pseudoCount, const vector< double >& _background ) : background( _background ), LLRMat( countMatrix.nRows(), 4 ), eth(0), infomat(countMatrix.nRows(),4), countmat( countMatrix )
{ 
    assert( background.size() == 4 );
    
    pwm = compWtmx( countMatrix, pseudoCount );
    init();
	
}



double Motif::LLR( const Sequence& elem ) const
{
    int l = pwm.nRows();
    if ( elem.size() != l ) return GSL_NEGINF;
    if ( elem.containsMissing() ) return GSL_NEGINF;
    
    double result = 0;
    for ( int i = 0; i < l; i++ ) {
        result += LLRMat( i, elem[ i ] ); 
    }
    
    return result;
}
double Motif::LLR2(  Sequence& elem ) 
{
    int l = pwm.nRows();
  
    if ( elem.containsMissing() ) return GSL_NEGINF;
    
    double result = 0;
if ( elem.size() != l ) {
    for ( int i = 0; i < l; i++ ) {
	 if ( i < elem.size() )
        result += LLRMat( i, elem[ i ] ); 
    }
    }
	else{ for ( int i = 0; i < l; i++ ) {
		 
		result += LLRMat( i, elem[ i ] ); 
	    } 
	}
    return result;
}
Matrix Motif::getInfoMat()
{
Matrix enepwm = getPwm();
Matrix e(enepwm.nRows(),4,0);
double Expected_info=0;
    
    for ( int i = 0; i < enepwm.nRows(); i++ ) {
        for ( int j = 0; j < enepwm.nCols(); j++ ) {
            e(i,j)=  log( 4*enepwm.getElement(i,j));	
        }	
    }

return  e;
}

Matrix Motif::getKLMat()
{
Matrix enepwm = getPwm();
Matrix e(enepwm.nRows(),4,0);
double Expected_info=0;
    
    for ( int i = 0; i < enepwm.nRows(); i++ ) {
        for ( int j = 0; j < enepwm.nCols(); j++ ) {
            e(i,j)= enepwm.getElement( i, j ) * log( enepwm.getElement(i,j)/.25 );	
        }	
    }

return  e;
}
double Motif::getInformation()
{
Matrix enepwm = getPwm();
double Expected_info=0;
    
    for ( int i = 0; i < enepwm.nRows(); i++ ) {
        for ( int j = 0; j < enepwm.nCols(); j++ ) {
            Expected_info = Expected_info  - enepwm.getElement( i, j ) * log( enepwm.getElement(i,j)/.25 );	
        }	
    }


return Expected_info;
}

double Motif::getInformation2()
{
Matrix _pwm = getPwm();
Matrix Ipwm = getInfoMat();
double Expected_info=0;
    
    for ( int i = 0; i < _pwm.nRows(); i++ ) {
        for ( int j = 0; j < _pwm.nCols(); j++ ) {
            Expected_info = Expected_info  - _pwm.getElement( i, j ) * Ipwm.getElement(i,j) ;	
        }	
    }


return Expected_info;
}

Matrix Motif::getLLRMatpos()
{
Matrix m = getLLRMat();
vector< int > maxSite22;
int l = pwm.nRows();
    for ( int i = 0; i < l; i++ ) {
        int b_max;
        max( pwm.getRow( i ), b_max );
        maxSite22.push_back( b_max );	
    }
    
    
    vector< double >  maxLLR2;
    for ( int i = 0; i < l; i++ ) {
        maxLLR2.push_back( LLRMat( i, maxSite22[ i ] ) );	
    }

for ( int i = 0; i < LLRMat.nRows(); i++ ) {
        for ( int j = 0; j < LLRMat.nCols(); j++ ) {
            m(i,j) = maxLLR2[i]  - LLRMat(i,j);
		
		
		
        }	
    }
return m;
}

double Motif::getRMSE()  
{
Matrix thisA = getLLRMatpos();
Matrix enepwmA = getPwm();
double RMSE=0;
double n =0;

vector< double >meancol;
    
    for ( int i = 0; i < thisA.nRows(); i++ ) {
		double m=0;
                for ( int j = 0; j < thisA.nCols(); j++ ) {
		m = m + thisA(i,j)*enepwmA(i,j);
		}
		meancol.push_back( m );
	    }
vector< double >scol;
   for ( int i = 0; i < enepwmA.nRows(); i++ ) {
	double s=0;
        for ( int j = 0; j < enepwmA.nCols(); j++ ) {
	n=n+1;
            RMSE = RMSE + enepwmA(i,j) * pow(thisA(i,j) - meancol[i],2);
		s= s +	enepwmA(i,j) * pow(thisA(i,j) - meancol[i],2);
        }	
	scol.push_back(s );
    }
ofstream j("ja");
j<< " mean first row, sigma next row" << endl;
j<< meancol << endl;
j<< scol << endl;
j.close();
	
return sqrt(RMSE );

}
double Motif::getRMSEdc( string& file )  
{
Matrix thisA = getLLRMatpos();
Matrix enepwmA = getPwm();
double RMSE=0;
double n =0;

vector< double >meancol;
    
    for ( int i = 0; i < thisA.nRows(); i++ ) {
		double m=0;
                for ( int j = 0; j < thisA.nCols(); j++ ) {
		m = m + thisA(i,j)*enepwmA(i,j);
		}
		meancol.push_back( m );
	    }
vector< double >scol;
   for ( int i = 0; i < enepwmA.nRows(); i++ ) {
	double s=0;
        for ( int j = 0; j < enepwmA.nCols(); j++ ) {
	n=n+1;
            RMSE = RMSE + enepwmA(i,j) * pow(thisA(i,j) - meancol[i],2);
		s= s +	enepwmA(i,j) * pow(thisA(i,j) - meancol[i],2);
        }	
	scol.push_back(s );
    }
ofstream j( file.c_str() );
j<< " mean first row, sigma next row" << endl;
j<< meancol << endl;
j<< scol << endl;
j.close();
	
return sqrt(RMSE );

}

double Motif::getMean()
{
Matrix eneLLR =getLLRMat();
Matrix enepwm = getPwm();
   
vector< int > maxSite22;
int l = enepwm.nRows();
    for ( int i = 0; i < l; i++ ) {
        int b_max;
        max( pwm.getRow( i ), b_max );
        maxSite22.push_back( b_max );	
    }
    
    
    vector< double >  maxLLR2;
    for ( int i = 0; i < l; i++ ) {
        maxLLR2.push_back( LLRMat( i, maxSite22[ i ] ) );	
    }

double Expected_energy=0;
    
    for ( int i = 0; i < eneLLR.nRows(); i++ ) {
        for ( int j = 0; j < eneLLR.nCols(); j++ ) {
            Expected_energy = Expected_energy + maxLLR2[i]*enepwm.getElement(i,j) - eneLLR.getElement( i, j ) * enepwm.getElement(i,j);
		
		
		
        }	
    }


return Expected_energy;
}
double Motif::getgenomicMean()
{
Matrix eneLLR =getLLRMat();
double Expected_energy=0;
    
    for ( int i = 0; i < eneLLR.nRows(); i++ ) {
        for ( int j = 0; j < eneLLR.nCols(); j++ ) {
            Expected_energy = Expected_energy  - eneLLR.getElement( i, j ) * background[j] ;	
        }	
    }

Expected_energy = Expected_energy + getMaxLLR();  

return Expected_energy;
}

int SeqAnnotator::annotscan3( const Sequence& seq,  vector< int >& en  ,  vector< int >& en2  ) 
{
    
    ofstream f("f.t",ios::app);
    ofstream ff("ff.t",ios::app);
    
	int i = 0;
    
	
    for ( int i = 0; i < seq.size(); i++ ) {
	
        
	int binnum = 0;
        for ( int k = 0; k < 2; k++ ) {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) continue;
            double energy;
            
            
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
	   
            
binnum = ceil( energy );
	if( k == 0 ){
         en[ binnum ]= en[ binnum ] + 1;
	f << energy << '\t' ;
}

	if(k == 1 ){
	  en2[binnum ] = en2[binnum ] +1;
	ff << energy << '\t' ;
}
            
           
           
               
	
         
        }	

    }
    f.close();
ff.close();
    return en.size();
}



double Motif::getRMSE_X_djor()  
{
double RMSE=0;
vector< double > ea ; 
double etemp=0;
   for ( int i = 0; i < pwm.nRows(); i++ ) {
	etemp=0;
        for ( int j = 0; j < pwm.nCols(); j++ ) {
	etemp = etemp + pwm(i,j)*LLRMat(i,j);	    
	}
	ea.push_back(etemp);
}
vector< double > ea2 ; 
double etemp2=0;
   for ( int i = 0; i < pwm.nRows(); i++ ) {
	etemp2=0;
        for ( int j = 0; j < pwm.nCols(); j++ ) {
	etemp2 = etemp2 - LLRMat(i,j);	    
	}
	ea2.push_back(etemp2/4);
}
double sum=0;
double sum2=0;
   for ( int i = 0; i < pwm.nRows(); i++ ) {
        for ( int j = 0; j < pwm.nCols(); j++ ) {
		
  
		RMSE = RMSE + pwm(i,j)*pow(LLRMat(i,j) - ea[i],2);
        }	
    }
	RMSE = sqrt(RMSE);
vector< int > maxSite2;
   for ( int i = 0; i < pwm.nRows(); i++ ) {
        int b_max;
        max( pwm.getRow( i ), b_max );
        maxSite2.push_back( b_max );	
    }
    
    
    vector< double >  maxLLR2;
    for ( int i = 0; i < pwm.nRows(); i++ ) {
        maxLLR2.push_back( LLRMat( i, maxSite2[ i ] ) );	
    }

for ( int i = 0; i < pwm.nRows(); i++ ) {
		sum = sum - maxLLR2[i] -ea2[i];
		sum2 = sum2  - maxLLR2[i] + ea[i];
}


return RMSE;


}

 void align2( string& file, string& file2, int nrandStarts, int nUpdates , int wid, int l   )
{
gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		

int  initp = 1; 
string  ipfile;


string seqFile, dcFile;
seqFile=file;
vector< Sequence > seqs;
    vector< string > seqNames;
 int rval;
    rval = readSequences( seqFile, seqs, seqNames );  
    assert( rval != RET_ERROR );
    int nSeqs = seqs.size();
SiteVec seqSitesP;
 double gcContent = .5;

vector< Motif > motifs2;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );
double energyThr =2;
vector< double > energyThrs( 1, energyThr );
	
double pwm_best = 1;
double pwm_model =2.1;
double pwmrmse = 1.8;
int pin;
vector< int > pos(nSeqs,initp);  
    if ( !ipfile.empty() ) {
        ifstream fcoop( ipfile.c_str() );
        if ( !fcoop ) {
            cerr << "Cannot open the ip file " << endl;
            exit( 1 );
        }  
       for(int i = 0 ; i < nSeqs ; i++){
         fcoop >> pin;
	 pos[i] = pin;
        }        
    } 


vector< int > posbest(nSeqs,initp);

	ofstream core("corefile.txt");
		for( int jj = 0; jj < nSeqs ; jj++) {
			
			if(pos[jj] + l > seqs[jj].size() ) { pos[jj] = 1; }   
			
			
			core << ">" << seqNames[jj] << endl;
		    
			for( int ii = pos[jj]; ii < pos[jj] +l ; ii++) {  
			core << ALPHABET[ seqs[jj][ii] ];
			}
			core << endl;
			
		}
	core.close();
	dcFile = "corefile.txt";


	Matrix dcm = countmatrix( dcFile );
	
	Motif dcmm( dcm, background);
	
	
	
for ( int k = 0; k < nUpdates; k++ ) {
 for ( int j = 0; j < seqs.size(); j++ ) {
	ofstream core("corefile.txt");
	for( int jj = 0; jj < nSeqs ; jj++) {
			if(j == jj) { continue;}  
			if(pos[jj] + l > seqs[jj].size() ) { pos[jj] = 1; }   
			
			
			core << ">" << seqNames[jj] << endl;
		    
			for( int ii = pos[jj]; ii < pos[jj] +l ; ii++) {  
			core << ALPHABET[ seqs[jj][ii] ];
			}
			core << endl;
			
	}
	core.close();
	dcFile = "corefile.txt";
	vector< Motif > motifs2;

	Matrix dcm = countmatrix( dcFile );
	Motif dcmm( dcm, background);
	motifs2.push_back( dcmm ); 
 
        SeqAnnotator ann( motifs2, energyThrs );
	SiteVec seqSites;
	 ann.annot( seqs[ j ], seqSites );
	
	 Site m = ann.siteMax( seqSites );
	
	int p = m.start;
	pos[j] = p ;
	
 

  }  
 
	ofstream core2("corefile2.txt");
		for( int jj = 0; jj < nSeqs ; jj++) {
			
			
			
			int start = pos[jj];
			if(pos[jj] + l > seqs[jj].size() ) { start = 1; } 
			
			core2 << ">" <<  seqNames[jj] << endl;
			for( int ii = 0; ii < l ; ii++) {  
			 int p = start + ii;
		
			core2 << ALPHABET[ seqs[jj][p] ];
			}
			core2 << endl;
			
		}
	core2.close();
	dcFile = "corefile2.txt";
	Matrix dcm = countmatrix( dcFile );
	
	Motif dcmm( dcm, background);
	
	
	
	pwmrmse = pwm_model;  
	
	pwm_model = dcmm.getInformation();  
	
		
		if ( pwm_model > pwm_best ) {
			pwm_best = pwm_model;	
			ofstream optcore("coreOPT.txt");	
			ofstream ipf("initialPOS.txt");
			for( int jj = 0; jj < nSeqs ; jj++) {
				int start = pos[jj];
				optcore << ">" << start << seqNames[jj] << endl;
				ipf << start << '\t' ;
				posbest[jj] = start;
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start + ii;
				optcore << ALPHABET[ seqs[jj][p] ];
				}
			optcore << endl;
			}
			ipf << endl;
			optcore.close();
			ipf.close();
			
		} 
	
	if( abs(pwm_model - pwmrmse) < .01 ) break;

}

string opt = "coreOPT.txt";
for ( int r = 0; r < nrandStarts; r++){
	for(int i = 0; i < nSeqs; i++){
	
	pos[i]=  gsl_rng_uniform_int(rng,wid); 
	}
	
	for ( int k = 0; k < nUpdates; k++ ) {
	
	 for ( int j = 0; j < seqs.size(); j++ ) {
		ofstream core("corefile.txt"); 
			for( int jj = 0; jj < nSeqs ; jj++) {
				if(j == jj) { continue;}
				if(pos[jj] + l > seqs[jj].size() ) { pos[jj] = 1; }   
				
				core << ">" << seqNames[jj] << endl;
				for( int ii = pos[jj]; ii < pos[jj] +l ; ii++) {  
				core << ALPHABET[ seqs[jj][ii] ];
				}
				core << endl;
				
			}
		core.close();
		dcFile = "corefile.txt";
		vector< Motif > motifs2; 
		Matrix dcm = countmatrix( dcFile );
		Motif dcmm( dcm, background);
		motifs2.push_back( dcmm ); 
	 
		SeqAnnotator ann( motifs2, energyThrs );
		SiteVec seqSites;
		 ann.annot( seqs[ j ], seqSites );
		 Site m = ann.siteMax( seqSites );
		int p = m.start;
		pos[j] = p ;
	  }  
	
	
	
	
	
		ofstream core2("corefile2.txt");
			for( int jj = 0; jj < nSeqs ; jj++) {
				
				
				
				int start = pos[jj];
				if(pos[jj] + l > seqs[jj].size() ) { start = 1; } 
				
				core2 << ">" <<  seqNames[jj] << endl;
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start + ii;
		
				core2 << ALPHABET[ seqs[jj][p] ];
				}
				core2 << endl;
				
			}
		core2.close();
		dcFile = "corefile2.txt";
		Matrix dcm = countmatrix( dcFile );
		
		Motif dcmm( dcm, background);
		
		
	
		pwmrmse = pwm_model;  
		
		pwm_model = dcmm.getInformation();  
		
			
			if ( pwm_model > pwm_best ) {
				pwm_best = pwm_model;	
				ofstream optcore("coreOPT.txt");
				ofstream ipf("initialPOS.txt");
				for( int jj = 0; jj < nSeqs ; jj++) {
					int start = pos[jj];
					optcore << ">" << start << seqNames[jj] << endl;
					ipf << start << '\t' ;
					for( int ii = 0; ii < l ; ii++) {  
					 int p = start + ii;
					optcore << ALPHABET[ seqs[jj][p] ];
					}
				optcore << endl;
				}
				optcore.close();
				ipf << endl;
				ipf.close();
			
			} 
		
		if( abs(pwm_model - pwmrmse) < .01 ) break;

	}
}
opt = "coreOPT.txt";
countmatrixflip(opt); 
string  cp = "rcRead.txt";
ifstream finp( cp.c_str() );
ofstream dest( file2.c_str() );
while( !finp.eof() )
{
string temp = "";
getline(finp,temp);
dest << temp << endl;
}
finp.close();
dest.close();
}
  void align( string& file,int nrandStarts, int nUpdates , int wid, int l   ){


	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		

int  initp = 1; 
string  ipfile;


string seqFile, dcFile;
seqFile=file;
vector< Sequence > seqs;
    vector< string > seqNames;
 int rval;
    rval = readSequences( seqFile, seqs, seqNames );  
    assert( rval != RET_ERROR );
    int nSeqs = seqs.size();
SiteVec seqSitesP;
 double gcContent = .5;

vector< Motif > motifs2;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );
double energyThr =2;
vector< double > energyThrs( 1, energyThr );
	
double pwm_best = 1;
double pwm_model =2.1;
double pwmrmse = 1.8;
int pin;
vector< int > pos(nSeqs,initp);  
    if ( !ipfile.empty() ) {
        ifstream fcoop( ipfile.c_str() );
        if ( !fcoop ) {
            cerr << "Cannot open the ip file " << endl;
            exit( 1 );
        }  
       for(int i = 0 ; i < nSeqs ; i++){
         fcoop >> pin;
	 pos[i] = pin;
        }        
    } 


vector< int > posbest(nSeqs,initp);

	ofstream core("corefile.txt");
		for( int jj = 0; jj < nSeqs ; jj++) {
			
			if(pos[jj] + l > seqs[jj].size() ) { pos[jj] = 1; }   
			
			
			core << ">" << seqNames[jj] << endl;
		    
			for( int ii = pos[jj]; ii < pos[jj] +l ; ii++) {  
			core << ALPHABET[ seqs[jj][ii] ];
			}
			core << endl;
			
		}
	core.close();
	dcFile = "corefile.txt";


	Matrix dcm = countmatrix( dcFile );
	
	Motif dcmm( dcm, background);
	
	
	
for ( int k = 0; k < nUpdates; k++ ) {
 for ( int j = 0; j < seqs.size(); j++ ) {
	ofstream core("corefile.txt");
	for( int jj = 0; jj < nSeqs ; jj++) {
			if(j == jj) { continue;}  
			if(pos[jj] + l > seqs[jj].size() ) { pos[jj] = 1; }   
			
			
			core << ">" << seqNames[jj] << endl;
		    
			for( int ii = pos[jj]; ii < pos[jj] +l ; ii++) {  
			core << ALPHABET[ seqs[jj][ii] ];
			}
			core << endl;
			
	}
	core.close();
	dcFile = "corefile.txt";
	vector< Motif > motifs2;

	Matrix dcm = countmatrix( dcFile );
	Motif dcmm( dcm, background);
	motifs2.push_back( dcmm ); 
 
        SeqAnnotator ann( motifs2, energyThrs );
	SiteVec seqSites;
	 ann.annot( seqs[ j ], seqSites );
	
	 Site m = ann.siteMax( seqSites );
	
	int p = m.start;
	pos[j] = p ;
	
 

  }  
 
	ofstream core2("corefile2.txt");
		for( int jj = 0; jj < nSeqs ; jj++) {
			
			
			
			int start = pos[jj];
			if(pos[jj] + l > seqs[jj].size() ) { start = 1; } 
			
			core2 << ">" <<  seqNames[jj] << endl;
			for( int ii = 0; ii < l ; ii++) {  
			 int p = start + ii;
		
			core2 << ALPHABET[ seqs[jj][p] ];
			}
			core2 << endl;
			
		}
	core2.close();
	dcFile = "corefile2.txt";
	Matrix dcm = countmatrix( dcFile );
	
	Motif dcmm( dcm, background);
	
	
	
	pwmrmse = pwm_model;  
	
	pwm_model = dcmm.getInformation();  
	
		
		if ( pwm_model > pwm_best ) {
			pwm_best = pwm_model;	
			ofstream optcore("coreOPT.txt");	
			ofstream ipf("initialPOS.txt");
			for( int jj = 0; jj < nSeqs ; jj++) {
				int start = pos[jj];
				optcore << ">" << start << seqNames[jj] << endl;
				ipf << start << '\t' ;
				posbest[jj] = start;
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start + ii;
				optcore << ALPHABET[ seqs[jj][p] ];
				}
			optcore << endl;
			}
			ipf << endl;
			optcore.close();
			ipf.close();
			
		} 
	
	if( abs(pwm_model - pwmrmse) < .01 ) break;

}

string opt = "coreOPT.txt";
for ( int r = 0; r < nrandStarts; r++){
	for(int i = 0; i < nSeqs; i++){
	
	pos[i]=  gsl_rng_uniform_int(rng,wid); 
	}
	
	for ( int k = 0; k < nUpdates; k++ ) {
	
	 for ( int j = 0; j < seqs.size(); j++ ) {
		ofstream core("corefile.txt"); 
			for( int jj = 0; jj < nSeqs ; jj++) {
				if(j == jj) { continue;}
				if(pos[jj] + l > seqs[jj].size() ) { pos[jj] = 1; }   
				
				core << ">" << seqNames[jj] << endl;
				for( int ii = pos[jj]; ii < pos[jj] +l ; ii++) {  
				core << ALPHABET[ seqs[jj][ii] ];
				}
				core << endl;
				
			}
		core.close();
		dcFile = "corefile.txt";
		vector< Motif > motifs2; 
		Matrix dcm = countmatrix( dcFile );
		Motif dcmm( dcm, background);
		motifs2.push_back( dcmm ); 
	 
		SeqAnnotator ann( motifs2, energyThrs );
		SiteVec seqSites;
		 ann.annot( seqs[ j ], seqSites );
		 Site m = ann.siteMax( seqSites );
		int p = m.start;
		pos[j] = p ;
	  }  
	
	
	
	
	
		ofstream core2("corefile2.txt");
			for( int jj = 0; jj < nSeqs ; jj++) {
				
				
				
				int start = pos[jj];
				if(pos[jj] + l > seqs[jj].size() ) { start = 1; } 
				
				core2 << ">" <<  seqNames[jj] << endl;
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start + ii;
		
				core2 << ALPHABET[ seqs[jj][p] ];
				}
				core2 << endl;
				
			}
		core2.close();
		dcFile = "corefile2.txt";
		Matrix dcm = countmatrix( dcFile );
		
		Motif dcmm( dcm, background);
		
		
	
		pwmrmse = pwm_model;  
		
		pwm_model = dcmm.getInformation();  
		
			
			if ( pwm_model > pwm_best ) {
				pwm_best = pwm_model;	
				ofstream optcore("coreOPT.txt");
				ofstream ipf("initialPOS.txt");
				for( int jj = 0; jj < nSeqs ; jj++) {
					int start = pos[jj];
					optcore << ">" << start << seqNames[jj] << endl;
					ipf << start << '\t' ;
					for( int ii = 0; ii < l ; ii++) {  
					 int p = start + ii;
					optcore << ALPHABET[ seqs[jj][p] ];
					}
				optcore << endl;
				}
				optcore.close();
				ipf << endl;
				ipf.close();
			
			} 
		
		if( abs(pwm_model - pwmrmse) < .01 ) break;

	}
}
opt = "coreOPT.txt";
countmatrixflip(opt); 

}


Matrix countmatrixS( const string& file )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  





	 Matrix m(4,length,0);
			
ofstream ja("jake.txt");
	for( int j = 0; j < seq.size(); j++ ){
		ja << "[ ";
		Sequence readseq(seq[j]);
			
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;  ja<< " 1 0 0 0 " ;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;  ja << " 0 1 0 0 ";}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;  ja<< " 0 0 1 0 ";}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1; ja << " 0 0 0 1 ";}
				
			}
	ja << " ]";
	ja<< endl;
	
	


	} 
	ja.close();
	Matrix countMatrix = m.transpose();
	

vector< double > _background(4,.25);

	 double pseudoCount =1;


	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;

	   return pwm;
}

void alignNEW(string& file,int nrandStarts, int nUpdates , int wid, int l   ){


int  initp = 1; 
string  ipfile;


string seqFile, dcFile;
seqFile=file;





    vector< Sequence > seqs;
    vector< string > seqNames;
    int rval;
    rval = readSequences( seqFile, seqs, seqNames );  
    assert( rval != RET_ERROR );
    int nSeqs = seqs.size();
    SiteVec seqSitesP;
    double gcContent = .5;
    vector< Motif > motifs2;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );
    double energyThr =2;
    vector< double > energyThrs( 1, energyThr );
	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );

double pwm_best = 1;
double pwm_model2 =1.1;
double pwm_model= 4;
double pwm_shift_best = 1;
double pwm_shift0 = 1;
double pwm_shift = 1;
int pin;
vector< int > stran(nSeqs,0);
vector< int > pos(nSeqs,initp);  
vector<double > sensitivSites(nSeqs,0);
    if ( !ipfile.empty() ) {
        ifstream fcoop( ipfile.c_str() );
        if ( !fcoop ) {
            cerr << "Cannot open the ip file " << endl;
            exit( 1 );
        }  
       for(int i = 0 ; i < nSeqs ; i++){
         fcoop >> pin;
	 pos[i] = pin;
        }        
    } 
	vector< int > posbest(nSeqs,initp);
	
	ofstream core("corefile.txt");
		for( int jj = 0; jj < nSeqs ; jj++) {
			if(pos[jj] + l > seqs[jj].size() ) { pos[jj] = 1; }   
			core << ">" << seqNames[jj] << endl;
			for( int ii = pos[jj]; ii < pos[jj] +l ; ii++) {  
			core << ALPHABET[ seqs[jj][ii] ];
			}
			core << endl;
		}
	core.close();
	dcFile = "corefile.txt";
	Matrix dcma = countmatrixS( dcFile );
	Motif dcmma( dcma, background);
	
	pwm_best = dcmma.getInformation();
	
for ( int k = 0; k < nUpdates; k++ ) {
	ofstream core("corefile.txt");
	
		for( int jj = 0; jj < nSeqs ; jj++) {
				if(pos[jj] + l > seqs[jj].size() ) { pos[jj] = 1; }  			
				int start = pos[jj];
				
				if(stran[jj]){
					core << ">" << seqNames[jj] << endl;
					if(  start+ l > seqs[jj].size() ){continue;}
					for( int ii = 0; ii < l ; ii++) {  
					 int p = start + ii;
					core << ALPHABET[ seqs[jj][p] ];
					}
				}
				else{  
					if( k == 0 ){                                         
					    for( int ii = 0; ii < l ; ii++) {  
						 int p = start + ii;
						core << ALPHABET[ seqs[jj][p] ];
						}
						core << endl;
						continue; 
					}
					if(  start  <= 0 ){continue;}
					core << ">" << seqNames[jj] << endl;
					for( int ii = 0; ii < l ; ii++) {  
					 int p = start ;
					core << ALPHABET[complement(symbolToInt(ALPHABET[ seqs[jj][p+l -1 - ii ] ]))];;
					}
				}
				core << endl;
		}
	
	core.close();
	dcFile = "corefile.txt";
	vector< Motif > motifs2;

	Matrix dcm = countmatrixS( dcFile );     

	Motif dcmm( dcm, background);
	pwm_model2 = dcmm.getInformation();
	motifs2.push_back( dcmm ); 
        SeqAnnotator ann( motifs2, energyThrs );
 for ( int j = 0; j < seqs.size(); j++ ) {       
	SiteVec seqSites;
	 ann.annot( seqs[ j ], seqSites );
	 Site m = ann.siteMax( seqSites );
	int p = m.start;
	int st= m.strand;
	stran[j] = st;
	pos[j] = p ;	
  }  
	ofstream core2("corefile2.txt");
		for( int jj = 0; jj < nSeqs ; jj++) {
			int start = pos[jj];
			if(pos[jj] + l > seqs[jj].size() ) { start = 1; }	
			core2 << ">" <<  seqNames[jj] << stran[jj] << endl;
			if(stran[jj]){
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start + ii;
				core2 << ALPHABET[ seqs[jj][p] ];
				}
			}
			else{  
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start ;
				core2 << ALPHABET[complement(symbolToInt(ALPHABET[ seqs[jj][p+ motifs2[ 0].length() -1 - ii ] ]))];;
				}
			}
			core2 << endl;
		}
	core2.close();
	dcFile = "corefile2.txt";

	Matrix dcm2 = countmatrixS( dcFile );
	Motif dcmm2( dcm2, background);
	pwm_model= pwm_model2;  
	pwm_model2 = dcmm2.getInformation();  
	
	
		if ( pwm_model2 > pwm_best ) {
			pwm_best = pwm_model2;	
			ofstream optcore("coreOPT.txt");	
			ofstream ipf("initialPOS.txt");
			for( int jj = 0; jj < nSeqs ; jj++) {
				int start = pos[jj];
				ipf << start << '\t' ;
					if(pos[jj] + l > seqs[jj].size() ) { start = 1; } 
			optcore << ">" <<  seqNames[jj] << stran[jj] << endl;
			if(stran[jj]){
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start + ii;
				optcore << ALPHABET[ seqs[jj][p] ];
				}
			}
			else{
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start ;
				optcore << ALPHABET[complement(symbolToInt(ALPHABET[ seqs[jj][p+ motifs2[ 0].length() -1 - ii ] ]))];;
				}
			}
			optcore << endl;
			}
			ipf << endl;
			optcore.close();
			ipf.close();
			
		} 
	
	
	pwm_shift_best = pwm_model2;
	pwm_shift0 = pwm_model2;

	for(int shif = 0; shif < 0 ; shif++ ){


		ostringstream numberm;
			numberm << shif;
	
				ofstream core3(("corefile3.txt" + numberm.str()).c_str());
			
			

			int miss = 0;
			for( int jj = 0; jj < nSeqs ; jj++) {
				int start = pos[jj];
				if(pos[jj] + l + shif > seqs[jj].size() ) { miss=miss +1; continue; }  	
				if(stran[jj]){
					core3 << ">" <<  seqNames[jj] << stran[jj] << endl;
					for( int ii = 0; ii < l ; ii++) {  
					 int p = start + ii +shif;
					core3 << ALPHABET[ seqs[jj][p] ];
					}
				}
				else{
				    if(  start -shif <= 0 ){miss=miss +1;   continue;}
					core3 << ">" <<  seqNames[jj] << stran[jj] << endl;
					for( int ii = 0; ii < l ; ii++) {  
					 int p = start - shif ;
					core3 << ALPHABET[complement(symbolToInt(ALPHABET[ seqs[jj][p+ motifs2[ 0].length() -1 - ii ] ]))];;
					}
				}
				core3 << endl;
			}
		core3.close();
	
		if( nSeqs - miss < nSeqs-2 ) { break;}  

		dcFile =("corefile3.txt" + numberm.str());
		Matrix dcm = countmatrixS( dcFile );
		Motif dcmm( dcm, background);
		
		
		pwm_shift = dcmm.getInformation();  
		
		
			if ( pwm_shift >  pwm_shift_best) {   
				pwm_shift_best = pwm_shift;
			}  
	  
				
							
			if( pwm_shift_best > pwm_shift0 ){     
				for( int jj = 0; jj < nSeqs ; jj++) {		
					if(stran[jj]){
						
						pos[jj] = pos[jj] + shif;
					} 
					else{  pos[jj] =pos[jj] - shif; }
				}
			} 
			if ( pwm_shift_best > pwm_best ) {
			copyf( dcFile, "coreOPT.txt");
			
			} 


	} 



}

for ( int r = 0; r < nrandStarts; r++){
	for(int i = 0; i < nSeqs; i++){
	pos[i]=  gsl_rng_uniform_int(rng,wid); 
	}
	for(int i = 0; i < nSeqs; i++){
	stran[i]=  gsl_rng_uniform_int(rng,2); 
	}
	for ( int k = 0; k < nUpdates; k++ ) {
	 for ( int j = 0; j < seqs.size(); j++ ) {
		ofstream core("corefile.txt"); 
			for( int jj = 0; jj < nSeqs ; jj++) {
				if(pos[jj] + l > seqs[jj].size() ) { pos[jj] = 1; }   
				core << ">" << seqNames[jj] << endl;
				for( int ii = pos[jj]; ii < pos[jj] +l ; ii++) {  
				core << ALPHABET[ seqs[jj][ii] ];
				}
				core << endl;
			}
		core.close();
		dcFile = "corefile.txt";
		vector< Motif > motifs2; 
		Matrix dcm = countmatrixS( dcFile );
		Motif dcmm( dcm, background);
		motifs2.push_back( dcmm ); 
	 
		SeqAnnotator ann( motifs2, energyThrs );
		SiteVec seqSites;
		 ann.annot( seqs[ j ], seqSites );
		 Site m = ann.siteMax( seqSites );
		int p = m.start;
		pos[j] = p ;
		int st= m.strand;
		stran[j] = st;
	  }  
	
	
	
	
	
		ofstream core2("corefile2.txt");
			for( int jj = 0; jj < nSeqs ; jj++) {
				int start = pos[jj];
				core2 << ">" << seqNames[jj] << endl;
			if(stran[jj]){
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start + ii;
				core2 << ALPHABET[ seqs[jj][p] ];
				}
			}
			else{
				if( k == 0 ){ 
				    for( int ii = 0; ii < l ; ii++) {  
					 int p = start + ii;
					core2 << ALPHABET[ seqs[jj][p] ];
					}
					core2 << endl;
					continue; 
				}
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start ;
				core2 << ALPHABET[complement(symbolToInt(ALPHABET[ seqs[jj][p+l -1 - ii ] ]))];;
				}
			}
			core2 << endl;
		} 
		core2.close();
		dcFile = "corefile2.txt";
		Matrix dcm = countmatrixS( dcFile );
		Motif dcmm( dcm, background);
		pwm_model= pwm_model2;  
		pwm_model2 = dcmm.getInformation();  
	ofstream ipf("initialPOS.txt");
		if ( pwm_model2 > pwm_best ) {
			pwm_best = pwm_model2;	
			string coreop = "coreOPT.txt";
			copyf( dcFile, coreop );   
			for( int jj = 0; jj < nSeqs ; jj++) {
				int start = pos[jj];
					if(pos[jj] + l > seqs[jj].size() ) { start = 1; } 
				ipf << start << '\t' ;
			} 
		ipf << endl;
	ipf.close();
		
	
	int miss=0;
	for(int shif = 0; shif < 0 ; shif++ ){
		ostringstream numberm;
		numberm << shif;
			ofstream core3(("coreOPT.txt" + numberm.str()).c_str());
		for( int jj = 0; jj < nSeqs ; jj++) {
			int start = pos[jj];
			if(pos[jj] + l + shif > seqs[jj].size() ) { miss=miss + 1; continue; } 
			
			if(stran[jj]){
				core3 << ">" <<  seqNames[jj] << stran[jj] << endl;
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start + ii +shif;
				core3 << ALPHABET[ seqs[jj][p] ];
				}
			}
			else{
			    if(  start -shif < 0 ){miss=miss + 1; continue;}
				core3 << ">" <<  seqNames[jj] << stran[jj] << endl;
				for( int ii = 0; ii < l ; ii++) {  
				 int p = start - shif ;
				core3 << ALPHABET[complement(symbolToInt(ALPHABET[ seqs[jj][p+ l -1 - ii ] ]))];;
				}
			}
			core3 << endl;
		}
	core3.close();
	if( nSeqs - miss < 10 ) { break;}  
	dcFile =("coreOPT.txt" + numberm.str());

	Matrix dcm = countmatrixS( dcFile );
		Motif dcmm( dcm, background);
		
		
		pwm_shift = dcmm.getInformation();  
		
		
			if ( pwm_shift >  pwm_shift_best) {   
				pwm_shift_best = pwm_shift;
			}  
	  
				
							
			if( pwm_shift_best > pwm_shift0 ){
				for( int jj = 0; jj < nSeqs ; jj++) {		
					if(stran[jj]){
						
						pos[jj] = pos[jj] + shif;
					} 
					else{  pos[jj] =pos[jj] - shif; }
				}
			} 
				
			if ( pwm_shift_best > pwm_best ) {
			copyf( dcFile, "coreOPT.txt");
			
			} 







} 
		} 
ipf.close();
string opt2 = "coreOPT.txt";  
string optr="rcRead.txt";
countmatrixflip(dcFile );
Matrix optm = countmatrixS( optr );
Motif dcmm2( optm, background);
if (  dcmm2.getInformation() > pwm_best ) {
		copyf( optr, opt2 );   		
}

		if( abs(pwm_model2 - pwm_model) < .01 ) break;

	}
}

}






double Motif::PvalueSum9( double alpha ) 
{
int l = pwm.nRows();
vector< double > freqs;
vector< double > eners;
freqs.clear();
eners.clear();
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = background[p0];
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*background[p1];
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*background[p2];
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*background[p3];
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*background[p4];
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*background[p5];
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*background[p6];
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*background[p7];
			double enerf7 = enerf6 + LLRMat(7,p7);
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*background[p8];
		double enerf8 = enerf7 + LLRMat(8,p8);
	
freqs.push_back( pvalf8);
eners.push_back( maxLLR - enerf8 );  


    }
}
}
}
}
}
}
}
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
	

ofstream pvalfile("pv.txt");
double ethres =0;
vector< int > ne;
ofstream pvalfile2("pv2.txt");
int ntemp =0;
for(int i = -1 ; i < 160 ; i++){  
ethres = .5*i ;
double pvalf=0;
int n =0;

	for ( int j = 0; j < eners.size() ; j++ ) {
		n = n + 1;
		pvalf = freqs[j] + pvalf;	
		if( ethres <= eners[j] ) { 	
		
		pvalfile << ethres << '\t' << pvalf << endl;
		pvalfile2<< ethres<< '\t' << abs(n - ntemp) << endl;
		ntemp = n;
		break;
		}
}  
}
pvalfile.close();
pvalfile2.close();

double pvalf2=0;

double sum=0;
	for ( int j = 0; j < eners.size() ; j++ ) {
		
		pvalf2 = freqs[j] + pvalf2;	
		if( alpha <= pvalf2 ) { 	
		
		sum = eners[j];
		break;
		}
	}  

	
    return  sum ; 

}

double Motif::PvalueSum10( double alpha ) 
{
int l = pwm.nRows();
vector< double > freqs;
vector< double > eners;
freqs.clear();
eners.clear();
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = background[p0];
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*background[p1];
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*background[p2];
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*background[p3];
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*background[p4];
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*background[p5];
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*background[p6];
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*background[p7];
			double enerf7 = enerf6 + LLRMat(7,p7);
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*background[p8];
		double enerf8 = enerf7 + LLRMat(8,p8);
	for( int p9 = 0 ; p9 < 4 ; p9++){
	double pvalf9 = pvalf8*background[p9];
	double enerf9 = enerf8 + LLRMat(9,p9);

freqs.push_back( pvalf9);
eners.push_back( maxLLR - enerf9 );  

  }
    }
}
}
}
}
}
}
}
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
	

ofstream pvalfile("pv.txt");
double ethres =0;
vector< int > ne;
ofstream pvalfile2("pv2.txt");
int ntemp =0;
for(int i = -1 ; i < 160 ; i++){  
ethres = .5*i ;
double pvalf=0;
int n =0;

	for ( int j = 0; j < eners.size() ; j++ ) {
		n = n + 1;
		pvalf = freqs[j] + pvalf;	
		if( ethres <= eners[j] ) { 	
		
		pvalfile << ethres << '\t' << pvalf << endl;
		pvalfile2<< ethres<< '\t' << abs(n - ntemp) << endl;
		ntemp = n;
		break;
		}
}  
}
pvalfile.close();
pvalfile2.close();

double pvalf2=0;

double sum=0;
	for ( int j = 0; j < eners.size() ; j++ ) {
		
		pvalf2 = freqs[j] + pvalf2;	
		if( alpha <= pvalf2 ) { 	
		
		sum = eners[j];
		break;
		}
	}  

	
    return  sum ; 

}

double mutualinformatio( const vector< Motif >& motifs, double alphadc )
{
vector< double > freqs;
vector< double > freqscb;
vector< double > freqsdc;
vector< double > freqsdu;

freqs.clear();

Matrix pwmdc =motifs[0].getPwm();
Matrix pwmdu =motifs[1].getPwm();
Matrix pwm =motifs[3].getPwm();



double normdc =0;
double normdu = 0;
double normcb = 0;
double normcbm =0;
int index =0;
double mi=0;
double mi2=0;
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double pvalf0dc=1.0;
double pvalf0du=1.0;
pvalf0 = pwm(0,p0);
pvalf0dc = pwmdc(0,p0);
pvalf0du = pwmdu(0,p0);



	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*pwm(1,p1);
	double pvalf1dc = pvalf0dc*pwmdc(1,p1);
	double pvalf1du = pvalf0du*pwmdu(1,p1);
	

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*pwm(2,p2);
		double pvalf2dc = pvalf1dc*pwmdc(2,p2);
		double pvalf2du = pvalf1du*pwmdu(2,p2);
	
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*pwm(3,p3);
			double pvalf3dc = pvalf2dc*pwmdc(3,p3);
			double pvalf3du = pvalf2du*pwmdu(3,p3);
		
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*pwm(4,p4);
				double  pvalf4dc = pvalf3dc*pwmdc(4,p4);
				double  pvalf4du = pvalf3du*pwmdu(4,p4);
			
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*pwm(5,p5);
					double pvalf5dc = pvalf4dc*pwmdc(5,p5);
					double pvalf5du = pvalf4du*pwmdu(5,p5);
				
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*pwm(6,p6);
				double pvalf6dc = pvalf5dc*pwmdc(6,p6);
				double pvalf6du = pvalf5du*pwmdu(6,p6);
		
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*pwm(7,p7);
			double pvalf7dc = pvalf6dc*pwmdc(7,p7);
			double pvalf7du = pvalf6du*pwmdu(7,p7);
	
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*pwm(8,p8);
		double pvalf8dc = pvalf7dc*pwmdc(8,p8);
		double pvalf8du = pvalf7du*pwmdu(8,p8);

	for( int p9 = 0 ; p9 < 4 ; p9++){
	double pvalf9 = pvalf8*pwm(9,p9);
	double pvalf9dc = pvalf8dc*pwmdc(9,p9);
	double pvalf9du = pvalf8du*pwmdu(9,p9);

for( int p10 = 0 ; p10 < 4 ; p10++){
double pvalf10 = pvalf9*pwm(10,p10);
double pvalf10dc = pvalf9dc*pwmdc(10,p10);
double pvalf10du = pvalf9du*pwmdu(10,p10);
mi = mi + alphadc*pvalf10dc * log2(pvalf10dc/pvalf10) +(1-alphadc)*pvalf10du* log2(pvalf10du/pvalf10);
mi2= mi2 + alphadc*pvalf10dc * log2(pvalf10dc/(alphadc*pvalf10dc+ (1-alphadc)*pvalf10du)) +(1-alphadc)*pvalf10du* log2(pvalf10du/(alphadc*pvalf10dc+ (1-alphadc)*pvalf10du));
index = index +1;
if(index == 10000) {  
}
if(index == 10001) { }
if(index == 4000000) { 
 }
if(index == 4000001) { 

}
normcb = normcb + pvalf10;
normdc =normdc + pvalf10dc;
normdu =normdu + pvalf10du;
normcbm = normcbm + alphadc*pvalf10dc + (1-alphadc)*pvalf10du;
 
}
  }
    }
}
}
}
}
}
}
}
}


	


return mi;
}
double Motif::PvalueSumf( double alpha ) 
{
int l = pwm.nRows();
vector< double > freqs;
vector< double > eners;
freqs.clear();
eners.clear();

int nindex = int(maxLLR-minLLR);
vector< vector<double> > bins(int(maxLLR-minLLR)+1, vector< double >(0));
int index =0;
double sumbin=0;
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = pwm(0,p0);
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*pwm(1,p1);
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*pwm(2,p2);
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*pwm(3,p3);
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*pwm(4,p4);
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*pwm(5,p5);
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*pwm(6,p6);
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*pwm(7,p7);
			double enerf7 = enerf6 + LLRMat(7,p7);
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*pwm(8,p8);
		double enerf8 = enerf7 + LLRMat(8,p8);
	for( int p9 = 0 ; p9 < 4 ; p9++){
	double pvalf9 = pvalf8*pwm(9,p9);
	double enerf9 = enerf8 + LLRMat(9,p9);
for( int p10 = 0 ; p10 < 4 ; p10++){
double pvalf10 = pvalf9*pwm(10,p10);
double enerf10 = enerf9 + LLRMat(10,p10);
index = int( maxLLR-enerf10);
bins[index].push_back( pvalf10 );
freqs.push_back( pvalf10 );
eners.push_back( maxLLR -enerf10 );  
}
  }
    }
}
}
}
}
}
}
}
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
	
vector< double >bins2;
double z=0;
for(int e=0;e< bins.size() ; e++ ){  
double num = 0;
	for(int i=0; i<bins[e].size();i++){
			num=num + bins[e][i];
			}
z=z+num;
bins2.push_back( num );
}



ofstream pvalfile("pvf.txt");
pvalfile  << bins2 << endl;
for(int i = 0 ; i < bins2.size() ; i++ ){
}
pvalfile.close();
double pvalf=0;
for(int i = 0 ; i < bins2.size() ; i++){  

		pvalf = bins2[i]/z + pvalf;	
			if( alpha <= pvalf ) { 	
			
			
			
			
			pvalf=i;  
			break;
			}  
}


double pvalf2=0;

double sum=0;
	for ( int j = 0; j < eners.size() ; j++ ) {
		
		pvalf2 = freqs[j] + pvalf2;	
		if( alpha <= pvalf2 ) { 	
		
		sum = eners[j];
		break;
		}
	}  
sum = pvalf;



	
    return  sum ; 

}


double Motif::PvalueSumCond( double alpha )
{
return 1.0;
}
double Motif::PvalueSum( double alpha ) 
{
int l = pwm.nRows();
vector< double > freqs;
vector< double > eners;
freqs.clear();
eners.clear();
int nindex = int(10*(maxLLR-minLLR));


vector< vector<double> > bins(int(10*(maxLLR-minLLR))+1, vector< double >(0));
int index =0;
double sumbin=0;
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = background[p0];
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*background[p1];
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*background[p2];
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*background[p3];
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*background[p4];
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*background[p5];
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*background[p6];
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*background[p7];
			double enerf7 = enerf6 + LLRMat(7,p7);
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*background[p8];
		double enerf8 = enerf7 + LLRMat(8,p8);
	
freqs.push_back( pvalf8 );

index = int( 10*(maxLLR-enerf8));
bins[index].push_back( pvalf8 );
eners.push_back( 10*(maxLLR - enerf8 ));  
  
    }
}
}
}
}
}
}
}
}
vector< double >bins2;
double z=0;
for(int e=0;e< bins.size() ; e++ ){  
double num = 0;
	for(int i=0; i<bins[e].size();i++){
			num=num + bins[e][i];
			}
z=z+num;
bins2.push_back( num );
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
	

double ethres =0;
vector< int > ne;
int ntemp =0;
double pvalf=0;
for(int i = 0 ; i < bins2.size() ; i++){  

		pvalf = bins2[i]/z + pvalf;	
			if( alpha <= pvalf ) { 	
			
		
			
			
			pvalf=i;  
			break;
			}  
}
ofstream cdfbackground("cdfbackground.txt");
double pvalfb=0;
for(int i = 0 ; i < bins2.size() ; i++){  

		pvalfb = bins2[i]/z + pvalfb;	
		cdfbackground << i << '\t' << pvalfb<< endl ;
			
			
		
			
			
			
			
		
}
cdfbackground << endl;
cdfbackground.close();
ofstream pvalfile("backgroundCounts.txt");
for(int i = 0 ; i < bins2.size() ; i++ ){
pvalfile << bins2[i] << endl;
}
pvalfile.close();

double pvalf2=0;

double sum=0;
	for ( int j = 0; j < eners.size() ; j++ ) {
		
		pvalf2 = freqs[j] + pvalf2;	
		if( alpha <= pvalf2 ) { 	
		
		sum = eners[j];
		break;
		}
	}  
sum = pvalf;

	
    return  sum ; 

}


double Motif::PvalueSumdc( double alpha ) 
{
int l = pwm.nRows();
vector< double > freqs;
vector< double > eners;
freqs.clear();
eners.clear();
int nindex = int(10*(maxLLR-minLLR));



vector< vector<double> > bins(int(10*(maxLLR-minLLR))+1, vector< double >(0));
int index =0;
double sumbin=0;
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = background[p0];
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*background[p1];
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*background[p2];
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*background[p3];
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*background[p4];
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*background[p5];
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*background[p6];
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*background[p7];
			double enerf7 = enerf6 + LLRMat(7,p7);
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*background[p8];
		double enerf8 = enerf7 + LLRMat(8,p8);
	for( int p9 = 0 ; p9 < 4 ; p9++){
	double pvalf9 = pvalf8*background[p9];
	double enerf9 = enerf8 + LLRMat(9,p9);
for( int p10 = 0 ; p10 < 4 ; p10++){
double pvalf10 = pvalf9*background[p10];
double enerf10 = enerf9 + LLRMat(10,p10);
freqs.push_back( pvalf10 );

index = int( 10*(maxLLR-enerf10));
bins[index].push_back( pvalf10 );
eners.push_back( 10*(maxLLR - enerf10) );  
}
  }
    }
}
}
}
}
}
}
}
}
vector< double >bins2;
double z=0;
for(int e=0;e< bins.size() ; e++ ){  
double num = 0;
	for(int i=0; i<bins[e].size();i++){
			num=num + bins[e][i];
			}
z=z+num;
bins2.push_back( num );
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
	

double ethres =0;
vector< int > ne;
int ntemp =0;
double pvalf=0;
for(int i = 0 ; i < bins2.size() ; i++){  

		pvalf = bins2[i]/z + pvalf;	
			if( alpha <= pvalf ) { 	
			
		
			
			
			pvalf=i;  
			break;
			}  
}
ofstream cdfbackgrounddc("cdfbackgrounddc.txt");
double pvalfb=0;
for(int i = 0 ; i < bins2.size() ; i++){  

		pvalfb = bins2[i]/z + pvalfb;	
		cdfbackgrounddc << i << '\t' << pvalfb<< endl ;
			
			
		
			
			
			
			
		
}
cdfbackgrounddc.close();
ofstream pvalfile("backgroundCountsdc.txt");
for(int i = 0 ; i < bins2.size() ; i++ ){
pvalfile << bins2[i] << endl;
}
pvalfile.close();

double pvalf2=0;

double sum=0;
	for ( int j = 0; j < eners.size() ; j++ ) {
		
		pvalf2 = freqs[j] + pvalf2;	
		if( alpha <= pvalf2 ) { 	
		
		sum = eners[j];
		break;
		}
	}  
sum = pvalf;

	
    return  sum ; 

}


double Motif::PvalueSumdu( double alpha ) 
{
int l = pwm.nRows();
vector< double > freqs;
vector< double > eners;
freqs.clear();
eners.clear();
int nindex = int(10*(maxLLR-minLLR));


vector< vector<double> > bins(int(10*(maxLLR-minLLR))+1, vector< double >(0));
int index =0;
double sumbin=0;
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = background[p0];
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*background[p1];
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*background[p2];
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*background[p3];
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*background[p4];
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*background[p5];
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*background[p6];
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*background[p7];
			double enerf7 = enerf6 + LLRMat(7,p7);
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*background[p8];
		double enerf8 = enerf7 + LLRMat(8,p8);
	for( int p9 = 0 ; p9 < 4 ; p9++){
	double pvalf9 = pvalf8*background[p9];
	double enerf9 = enerf8 + LLRMat(9,p9);
for( int p10 = 0 ; p10 < 4 ; p10++){
double pvalf10 = pvalf9*background[p10];
double enerf10 = enerf9 + LLRMat(10,p10);
freqs.push_back( pvalf10 );

index = int( 10*(maxLLR-enerf10));
bins[index].push_back( pvalf10 );
eners.push_back( 10*(maxLLR - enerf10) );  
}
  }
    }
}
}
}
}
}
}
}
}
vector< double >bins2;
double z=0;
for(int e=0;e< bins.size() ; e++ ){  
double num = 0;
	for(int i=0; i<bins[e].size();i++){
			num=num + bins[e][i];
			}
z=z+num;
bins2.push_back( num );
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
	

double ethres =0;
vector< int > ne;
int ntemp =0;
double pvalf=0;
for(int i = 0 ; i < bins2.size() ; i++){  

		pvalf = bins2[i]/z + pvalf;	
			if( alpha <= pvalf ) { 	
			
		
			
			
			pvalf=i;  
			break;
			}  
}
ofstream cdfbackground("cdfbackgrounddu.txt");
double pvalfb=0;
for(int i = 0 ; i < bins2.size() ; i++){  

		pvalfb = bins2[i]/z + pvalfb;	
		cdfbackground << i << '\t' << pvalfb<< endl ;
			
			
		
			
			
			
			
		
}
cdfbackground.close();
ofstream pvalfile("backgroundCountsdu.txt");
for(int i = 0 ; i < bins2.size() ; i++ ){
pvalfile << bins2[i] << endl;
}
pvalfile.close();

double pvalf2=0;

double sum=0;
	for ( int j = 0; j < eners.size() ; j++ ) {
		
		pvalf2 = freqs[j] + pvalf2;	
		if( alpha <= pvalf2 ) { 	
		
		sum = eners[j];
		break;
		}
	}  
sum = pvalf;

	
    return  sum ; 

}


double Motif::PvalueSumcb( double alpha ) 
{
int l = pwm.nRows();
vector< double > freqs;
vector< double > eners;
freqs.clear();
eners.clear();
int nindex = int(10*(maxLLR-minLLR));


vector< vector<double> > bins(int(10*(maxLLR-minLLR))+1, vector< double >(0));
int index =0;
double sumbin=0;
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = background[p0];
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*background[p1];
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*background[p2];
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*background[p3];
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*background[p4];
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*background[p5];
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*background[p6];
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*background[p7];
			double enerf7 = enerf6 + LLRMat(7,p7);
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*background[p8];
		double enerf8 = enerf7 + LLRMat(8,p8);
	for( int p9 = 0 ; p9 < 4 ; p9++){
	double pvalf9 = pvalf8*background[p9];
	double enerf9 = enerf8 + LLRMat(9,p9);
for( int p10 = 0 ; p10 < 4 ; p10++){
double pvalf10 = pvalf9*background[p10];
double enerf10 = enerf9 + LLRMat(10,p10);
freqs.push_back( pvalf10 );

index = int( 10*(maxLLR-enerf10));
bins[index].push_back( pvalf10 );
eners.push_back( 10*(maxLLR - enerf10) );  
}
  }
    }
}
}
}
}
}
}
}
}
vector< double >bins2;
double z=0;
for(int e=0;e< bins.size() ; e++ ){  
double num = 0;
	for(int i=0; i<bins[e].size();i++){
			num=num + bins[e][i];
			}
z=z+num;
bins2.push_back( num );
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
	

double ethres =0;
vector< int > ne;
int ntemp =0;
double pvalf=0;
for(int i = 0 ; i < bins2.size() ; i++){  

		pvalf = bins2[i]/z + pvalf;	
			if( alpha <= pvalf ) { 	
			
		
			
			
			pvalf=i;  
			break;
			}  
}
ofstream cdfbackground("cdfbackgroundcb.txt");
double pvalfb=0;
for(int i = 0 ; i < bins2.size() ; i++){  

		pvalfb = bins2[i]/z + pvalfb;	
		cdfbackground << i << '\t' << pvalfb<< endl ;
			
			
		
			
			
			
			
		
}
cdfbackground << endl;
cdfbackground.close();
ofstream pvalfile("backgroundCountscb.txt");
for(int i = 0 ; i < bins2.size() ; i++ ){
pvalfile << bins2[i] << endl;
}
pvalfile.close();

double pvalf2=0;

double sum=0;
	for ( int j = 0; j < eners.size() ; j++ ) {
		
		pvalf2 = freqs[j] + pvalf2;	
		if( alpha <= pvalf2 ) { 	
		
		sum = eners[j];
		break;
		}
	}  
sum = pvalf;

	
    return  sum ; 

}


double Motif::PvalueSumt( double alpha ) 
{
int l = pwm.nRows();


vector< double > freqs;
vector< double > eners;
freqs.clear();
eners.clear();
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = background[p0];
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*background[p1];
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*background[p2];
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*background[p3];
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*background[p4];
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*background[p5];
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*background[p6];
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*background[p7];
			double enerf7 = enerf6 + LLRMat(7,p7);
	
freqs.push_back( pvalf7 );
eners.push_back( maxLLR - enerf7);  
}
}
}
}
}
}
}
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
ofstream cumdistf("fdis.txt");
for ( int i = 0; i < eners.size() ; i++ ) {
	cumdistf << freqs[i] << '\t' << eners[i] << endl;
	
}
cumdistf.close();

ofstream cumdist("cdis.txt");	
double sum=0;
	for ( int i = 0; i < eners.size() ; i++ ) {
		sum = sum + freqs[i];
		if( sum > alpha ){cumdist << sum << '\t' <<eners[i] << endl; break;}
	}
	
cumdist.close();

ofstream pvalfile("pv.txt");
double ethres =0;
vector< int > ne;
ofstream pvalfile2("pv2.txt");
int ntemp =0;
for(int i = -1 ; i <60; i++){  
ethres = .5*i ;
double pvalf=0;
int n =0;

	for ( int j = 0; j < eners.size() ; j++ ) {
		n = n + 1;
		pvalf = freqs[j] + pvalf;	
		if( ethres <= eners[j] ) { 	
		
		pvalfile << ethres << '\t' << pvalf << endl;
		pvalfile2<< ethres<< '\t' << abs(n - ntemp) << endl;
		ntemp = n;
		break;
		}
}  
}
pvalfile.close();
pvalfile2.close();

double pvalf2=0;


	for ( int j = 0; j < eners.size() ; j++ ) {
		
		pvalf2 = freqs[j] + pvalf2;	
		if( alpha <= pvalf2 ) { 	
		
		sum = eners[j];
		break;
		}
	}  

	
    return  sum ; 

}
double Motif::PvalueSumt6( double alpha ) 
{
int l = pwm.nRows();


vector< double > freqs;
vector< double > eners;
freqs.clear();
eners.clear();
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = background[p0];
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*background[p1];
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*background[p2];
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*background[p3];
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*background[p4];
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*background[p5];
					double enerf5 = enerf4 + LLRMat(5,p5);
		
	
freqs.push_back( pvalf5 );
eners.push_back( maxLLR - enerf5);  
}
}
}
}
}
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
ofstream cumdistf("fdis6.txt");
for ( int i = 0; i < eners.size() ; i++ ) {
	cumdistf << freqs[i] << '\t' << eners[i] << endl;
	
}
cumdistf.close();

ofstream cumdist("cdis6.txt");	
double sum=0;
	for ( int i = 0; i < eners.size() ; i++ ) {
		sum = sum + freqs[i];
		if( sum > alpha ){cumdist << sum << '\t' <<eners[i] << endl; break;}
	}
	
cumdist.close();

ofstream pvalfile("pv6.txt");
double ethres =0;
vector< int > ne;
ofstream pvalfile2("pv26.txt");
int ntemp =0;
for(int i = -1 ; i <60; i++){  
ethres = .5*i ;
double pvalf=0;
int n =0;

	for ( int j = 0; j < eners.size() ; j++ ) {
		n = n + 1;
		pvalf = freqs[j] + pvalf;	
		if( ethres <= eners[j] ) { 	
		
		pvalfile << ethres << '\t' << pvalf << endl;
		pvalfile2<< ethres<< '\t' << abs(n - ntemp) << endl;
		ntemp = n;
		break;
		}
}  
}
pvalfile.close();
pvalfile2.close();

double pvalf2=0;


	for ( int j = 0; j < eners.size() ; j++ ) {
		
		pvalf2 = freqs[j] + pvalf2;	
		if( alpha <= pvalf2 ) { 	
		
		sum = eners[j];
		break;
		}
	}  

	
    return  sum ; 

}

double Motif::Pvalue( double alpha ) 
{
int l = pwm.nRows();

vector< double > freqs;
vector< double > eners;
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf=1;
double enerf=0;
pvalf = pvalf*pwm(0,p0);
enerf = enerf + LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf*pwm(1,p1);
	double enerf1 = enerf + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*pwm(2,p2);
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*pwm(3,p3);
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*pwm(4,p4);
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*pwm(5,p5);
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*pwm(6,p6);
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*pwm(7,p7);
			double enerf7 = enerf6 + LLRMat(7,p7);
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*pwm(8,p8);
		double enerf8 = enerf7 + LLRMat(8,p8);
	for( int p9 = 0 ; p9 < 4 ; p9++){
	double pvalf9 = pvalf8*pwm(9,p9);
	double enerf9 = enerf8 + LLRMat(9,p9);
for( int p10 = 0 ; p10 < 4 ; p10++){
double pvalf10 = pvalf9*pwm(10,p10);
double enerf10 = enerf9 + LLRMat(10,p10);
freqs.push_back( pvalf10 );
eners.push_back( maxLLR - enerf10 );  
}
  }
    }
}
}
}
}
}
}
}
}
	
	std::sort(freqs.begin(), freqs.end(), siteSortPredicate2);
	std::sort(eners.begin(), eners.end(), siteSortPredicate2);
	
ofstream cumdistf("fdis.txt");
for ( int i = 0; i < eners.size() ; i++ ) {
	cumdistf << freqs[i] << '\t' << eners[i] << endl;
}
cumdistf.close();

ofstream cumdist("cdis.txt", ios::app);	
double sum=0;
	for ( int i = 0; i < eners.size() ; i++ ) {
		sum = sum + freqs[i];
		if( sum > alpha ){cumdist << sum << '\t' <<eners[i] << endl; break;}
	}
	
cumdist.close();

	
    return  sum ; 
}

double Motif::energy( const Sequence& elem ) const
{
	
	
    return ( -LLR( elem ) +   maxLLR );	  
}

void Motif::sample( const gsl_rng* rng, Sequence& elem, bool strand ) const
{
    assert( rng != NULL );
    
    int l = pwm.nRows();
    Sequence sampleElem;
    for ( int i = 0; i < l; i++ ) {
        
        vector< double > distr = pwm.getRow( i );
        
        
        int nt = sampleMul( rng, distr );
        sampleElem.push_back( nt );
    }		
    
    if ( strand == 0 ) elem = sampleElem.compRevCompl();
    else elem = sampleElem;
}

int Motif::load( const string& file, const vector< double >& background, string& name )
{
    vector< Motif > motifs;
    vector< string > names;
    int rval = readMotifs( file, background, motifs, names );
    if ( rval == RET_ERROR ) return RET_ERROR;
    
    copy( motifs[ 0 ] );
    name = names[ 0 ];
    return rval;				
}

int Motif::load( const string& file, const vector< double >& background )
{
    string name;
    int rval = load( file, background, name );
    
    return rval;	
}

ostream& operator<<( ostream& os, const Motif& motif )
{
	Sequence s = motif.getMaxSite() ;
	os << "DE"<<'\t'<<rand()<<'\t'<<"REL"<<endl;
	for(int i =0; i< motif.pwm.nRows() ; i++){
		os << i << '\t';
		for ( int j = 0; j < 4; j++ ) {			
        	   os << motif.countmat(i,j) << '\t' ;
        	}
	
		os<< ALPHABET2[s[i]] << endl;
	
	}
	os << "XX" << endl;

    
    return os;
}
double Motif::information( const Sequence& elem ) const
{

    return ( INFO( elem ) ); 
}

double Motif::INFO( const Sequence& elem ) const
{          
    int l = pwm.nRows();
    if ( elem.size() != l ) return GSL_NEGINF;
    if ( elem.containsMissing() ) return GSL_NEGINF;
   
    double result = 0;
    for ( int i = 0; i < l; i++ ) {
        result += infomat( i, elem[ i ] ); 
    }

    return result;
}
void Motif::init()
{
    int l = pwm.nRows();
    
    
    for ( int i = 0; i < l; i++ ) {
        for ( int j = 0; j < 4; j++ ) {			
            LLRMat( i, j ) =  log( pwm( i, j ) );  
        }
    }
    
    
    for ( int i = 0; i < l; i++ ) {
        int b_max;
        max( pwm.getRow( i ), b_max );
        maxSite.push_back( b_max );	
    }
    
    
    maxLLR = 0;
    for ( int i = 0; i < l; i++ ) {
        maxLLR += LLRMat( i, maxSite[ i ] );	
    }


    
    for ( int i = 0; i < l; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
            infomat(i,j)=  log( 4*pwm(i,j));	
        }	
    }

	Matrix a =  pwm*(-1);
    
    for ( int i = 0; i < l; i++ ) {
        int b_min;
        max( a.getRow( i ), b_min );
        minSite.push_back( b_min );	
    }
    
    
    minLLR = 0;
    for ( int i = 0; i < l; i++ ) {
        minLLR += LLRMat( i, minSite[ i ] );	
    }
   
    

}
int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs, vector< string >& names )
{
    
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }
    motifs.clear(); 
    names.clear();

    string line;
    
    
    do {
        getline( fin, line );
        
        if ( line[ 0 ] != '>' ) continue;
        
        
        int MAX_SIZE = 100;
        char lineStr[ MAX_SIZE ];
        strcpy( lineStr, ( line.substr( 1 ) ).c_str() );
        char *name, *lengthStr, *pseudoCountStr;
        name = strtok( lineStr, " \t" );
        lengthStr = strtok( NULL, " \t" );
        pseudoCountStr = strtok( NULL, " \t" );
        int length;
        double pseudoCount;
        if ( lengthStr ) length = atoi( lengthStr );
        else { return RET_ERROR; }
        if ( pseudoCountStr ) pseudoCount = atof( pseudoCountStr );
        else pseudoCount = PSEUDO_COUNT;
        
        
        Matrix countMat( length, NBASES );
        for ( int i = 0; i < length; ++i ) {
            for ( int j = 0; j < NBASES; ++j ) {
                fin >> countMat( i, j );
            }	
        }
        
        
        names.push_back( string( name ) );
        motifs.push_back( Motif( countMat, pseudoCount, background ) );	
    } while ( !fin.eof() );
                                    
    return 0;
}

int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs )
{
    vector< string > names;
    return readMotifs( file, background, motifs, names );	
}

ostream& operator<<( ostream& os, const Site& site )
{
    char strandChar = site.strand ? '+' : '-';
    os <<"\n"<< site.start + 1 << "\t" << strandChar << "\t" << site.factorIdx << "\t" << site.energy << "\t" << site.wtRatio;
    
    return os;
}

bool siteOverlap( const Site& a, const Site& b, const vector< Motif >& motifs )
{
    if ( a.start + motifs[ a.factorIdx ].length() <= b.start ) return false;
    if ( b.start + motifs[ b.factorIdx ].length() <= a.start ) return false;
	
    
    return true;	
}
bool siteOverlap2( const Site& a, const Site& b, const vector< Motif >& motifs )
{
    if ( a.start + motifs[ a.factorIdx ].length() <= b.start ) return false;
    if ( b.start + motifs[ b.factorIdx ].length() <= a.start ) return false;
	if( a.start == b.start ) return false;
    
    return true;	
}

int maxstring( string& file )
{

 ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
			
		
	}
string templ;
int tem = 0;
int tempmax = 0;
	for( int j = 0; j < seq.size(); j++ ){

		for( int i = 0;  i< seq[j].size(); i++){
			templ[i] = seq[j][i];
			tem = templ.length();
			if( tem > tempmax) tempmax = tem;
		
		}
	}
	
return tempmax;

}

void SensitivityFN( const string& core , const string& dc2 ){  
ofstream TPcor( "TPdc.txt" );
ofstream TNcor("TNdc.txt" );
 ifstream dca( core.c_str() );
    if ( !dca ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq ;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!dca.eof()){
		temp = "";

		getline(dca, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(dca, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
ifstream dcb( dc2.c_str() );
    if ( !dcb ) {
      
    }

	string tempb;
	vector <string> seq_nameb;
	vector <string> seqb ;
	
	string templb;
int temb = 0;
int tempmaxb = 0;
	seqb.clear();
	while(!dcb.eof()){
		tempb= "";

		getline(dcb, tempb);
		if(tempb.length() == 0){
			break;
		}
	
		string nameb (tempb, 1, tempb.length() - 1);
		seq_nameb.push_back(nameb);
		
		getline(dcb, tempb);
		seqb.push_back(tempb);
		temb = tempb.length();
				
		
	}
int p =0;
int n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="aaatt";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}
TPcor.close();
TNcor.close();

p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="aatt";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="aaa";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="aa";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="ttt";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="tt";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="gg";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}



p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="ggg";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="aa";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="ga";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="at";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="ta";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


p =0;
n=0;

	for(int j=0; j< seqb.size(); j++ ){
	const string ast ="tc";  
	const string bst = seqb[j] ;
	size_t  found;
	 
	 found=bst.find(ast);
 	 if (found!=string::npos){ p=p+1;  continue;}
	
	
	else{  n = n+1;  }
	
	}


}
void SensitivityFN2( const string& core , const string& dc2, const string& tpout , const string& tnout  ){  
 ofstream tpou( tpout.c_str() );
 ofstream tnou( tnout.c_str() );
 ifstream dcore( core.c_str() );
    if ( !dcore ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq ;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!dcore.eof()){
		temp = "";

		getline(dcore, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(dcore, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
ifstream dcb( dc2.c_str() );
    if ( !dcb ) {
      
    }

	string tempb;
	vector <string> seq_nameb;
	vector <string> seqb ;
	
	string templb;
int temb = 0;
int tempmaxb = 0;
	seqb.clear();
	while(!dcb.eof()){
		tempb= "";

		getline(dcb, tempb);
		if(tempb.length() == 0){
			break;
		}
	
		string nameb (tempb, 1, tempb.length() - 1);
		seq_nameb.push_back(nameb);
		
		getline(dcb, tempb);
		seqb.push_back(tempb);
		temb = tempb.length();
				
		
	}
int p =0;
int n=0;

ofstream negs("fpsites.txt");

for(int j=0; j< seqb.size(); j++ ){
	for(int i=0; i< seq.size(); i++ ){
	
	const string ast = seq[i] ;
	const string bst = seqb[j] ;
	if (ast  == bst ) {p=p+1; break; }
	
	if (i == seq.size() -1 ) { n = n +1; negs << ">" << seq_nameb[j] << endl << seqb[j] << endl; }
	}
}
negs.close();
tpou.close();
tnou.close();

}
void Sensitivity( string& dc , string& dc2 , string& tar){
 ifstream dca( dc.c_str() );
    if ( !dca ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq ;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!dca.eof()){
		temp = "";

		getline(dca, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(dca, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
ifstream dcb( dc2.c_str() );
    if ( !dcb ) {
      
    }

	string tempb;
	vector <string> seq_nameb;
	vector <string> seqb ;
	
	string templb;
int temb = 0;
int tempmaxb = 0;
	seqb.clear();
	while(!dcb.eof()){
		tempb= "";

		getline(dcb, tempb);
		if(tempb.length() == 0){
			break;
		}
	
		string nameb (tempb, 1, tempb.length() - 1);
		seq_nameb.push_back(nameb);
		
		getline(dcb, tempb);
		seqb.push_back(tempb);
		temb = tempb.length();
				
		
	}
	int p =0;
	int n=0;
	vector <string> seq_namet;
	vector <string> seqt ;
	int flag =0;
	for(int i=0; i< seq.size(); i++ ){
		for(int j=0; j< seqb.size(); j++ ){
			  string ast = seq[i] ;                
			  string bstname=seq_nameb[j];         
		          string bsn= bstname.substr(bstname.size()-3,bstname.size()-1);    
			  string snam = seq_name[i].substr(1,3); 
			
		 if( bstname.find( snam )!=string::npos ){ 
			  
			if(  bsn =="wil" || bsn =="ana"){  
			
			 string bst = seqb[j] ;         
			size_t  found;
			 found=bst.find(ast);           
			if (found!=string::npos) {
				if(found > 4){seq_namet.push_back(seq_nameb[j]); seqt.push_back( seqb[j].substr(found-5,20) ); continue;}
				else { seq_namet.push_back(seq_nameb[j]); seqt.push_back( seqb[j].substr(0,20) ); continue; }
			};   
		 	int n;
		       char c;
			string ast2;
			  
		       for ( n=0 ; ast[n]!='\0' ; n++)
			 {
			    c=ast[n];
			 ast2+=toupper(c);               
			 }
			
		   if(!isupper(ast[1])){
			  
			found=bst.find(ast2);
		 	 if (found!=string::npos) { 
				if(found > 4){seq_namet.push_back(seq_nameb[j]); seqt.push_back( seqb[j].substr(found-5,20) ); continue;}
				else { seq_namet.push_back(seq_nameb[j]); seqt.push_back( seqb[j].substr(0,20) ); continue; }
			};
		   }  
		   Sequence s(ast2);
		
		   Sequence r = s.compRevCompl() ;	
			
		   ostringstream sst;              
			sst << r ;
			
		       string rr = sst.str();   
	
string ast3;
			
		       for ( n=0 ; rr[n]!='\0' ; n++)
			 {
			    c=rr[n];
			 ast3+=toupper(c);               
			 }

			found=bst.find(ast3);
		 	 if (found!=string::npos) {
				
				if(found > 4){seq_namet.push_back(seq_nameb[j]); seqt.push_back( seqb[j].substr(found-5,20) ); continue;}
				else { seq_namet.push_back(seq_nameb[j]); seqt.push_back( seqb[j].substr(0,20) ); continue; }
			};
		  } 
		 } 
		} 
	} 
ofstream fout( tar.c_str() );
    
  
        for ( int i = 0; i < seqt.size(); i++ ) {
            fout << ">" << seq_namet[ i ] << endl;
            fout << seqt[ i ] << endl;
        }
  fout.close();

}

void scorematrix(Matrix e , string& file, double maxLLR, string& file2 ) 
{
    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int l = e.nRows() ;  

ofstream ener( file2.c_str() );
		  for( int j = 0; j < seq.size() ; j++ ){
			
					
					
					Sequence readseq(seq[j]);
				
			
				    double result = 0;
				    for ( int i = 0; i < l; i++ ) {
					if( i < readseq.size() ){if( readseq[i] < 4 && readseq[i] >= 0 ){
				       result += e( i, readseq[ i ] ); 
				    	}}
				   }
				double energy = 0;
				energy = maxLLR - result;
				
			ener<< seq_name[j] << " & "<< seq[j] << " & " << fixed << setprecision(2) << energy <<"\\\\" << endl;  
		    }

		
		
ener.close();


}

void scorematrixw(Matrix e , string& file, double maxLLR, string& file2 ) 
{
    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int l = e.nRows() ;  

ofstream ener( file2.c_str() );
		  for( int j = 0; j < seq.size() ; j++ ){
			
					
					
					Sequence readseq(seq[j]);
				
			
				    double result = 0;
				    for ( int i = 0; i < l; i++ ) {
					if( i < readseq.size() ){if( readseq[i] < 4 && readseq[i] >= 0 ){
				       result += e( i, readseq[ i ] ); 
				    	}}
				   }
				double energy = 0;
				energy = result;
				
			ener<< seq_name[j] << " & "<< seq[j] << " & " << fixed << setprecision(2) << energy <<"\\\\" << endl;  
			
		    }

		
		
ener.close();


}



Matrix countmatrixflip( const string& file )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  



ofstream otes("rcRead.txt");
	 Matrix m(4,length,0);
			

	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[j]);
		Sequence RCreadseq( readseq.compRevCompl() );  
		int na = 0;
		for(int i = 2;  i< readseq.size()-2; i++){
			if(readseq[i] == 0)            
				{ na = na + 1;}
		}
		int naRC = 0;
		for(int i = 2;  i< RCreadseq.size()-2; i++){
			if(RCreadseq[i] == 0)            
				{ naRC = naRC + 1;}
		}
		if( na > naRC ){
			otes << ">" << seq_name[j] << endl << readseq <<endl;
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		}  
		else {
			otes << ">" << seq_name[j] << endl << RCreadseq <<endl;
			for( int i = 0;  i< RCreadseq.size(); i++)  
			{
				if(RCreadseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(RCreadseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(RCreadseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(RCreadseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		} 


	} 
	Matrix countMatrix = m.transpose();
	

	 double pseudoCount =1;

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	otes.close();
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	   return pwm;
}

void countmatrix( const string& file, vector< string >& names, vector< Motif >& motifs,vector< double >& background)
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
	
		
		getline(seq_file, temp);
		seq.push_back(temp);
			seq_name.push_back(name + " " + temp);

				 Sequence s(temp);
		
		   Sequence r = s.compRevCompl() ;	
		
		   ostringstream sst;              
			sst << r ;
			
		       string rr = sst.str();   
			
		seq_name.push_back(name+ " " + rr);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  





	for( int j = 0; j < seq.size(); j++ ){
		Matrix m(4,length,0);
	 	Matrix m2(4,length,0);	
		Sequence readseq(seq[j]);
		Sequence RCreadseq( readseq.compRevCompl() );  
		if( 1){

			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		} 
		if( 1){
			
			for( int i = 0;  i< RCreadseq.size(); i++) 
			{
				if(RCreadseq[i] == 0)            
				{ m2(0,i) = m2(0,i) + 1;}
				if(RCreadseq[i] == 1)
				{ m2(1,i) = m2(1,i) + 1;}
				if(RCreadseq[i] == 2)
				{ m2(2,i) = m2(2,i) + 1;}
				if(RCreadseq[i] == 3)
				{ m2(3,i) = m2(3,i) + 1;}
				
			}
		} 
			Matrix countMatrix = m.transpose();
			Matrix countMatrix2 = m2.transpose();

	 double pseudoCount =1;
	
	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	  Matrix   pwm2 = compWtmx( countMatrix2, pseudoCount ) ;
		Motif dcmm( pwm, background);
		Motif dcmm2( pwm2, background);
		motifs.push_back( dcmm ); 
		motifs.push_back( dcmm2 ); 
	} 
names = seq_name;
}

double Motif::Dotp( Motif m1, Motif m2)
{
	double norm_info=0;
	double Expected_info=0;
	Matrix ma =	m1.getPwm();
	Matrix mb = 	m2.getPwm();
	 for ( int i = 0; i < ma.nRows(); i++ ) {
        for ( int j = 0; j < ma.nCols(); j++ ) {
            Expected_info = Expected_info  + ma.getElement( i, j ) * mb.getElement(i,j) ;	
		norm_info= norm_info + ma.getElement( i, j ) * ma.getElement(i,j) ;	
        }	
    	}
return double(Expected_info/norm_info);	
}


void countmatrix_unique( const string& file, vector< string >& names, vector< Motif >& motifs,vector< double >& background)
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
	
		
		getline(seq_file, temp);
		seq.push_back(temp);
			seq_name.push_back(name + " " + temp);

				 Sequence s(temp);
		
		   Sequence r = s.compRevCompl() ;	
		
		   ostringstream sst;              
			sst << r ;
			
		       string rr = sst.str();   
			
		seq_name.push_back(name+ " " + rr);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  





	for( int j = 0; j < seq.size(); j++ ){
		Matrix m(4,length,0);
	 	Matrix m2(4,length,0);	
		Sequence readseq(seq[j]);
		Sequence RCreadseq( readseq.compRevCompl() );  
		if( 1){

			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		} 
		if( 1){
			
			for( int i = 0;  i< RCreadseq.size(); i++) 
			{
				if(RCreadseq[i] == 0)            
				{ m2(0,i) = m2(0,i) + 1;}
				if(RCreadseq[i] == 1)
				{ m2(1,i) = m2(1,i) + 1;}
				if(RCreadseq[i] == 2)
				{ m2(2,i) = m2(2,i) + 1;}
				if(RCreadseq[i] == 3)
				{ m2(3,i) = m2(3,i) + 1;}
				
			}
		} 
			Matrix countMatrix = m.transpose();
			Matrix countMatrix2 = m2.transpose();

	 double pseudoCount =1;
	
	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	  Matrix   pwm2 = compWtmx( countMatrix2, pseudoCount ) ;
		Motif dcmm( pwm, background);
		Motif dcmm2( pwm2, background);
		motifs.push_back( dcmm ); 
		motifs.push_back( dcmm2 ); 
	} 
double eucdist =0;
vector< double > euc( motifs.size(),0 );
 vector< Motif > motifsexui;
 vector< Motif > motifsexuid;
 vector< string > motifsexuin;
 vector< string > motifsexuidn;
	int flag = 0;
	for( int i = 0 ; i < motifs.size() ; i++){
	
		flag = 0;
		for(;;){
			vector< double > euc( motifs.size(),0 );
			for( int j=0; j< motifs.size(); j++ ){
				Motif t = motifs[j];	
				eucdist=t.Dotp( motifs[j], motifs[i]);
				euc.push_back( eucdist );
			}
			double gass = 1.0;
			if(  count(euc.begin(),euc.end(), gass ) > 1   ){ flag =1; break; }
			break;
		}
		if( flag == 1 ){ motifsexuid.push_back( motifs[i] );  motifsexuidn.push_back( seq_name[i] ); continue; }  
		if( flag ==0 ){ motifsexui.push_back( motifs[i] ) ;  motifsexuin.push_back( seq_name[i] ); }
	}


      
   
 vector< Motif > motifsexui2 = motifsexui;

	
		for( int i = 0 ; i < motifsexuid.size() ; i++){  
	
		
		
			vector< double > euc( motifs.size(),0 );
			for( int j=0; j< motifsexui.size(); j++ ){      
				Motif t = motifs[j];	
				eucdist=t.Dotp( motifsexuid[i], motifsexui[j]);
				euc.push_back( eucdist );
				
			}
			double gass = 1.0;
			int ct = count(euc.begin(),euc.end(), gass ) ;
			if(  ct > 0   ){ continue; }
			if( ct ==0 ){ motifsexui.push_back( motifsexuid[i] ); motifsexuin.push_back( motifsexuidn[i] ); continue;}
			
			
		
		
		
	}

   

motifs.clear();
for( int j=0; j< motifsexui.size(); j++ ){  
motifs.push_back( motifsexui[j] );
names.push_back( motifsexuin[j] );
}


   


}



Matrix countmatrixfull( const string& file )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  



	 Matrix m(4,length,0);
			

	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[j]);
		Sequence RCreadseq( readseq.compRevCompl() );  
		int na = 0;
		int ng =0;
		for(int i = 0;  i< 2; i++){
			if(readseq[i] == 2)            
				{ ng = ng + 1;}
		}
		for(int i = 2;  i< readseq.size()-2; i++){
			if(readseq[i] == 0)            
				{ na = na + 1;}
		}
		int naRC = 0;
		int ngRC =0;
		for(int i = 0;  i< 2; i++){
			if(RCreadseq[i] == 2)            
				{ ngRC = ngRC + 1;}
		}
		for(int i = 2;  i< RCreadseq.size()-2; i++){
			if(RCreadseq[i] == 0)            
				{ naRC = naRC + 1;}
		}
		if( na >= naRC|| ng >=ngRC  ) {   
					
				int pos = 0;
				int flag = 0;
				for( int p = 0;  p< readseq.size(); p++)  
					{
						if(readseq[p] == 2){ pos =p; flag=1; break;}
					}  
				if( flag==0){ continue;} 
				
					
					
					
					for( int i = pos;  i< readseq.size(); i++)  
					{
						if(readseq[i] == 0)            
						{ m(0,i) = m(0,i) + 1;}
						if(readseq[i] == 1)
						{ m(1,i) = m(1,i) + 1;}
						if(readseq[i] == 2)
						{ m(2,i) = m(2,i) + 1;}
						if(readseq[i] == 3)
						{ m(3,i) = m(3,i) + 1;}
						
					}
				}  
		else {
			
			for( int i = 0;  i< RCreadseq.size(); i++)  
			{
				if(RCreadseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(RCreadseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(RCreadseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(RCreadseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		} 


	} 
	Matrix countMatrix = m.transpose();
	
	 double pseudoCount =1;
	
	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	   return pwm;
}



Matrix countmatrix( const string& file )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  
ofstream otes("rcRead.txt");



	 Matrix m(4,length,0);
	

	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[j]);

		Sequence RCreadseq( readseq.compRevCompl() );  

		int na = 0;
		for(int i = 3;  i< readseq.size()-3; i++){
			if(readseq[i] == 0)            
				{ na = na + 1;}
		}
		int naRC = 0;
		for(int i = 3;  i< RCreadseq.size()-3; i++){
			if(RCreadseq[i] == 0)            
				{ naRC = naRC + 1;}
		}
		

		if( 1 ){
			
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		}  

	} 

	Matrix countMatrix = m.transpose();
	
otes.close();
	 double pseudoCount =1;
	
	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	   return pwm;
}


Matrix countmatrixSkip( const string& file )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  
ofstream otes("rcRead.txt");



	 Matrix m(4,length,0);
			

	for( int j = 0; j < seq.size(); j++ ){
		int flag=0;
		Sequence readseq(seq[j]);
		Sequence RCreadseq( readseq.compRevCompl() );  
		int na = 0;
		int ng =0;
		for(int i = 0;  i< 2; i++){
			if(readseq[i] == 2)            
				{ ng = ng + 1;}
		}
		for(int i = 2;  i< readseq.size()-2; i++){
			if(readseq[i] == 0)            
				{ na = na + 1;}
		}
		int naRC = 0;
		int ngRC =0;
		for(int i = 0;  i< 3; i++){
			if(RCreadseq[i] == 2)            
				{ ngRC = ngRC + 1;}
		}
		for(int i = 2;  i< RCreadseq.size()-2; i++){
			if(RCreadseq[i] == 0)            
				{ naRC = naRC + 1;}
		}
		if( na >= naRC && ng >=ngRC  ) {   
			otes << ">" << seq_name[j] << endl << readseq <<endl;
			int posG=0;
			for( int i = 0;  i< readseq.size(); i++)  
			{
			  if(readseq[i] == 2){ posG=i; break;}
			}
			for( int i = posG;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
			flag=1;
		}  
		if( flag==0 && na >= naRC ){
			otes << ">" << seq_name[j] << endl << readseq <<endl;
			int posG=0;
			for( int i = 0;  i< readseq.size(); i++)  
			{
			  if(readseq[i] == 2){ posG=i; break;}
			}
			for( int i = posG;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
			flag=1;
		}
		if( flag==0 && ng >=ngRC  ) {
			otes << ">" << seq_name[j] << endl << RCreadseq <<endl;
			int posG=0;
			for( int i = 0;  i< readseq.size(); i++)  
			{
			  if(readseq[i] == 2){ posG=i; break;}
			}
			for( int i = posG;  i< RCreadseq.size(); i++)  
			{
				if(RCreadseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(RCreadseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(RCreadseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(RCreadseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
			flag=1;	
		} 
		if( flag==0  ) {
			
		} 


	} 
	Matrix countMatrix = m.transpose();
	
otes.close();
	 double pseudoCount =1;
	
	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	   return pwm;
}

Matrix countmatrixtw( const string& file )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  




	 Matrix m(4,length,0);
			

	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[j]);
		
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
	
		


	} 
	Matrix countMatrix = m.transpose();
	
	 double pseudoCount =1;
	
	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	   return pwm;
}


Matrix countmat_int( const string& file )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  



	 Matrix m(4,length,0);
			

	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[j]);
		Sequence RCreadseq( readseq.compRevCompl() );  
		int na = 0;
		for(int i = 2;  i< readseq.size()-2; i++){
			if(readseq[i] == 0)            
				{ na = na + 1;}
		}
		int naRC = 0;
		for(int i = 2;  i< RCreadseq.size()-2; i++){
			if(RCreadseq[i] == 0)            
				{ naRC = naRC + 1;}
		}
		if( na >= naRC|| length < 10){
			
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		}  
		else {
			
			for( int i = 0;  i< RCreadseq.size(); i++)  
			{
				if(RCreadseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(RCreadseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(RCreadseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(RCreadseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		} 


	} 
	Matrix countMatrix = m.transpose();
	
	 double pseudoCount =1;
	
	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	   return countMatrix ;
}

Matrix randSites2Matrix( const string& file )
{
std::srand(std::time(0));
    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
        
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
		
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  
	
	int size = seq.size();
	ofstream ener( "Data\energydc.txt" );

	vector <int> indices(size, 0);
	for(int i = 0; i < indices.size(); i++)
		indices[i] = i;
	
	

	random_shuffle(indices.begin(), indices.end(), Random);
ofstream otes("Data\rcReaddc.txt");



	 Matrix m(4,length,0);
			

	for( int j = 0; j < seq.size()/2; j++ ){

		Sequence readseq(seq[indices[j]]);
		Sequence RCreadseq( readseq.compRevCompl() );  
		int na = 0;
		for(int i = 2;  i< readseq.size()-2; i++){
			if(readseq[i] == 0)            
				{ na = na + 1;}
		}
		int naRC = 0;
		for(int i = 2;  i< RCreadseq.size()-2; i++){
			if(RCreadseq[i] == 0)            
				{ naRC = naRC + 1;}
		}
		if( na >= naRC|| length < 10){
			otes << ">" << seq_name[indices[j]] << endl << readseq <<endl;
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		}  
		else {
			otes << ">" << seq_name[indices[j]] << endl << RCreadseq <<endl;
			for( int i = 0;  i< RCreadseq.size(); i++)  
			{
				if(RCreadseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(RCreadseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(RCreadseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(RCreadseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		} 


	} 


	Matrix countMatrix = m.transpose();
	
	 double pseudoCount =.0001;

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	    Matrix pwm( l, 4 );

	    
	    for ( int i = 0; i < l; i++ ) {
		double n = 0;       
		for ( int j = 0; j < 4; j++ ) {
		    n += countMatrix( i, j );
		}
		for ( int j = 0; j < 4; j++ ) {
		    pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
		}	
	    }

	vector< double > backgro(4,.25);
		 Matrix LLRMat(l,4,0);	
	    	Sequence maxSite4;	
	   	 double maxLLR;	
	    
	    for ( int i = 0; i < l; i++ ) {
		for ( int j = 0; j < 4; j++ ) {			
		    LLRMat( i, j ) = log( pwm( i, j ) ) ; 
		}
	    }
	     
	    
	    for ( int i = 0; i < l; i++ ) {
		int b_max;
		max( pwm.getRow( i ), b_max );
		maxSite4.push_back( b_max );	
	    }
	  
	    
	    maxLLR = 0;
	    for ( int i = 0; i < l; i++ ) {
		maxLLR += LLRMat( i, maxSite4[ i ] );	
	    }

	vector< double > e;
	
		  for( int j = 0; j < seq.size()/2 ; j++ ){
			
					
					
					Sequence readseq(seq[indices[j]]);
				
			
				    double result = 0;
				    for ( int i = 0; i < l; i++ ) {
					if( i < readseq.size() )
				       result += LLRMat( i, readseq[ i ] ); 
				    }
			
				double energy = maxLLR - result;
				e.push_back( energy );
		    }
		ener <<  fixed << setprecision(2)<< e <<  endl;


ener.close();
return pwm;
}


Matrix countmatrixerr( const string& file )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  




	 Matrix m(4,length,1);
			

	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[j]);

		for( int i = 0;  i< readseq.size(); i++)  
		{
			if(readseq[i] == 0)            
			{ m(0,i) = m(0,i) + 1;}
			if(readseq[i] == 1)
			{ m(1,i) = m(1,i) + 1;}
			if(readseq[i] == 2)
			{ m(2,i) = m(2,i) + 1;}
			if(readseq[i] == 3)
			{ m(3,i) = m(3,i) + 1;}
			
		}

	} 
	Matrix countMatrix = m.transpose();
	

	 double pseudoCount =1;

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	   return pwm;
}
void Motif::countmatrix3( string& file, string& file2, string& file3  )
{
int l = pwm.nRows();
vector< double > freqs;
vector< double > enersa;
freqs.clear();
enersa.clear();
int nindex = int(maxLLR-minLLR);


vector< vector<double> > bins(int(maxLLR-minLLR)+1, vector< double >(0));
int index =0;
double sumbin=0;
int kk = 0;
for( int p0 = 0 ; p0 < 4 ; p0++){
double pvalf0=1.0;
double enerf0=0;
pvalf0 = background[p0];
enerf0 = LLRMat(0,p0);

	for( int p1 = 0 ; p1 < 4 ; p1++){
	double pvalf1 = pvalf0*background[p1];
	double enerf1 = enerf0 + LLRMat(1,p1);

		for( int p2 = 0 ; p2 < 4 ; p2++){
		double pvalf2 = pvalf1*background[p2];
		double enerf2 = enerf1 + LLRMat(2,p2);
			for( int p3 = 0 ; p3 < 4 ; p3++){
			double pvalf3 = pvalf2*background[p3];
			double enerf3 = enerf2 + LLRMat(3,p3);
				for( int p4 = 0 ; p4 < 4 ; p4++){
			       double  pvalf4 = pvalf3*background[p4];
				double enerf4 = enerf3 + LLRMat(4,p4);
					for( int p5 = 0 ; p5 < 4 ; p5++){
					double pvalf5 = pvalf4*background[p5];
					double enerf5 = enerf4 + LLRMat(5,p5);
				for( int p6 = 0 ; p6 < 4 ; p6++){
				double pvalf6 = pvalf5*background[p6];
				double enerf6 = enerf5 + LLRMat(6,p6);
			for( int p7 = 0 ; p7 < 4 ; p7++){
			double pvalf7 = pvalf6*background[p7];
			double enerf7 = enerf6 + LLRMat(7,p7);
		for( int p8 = 0 ; p8 < 4 ; p8++){
		double pvalf8 = pvalf7*background[p8];
		double enerf8 = enerf7 + LLRMat(8,p8);
	for( int p9 = 0 ; p9 < 4 ; p9++){
	double pvalf9 = pvalf8*background[p9];
	double enerf9 = enerf8 + LLRMat(9,p9);
for( int p10 = 0 ; p10 < 4 ; p10++){
double pvalf10 = pvalf9*background[p10];
double enerf10 = enerf9 + LLRMat(10,p10);
freqs.push_back( pvalf10 );

index = int( maxLLR-enerf10);
bins[index].push_back( pvalf10 );
enersa.push_back( maxLLR - enerf10 );  
}
  }
    }
}
}
}
}
}
}
}
}
vector< double >bins2;
double z=0;
for(int e=0;e< bins.size() ; e++ ){  
double num = 0;
	for(int i=0; i<bins[e].size();i++){
			num=num + bins[e][i];
			}
z=z+num;
bins2.push_back( num );
}

double ethres =0;
vector< double > ne(bins2.size());
int ntemp =0;

	double pvalf=0;

	
	ofstream cdfbackgrounddc("Data/cdfbackgrounddca.txt");
	double pvalfb=0;
	for(int i = 0 ; i < bins2.size() ; i++){  

			pvalfb = bins2[i]/z + pvalfb;	
			cdfbackgrounddc << i << '\t' << pvalfb<< endl ;
			ne[i]=pvalfb;
				
				
			
				
				
				
				
			
	}
	
	cdfbackgrounddc.close();
	ofstream pvalfile("Data/backgroundCountsdca.txt");
	
	
	for(int i = 0 ; i < bins2.size() ; i++ ){
	
	pvalfile << bins2[i] << endl;
	
	}
	pvalfile.close();

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
        
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem >= tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  


	    	Sequence maxSite3;	
	   	
	    for ( int i = 0; i < l; i++ ) {
		int b_max;
		max( pwm.getRow( i ), b_max );
		maxSite3.push_back( b_max );	
	    }
	


	vector< double > e;
	
ofstream ener( file2.c_str() );
ofstream pvala( file3.c_str() );
		  for( int j = 0; j < seq.size() ; j++ ){
			
					
					
					Sequence readseq(seq[j]);
				
			
				    double result = 0;
				    for ( int i = 0; i < l; i++ ) {
					if( i < readseq.size() ){
		
						if( readseq[i] < 4 && readseq[i] >= 0 ){
					       		result += LLRMat( i, readseq[ i ] ); 
						}
						
						if( readseq[i] ==4 ){
					       		result += LLRMat( i, maxSite3[ i ] ); 
						}
							
					 } 
				        } 
					
				double energy = 0;
				energy = maxLLR - result;
				e.push_back( energy );
			
			ener << fixed << setprecision(5) << energy << endl; 
			pvala  << int(energy) << '\t' << ne[int(energy)] << endl;
			
		    }

pvala.close();
ener.close();

}


Matrix countmatrix2(string& file, string& file2 )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
        
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem >= tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  


int size = seq.size();
std::srand(std::time(0));


	vector <int> indices(size, 0);
	for(int i = 0; i < indices.size(); i++)
		indices[i] = i;



	


	 Matrix m(4,length,0);
			
	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[indices[j]]);
			
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
	} 
	Matrix countMatrix = m.transpose();
	
	 
 double pseudoCount =1;

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	    Matrix pwm( l, 4,0 );
	    
	    for ( int i = 0; i < l; i++ ) {
		double n = 0;       
		for ( int j = 0; j < 4; j++ ) {
		    n += countMatrix( i, j );
			
		}
		for ( int j = 0; j < 4; j++ ) {
		    pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
		}	
	    }
	vector< double > backgro(4,.25);
		 Matrix LLRMat(l,4,0);	
	    	Sequence maxSite3;	
	   	 double maxLLR;	
	    
	    for ( int i = 0; i < l; i++ ) {
		for ( int j = 0; j < 4; j++ ) {			
		    LLRMat( i, j ) = log( pwm( i, j ) );
			
		}
	    }
	   
	    
	    for ( int i = 0; i < l; i++ ) {
		int b_max;
		max( pwm.getRow( i ), b_max );
		maxSite3.push_back( b_max );	
	    }
	  
	    
	    maxLLR = 0;
	    for ( int i = 0; i < l; i++ ) {
		
		maxLLR += LLRMat( i, maxSite3[ i ] );	
		
		
	    }
	










	
	
	



	


	vector< double > e;
	
ofstream ener( file2.c_str() );
		  for( int j = 0; j < seq.size() ; j++ ){  
			
					
					Sequence readseq(seq[indices[j]]);
				    double result = 0;
				    for ( int i = 0; i < l; i++ ) {
					if( i < readseq.size() ){
		
						if( readseq[i] < 4 && readseq[i] >= 0 ){
					       		result += LLRMat( i, readseq[ i ] ); 
						}
						
						if( readseq[i] ==4 ){
					       		result += LLRMat( i, maxSite3[ i ] ); 
						}
							
					 } 
				        } 
					
				double energy = 0;
				energy = maxLLR - result;
				e.push_back( energy );
			
			ener << fixed << setprecision(2) << energy << endl; 
		    }

		
		
ener.close();
return pwm;
}


Matrix countmatrix2w(string& file, string& file2 , Motif _motif)
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
        
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem >= tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  


int size = seq.size();
std::srand(std::time(0));


	vector <int> indices(size, 0);
	for(int i = 0; i < indices.size(); i++)
		indices[i] = i;

 
    Matrix LLRcb=_motif.getLLRMatpos();
	Matrix pwmCB =_motif.getPwm();
	double minLLRcb =_motif.getMinLLR();
	double maxLLRcb =_motif.getMaxLLR();
    
				
   


	


	 Matrix m(4,length,0);
			
	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[indices[j]]);
			
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
	} 
	Matrix countMatrix = m.transpose();
	
	 
 double pseudoCount =1;

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	    Matrix pwm( l, 4,0 );
	    
	    for ( int i = 0; i < l; i++ ) {
		double n = 0;       
		for ( int j = 0; j < 4; j++ ) {
		    n += countMatrix( i, j );
			
		}
		for ( int j = 0; j < 4; j++ ) {
		    pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
		}	
	    }
	vector< double > backgro(4,.25);
		 Matrix LLRMat(l,4,0);	
	    	Sequence maxSite3;	
	   	 double maxLLR;	
	    
	    for ( int i = 0; i < l; i++ ) {
		for ( int j = 0; j < 4; j++ ) {			
		    LLRMat( i, j ) = -(-maxLLRcb/double(l)  +LLRcb(i,j)+  log( pwmCB( i, j )/pwm(i,j) )  ) ;   
			
		}
	    }
	   
	    
	    for ( int i = 0; i < l; i++ ) {
		int b_max;
		max( pwm.getRow( i ), b_max );
		maxSite3.push_back( b_max );	
	    }
	  
	    
	    maxLLR = 0;
	    for ( int i = 0; i < l; i++ ) {
		
		maxLLR += LLRMat( i, maxSite3[ i ] );	
		
		
	    }
	










	
	
	



	


	vector< double > e;
	
ofstream ener( file2.c_str() );
		  for( int j = 0; j < seq.size() ; j++ ){  
			
					
					Sequence readseq(seq[indices[j]]);
				    double result = 0;
				    for ( int i = 0; i < l; i++ ) {
					if( i < readseq.size() ){
		
						if( readseq[i] < 4 && readseq[i] >= 0 ){
					       		result += LLRMat( i, readseq[ i ] ); 
						}
						
						if( readseq[i] ==4 ){
					       		result += LLRMat( i, maxSite3[ i ] ); 
						}
							
					 } 
				        } 
					
				double energy = 0;
				energy = maxLLRcb - result;
				e.push_back( energy );
			
			ener << fixed << setprecision(2) << energy << endl; 
		    }

		
		
ener.close();
return pwm;
}


Matrix readDUSites2Matrix( const string& file)
{
std::srand(std::time(0));
    ifstream seq_file( file.c_str() );
 

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	int tem = 0;
	int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
		
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
			

				tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}


int l = tempmax ;  

	 Matrix m(4,l,0);	

	
		
	
	
	int size = seq.size();

	vector <int> indices(size, 0);
	for(int i = 0; i < indices.size(); i++)
		indices[i] = i;
	
	


	random_shuffle(indices.begin(), indices.end(), Random);

	for( int j = 0; j < seq.size()/3; j++ ){

		Sequence readseq(seq[indices[j]]);

		for( int i = 0;  i< readseq.size(); i++)  
		{
			if(readseq[i] == 0)            
			{ m(0,i) = m(0,i) + 1;}
			if(readseq[i] == 1)
			{ m(1,i) = m(1,i) + 1;}
			if(readseq[i] == 2)
			{ m(2,i) = m(2,i) + 1;}
			if(readseq[i] == 3)
			{ m(3,i) = m(3,i) + 1;}
		}

	} 

	 double pseudoCount =1;

	Matrix countMatrix = m.transpose();

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	   return pwm;

}

int readSites2Matrix( const string& file )
{
std::srand(std::time(0));
    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
        return RET_ERROR;
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	
	seq.clear();
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
		
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
			
		
	}
	
	int size = seq.size();
	ofstream ener( "Data/energy.txt" );
ofstream she("Data/beauty.txt");  
	vector <int> indices(size, 0);
	for(int i = 0; i < indices.size(); i++)
		indices[i] = i;
	
	

for( int s =0; s < 100; s++){
	random_shuffle(indices.begin(), indices.end(), Random);
	 Matrix m(4,11,0);
			
	for( int j = 0; j < seq.size()/3; j++ ){

		Sequence readseq(seq[indices[j]]);

		for( int i = 0;  i< readseq.size(); i++)  
		{
			if(readseq[i] == 0)            
			{ m(0,i) = m(0,i) + 1;}
			if(readseq[i] == 1)
			{ m(1,i) = m(1,i) + 1;}
			if(readseq[i] == 2)
			{ m(2,i) = m(2,i) + 1;}
			if(readseq[i] == 3)
			{ m(3,i) = m(3,i) + 1;}
			
		}
	
	
	} 
	Matrix countMatrix = m.transpose();
	
	 double pseudoCount =1;

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	    Matrix pwm( l, 4 );

	    
	    for ( int i = 0; i < l; i++ ) {
		double n = 0;       
		for ( int j = 0; j < 4; j++ ) {
		    n += countMatrix( i, j );
		}
		for ( int j = 0; j < 4; j++ ) {
		    pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
		}	
	    }
	
	vector< double > backgro(4,.25);
	Motif c(countMatrix, pseudoCount,backgro);
	she << c ;
		 Matrix LLRMat(l,4,0);	
	    	Sequence maxSite4;	
	   	 double maxLLR;	
	    
	    for ( int i = 0; i < l; i++ ) {
		for ( int j = 0; j < 4; j++ ) {			
		    LLRMat( i, j ) = log( pwm( i, j ) / backgro[ j ] );
		}
	    }
	     
	    
	    for ( int i = 0; i < l; i++ ) {
		int b_max;
		max( pwm.getRow( i ), b_max );
		maxSite4.push_back( b_max );	
	    }
	  
	    
	    maxLLR = 0;
	    for ( int i = 0; i < l; i++ ) {
		maxLLR += LLRMat( i, maxSite4[ i ] );	
	    }

	vector< double > e;
	
		  for( int j = 0; j < seq.size()/2 ; j++ ){
			
					
					
					Sequence readseq(seq[indices[j]]);
				
			
				    double result = 0;
				    for ( int i = 0; i < l; i++ ) {
					if( i < readseq.size() )
				       result += LLRMat( i, readseq[ i ] ); 
				    }
			
				double energy = maxLLR - result;
				e.push_back( energy );
		    }
		ener <<  fixed << setprecision(2)<< e <<  endl;


Matrix m2(4,11,0);

for( int j = seq.size()/3; j < seq.size(); j++ ){
        Sequence readseq(seq[indices[j]]);

	for( int i = 0;  i< readseq.size(); i++)  
	{
		if(readseq[i] == 0)            
		{ m2(0,i) = m2(0,i) + 1;}
		if(readseq[i] == 1)
		{ m2(1,i) = m2(1,i) + 1;}
		if(readseq[i] == 2)
		{ m2(2,i) = m2(2,i) + 1;}
		if(readseq[i] == 3)
		{ m2(3,i) = m2(3,i) + 1;}
		
	}
	
} 
Matrix countMatrix2 = m2.transpose();


    assert( countMatrix2.nCols() == 4 && pseudoCount >= 0 );
    
    int l2 = countMatrix2.nRows();		
    Matrix pwm2( l2, 4 );

    
    for ( int i = 0; i < l2; i++ ) {
        double n = 0;       
        for ( int j = 0; j < 4; j++ ) {
            n += countMatrix2( i, j );
        }
        for ( int j = 0; j < 4; j++ ) {
            pwm2( i, j ) = ( countMatrix2( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
        }	
    }

Motif cc(countMatrix2, pseudoCount,backgro);
	she << cc ;
	 Matrix LLRMat2(l2,4,0);	
    	Sequence maxSite2;	
   	 double maxLLR2;	
    
    for ( int i = 0; i < l2; i++ ) {
        for ( int j = 0; j < 4; j++ ) {			
            LLRMat2( i, j ) = log( pwm2( i, j ) / backgro[ j ] );
        }
    }
     
    
    for ( int i = 0; i < l2; i++ ) {
        int b_max2;
        max( pwm2.getRow( i ), b_max2 );
        maxSite2.push_back( b_max2 );	
    }
  
    
    maxLLR2 = 0;
    for ( int i = 0; i < l2; i++ ) {
        maxLLR2 += LLRMat2( i, maxSite2[ i ] );	
    }

vector< double > e2;

	  for( int j = seq.size()/2; j < seq.size(); j++ ){
		
				
				
				Sequence readseq2(seq[indices[j]]);
			if ( seq[j].empty() ) continue;
			
			    double result2 = 0;
			    for ( int i = 0; i < l2; i++ ) {
				if( i < readseq2.size() )
			       result2 += LLRMat2( i, readseq2[ i ] ); 
			    }
			
			double energy2 = maxLLR2 - result2;
			e2.push_back( energy2 );
	    }
	  ener <<  e2 <<  endl;
}
ener.close();
she.close();
return 0;
}


int readSites2Matrixw( const string& file, int len, Motif _motif )
{
std::srand(std::time(0));
    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
        return RET_ERROR;
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
		string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();


while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

map< string, int > speciesIdxMap2seq;
	   
speciesIdxMap2seq["mel"] =   1;
speciesIdxMap2seq["sim"] =    1;
speciesIdxMap2seq["sec"] =   1;
speciesIdxMap2seq["ere"] =  1;
speciesIdxMap2seq["yak"] =    1;
speciesIdxMap2seq["ana"] =  1;
speciesIdxMap2seq["pse"] =   1;
speciesIdxMap2seq["per"] =  1;
speciesIdxMap2seq["vir"] =  1;
speciesIdxMap2seq["moj"] =  1;
speciesIdxMap2seq["gri"] = 1;
speciesIdxMap2seq["wil"] =    1;

	   
int length = tempmax ;  

	int size = seq.size();
	ofstream ener( "energy.txt" );

ofstream ener2( "energy2.txt" );
	vector <int> indices(size, 0);
	for(int i = 0; i < indices.size(); i++)
		indices[i] = i;





 
    Matrix LLRcb=_motif.getLLRMatpos();
	Matrix pwmCB =_motif.getPwm();
	double minLLRcb =_motif.getMinLLR();
	double maxLLRcb =_motif.getMaxLLR();
    
				
   






	
	
for( int s =0; s < 1000; s++){
	random_shuffle(indices.begin(), indices.end(), Random);

	 Matrix m(4,length,0);
			
	
	for( int j = 0; j <  len; j++ ){
		 string bstname=seq_name[j];         
			
		  string bsn= bstname.substr(bstname.size()-3,bstname.size()-1);    
		
			
		Sequence readseq(seq[indices[j]]);
			
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;} 
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;} 
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1; }
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1; }
				
			}
	
	


	} 


	Matrix countMatrix = m.transpose();
	
	 double pseudoCount =.1;

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	    Matrix pwm( l, 4 );

	    
	    for ( int i = 0; i < l; i++ ) {
		double n = 0;       
		for ( int j = 0; j < 4; j++ ) {
		    n += countMatrix( i, j );
		}
		for ( int j = 0; j < 4; j++ ) {
		    pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
		}	
	    }

	vector< double > backgro(4,.25);
		 Matrix LLRMat(l,4,0);	
	    	Sequence maxSite;	
	   	 double maxLLR;	
	    
	    for ( int i = 0; i < l; i++ ) {
		for ( int j = 0; j < 4; j++ ) {			
		    LLRMat( i, j ) =-(-maxLLRcb/double(l)  +LLRcb(i,j)+  log( pwmCB( i, j )/pwm(i,j) )  ) ;
		}
	    }
	     
	    
	    for ( int i = 0; i < l; i++ ) {
		int b_max;
		max( pwm.getRow( i ), b_max );
		maxSite.push_back( b_max );	
	    }
	  
	    
	    maxLLR = 0;
	    for ( int i = 0; i < l; i++ ) {
		maxLLR += LLRMat( i, maxSite[ i ] );	
	    }

	vector< double > e;
	
		  for( int j = 0; j < len ; j++ ){
			
					
					
					Sequence readseq(seq[indices[j]]);
				
			
				    double result = 0;
				    for ( int i = 0; i < l; i++ ) {
					if( i < readseq.size() )
				       result += LLRMat( i, readseq[ i ] ); 
				    }
			
				double energy = maxLLRcb - result;
				e.push_back( energy );
		    }
		ener << e <<  endl;

Matrix m2(4,l,0);

for( int j = len; j < seq.size(); j++ ){
        Sequence readseq(seq[indices[j]]);

	for( int i = 0;  i< readseq.size(); i++)  
	{
		if(readseq[i] == 0)            
		{ m2(0,i) = m2(0,i) + 1;}
		if(readseq[i] == 1)
		{ m2(1,i) = m2(1,i) + 1;}
		if(readseq[i] == 2)
		{ m2(2,i) = m2(2,i) + 1;}
		if(readseq[i] == 3)
		{ m2(3,i) = m2(3,i) + 1;}
		
	}
	
} 
Matrix countMatrix2 = m2.transpose();


    assert( countMatrix2.nCols() == 4 && pseudoCount >= 0 );
    
    int l2 = countMatrix2.nRows();		
    Matrix pwm2( l2, 4 );

    
    for ( int i = 0; i < l2; i++ ) {
        double n = 0;       
        for ( int j = 0; j < 4; j++ ) {
            n += countMatrix2( i, j );
        }
        for ( int j = 0; j < 4; j++ ) {
            pwm2( i, j ) = ( countMatrix2( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
        }	
    }

	 Matrix LLRMat2(l2,4,0);	
    	Sequence maxSite2;	
   	 double maxLLR2;	
    
    for ( int i = 0; i < l2; i++ ) {
        for ( int j = 0; j < 4; j++ ) {			
            LLRMat2( i, j ) = -(-maxLLRcb/double(l)  +LLRcb(i,j)+  log( pwmCB( i, j )/pwm(i,j) )  );
        }
    }
     
    
    for ( int i = 0; i < l2; i++ ) {
        int b_max2;
        max( pwm2.getRow( i ), b_max2 );
        maxSite2.push_back( b_max2 );	
    }
  
    
    maxLLR2 = 0;
    for ( int i = 0; i < l2; i++ ) {
        maxLLR2 += LLRMat2( i, maxSite2[ i ] );	
    }

vector< double > e2;

	  for( int j = len; j < seq.size(); j++ ){
		
				
				
				Sequence readseq2(seq[indices[j]]);
			if ( seq[j].empty() ) continue;
			
			    double result2 = 0;
			    for ( int i = 0; i < l2; i++ ) {
				if( i < readseq2.size() )
			       result2 += LLRMat2( i, readseq2[ i ] ); 
			    }
			
			double energy2 = maxLLRcb - result2;
			e2.push_back( energy2 );
	    }
	  ener2 <<  e2 <<  endl;

}















ofstream ener3( "energy3.txt" );
 Matrix m(4,length,0);
			

	for( int j = 0; j <  len; j++ ){
		 string bstname=seq_name[j];         
			
		  string bsn= bstname.substr(bstname.size()-3,bstname.size()-1);    
		
			
		Sequence readseq(seq[j]);
			
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1 ; }
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i)  + 1 ; }
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1 ; }
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i)  + 1 ; }
				
			}
	
	


	} 


	Matrix countMatrix = m.transpose();
	
	 double pseudoCount =.1;

	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	    Matrix pwm( l, 4 );

	    
	    for ( int i = 0; i < l; i++ ) {
		double n = 0;       
		for ( int j = 0; j < 4; j++ ) {
		    n += countMatrix( i, j );
		}
		for ( int j = 0; j < 4; j++ ) {
		    pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
		}	
	    }

	vector< double > backgro(4,.25);
		 Matrix LLRMat(l,4,0);	
	    	Sequence maxSite;	
	   	 double maxLLR;	
	    
	    for ( int i = 0; i < l; i++ ) {
		for ( int j = 0; j < 4; j++ ) {			
		    LLRMat( i, j ) = log( pwm( i, j ) / backgro[ j ] );
		}
	    }
	     
	    
	    for ( int i = 0; i < l; i++ ) {
		int b_max;
		max( pwm.getRow( i ), b_max );
		maxSite.push_back( b_max );	
	    }
	  
	    
	    maxLLR = 0;
	    for ( int i = 0; i < l; i++ ) {
		maxLLR += LLRMat( i, maxSite[ i ] );	
	    }

	vector< double > e;
	
		  for( int j = 0; j < seq.size() ; j++ ){
			
					
					
					Sequence readseq(seq[j]);
				
			
				    double result = 0;
				    for ( int i = 0; i < l; i++ ) {
					if( i < readseq.size() )
				       result += LLRMat( i, readseq[ i ] ); 
				    }
			
				double energy = maxLLR - result;
				e.push_back( energy );
		    }
		ener3 <<  fixed << setprecision(2)<< e <<  endl;








ener3.close();

ener2.close();
ener.close();
return 0;
}
int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, vector< string >& names, bool readEnergy )
{
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    sites.clear();
    names.clear();

    SiteVec currVec;
    int nrecords = 0;       
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;

        if ( line.substr( 0, 1 ) == ">" ) {
            stringstream ss( line.substr( 1 ) );
            string name; 
            ss >> name;
            names.push_back( name );
            nrecords++;
            if ( nrecords > 1 ) {
                sites.push_back( currVec );
                currVec.clear();
            }
        } else {
            int start;
            char strandChar;
            string factor;
            double energy = 0;
            stringstream ss( line );
            ss >> start >> strandChar >> factor;
            if ( readEnergy ) ss >> energy; 
            bool strand = strandChar == '+' ? 1 : 0;
            map<string, int>::const_iterator iter = factorIdxMap.find( factor );
            currVec.push_back( Site( start - 1, strand, iter->second , energy ) );
        }
    }

    sites.push_back( currVec );

    return 0;
}

int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, bool readEnergy )
{
    vector< string > names;
    return readSites( file, factorIdxMap, sites, names, readEnergy );
}


int SeqAnnotator::annot( const Sequence& seq, SiteVec& sites ) const
{
    sites.clear();
    
    
    for ( int i = 0; i < seq.size(); i++ ) {
        
        for ( int k = 0; k < motifs.size(); k++ ) {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) continue;
            double energy;
            
            
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
            if ( energy <= energyThrs[ k ] ) {  
                sites.push_back( Site( i, 1, k, energy ) );
            }	
            
            
            Sequence rcElem( seq, i, l, 0 );
            energy = motifs[ k ].energy( rcElem );
            if ( energy <= energyThrs[ k ] ) {
                sites.push_back( Site( i, 0, k, energy ) );
            }				
        }	
    }
    
    return sites.size();
}
bool siteSortPredicate(const Site& d1, const Site& d2)
{
  return d1.start < d2.start;
}
bool siteSortPredicate2(double d1, double d2)
{
  return d1 < d2;
}
    
int SeqAnnotator::annoty4( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) 
{
vector< double > cons = f.getCol( column );
tsites.clear();
Site t;
vector< double > ff(3);  
       typedef  vector< Site > sitestt;
	sitestt te;	
       vector< sitestt > sitest( motifs.size() );    
       vector< sitestt > sitesp( motifs.size() );    
	  for ( int k = 0; k < 2; k++ ) {
		
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }

	          }
		   
		  if( k == 0 ) {
	double mine =100;
	int gbest;
	int counter=1;
	sitestt  s = sitest[0];
	tsites.clear();
	
	while( counter > 0 && s.size() > 0 ){
			
			counter=0;
			Site tempsite;
			 
		
 			for ( int i = 0; i < s.size(); i++ ) {
	
				if ( mine >= s[i].energy ){  
				mine = s[i].energy;
				tempsite = s[i];
				gbest = i;
				counter++;
 				}
			}
	
	
	
			tsites.push_back (tempsite);
			sitestoverlap( s, tsites ); 
			
			
          } 
	
	std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
	  p.predictOcc( tsites , 5, cons, ff);
			if ( ff[0] >.5) {  
			
			
			  		
			   return 1; 
			}
			else{ tsites.clear();}
	  }  
	  }  


SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			
			kbest = ktempindex;
		}
		}
}
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 

vector< double  > sitest_tw_Z( 0 );
 
int N=1;
double previous = 0;
double current = 0;
do {

tsites.clear();
for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  
p.predictOcc( tsites, 5, cons, ff);
previous = ff[0];  
if( previous > .5 )
{ break ; }
N++;
tsites.push_back( siteMax( sitest[0] ) );  

	for( int i = 0; i < sitest[1].size(); i++){  
		
		sitesorderedtemp.push_back( sitest[1][i] )  ;  
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   
int ind = maxZ( sitest_tw_Z );	
double maxZscore = maxZs( sitest_tw_Z );
	sitesp[1].push_back( sitest[1][ ind ] ); 
 	sitestoverlap( sitest[1], sitesp[1] );
if (N > 6) break;

tsites.clear();

for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0]; 
if( current > .5 ) break;
} while ( abs(previous - current ) > .25 );
return tsites.size();

}








int SeqAnnotator::annotydorsalold(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) 
{
vector< double > cons = f.getCol( column );
tsites.clear();
Site t;
vector< double > ff(3);  
if( cons[2] ==0 ) {   
       typedef  vector< Site > sitestt;
	sitestt te;	
       vector< sitestt > sitest( motifs.size() );    
       vector< sitestt > sitesp( motifs.size() );    
	  for ( int k = 0; k < 2; k++ ) {
		
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }

	          }
		   
		  
	  }  

SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			
			kbest = ktempindex;
		}
		}
}
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 

vector< double  > sitest_tw_Z( 0 );

int N=1;
double previous = 0;
double current = 0;
do {

tsites.clear();
for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  
p.predictOcc( tsites, 5, cons, ff);
previous = ff[0];  
if( previous > .5 )
{ break ; }
N++;
tsites.push_back( siteMax( sitest[0] ) );  

	for( int i = 0; i < sitest[1].size(); i++){  
		
		sitesorderedtemp.push_back( sitest[1][i] )  ;  
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   
int ind = maxZ( sitest_tw_Z );	
double maxZscore = maxZs( sitest_tw_Z );
	sitesp[1].push_back( sitest[1][ ind ] ); 
 	sitestoverlap( sitest[1], sitesp[1] );
if (N > 6) break;

tsites.clear();

for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0]; 
if( current > .5 ) break;
} while ( abs(previous - current ) > .25 );

return tsites.size();


} 


else {

 typedef  vector< Site > sitestt;	
       vector< sitestt > sitest( motifs.size() );    
       vector< sitestt > sitesp( motifs.size() );    
	  for ( int k = 0; k < motifs.size(); k++ ) {
		
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		 if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		}
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		 if ( energy <= energyThrs[ k ] ) {	
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		}

	          }
	

 
		 if( k == 0 ) {
	double mine =100;
	int gbest;
	int counter=1;
	sitestt  s = sitest[0];
	tsites.clear();
	
	
	
	while( counter > 0 && s.size() > 0 ){
			
			counter=0;
			Site tempsite;
			 
	
	
 			for ( int i = 0; i < s.size(); i++ ) {
	
				if ( mine >= s[i].energy ){  
				mine = s[i].energy;
				tempsite = s[i];
				gbest = i;
				counter++;
	
 				}
			}
	
	
			if( counter != 0) {
			tsites.push_back (tempsite);
			sitestoverlap( s, tsites );
			} 
	
			
          } 
	 
	std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
	  p.predictOcc( tsites , 5, cons, ff);
			if ( ff[0] >.5) {  
			   
			  s.clear();		
			   return 1; 
			}
			else{ s.clear(); tsites.clear();}
	  }  
	  } 


SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			
			kbest = ktempindex;
		}
		}
}
double Zbestdd=0;
int ibestdd =0;
int kbestdd = 0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = k+1; i < sitest[0].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
		if( Zbest < Ztemp ) {
			Zbestdd = Ztemp;
			ibestdd = itempindex;			
			kbestdd = ktempindex;
		}
		}
}
if( Zbestdd < Zbest ) {
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) {  
			  tsites= sitesp[0];
			   return 1; }
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 
}
else{                                     
sitesp[0].push_back( sitest[0][kbestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
sitesp[0].push_back( sitest[0][ibestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) {  
			  tsites= sitesp[0];
			   return 1; }
}
sitesp[2].push_back (siteMax( sitest[2] ));
sitestoverlap( sitest[2], sitesp[2] );

int N = 1; 
double previous;
double current;

do {
vector< double > ff(3); 
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}

N++;
vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); 
SiteVec sitesorderedtemp = tsites;
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitest[k].size(); i++){
		sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[k][i] )  ;  
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   
vector< double > uniq = sitest_allk_alli_Z[k];  
std::unique(uniq.begin(),uniq.end() );
if ( uniq.size() == 1 ) {
 continue; 
}
int ind = maxZ( sitest_allk_alli_Z[k] );	
if ( cons[k] != 0 ) {
	sitesp[k].push_back( sitest[k][ ind ] ); 
 	sitestoverlap( sitest[k], sitesp[k] );
}
}  
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0];  
if( current > .5 )
{ break ; }
  
if(N ==6) { 
break; 
} 


} while( N < 3 ) ;

tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);

return sitest.size();

}


}




int SeqAnnotator::annotydorsal(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) 
{
vector< double > cons = f.getCol( column );
tsites.clear();
Site t;
vector< double > ff(3);  
if( 1 ) {  
    typedef  vector< Site > sitestt;
    sitestt te;	
    vector< sitestt > sitest( motifs.size() );    
    vector< sitestt > sitesp( motifs.size() );    
    for ( int k = 0; k < 2; k++ ) {
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	    int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }
                   } 

		     if( k == 0 ) {
			for ( int i = 0; i < seq.size(); i++ ) {
			      	   int l = d2.length();
			 	    if ( i + l > seq.size() ) break ;
			 	    double energy;
				    Sequence elem( seq, i, l, 1 );
			 	    energy = d2.energy( elem );
				    if ( energy <= energyThrs[ 0 ] ) {	
			 	    sitest[0].push_back( Site( i, 1, 0, energy ) );
				    }
		     		    Sequence rcElem( seq, i, l, 0 );
		      	      	    energy = d2.energy( rcElem );
				    if ( energy <= energyThrs[ 0] ) {
				    sitest[0].push_back( Site( i, 0, 0, energy ) );
				    }
			  }  
			  double mine =100;
			  int gbest;
			  int counter=1;
			  sitestt  s = sitest[k];
			  tsites.clear();
			  while( counter > 0 && s.size() > 0 ){
					counter=0;
					Site tempsite;
		 			for ( int i = 0; i < s.size(); i++ ) {
						if ( mine >= s[i].energy ){  
						mine = s[i].energy;
						tempsite = s[i];
						gbest = i;
						counter++;
		 				}
					}
					if( counter != 0) {
					tsites.push_back (tempsite);
					sitestoverlap( s, tsites );
					}
			  } 
			  mine =100;
			  gbest = 100;
			  counter=1;
			  while( counter > 0 && s.size() > 0 ){  
					counter=0;
					Site tempsite;
		 			for ( int i = 0; i < s.size(); i++ ) {
						if ( mine >= s[i].energy ){  
						mine = s[i].energy;
						tempsite = s[i];
						gbest = i;
						counter++;
		 				}
					}
					if( counter != 0) {
					tsites.push_back (tempsite);
					sitestoverlap( s, tsites );
					} 
			  } 
			  s.clear();
			} 
			else{
				double mine =100;
				int gbest;
				int counter=1;
				sitestt  s = sitest[k];
				while( counter > 0 && s.size() > 0 ){
					counter=0;
					Site tempsite;
		 			for ( int i = 0; i < s.size(); i++ ) {
						if ( mine >= s[i].energy ){  
						mine = s[i].energy;
						tempsite = s[i];
						gbest = i;
						counter++;
		 				}
					}
					if( counter != 0) {
					tsites.push_back (tempsite);
					sitestoverlap( s, tsites );
					} 
			  	} 
				s.clear();
			} 

	} 

	std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
	p.predictOcc( tsites , 5, cons, ff);

	if ( ff[0] >.5) {  
		return 1; 
	}
	else{ tsites.clear();}
	SiteVec sitesorderedtemp(0);
	double Ztemp=0;
	int ktempindex=0;
	int itempindex=0;
	double Zbest=0;
	int ibest =0;
	int  kbest=0;
	for ( int k = 0; k < sitest[0].size(); k++ ) {
			ktempindex = k;
	     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
			sitesorderedtemp.push_back( sitest[0][k] );
			sitesorderedtemp.push_back( sitest[1][i] );
			std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
			Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
			sitesorderedtemp.clear();
			itempindex=i;
			if( Zbest < Ztemp ) {
				Zbest = Ztemp;
				ibest = itempindex;			
				kbest = ktempindex;
			}
			}
	} 
	sitesp[0].push_back( sitest[0][kbest] );
	sitestoverlap( sitest[0], sitesp[0] ); 
	
	sitesp[1].push_back( sitest[1][ibest] );
	sitestoverlap( sitest[1], sitesp[1] ); 
	vector< double  > sitest_tw_Z( 0 );
	int N=1;
	double previous = 0;
	double current = 0;
	do {
		tsites.clear();
		for ( int k = 0; k < 2; k++ ) {
			for( int i = 0; i < sitesp[k].size(); i++){
				tsites.push_back( sitesp[k][i] );
			}
		}
		std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  
		p.predictOcc( tsites, 5, cons, ff);
		previous = ff[0];  
		if( previous > .5 )
		{ break ; }
		N++;
		tsites.push_back( siteMax( sitest[0] ) );  
			for( int i = 0; i < sitest[1].size(); i++){  
				
				sitesorderedtemp.push_back( sitest[1][i] )  ;  
				std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
				sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
				sitesorderedtemp.clear();
			}   
		int ind = maxZ( sitest_tw_Z );	
		
		double maxZscore = maxZs( sitest_tw_Z );
		
			sitesp[1].push_back( sitest[1][ ind ] ); 
		 	sitestoverlap( sitest[1], sitesp[1] );
		if (N > 6) break;
		tsites.clear();
		
		for ( int k = 0; k < 2; k++ ) {
			for( int i = 0; i < sitesp[k].size(); i++){
				tsites.push_back( sitesp[k][i] );
			}
		}
		std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
		p.predictOcc( tsites, 5, cons, ff);
		current = ff[0]; 
		if( current > .5 ) break;
		} 
	while ( abs(previous - current ) > .25 );

	return tsites.size();
} 
else {
typedef  vector< Site > sitestt;	
       vector< sitestt > sitest( motifs.size() );    
       vector< sitestt > sitesp( motifs.size() );    
	  for ( int k = 0; k < motifs.size(); k++ ) {
     		  for ( int i = 0; i < seq.size(); i++ ) {
		      	   int l = motifs[ k ].length();
		 	    if ( i + l > seq.size() ) break ;
		 	    double energy;
			    Sequence elem( seq, i, l, 1 );
		 	    energy = motifs[ k ].energy( elem );
			 if ( energy <= energyThrs[ k ] ) {	
		 	    sitest[k].push_back( Site( i, 1, k, energy ) );
			}
	     		    Sequence rcElem( seq, i, l, 0 );
	      	      	    energy = motifs[ k ].energy( rcElem );
			 if ( energy <= energyThrs[ k ] ) {	
			    sitest[k].push_back( Site( i, 0, k, energy ) );
			}
	          }

	  } 

std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
	  p.predictOcc( tsites , 5, cons, ff);
			if ( ff[0] >.5) {  
			  
			 
			   return 1; 
			}
			else{ tsites.clear();}

SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			
			kbest = ktempindex;
		}
		}
} 
double Zbestdd=0;
int ibestdd =0;
int kbestdd = 0;
if( sitest[0].size() > 1 ) {
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = k+1; i < sitest[0].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
		if( Zbest < Ztemp ) {
			Zbestdd = Ztemp;
			ibestdd = itempindex;			
			kbestdd = ktempindex;
		}
		}
} 
if( Zbestdd < Zbest ) {
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) { 
			  tsites= sitesp[0];
			   return 1; }
} 
else{                                     

sitesp[0].push_back( sitest[0][kbestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
sitesp[0].push_back( sitest[0][ibestdd] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) { 
			  tsites= sitesp[0];
			   return 1; }
}
}

int N = 1; 
double previous;
double current;

do {
vector< double > ff(3); 
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
N++;
vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); 
SiteVec sitesorderedtemp = tsites;
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitest[k].size(); i++){
		sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[k][i] )  ;  
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   
if(sitest[k].size() > 0  ) {
int ind = maxZ( sitest_allk_alli_Z[k] );	
if ( cons[k] != 0 ) {
	sitesp[k].push_back( sitest[k][ ind ] ); 
 	sitestoverlap( sitest[k], sitesp[k] );
}
}
}  
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
vector< double > ff2(3);
p.predictOcc( tsites, 5, cons, ff2);

current = ff[0];  
if( current > .5 )
{ break ; }
  
if(N ==6) { 
break; 
} 
 

} while( N < 3 ) ;

tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
return sitest.size();

}


}





int SeqAnnotator::annotydorsalZ(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) 
{
vector< double > cons = f.getCol( column );
tsites.clear();
Site t;
vector< double > ff(3);  


typedef  vector< Site > sitestt;	
       vector< sitestt > sitest( motifs.size() );    
       vector< sitestt > sitesp( motifs.size() );    
	  for ( int k = 0; k < 2; k++ ) {
     		 for ( int i = 0; i < seq.size(); i++ ) {
		      	   int l = motifs[ k ].length();
		 	    if ( i + l > seq.size() ) break ;
		 	    double energy;
			    Sequence elem( seq, i, l, 1 );
		 	    energy = motifs[ k ].energy( elem );
			 if ( energy <= energyThrs[ k ] ) {	
		 	    sitest[k].push_back( Site( i, 1, k, energy ) );
			}
	     		    Sequence rcElem( seq, i, l, 0 );
	      	      	    energy = motifs[ k ].energy( rcElem );
			 if ( energy <= energyThrs[ k ] ) {	
			    sitest[k].push_back( Site( i, 0, k, energy ) );
			}
	         }
		
	  } 


SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			
			kbest = ktempindex;
		}
		}
} 
double Zbestdd=0;
int ibestdd =0;
int kbestdd = 0;
if( sitest[0].size() > 1 ) {
for ( int k = 0; k < sitest[0].size(); k++ ) {
		ktempindex = k;
     		  for ( int i = k+1; i < sitest[0].size(); i++ ) {
		  	if( !siteOverlap( sitest[0][k] ,sitest[0][i] , this->motifs) )  { 
				sitesorderedtemp.push_back( sitest[0][k] );
				sitesorderedtemp.push_back( sitest[0][i] );
				std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
				Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
				sitesorderedtemp.clear();
				itempindex=i;

				if( Zbest < Ztemp ) {
					Zbestdd = Ztemp;
					ibestdd = itempindex;			
					kbestdd = ktempindex;
				}
			} 
		}

} 
	if( Zbestdd < Zbest && Zbest != 0) {
		sitesp[0].push_back( sitest[0][kbest] );
		sitestoverlap( sitest[0], sitesp[0] ); 
		p.predictOcc( sitesp[0] , 5, cons, ff);
					if ( ff[0] >.5) { 
					  tsites= sitesp[0];
					   return 1; }
		sitesp[1].push_back( sitest[1][ibest] );
		sitestoverlap( sitest[1], sitesp[1] ); 

		} 
	else{                                     
		

		sitesp[0].push_back( sitest[0][kbestdd] );
		sitestoverlap( sitest[0], sitesp[0] ); 
		
			if( ibestdd != 0 ) {
			sitesp[0].push_back( sitest[0][ibestdd] );
	
			}
	 }
	}

tsites.clear();

		for ( int k = 0; k < 2; k++ ) {
			for( int i = 0; i < sitesp[k].size(); i++){
				tsites.push_back( sitesp[k][i] );
			}
		}
		std::sort(tsites.begin(), tsites.end(), siteSortPredicate); 


return tsites.size();





std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
return sitest.size();




}


int SeqAnnotator::annotydorsalZ2( vector< Sequence > seq, vector< vector< vector< Site > > >& sitespp , Matrix f, Matrix e, ExprFunc& p, int column, vector< string >  names) 
{
double ethr = .1*column + .0001 ;
Site t;
vector< double > ff(3);  
 vector< vector< vector< Site > > > sitespt(seq.size(), vector< vector< Site > >(4));
typedef  vector< Site > sitestt;	
 
 
	for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
     		 for ( int i = 0; i < seq[e].size(); i++ ) {
		      	   int l = motifs[ k ].length();
		 	    if ( i + l > seq[e].size() ) break ;
		 	    double energy;
			    Sequence elem( seq[e], i, l, 1 );
		 	    energy = motifs[ k ].energy( elem );
				
			 
			 if ( energy <= energyThrs[k] ) {	
				
		 	    sitespt[e][k].push_back( Site( i, 1, k, energy ) );
			}
	     		    Sequence rcElem( seq[e], i, l, 0 );
	      	      	    energy = motifs[ k ].energy( rcElem );
				
			 if ( energy <= energyThrs[k] ) {	
				
			    sitespt[e][k].push_back( Site( i, 0, k, energy ) );
			}
	         }
		
	  } 
	} 
 
	for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		if( sitespt[e][k].size() > 0 ) {
		sitestoverlap( sitespt[e][k], sitespp[e][k] );
		}
	  }
	}
	for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		if( sitespt[e][k].size() > 0 ) {
		sitespp[e][k].push_back( siteMax(  sitespt[e][k] ) );
	
		}
	  }
	}

	p.predictTwist( sitespp );
	


return 1;





}

void SeqAnnotator::modifyEthres( vector< double >&t ) { 
energyThrs=t;
} 
int SeqAnnotator::annotydorsalZarend( vector< Sequence > seq, vector< vector< vector< Site > > >& sitespp , Matrix f, Matrix e, ExprFunc& p, int column, vector< string >  names) 
{
Site t;
vector< double > ff(3);  
 vector< vector< vector< Site > > > sitespt(seq.size(), vector< vector< Site > >(4));
typedef  vector< Site > sitestt;	

 
 
for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
		sitespp[e][k].clear() ;  
		sitespt[e][k].clear();
		
	  }
	}

  for ( int k = 0; k < motifs.size(); k++ ) {  
		
}	
	
	for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
     		 for ( int i = 0; i < seq[e].size(); i++ ) {
		      	   int l = motifs[ k ].length();
		 	    if ( i + l > seq[e].size() ) break ;
		 	    double energy;
			    Sequence elem( seq[e], i, l, 1 );
		 	    energy = motifs[ k ].energy( elem );
			 if ( energy <= energyThrs[ k ] ) {
			
			
		 	    sitespt[e][k].push_back( Site( i, 1, k, energy ) );  
			}
	     		    Sequence rcElem( seq[e], i, l, 0 );
	      	      	    energy = motifs[ k ].energy( rcElem );
			if ( energy <= energyThrs[k] ) {	
			
			    sitespt[e][k].push_back( Site( i, 0, k, energy ) ); 
			}
	         }
		
	  } 
	} 
	 
ofstream ecounts("Data/enum.txt");
for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
		int num = 0;
		while( !sitespt[e][k].empty()){
			if( sitespt[e][k].size() > 0 ) {
			sitespp[e][k].push_back(  siteMax(  sitespt[e][k]  ));  
			num = num +1;
			sitestoverlap( sitespt[e][k], sitespp[e][k] );  
			
			}
		}
		ecounts<< e << '\t' << k << '\t' << num <<endl;;
		
		
		
	  }
	}


ecounts.close();
	ofstream dusty("Data/stranddu.txt");
	ofstream dcsty("Data/stranddc.txt");
	ofstream cbstr("Data/strandcb.txt");
	ofstream twstr("Data/strandtw.txt");
for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
		for( int i=0; i < sitespp[e][k].size() ; i++ ){
		Site a = sitespp[e][k][i];
		if( k ==1){ dusty <<a.strand << '\t';}
		if( k ==0){ dcsty  << a.strand << '\t';}
		if( k ==2){ twstr << a.strand << '\t';}
		if( k ==3){ cbstr << a.strand << '\t';}
		
		}
	  }
	}
dusty.close();
dcsty.close();
cbstr.close();
twstr.close();
p.predictTwist( sitespp );
p.predictTwist2( sitespp );
return 1;




}
int SeqAnnotator::annotydorsalZarend2u( vector< Sequence > seq, vector< vector< vector< Site > > >& sitespp , Matrix f, Matrix e, ExprFunc& p, int column, vector< string >  names, vector< vector< Site > >& seqSitesbot, vector< vector< Site > >& seqSitesm1)
{
Site t;
vector< double > ff(3);  
 vector< vector< vector< Site > > > sitespt(seq.size(), vector< vector< Site > >(4));
typedef  vector< Site > sitestt;	

 
 
for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
		sitespp[e][k].clear() ;  
		sitespt[e][k].clear();
		
	  }
	}

  for ( int k = 0; k < motifs.size(); k++ ) {  
		
}	


Sequence aaaa("GGAATTTCC");

Sequence bbbb("AAAAAAAAA");







	
ofstream cou("Data/couMI.txt");
int co = 0; 
	for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
     		 for ( int i = 0; i < seq[e].size(); i++ ) {
		      	   int l = motifs[ k ].length();
			    if( k==0 ){ co = co + 1; }
		 	    if ( i + l > seq[e].size() ) break ;
		 	    double energy;
			    Sequence elem( seq[e], i, l, 1 );
		 	    energy = motifs[ k ].energy( elem );
				 
			 if ( energy <= energyThrs[ k ] ) {
			
			
		 	    sitespt[e][k].push_back( Site( i, 1, k, energy ) );  
			}
	     		    Sequence rcElem( seq[e], i, l, 0 );
	      	      	    energy = motifs[ k ].energy( rcElem );
			if ( energy <= energyThrs[k] ) {	
			
			    sitespt[e][k].push_back( Site( i, 0, k, energy ) ); 
			}
	         }
		
	  } 
	} 
	
	 
for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
		
		while( !sitespt[e][k].empty()){
			if( sitespt[e][k].size() > 0 ) {
			sitespp[e][k].push_back(  siteMax(  sitespt[e][k]  ));  
			
			sitestoverlap( sitespt[e][k], sitespp[e][k] );  
			
			}
		}
		
		
		
		
	  }
	}
cou.close();

p.predictTwist( sitespp );
p.predictTwist2( sitespp );
return 1;

}

int SeqAnnotator::annotydorsalZarend2( vector< Sequence > seq, vector< vector< vector< Site > > >& sitespp , Matrix f, Matrix e, ExprFunc& p, int column, vector< string >  names, vector< vector< Site > >& seqSitesbot) 
{
Site t;
vector< double > ff(3);  
 vector< vector< vector< Site > > > sitespt(seq.size(), vector< vector< Site > >(4));
typedef  vector< Site > sitestt;	

 
 
for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
		sitespp[e][k].clear() ;  
		sitespt[e][k].clear();
		
	  }
	}

  for ( int k = 0; k < motifs.size(); k++ ) {  
		
}	
	
	for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
     		 for ( int i = 0; i < seq[e].size(); i++ ) {
		      	   int l = motifs[ k ].length();
		 	    if ( i + l > seq[e].size() ) break ;
		 	    double energy;
			    Sequence elem( seq[e], i, l, 1 );
		 	    energy = motifs[ k ].energy( elem );
			 if ( energy <= energyThrs[ k ] ) {
			
			
		 	    sitespt[e][k].push_back( Site( i, 1, k, energy ) );  
			}
	     		    Sequence rcElem( seq[e], i, l, 0 );
	      	      	    energy = motifs[ k ].energy( rcElem );
			if ( energy <= energyThrs[k] ) {	
			
			    sitespt[e][k].push_back( Site( i, 0, k, energy ) ); 
			}
	         }
		
	  } 
	} 
	 
ofstream ecounts("enum.txt");
for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		
		int num = 0;
		while( !sitespt[e][k].empty()){
			if( sitespt[e][k].size() > 0 ) {
			sitespp[e][k].push_back(  siteMax(  sitespt[e][k]  ));  
			num = num +1;
			sitestoverlap( sitespt[e][k], sitespp[e][k] );  
			
			}
		}
		ecounts<< e << '\t' << k << '\t' << num <<endl;;
		
		
		
	  }
	}

vector< vector< vector< Site > > > sitesptex(seq.size(), vector< vector< Site > >(4));
for( int e = 0; e < seq.size(); e++ ) {
	sitesptex[e][0] = seqSitesbot[e];  
}

for( int e = 0; e < seq.size(); e++ ) {
sitestoverlap( sitesptex[e][0], sitespp[e][0] );  
}
for ( int e = 0; e < seq.size(); e++ ) {  
	  for ( int k = 0; k < motifs.size(); k++ ) {  
		if( k == 0 ) {
				
		
					
				
				
			while( !sitesptex[e][k].empty()){
				if( sitesptex[e][k].size() > 0 ) {
				sitespp[e][k].push_back(  siteMax(   sitesptex[e][0]   ));  
				
				sitestoverlap( sitesptex[e][k], sitespp[e][k] );  
				
				}
			}
		
			
			
			
		}
	  }
	}


ecounts.close();

p.predictTwist2( sitespp );
return 1;




}


 int SeqAnnotator::annoty2( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names) 
{

vector< double > cons = f.getCol( column );
tsites.clear();
Site t;
vector< double > ff(3);  
if(cons[2] ==0 ) {   
       typedef  vector< Site > sitestt;
	sitestt te;	
       vector< sitestt > sitest( motifs.size() );    
       vector< sitestt > sitesp( motifs.size() );    
	  for ( int k = 0; k < 2; k++ ) {
		
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }

	          }
		   
		  if( k == 0 ) {
			tsites.clear();
			
			tsites.push_back (siteMax( sitest[0] ));
			p.predictOcc( tsites , 5, cons, ff);
			if ( ff[0] >.5) {  
			   return 1; }
			else{ tsites.clear();}
		  }
	  }


SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		

		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			
			kbest = ktempindex;
		}
		}
}
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 


vector< double  > sitest_tw_Z( 0 );
 
int N=1;
double previous = 0;
double current = 0;
do {

tsites.clear();
for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}

std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  
p.predictOcc( tsites, 5, cons, ff);
previous = ff[0];  
if( previous > .5 )
{ break ; }
N++;
tsites.push_back( siteMax( sitest[0] ) );  

	for( int i = 0; i < sitest[1].size(); i++){  
		
		sitesorderedtemp.push_back( sitest[1][i] )  ;  
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   

int ind = maxZ( sitest_tw_Z );	
double maxZscore = maxZs( sitest_tw_Z );
	sitesp[1].push_back( sitest[1][ ind ] ); 
 	sitestoverlap( sitest[1], sitesp[1] );
if (N > 6) break;

tsites.clear();

for ( int k = 0; k < 2; k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0]; 
if( current > .5 ) break;
} while ( abs(previous - current ) > .25 );

return tsites.size();


}










else {

 typedef  vector< Site > sitestt;	
       vector< sitestt > sitest( motifs.size() );    
       vector< sitestt > sitesp( motifs.size() );    
	  for ( int k = 0; k < motifs.size(); k++ ) {
		
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		 if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		}
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		 if ( energy <= energyThrs[ k ] ) {	
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		}

	          }
	

 
		  if( k == 0 ) {

tsites.clear();
			
tsites.push_back (siteMax( sitest[0] ));


			
			p.predictOcc( tsites , 5, cons, ff);
			if ( ff[0] >.5) {  
			  
			   return 1; }
			else{ tsites.clear();}
		  }  
	  } 

SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
		
		ktempindex = k;
     		  for ( int i = 0; i < sitest[1].size(); i++ ) {
		sitesorderedtemp.push_back( sitest[0][k] );
		sitesorderedtemp.push_back( sitest[1][i] );
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
		sitesorderedtemp.clear();
		itempindex=i;
		
		if( Zbest < Ztemp ) {
			Zbest = Ztemp;
			ibest = itempindex;			
			kbest = ktempindex;
		}
		}

}

sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] ); 
p.predictOcc( sitesp[0] , 5, cons, ff);
			if ( ff[0] >.5) {  
			  tsites= sitesp[0];
			   return 1; }

sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] ); 
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
return sitest.size();

int N = 1; 
double previous;
double current;

do {
vector< double > ff(3); 
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
N++;
vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); 
SiteVec sitesorderedtemp = tsites;
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitest[k].size(); i++){
		sitesorderedtemp = tsites;
		sitesorderedtemp.push_back( sitest[k][i] )  ;  
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
	        sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
		sitesorderedtemp.clear();
	}   
vector< double > uniq = sitest_allk_alli_Z[k];  
std::unique(uniq.begin(),uniq.end() );
if ( uniq.size() == 1 ) {
 continue; 
}
int ind = maxZ( sitest_allk_alli_Z[k] );	
if ( cons[k] != 0 ) {
	sitesp[k].push_back( sitest[k][ ind ] ); 
 	sitestoverlap( sitest[k], sitesp[k] );
}
}  
tsites.clear();
for( int k = 0; k < motifs.size(); k++ ){
sitestoverlap( sitest[k], sitesp[k] );      
}
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){
		tsites.push_back( sitesp[k][i] );
	}
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0];  
if( current > .5 )
{ break ; }
  
if(N ==6) { 
break; 
} 


} while( N < 3 ) ;

tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);

return sitest.size();

}
}


int SeqAnnotator::annoty3( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites) 
{


vector< double > cons = f.getCol( column );
tsites.clear();
vector< double > ff(3);  

       typedef  vector< Site > sitestt;
       vector< sitestt > sitest( motifs.size() );    
       vector< sitestt > sitesp( motifs.size() );    
	{ int k =2;
		
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		    if ( energy <= energyThrs[ k ] ) {	
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		    }
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		    if ( energy <= energyThrs[ k ] ) {
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		    }
		  }
	  
double mine =100;
	int gbest;
	int counter=1;
	sitestt  s = sitest[k];
sitestt  sp;
sp.clear();
Site sns;
for( int i =0; i < seqSites.size() ; i++ ) {
	sns = seqSites[i];
	if( sns.factorIdx == 2 ) {
	  sp.push_back( sns );
	}
}
sitestoverlap( sitest[2], sp );
sitestoverlap( s, sp );
sp.clear();
	
	
	
	
	while( counter > 0 && s.size() > 0 ){
			
			counter=0;
			Site tempsite;
			 
	
	
 			for ( int i = 0; i < s.size(); i++ ) {
	
				if ( mine >= s[i].energy ){  
				mine = s[i].energy;
				tempsite = s[i];
				gbest = i;
				counter++;
				
 				}
			}
	
	
			if( counter != 0) {
			
			sp.push_back (tempsite);
			sitestoverlap( s, sp );
			} 
	
			
          } 
		s.clear();
	 
for( int i=0; i< seqSites.size(); i++){
tsites.push_back( seqSites[i] );
}
for( int i=0; i< sp.size(); i++){
tsites.push_back( sp[i] );
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);


double	pre=p.predictExpr( tsites, 5, cons) ;

if( abs(pre - e(ti,column) ) < .2 )
{sp.clear(); return 1; }
		
			else{   sitestoverlap( sitest[2],sp ) ; sp.clear(); } 



	} 

int ktempindex=0;
double Etemp=0;
double diffE= 3;
double diffbest=3;
int  kbest=0;
for ( int k = 0; k < sitest[2].size(); k++ ) {
		SiteVec sitesorderedtemp = tsites;  
		ktempindex = k;
     		  
		sitesorderedtemp.push_back( sitest[2][k] );
		
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
		sitesorderedtemp.clear();
		ktempindex=k;
		diffE = abs( Etemp- 0 ) ;    

		if(diffbest > diffE ) {
			diffbest = diffE;
			kbest = ktempindex;			
			
		}
		
}
sitesp[2].clear();
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] ); 

int N=1;
double previous = 0;
double current = 0;
if( !sitest[2].empty())  {
do {

if( N == 1 ) {
 { int k = 2;
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
}

N++;
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  
previous=p.predictExpr( tsites, 5, cons) ;
if( abs(previous - e(ti,column) ) < .2 )
{ break ; }
if( previous < .2 )
{ break ; }



if( !sitest[2].empty() ){
for ( int k = 0; k < sitest[2].size(); k++ ) {
		SiteVec sitesorderedtemp = tsites;
		ktempindex = k;
     		  
		sitesorderedtemp.push_back( sitest[2][k] );
		
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
		sitesorderedtemp.clear();
		ktempindex=k;
		diffE = abs( Etemp-  0 );

		if(diffbest > diffE ) {
			diffbest = diffE;
			kbest = ktempindex;			
			
		}
		
}
tsites.push_back( sitest[2][kbest] ) ;

std::sort(tsites.begin(), tsites.end(), siteSortPredicate);


current = p.predictExpr(tsites,5,  cons); 
if( abs(current - e(ti,column) ) < .2 ) break;
if( current  < .2 ) break;
sitesp[2].clear();
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] );


} 
else break;
} while( N < 3  ); 
} 



return tsites.size();
}

int SeqAnnotator::annotyd( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites, vector< SiteVec >& dd) 
{

vector< double > cons = f.getCol( column );
tsites.clear();
vector< double > ff(3);  

       typedef  vector< Site > sitestt;
       vector< sitestt > sitest( motifs.size() );    
       vector< sitestt > sitesp( motifs.size() );    
	{ int k =2;
		
     		  for ( int i = 0; i < seq.size(); i++ ) {
              	   int l = motifs[ k ].length();
         	    if ( i + l > seq.size() ) break ;
         	    double energy;
        	    Sequence elem( seq, i, l, 1 );
         	    energy = motifs[ k ].energy( elem );
		   
         	    sitest[k].push_back( Site( i, 1, k, energy ) );
		   
     		    Sequence rcElem( seq, i, l, 0 );
      	      	    energy = motifs[ k ].energy( rcElem );
		  
	            sitest[k].push_back( Site( i, 0, k, energy ) );
		  
		  }
	  }

int ktempindex=0;
double Etemp=0;
double diffE= 3;
double diffbest=3;
int  kbest=0;
for ( int k = 0; k < sitest[2].size(); k++ ) {
		SiteVec sitesorderedtemp = seqSites;
		ktempindex = k;
     		  
		sitesorderedtemp.push_back( sitest[2][k] );
		
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
		sitesorderedtemp.clear();
		ktempindex=k;
		diffE = abs( Etemp- 0 ) ;    

		if(diffbest > diffE ) {
			diffbest = diffE;
			kbest = ktempindex;			
			
		}
		
}
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] ); 


int N=1;
double previous = 0;
double current = 0;

do {

tsites.clear();
tsites = seqSites;

 { int k = 2;
	for( int i = 0; i < sitesp[k].size(); i++){

		tsites.push_back( sitesp[k][i] );
		
	}

}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  
previous=p.predictExpr( tsites, 5, cons) ;

if( abs(previous - e(ti,column) ) < .2 )
{ break ; }
if( previous < .2 )
{ break ; }
N++;




for ( int k = 0; k < sitest[2].size(); k++ ) {
		SiteVec sitesorderedtemp = tsites;
		ktempindex = k;
     		  
		sitesorderedtemp.push_back( sitest[2][k] );
		
		std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);	
		Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
		sitesorderedtemp.clear();
		ktempindex=k;
		diffE = abs( Etemp- 0 ) ;

		if(diffbest > diffE ) {
			diffbest = diffE;
			kbest = ktempindex;			
			
		}
		
}
tsites.push_back( sitest[2][kbest] ) ;

std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
current = p.predictExpr(tsites,5,  cons); 
if( abs(current - e(ti,column) ) < .2 ) break;
if( current  < .2 ) break;
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] ); 

} while( N < 3 || abs(previous - current ) > .25 );
Delete1( tsites,dd );  
return tsites.size();
}

bool SeqAnnotator::DeltaZ( vector< vector< Site > >& sitespa , vector< double >& cons , double Zbefore, ExprFunc& p, double Obefore)
{
vector< double > f(0);
vector< Site > tsitess(0);
tsitess.clear();
sitespa[1].pop_back();  
for ( int k = 0; k < motifs.size(); k++ ) {
	for( int i = 0; i < sitespa[k].size(); i++){  
		tsitess.push_back( sitespa[k][i] );
	}
}
std::sort(tsitess.begin(), tsitess.end(), siteSortPredicate);
double Zafter = p.predictZ( tsitess, cons);
p.predictOcc(tsitess,5,cons,f);
double Oafter = f[0];
if ( (Zbefore - Zafter)/Zbefore < .1  && (Obefore - Oafter)/Obefore < .1) { sitespa.clear(); return true; }
sitespa.clear();
return false; 
}



int SeqAnnotator::Delete1( vector< Site >& t,vector< vector< Site > >& tt )
{



int size = t.size();

	
	
	
int count=0;
for( int i = 0; i < size; i++){  
		SiteVec temp1;
		temp1=t;
		vector< Site >::iterator ptr2 = temp1.begin();
		ptr2 += i;
		
		
		
		
		temp1.erase(ptr2  );
		tt.push_back( temp1 );
		
	 	
		
}
	
 


return 1;
}
int SeqAnnotator::maxZ( vector< double >& Z )
{
 int tempindex = 0 ;
 double maxZ=(*min_element(Z.begin(), Z.end() ) );
 for ( int i = 0; i < Z.size(); i++ ) {
	
	if ( maxZ < Z[i] ){
		maxZ = Z[i];
		tempindex = i;
 	}
  }
return tempindex;
}
double SeqAnnotator::maxZs( vector< double >& Z )
{
 int tempindex = 0 ;
 double maxZ=(*min_element(Z.begin(), Z.end() ) );
 for ( int i = 0; i < Z.size(); i++ ) {
	
	if ( maxZ < Z[i] ){
		maxZ = Z[i];
		tempindex = i;
 	}
  }
return maxZ;
}
void SeqAnnotator::sitestoverlap( vector< Site >& sitest, vector< Site >& sitesp )
{
for ( int i = 0; i < sitesp.size(); i++ ) {
	
        vector< Site >::iterator ptr2 = sitest.begin();
	int j = 0;
	
	
	while( ptr2 != sitest.end() ){
		if (siteOverlap( sitesp[i], sitest[j], this->motifs) )  { 
		
		sitest.erase(ptr2);
	
	
		continue;  
		}  
		j++;
		ptr2++;
	}
	
 }

}


Site SeqAnnotator::siteMax( vector< Site >& sitess)
{
 Site tempsite;
 double mine=100000;
 for ( int i = 0; i < sitess.size(); i++ ) {
	
	if ( mine > sitess[i].energy ){  
		mine = sitess[i].energy;
		tempsite = sitess[i];
 	}
  }
 return tempsite;	
}
int SeqAnnotator::compEnergy( const Sequence& seq, SiteVec& sites ) const
{
    for ( int i = 0; i < sites.size(); i++ ) {
        Sequence elem( seq, sites[i].start, motifs[sites[i].factorIdx].length(), sites[i].strand );
        sites[i].energy = motifs[sites[i].factorIdx].energy( elem );
        sites[i].wtRatio = exp( sites[i].energy-10 );  
    }

    return 0;
}


ModelType getModelOption( const string& modelOptionStr )
{
    if ( toupperStr( modelOptionStr ) == "LOGISTIC" ) return LOGISTIC;
    if ( toupperStr( modelOptionStr ) == "DIRECT" ) return DIRECT;
    if ( toupperStr( modelOptionStr ) == "QUENCHING" ) return QUENCHING;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_UNLIMITED" ) return CHRMOD_UNLIMITED;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_LIMITED" ) return CHRMOD_LIMITED;
    if ( toupperStr ( modelOptionStr ) == "BINS" ) return BINS;
    cerr << "modelOptionStr is not a valid model option" << endl; 
    exit(1);
}

string getModelOptionStr( ModelType modelOption )
{
    if ( modelOption == LOGISTIC ) return "Logisitic";
    if ( modelOption == DIRECT ) return "Direct";
    if ( modelOption == QUENCHING ) return "Quenching";
    if ( modelOption == CHRMOD_UNLIMITED ) return "ChrMod_Unlimited";
    if ( modelOption == CHRMOD_LIMITED ) return "ChrMod_Limited";
    if ( modelOption == BINS ) return "Bins";
    return "Invalid";
}

string getIntOptionStr( FactorIntType intOption )
{
    if ( intOption == BINARY ) return "Binary";
    if ( intOption == GAUSSIAN ) return "Gaussian";
    if ( intOption == BINSF ) return "Binsf";
    return "Invalid";
}

ObjType getObjOption( const string& objOptionStr )
{
    if ( toupperStr( objOptionStr ) == "SSE" ) return SSE;
    if ( toupperStr( objOptionStr ) == "CORR" ) return CORR;
    if ( toupperStr( objOptionStr ) == "CROSS_CORR" ) return CROSS_CORR;

    cerr << "objOptionStr is not a valid option of objective function" << endl; 
    exit(1);
}

string getObjOptionStr( ObjType objOption )
{
    if ( objOption == SSE ) return "SSE";
    if ( objOption == CORR ) return "Corr";
    if ( objOption == CROSS_CORR ) return "Cross_Corr";

    return "Invalid";
}

string getSearchOptionStr( SearchType searchOption )
{
    if ( searchOption == UNCONSTRAINED ) return "Unconstrained";
    if ( searchOption == CONSTRAINED ) return "Constrained";

    return "Invalid";
}
double FactorIntFuncBinsf::compFactorInt( double normalInt, double dist, bool orientation )  const
{
     assert( dist >= 0 );
     return normalInt;
}

double FactorIntFuncBinary::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );
	
    double spacingTerm = ( dist < distThr ? normalInt : 1.0 );
    double orientationTerm = orientation ? 1.0 : orientationEffect;	
    return spacingTerm * orientationTerm;	
}

double FactorIntFuncGaussian::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );

    double GaussianInt = dist < distThr ? normalInt * exp( - ( dist * dist ) / ( 2.0 * sigma * sigma ) ) : 1.0;
    return max( 1.0, GaussianInt );    
}

double FactorIntFuncGeometric::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );
	
    double spacingTerm = max( 1.0, dist <= distThr ? normalInt : normalInt * pow( spacingEffect, dist - distThr ) ); 
    double orientationTerm = orientation ? 1.0 : orientationEffect;	
    return spacingTerm * orientationTerm;
}
 
ExprPar::ExprPar( int _nFactors, vector< vector< vector<double> > > _theV , vector< vector< vector<double> > > _theVr): factorIntMat() , theV(_theV),theVr(_theVr)  
{	
    assert( _nFactors > 0 );
	
    for ( int i = 0; i < _nFactors; i++ ) {
        maxBindingWts.push_back( ExprPar::default_weight );	
    }	

    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );       
}
ExprPar::ExprPar( int _nFactors ) : factorIntMat()
{	
    assert( _nFactors > 0 );
	
    for ( int i = 0; i < _nFactors; i++ ) {
        maxBindingWts.push_back( ExprPar::default_weight );	
    }	

    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );       

   if (modelOption == BINS ) {theV = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>(ExprPar::nbins,1)));
		theVr = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>(ExprPar::nbins,1)));}

   for ( int i = 0; i < _nFactors; i++ ) {
        double defaultEffect = modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : ExprPar::default_effect_Thermo;
        txpEffects.push_back( defaultEffect );
       
    }

}
	
ExprPar::ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, double _basalTxp ) : maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat ), txpEffects( _txpEffects ), repEffects( _repEffects ), basalTxp( _basalTxp )
{
    if ( !factorIntMat.isEmpty() ) assert( factorIntMat.nRows() == maxBindingWts.size() && factorIntMat.isSquare() ); 	
   
  
}

ExprPar::ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) : factorIntMat()
{	
    int _nFactors = actIndicators.size();
    assert( coopMat.isSquare() && coopMat.nRows() == _nFactors );
    assert( repIndicators.size() == _nFactors );
    int counter = 0;
	
    
    if ( estBindingOption ) {
        for ( int i = 0; i < _nFactors; i++ ) {
            double weight = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_weight ), log( max_weight ) ) ) : exp( pars[counter++] );
            maxBindingWts.push_back( weight );
        }
    } else {
        for ( int i = 0; i < _nFactors; i++ ) maxBindingWts.push_back( ExprPar::default_weight );
    }
    
   
if (modelOption == BINS ) {

theV = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));

for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = 0; j < i; j++ ) {         
            if ( coopMat( i, j ) ) {
		for(int k=0; k< ExprPar::nbins; k++){  
			  double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interaction ), log( max_interaction ) ) ) : exp( pars[counter++] );  
			  theV[i][j][k] = interaction;              
                 }
            }  
            else {
                 for(int k=0; k< ExprPar::nbins; k++){
			  theV[i][j][k] = 1;    
		     }
            }      
        }
    }
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = i ; j < _nFactors; j++ ) {            
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {
            theV[i][j][k] =  theV[j][i][k] ;
	     
            }       
        }
    }    
theVr = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));

for ( int i = 0; i < _nFactors; i++ ) {
      
 	for ( int j = 0; j < i; j++ ) { 
            if ( repIndicators[ i] ) {
		for(int k=0; k< ExprPar::nbins; k++){  
			 double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interactionr ), log( max_interactionr ) ) ) : exp( pars[counter++] );  
			  theVr[i][j][k] = interaction;              
                 }
            }  
            else {
                 for(int k=0; k< ExprPar::nbins; k++){
			  theVr[i][j][k] = 1;    
		     }
            }      
        }
    }
    for ( int i = 0; i < _nFactors; i++ ) {
        for ( int j = i + 1; j < _nFactors; j++ ) {                      
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {		
            theVr[i][j][k] =  theVr[j][i][k] ;
	     
            }       
        }
    }       
     
    for ( int i = 0; i < _nFactors; i++ ) {
        
            double effect = searchOption == CONSTRAINED ?  inverse_infty_transform( pars[counter++],  min_effect_Thermo ,  max_effect_Thermo   ) :  pars[counter++] ;
            txpEffects.push_back( effect ); 
        
    }
 
   
}

}

void ExprPar::getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );  
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();

   		
    
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : log( maxBindingWts[i] );
            pars.push_back( weight );
        }
    }


	for ( int i = 0; i < nFactors(); i++ ) {
       	 for ( int j = 0; j <= i; j++ ) {
            if (coopMat(i,j) ) {
		for(int k=0; k< ExprPar::nbins ; k++){
			
               double interaction = searchOption == CONSTRAINED ? infty_transform( log( theV[i][j][k] ), log( min_interaction ), log( max_interaction ) ) : log( theV[i][j][k] ); 
	       pars.push_back( interaction);
		             
                 }

            }  
          }
        }
        for ( int i = 0; i < nFactors(); i++ ) {
          for ( int j = 0; j < i; j++ ) {   
            if ( repIndicators[i] ) {
		for(int k=0; k< ExprPar::nbins ; k++){
         		 double interaction = searchOption == CONSTRAINED ? infty_transform( log( theVr[i][j][k] ), log( min_interactionr ), log( max_interactionr ) ) : log( theVr[i][j][k] );
	  	
			  pars.push_back( interaction);
		     
                 }

             }  
           }
         }    
	 for ( int i = 0; i < nFactors(); i++ ) {
 double effect = searchOption == CONSTRAINED ? infty_transform(  txpEffects[i] , min_effect_Thermo ,  max_effect_Thermo  ) :  txpEffects[i] ;
            pars.push_back( effect );
	  }

}

void ExprPar::getFreePars2( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) 
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );  
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();

   	
    
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : log( maxBindingWts[i] );
            pars.push_back( weight );
        }
    }


	for ( int i = 0; i < nFactors(); i++ ) {
       	 for ( int j = 0; j <= i; j++ ) {
            if (coopMat(i,j) ) {
		for(int k=0; k< ExprPar::nbins ; k++){
			
               double interaction = searchOption == CONSTRAINED ? infty_transform( log( theV[i][j][k] ), log( min_interaction ), log( max_interaction ) ) : log( theV[i][j][k] ); 
	       pars.push_back( interaction);
		             
                 }

            }  
          }
        }
        for ( int i = 0; i < nFactors(); i++ ) {
          for ( int j = 0; j < i; j++ ) {   
            if ( repIndicators[i] ) {
		for(int k=0; k< ExprPar::nbins ; k++){
         		 double interaction = searchOption == CONSTRAINED ? infty_transform( log( theVr[i][j][k] ), log( min_interactionr ), log( max_interactionr ) ) : log( theVr[i][j][k] );
	  	
			  pars.push_back( interaction);
		     
                 }

             }  
           }
         }    
	 for ( int i = 0; i < nFactors(); i++ ) {
 double effect = searchOption == CONSTRAINED ? infty_transform(  txpEffects[i] , min_effect_Thermo ,  max_effect_Thermo  ) :  txpEffects[i] ;
            pars.push_back( effect );
	  }

}
void ExprPar::print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const 
{
    
    
    for ( int i = 0; i < nFactors(); i++ ) {
        os << "motif"<<motifNames[i] << "\t" << "mBwt" << maxBindingWts[i] << "\t" <<"txpE"<< txpEffects[i];
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) os << "\t" << repEffects[i];
        os << endl;
    }
if(modelOption== BINS){
     
for(int i =0; i < theV.size(); i++){
	for(int j =0; j < theV[i].size(); j++){
		for(int k =0; k < theV[i][j].size(); k++)  
		{  if ( coopMat(i,j)){
			os << "theV" << "\t" << theV[i][j][k] ;
			}	
			
  		}
	}
}
}

}

int ExprPar::load( const string& file, IntMatrix& coopMat, const vector< bool >& repIndicators )
{
    
    ifstream fin( file.c_str() );
    if ( !fin ){ cerr << "Cannot open parameter file " << file << endl;	exit( 1 ); } 
    
    vector< string > motifNames( nFactors() );
    for ( int i = 0; i < nFactors(); i++ ) {
        fin >> motifNames[i] >> maxBindingWts[i] >> txpEffects[i]; 
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) fin >> repEffects[i];
    }

    
    map< string, int > factorIdxMap;
    for ( int i = 0; i < nFactors(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
    
theV = vector< vector< vector<double> > >(nFactors(),vector< vector<double> >(nFactors(), vector< double>( ExprPar::nbins,1)));
theVr = vector< vector< vector<double> > >(nFactors(),vector< vector<double> >(nFactors(), vector< double>( ExprPar::nbins,1)));
    string value, sym;
	fin >> sym;
	fin >> sym;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            	if ( coopMat( i, j ) ) {
			for(int k=0; k< ExprPar::nbins; k++){
			   fin >> value;  
			  double interaction = atof( value.c_str() ) ; 
			  theV[i][j][k] = interaction;              
                	 }
            	}  
                else {
                     for(int k=0; k< ExprPar::nbins; k++){
			  theV[i][j][k] = 1;    
		     }
                 }      
         }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = i + 1; j < nFactors(); j++ ) {
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {
            theV[i][j][k] =  theV[j][i][k] ;
            }       
        }
    }    

for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( repIndicators[ i] ) {
		fin >> sym;
	fin >> sym;
		for(int k=0; k< ExprPar::nbins; k++){  
			fin >> value;  
			  double interaction = atof( value.c_str() ) ; 
			  theVr[i][j][k] = interaction;              
                 }
            }  
            else {
                 for(int k=0; k< ExprPar::nbins; k++){
			  theVr[i][j][k] = 1;    
		     }
            }      
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = i + 1; j < nFactors(); j++ ) {                      
	    for ( int k = 0; k <  ExprPar::nbins; k++ ) {		
            theVr[i][j][k] =  theVr[j][i][k] ;
            }       
        }
    }       
   

    
  return 0;
}

void ExprPar::adjust()
{
    
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) ) maxBindingWts[i] *= 2.0;
        if ( maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) ) maxBindingWts[i] /= 2.0;
    }
    
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
           
	   
             for(int k =0; k < ExprPar::nbins; k++){
                    if ( theV[i][j][k] < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) ) { 
                        theV[i][j][k] *= 10.0; 
                       
                         
	                   
                                           
                                           theV[i][j][k] = theV[j][i][k] ;
                         
                     }
                                         
	                                  
		                         
                     if ( theV[i][j][k] > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) { 
                                          
                                          
                          theV[i][j][k] /= 2.0; 
                         
                           
	                     
                                            theV[i][j][k] = theV[j][i][k] ; 
                         
                     }
               }
            
        }
    }

for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
           
	   
             for(int k =0; k < ExprPar::nbins; k++){
                    if ( theVr[i][j][k] < ExprPar::min_interactionr * ( 1.0 + ExprPar::delta ) ) { 
                        theVr[i][j][k] *= 10.0; 
                       
                         
	                   
                                           
                                           theVr[i][j][k] = theVr[j][i][k] ;
                         
                     }
                                         
	                                  
		                         
                     if ( theVr[i][j][k] > ExprPar::max_interactionr * ( 1.0 - ExprPar::delta ) ) { 
                                          
                                          
                          theVr[i][j][k] /= 2.0; 
                         
                           
	                     
                                            theVr[i][j][k] = theVr[j][i][k] ; 
                         
                     }
               }
            
        }
    }

    for ( int i = 0; i < nFactors(); i++ ) {
        
            if ( txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) ) txpEffects[i] *= 2.0;
            if ( txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) ) txpEffects[i] /= 2.0;
        
        
    }

}

ModelType ExprPar::modelOption = CHRMOD_UNLIMITED;  
SearchType ExprPar::searchOption = UNCONSTRAINED;   
int ExprPar::estBindingOption = 1;                  
vector< int > ExprFunc::Bc(5,1);  
vector< int > ExprFunc::Bct(5,1);  
vector< int > ExprFunc::B12(5,1);
vector< int > ExprFunc::B03(5,1);
vector< int > ExprFunc::B13(5,1);
vector< int > ExprFunc::B23(5,1);
vector< Sequence > ExprPredictor::seqsy;
vector< string > ExprPredictor::seqNmes;
vector< Sequence > ExprPredictor::seqsya;
vector< string > ExprPredictor::seqNmesa;
Matrix ExprPredictor::factorExprData2;
Matrix ExprPredictor::exprData2;
double ExprPar::default_weight = 10;
double ExprPar::default_interaction = .95;
double ExprPar::default_effect_Logistic = 5;
double ExprPar::default_effect_Thermo = 5;
double ExprPar::default_repression = 1.0E-2;
double ExprPar::default_basal_Logistic = -9.0;
double ExprPar::default_basal_Thermo = 0.001;
double ExprPar::min_weight = 1;		
double ExprPar::max_weight = 100;		
double ExprPar::min_interaction = .95;	
double ExprPar::max_interaction = 100;
double ExprPar::min_interactionr = .001;	
double ExprPar::max_interactionr = 1.1;
double ExprPar::min_effect_Logistic = 5;	
double ExprPar::max_effect_Logistic = 5;
double ExprPar::min_effect_Thermo = 1; 
double ExprPar::max_effect_Thermo = 10;
double ExprPar::min_repression = 1.0E-3;
double ExprPar::max_repression = 500; 
double ExprPar::min_basal_Logistic = -9.0;	
double ExprPar::max_basal_Logistic = -9.0;
double ExprPar::min_basal_Thermo =.001;	
double ExprPar::max_basal_Thermo = 0.001;
double ExprPar::delta = 0.001;  
int ExprPar::nbins = 5;  
 

ExprFunc::ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par,  const vector < int >& _B, const vector < int >& _Br  ) : motifs( _motifs ), intFunc( _intFunc ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), par( _par ),   B( _B ), Br( _Br )
{ 
    int nFactors = par.nFactors();
    
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );
    Bc = vector< int >(5,0);
    Bct = vector< int >(5,0);
    B12 = vector< int >(5,0);
    B03 = vector< int >(5,0);
    B13= vector< int >(5,0);
    B23= vector< int >(5,0);
}

void ExprFunc::predictTwist( vector< vector< vector< Site > > >& _d )
{
	ofstream dorsla("dorsla.txt");
	
		
     		
for ( int e = 0; e < _d.size(); e++ ) {  
	
		for(int j =0; j < _d[e][0].size(); j++ ) {
		


		for(int t=0; t< _d[e][1].size();t++ ) {
Site a =_d[e][0][j];
Site b= _d[e][1][t];
double dist =0 ;
if( !siteOverlap( a,b,motifs) ) {

if( a.start > b.start ) {  dist = abs( a.start - (b.start + motifs[1].length() ) );  }  
if( b.start > a.start ) { dist =  abs( b.start - (a.start + motifs[0].length() ) );  }
}

double maxInt=1;

if (      dist <B[0] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 0 ]; ExprFunc::Bc[0]++; } 
	if (      dist >=B[0] && dist <= B[1] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 1 ]; ExprFunc::Bc[1]++;}
	if (      dist > B[1] && dist <= B[2] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 2 ]; ExprFunc::Bc[2]++;}

	if (      dist > B[2] && dist <= B[3] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 3 ]; ExprFunc::Bc[3]++;}
	if (      dist > B[3] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 4 ]; ExprFunc::Bc[4]++;}
	
		}
	}
}

ofstream dc("dctwsites.txt");
ofstream dust2("stranddctw.txt");
ofstream dcsp("dctwsp.txt");
int countsdc =0;
for ( int e = 0; e < _d.size(); e++ ) {  
	
		for(int j =0; j < _d[e][0].size(); j++ ) {
			bool count =0;
			countsdc++;
			vector< int > dt;
			if( _d[e][2].size() ==0){ dorsla <<  999 << '\t' ;  ExprFunc::Bct[4]++;}
			
			for(int t=0; t< _d[e][2].size();t++ ) {
				
				Site a =_d[e][0][j];
				Site b= _d[e][2][t];
				double dist =0 ;

				dust2 << a.strand << '\t';
				if( !siteOverlap( a,b,motifs) ) {
				if( a.start > b.start ) {  dist = abs( a.start - (b.start + motifs[2].length() ) );  }
				if( b.start > a.start ) { dist =  abs( b.start - (a.start + motifs[0].length() ) );  }
				}
				
				dt.push_back( dist );
				double maxInt=1;
				
				dorsla << dist << '\t' ; 
					
			} 
	
	int dist_opt = 50000;
	for(int t=0; t< dt.size();t++ ) {
		int dist_temp =dt[t];
		if( dist_temp < dist_opt ) {
		dist_opt= dist_temp;
		} 		
	}
	
	
	
	if (      dist_opt < B[0] )                 {  ExprFunc::Bct[0]++; } 
	if (      dist_opt >=B[0] && dist_opt <= B[1] ) {  ExprFunc::Bct[1]++; }
	if (      dist_opt > B[1] && dist_opt <= B[2] ) {  ExprFunc::Bct[2]++;}
	if (      dist_opt > B[2] && dist_opt <= B[3] ) { ExprFunc::Bct[3]++; }
	if (      dist_opt > B[3] )                 {  ExprFunc::Bct[4]++; }
	
	if( _d[e][2].size() > 0 ) { 

		dc << dist_opt << '\t'; 



		
				
		Site aa =_d[e][0][j];
				
				  
							int jj = 0;
				 	if(  aa.start+motifs[ aa.factorIdx ].length()+5 > ExprPredictor::seqsya[e].size() ){continue;}
					 if(  aa.start-5  <= 0 ){continue;}
					dcsp << ">" << ExprPredictor::seqNmesa[e] << "_" << '\t' ;
					int shift = motifs[ aa.factorIdx ].length()+9;
					for( int p =0 ; p < shift; p++ ){
						jj = aa.start-4 ;
						if(aa.strand){
							dcsp <<  ALPHABET[ ExprPredictor::seqsya[e][ jj+p ]  ] ;
						}
						else{ dcsp << ALPHABET[complement(symbolToInt(ALPHABET[ ExprPredictor::seqsya[e][jj+ shift -1 - p ] ]))];} 
					} 
					dcsp << '\t';
		

		dcsp << _d[e][0][j] << '\t' << dist_opt << endl;

	}
     }  
} 

dust2.close();
dorsla.close();
dc << endl;
dc.close();
dcsp.close();
ofstream dut("dutwsites.txt");
ofstream dct("dutwdist.txt");
ofstream dust("stranddutw.txt");
int countsdu =0;
for ( int e = 0; e < _d.size(); e++ ) {  
	

		for(int j =0; j < _d[e][1].size(); j++ ) {
		bool count =0;
		countsdu++;
vector< int > dt;
if( _d[e][2].size() ==0){  ExprFunc::B12[4]++;}
		for(int t=0; t< _d[e][2].size();t++ ) {
Site a =_d[e][1][j];
Site b= _d[e][2][t];
double dist =0 ;
dust << a.strand << '\t';
if( !siteOverlap( a,b,motifs) ) {
if( a.start > b.start ) {  dist = abs( a.start - (b.start + motifs[2].length() ) );  }
if( b.start > a.start ) { dist =  abs( b.start - (a.start + motifs[1].length() ) );  }
}
dt.push_back( dist );

		} 
			int dist_opt = 50000;
	for(int t=0; t< dt.size();t++ ) {
		int dist_temp =dt[t];
		if( dist_temp < dist_opt ) {
		dist_opt= dist_temp;
		} 		
	}  
	
	double maxInt=1;
	
	dct << dist_opt << '\t' ;
	if (      dist_opt < B[0] )                 {  ExprFunc::B12[0]++;} 
		if (      dist_opt >=B[0] && dist_opt <= B[1] ) {  ExprFunc::B12[1]++; }
		if (      dist_opt > B[1] && dist_opt <= B[2] ) { ExprFunc::B12[2]++; }
		if (      dist_opt > B[2] && dist_opt <= B[3] ) {  ExprFunc::B12[3]++; }
		if (      dist_opt > B[3] )                 {  ExprFunc::B12[4]++; }
		

		if( _d[e][2].size() > 0 ) { 

			dut << dist_opt << '\t'; 
		} 
	}  
}  

dut << endl;
dut.close();
dct << endl;
dct.close();
dust.close();
for ( int e = 0; e < _d.size(); e++ ) {  
	
		for(int j =0; j < _d[e][0].size(); j++ ) {

		for(int t=0; t< _d[e][3].size();t++ ) {
Site a =_d[e][0][j];
Site b= _d[e][3][t];
double dist =0 ;
if( !siteOverlap( a,b,motifs) ) {
if( a.start > b.start ) {  dist = abs( a.start - (b.start + motifs[3].length() ) );  }
if( b.start > a.start ) { dist =  abs( b.start - (a.start + motifs[0].length() ) );  }
}
double maxInt=1;
if (      dist <B[0] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 0 ]; ExprFunc::B03[0]++; } 
	if (      dist >=B[0] && dist <= B[1] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 1 ]; ExprFunc::B03[1]++;}
	if (      dist > B[1] && dist <= B[2] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 2 ]; ExprFunc::B03[2]++;}
	if (      dist > B[2] && dist <= B[3] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 3 ]; ExprFunc::B03[3]++;}
	if (      dist > B[3] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 4 ]; ExprFunc::B03[4]++;}
	
		}
	}
}

for ( int e = 0; e < _d.size(); e++ ) {  
	
		for(int j =0; j < _d[e][1].size(); j++ ) {

		for(int t=0; t< _d[e][3].size();t++ ) {
Site a =_d[e][1][j];
Site b= _d[e][3][t];
double dist =0 ;
if( !siteOverlap( a,b,motifs) ) {
if( a.start > b.start ) {  dist = abs( a.start - (b.start + motifs[3].length() ) );  }
if( b.start > a.start ) { dist =  abs( b.start - (a.start + motifs[1].length() ) );  }
}
double maxInt=1;
if (      dist <B[0] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 0 ]; ExprFunc::B13[0]++; } 
	if (      dist >=B[0] && dist <= B[1] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 1 ]; ExprFunc::B13[1]++;}
	if (      dist > B[1] && dist <= B[2] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 2 ]; ExprFunc::B13[2]++;}
	if (      dist > B[2] && dist <= B[3] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 3 ]; ExprFunc::B13[3]++;}
	if (      dist > B[3] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 4 ]; ExprFunc::B13[4]++;}
	
		}
	}
}
ofstream cbt("cbtwdist.txt");
for ( int e = 0; e < _d.size(); e++ ) {  
	
if( _d[e][2].size() ==0){ExprFunc::B23[4]++;}
		for(int j =0; j < _d[e][2].size(); j++ ) {

		for(int t=0; t< _d[e][3].size();t++ ) {
Site a =_d[e][2][j];
Site b= _d[e][3][t];
double dist =0 ;
if( !siteOverlap( a,b,motifs) ) {
if( a.start > b.start ) {  dist = abs( a.start - (b.start + motifs[3].length() ) );  }
if( b.start > a.start ) { dist =  abs( b.start - (a.start + motifs[2].length() ) );  }
}
double maxInt=1;
cbt << dist << '\t' ;
if (      dist <B[0] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 0 ]; ExprFunc::B23[0]++; } 
	if (      dist >=B[0] && dist <= B[1] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 1 ]; ExprFunc::B23[1]++;}
	if (      dist > B[1] && dist <= B[2] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 2 ]; ExprFunc::B23[2]++;}
	if (      dist > B[2] && dist <= B[3] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 3 ]; ExprFunc::B23[3]++;}
	if (      dist > B[3] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 4 ]; ExprFunc::B23[4]++;}
	
		}
	}
}
cbt << endl;
cbt.close();
}



void ExprFunc::predictTwist2( vector< vector< vector< Site > > >& _d )
{

ofstream du( "Data/DorsalUnCondit2.txt" );
ofstream dc( "Data/DorsalCondit2.txt" );
for ( int e = 0; e < _d.size(); e++ ) {  

		for(int j =0; j < _d[e][0].size(); j++ ) {
			int flagd = 0;
			int flagt = 0;
			if(_d[e][2].size() ==0 ){
			vector< Site > tsites;
			tsites = _d[e][0];
			Site a =_d[e][0][j];
			for( int i = 0; i < tsites.size() ; i++ ) {	
				if( tsites[i].start == a.start ) { 
				int jj = 0;
				if(  a.start+motifs[ tsites[i].factorIdx ].length()+2 > ExprPredictor::seqsya[e].size() ){continue;}
	 			if(  tsites[i].start-2  <= 0 ){continue;}
				du << ">" << ExprPredictor::seqNmesa[e] << "_" << tsites[i].start <<endl ;
				for( int p =0 ; p < motifs[ tsites[i].factorIdx ].length()+2; p++ ){
					jj = tsites[i].start-1 ;
					if(tsites[i].strand){
						du <<  ALPHABET[ ExprPredictor::seqsya[e][ jj+p ]  ] ;
					}
					else{ du << ALPHABET[complement(symbolToInt(ALPHABET[ ExprPredictor::seqsya[e][jj+ motifs[ tsites[i].factorIdx ].length() -1 - p ] ]))];} 
				} 
				du << endl;
				break;  
				} 
			}  
		        } 
		for(int t=0; t< _d[e][2].size();t++ ) {

			Site a =_d[e][0][j];
			Site b= _d[e][2][t];
			double dist =0 ;
		
			if( !siteOverlap( a,b,motifs) ) {
			if( a.start > b.start ) {  dist = abs( a.start - (b.start + motifs[2].length() ) );  }
			if( b.start > a.start ) { dist =  abs( b.start - (a.start + motifs[0].length() ) );  }
			}
			double maxInt=1;

			if( dist <B[2] ){ 
				
				vector< Site > tsites;
				tsites = _d[e][0];
		
				for( int i = 0; i < tsites.size() ; i++ ) {	
				   if( tsites[i].start == a.start ) { 
							int jj = 0;
				 	if(  a.start+motifs[ tsites[i].factorIdx ].length()+2 > ExprPredictor::seqsya[e].size() ){continue;}
					 if(  tsites[i].start-2  <= 0 ){continue;}
					dc << ">" << ExprPredictor::seqNmesa[e] << "_" << tsites[i].start <<endl ;
					for( int p =0 ; p < motifs[ tsites[i].factorIdx ].length()+2; p++ ){
						jj = tsites[i].start-1 ;
						if(tsites[i].strand){
							dc <<  ALPHABET[ ExprPredictor::seqsya[e][ jj+p ]  ] ;
						}
						else{ dc << ALPHABET[complement(symbolToInt(ALPHABET[ ExprPredictor::seqsya[e][jj+ motifs[ tsites[i].factorIdx ].length() -1 - p ] ]))];} 
					} 
					dc << endl;
					flagd=1;
					flagt=1;
					break;
				   } 
				}  
			
			 } 
			if( flagd==1){ break;}
			
		}
	}
}

for ( int e = 0; e < _d.size(); e++ ) {  
	


		for(int j =0; j < _d[e][1].size(); j++ ) {
			int flag=0;
			int flagd = 0; 
			
			
				for(int ii =0; ii < _d[e][0].size(); ii++ ) {  
				

					Site a =_d[e][1][j];
					Site b= _d[e][0][ii];
					double dist =0 ;
					
					if( !siteOverlap( a,b,motifs) ) {
					if( a.start > b.start ) {  dist = abs( a.start - (b.start + motifs[0].length() ) );  }
					if( b.start > a.start ) { dist =  abs( b.start - (a.start + motifs[1].length() ) );  }
					}

					if ( dist <=B[1] ){ flag= 1; break;}



				} 
			
			if( flag==1 ){ continue;}
					if(_d[e][2].size() ==0 ){
					vector< Site > tsites;
							tsites = _d[e][1];
						Site a =_d[e][1][j];
							for( int i = 0; i < tsites.size() ; i++ ) {	
								if( tsites[i].start == a.start ) { 
								int jj = 0;
					if( a.start+ motifs[ tsites[i].factorIdx ].length()+2 > ExprPredictor::seqsya[e].size() ){continue;}
					 if(  tsites[i].start- 2 <= 0 ){continue;}
								du << ">" << ExprPredictor::seqNmesa[e] << "_" << tsites[i].start <<endl ;
								for( int p =0 ; p < motifs[ tsites[i].factorIdx ].length()+2; p++ ){
									jj = tsites[i].start-1 ;
									if(tsites[i].strand){
										du <<  ALPHABET[ ExprPredictor::seqsya[e][ jj+p ]  ] ;
									}
									else{ du << ALPHABET[complement(symbolToInt(ALPHABET[ ExprPredictor::seqsya[e][jj+ motifs[ tsites[i].factorIdx ].length() -1 - p ] ]))];} 
								} 
								du << endl;
								} 
							}  
						} 
		for(int t=0; t< _d[e][2].size();t++ ) {
				
				Site a =_d[e][1][j];
				Site b= _d[e][2][t];
				double dist =0 ;
				
				if( !siteOverlap( a,b,motifs) ) {
				if( a.start > b.start ) {  dist = abs( a.start - (b.start + motifs[2].length() ) );  }
				if( b.start > a.start ) { dist =  abs( b.start - (a.start + motifs[1].length() ) );  }
				}
				double maxInt=1;
				
				
				
			if (      dist <B[2] ){
				vector< Site > tsites;
					tsites = _d[e][1];
					for( int i = 0; i < tsites.size() ; i++ ) {	
						if( tsites[i].start == a.start ) { 
						int jj = 0;
					if( a.start+ motifs[ tsites[i].factorIdx ].length()+2 > ExprPredictor::seqsya[e].size() ){continue;}
					 if(  tsites[i].start-2  <= 0 ){continue;}
						dc << ">" << ExprPredictor::seqNmesa[e] << "_" << tsites[i].start <<endl ;
						for( int p =0 ; p < motifs[ tsites[i].factorIdx ].length()+2; p++ ){
							jj = tsites[i].start-1 ;
							if(tsites[i].strand){
								dc <<  ALPHABET[ ExprPredictor::seqsya[e][ jj+p ]  ] ;
							}
							else{ dc << ALPHABET[complement(symbolToInt(ALPHABET[ ExprPredictor::seqsya[e][jj+ motifs[ tsites[i].factorIdx ].length() -1 - p ] ]))];} 
						} 
						dc << endl;
						flagd=1;
						break;
						} 
					}  



			  } 
			if( flagd==1){ break;}
			else{ 
				vector< Site > tsites;
					tsites = _d[e][1];
					for( int i = 0; i < tsites.size() ; i++ ) {	
						if( tsites[i].start == a.start ) { 
								int jj = 0;
					if(  a.start+ motifs[ tsites[i].factorIdx ].length()+2 > ExprPredictor::seqsya[e].size() ){continue;}
					 if(  tsites[i].start-2 <= 0 ){continue;}
						du << ">" << ExprPredictor::seqNmesa[e] << "_" << tsites[i].start <<endl ;
						for( int p =0 ; p < motifs[ tsites[i].factorIdx ].length()+2; p++ ){
							jj = tsites[i].start-1 ;
							if(tsites[i].strand){
								du <<  ALPHABET[ ExprPredictor::seqsya[e][ jj+p ]  ] ;
							}
							else{ du << ALPHABET[complement(symbolToInt(ALPHABET[ ExprPredictor::seqsya[e][jj+ motifs[ tsites[i].factorIdx ].length() -1 - p ] ]))];} 
						} 
						du << endl;
						flagd=1;
						break;
						} 
					}  
			}
			if( flagd==1){ break;}  
		}
	}
}



du.close();
dc << endl;
dc.close();
}


Matrix countmatrix( const string& file, int & nc )
{

    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
       
    }

	string temp;
	vector <string> seq_name;
	vector <string> seq;
	
	string templ;
int tem = 0;
int tempmax = 0;
	seq.clear();
	while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
	
		string name (temp, 1, temp.length() - 1);
		seq_name.push_back(name);
		
		getline(seq_file, temp);
		seq.push_back(temp);
		tem = temp.length();
			if( tem > tempmax) tempmax = tem;	
		
	}
	

int length = tempmax ;  



	 Matrix m(4,length,0);
		nc =seq.size();	

	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[j]);
		Sequence RCreadseq( readseq.compRevCompl() );  
		int na = 0;
		for(int i = 2;  i< readseq.size()-2; i++){
			if(readseq[i] == 0)            
				{ na = na + 1;}
		}
		int naRC = 0;
		for(int i = 2;  i< RCreadseq.size()-2; i++){
			if(RCreadseq[i] == 0)            
				{ naRC = naRC + 1;}
		}
		if( na >= naRC|| length < 10){
			
			for( int i = 0;  i< readseq.size(); i++)  
			{
				if(readseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		}  
		else {
			
			for( int i = 0;  i< RCreadseq.size(); i++)  
			{
				if(RCreadseq[i] == 0)            
				{ m(0,i) = m(0,i) + 1;}
				if(RCreadseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(RCreadseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(RCreadseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				
			}
		} 


	} 
	Matrix countMatrix = m.transpose();
	
	 double pseudoCount =1;
	
	    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	    
	    int l = countMatrix.nRows();		
	
	  Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
	   return pwm;
}



void ExprFunc::phaseShift( vector< vector< vector< Site > > >& _d )
{
}

void ExprFunc::predictOcc( const SiteVec& _sites, int length, const vector< double >& factorConcs , vector< double >& fOcc)  

{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();	
    
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[ sites[i].factorIdx ] * sites[i].wtRatio );	
    }
    


if ( modelOption == BINS ) {


	vector< double > Y( n + 1 );
    Y[0] = 0;
	vector< double > Yt( n + 1 );
    Yt[0] = 0;
    
	
	double Z_bind = compPartFunc();
	for ( int k = 0; k < motifs.size(); k++ ) {  
	     for ( int i = 1; i <= n; i++ ) {
		int factorIdx;
		factorIdx = k;
		int match = sites[ i ].factorIdx == factorIdx ? 1 : 0; 
		double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
	
		}
		Y[ i ] = bindingWts[ i ] * sum;    
		Yt[i] = Y[i] + Yt[i - 1];
              }
	      double Y_total = Yt[n];
              factorOcc.push_back( Y_total / Z_bind );
	  }
        double totalEffect = 0;
fOcc = factorOcc;

      for ( int i = 0; i < motifs.size(); i++ ) {
            double effect = factorOcc[i];  
	
            totalEffect += effect;
	}
            
          
        
        
        
   
}
}

double ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();	
    
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }	
    
    
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );	
    }

	vector< double > Y( n + 1 );
    Y[0] = 0;
	vector< double > Yt( n + 1 );
    Yt[0] = 0;
    
	
	double Z_bind = compPartFunc();
	for ( int k = 0; k < motifs.size(); k++ ) {  
	     for ( int i = 1; i <= n; i++ ) {
		int factorIdx;
		factorIdx = k;
		int match = sites[ i ].factorIdx == factorIdx ? 1 : 0; 
		
		double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
	
		}
		Y[ i ] = bindingWts[ i ] * sum;    
		Yt[i] = Y[i] + Yt[i - 1];
              }
	      double Y_total = Yt[n];
	
              factorOcc.push_back( Y_total / Z_bind );
	  }

        double totalEffect = 0;

     
           
	
            
	
	 totalEffect =par.txpEffects[0] * factorOcc[0]   -5;
            
          
        
	
	
	 
       return      1/(1+exp(-totalEffect) );                    
        
   
}

double ExprFunc::predictZ( const SiteVec& _sites, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();	

    
    int n = _sites.size();
    sites = _sites;

    sites.insert( sites.begin(), Site() );  

    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );

    for ( int i = 1; i <= n; i++ ) {
        int j; 
        for ( j = i - 1; j >= 1; j-- ) {

            if ( ( sites[i].start - sites[j].start ) > range ) break; 
        }

        int boundary = j;
        boundaries.push_back( boundary );
    }	
   
    
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
	

        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );


    }

	
    
	
	double Z_bind = compPartFunc();


       return Z_bind;
        
   
}

    
double ExprFunc::compPartFunc()
{
    int n = sites.size() - 1;
    
	
    Z = vector< double >( n + 1 );
    Z[0] = 1.0;
    Zt = vector< double >( n + 1 );
    Zt[0] = 1.0;
	
	
double sum=0;  
	for ( int i = 1; i <= n; i++ ) {
		sum = Zt[boundaries[i]];  
	
		for ( int j = boundaries[i] + 1; j < i; j++ ) {
			if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
			sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
			
			
		}
		Z[i] = bindingWts[ i ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
	}
	
	
    double Z_bind = Zt[n];
	return Z_bind;
}

ModelType ExprFunc::modelOption = QUENCHING;  





double ExprFunc::compFactorInt( const Site& a, const Site& b ) const
{
double dist = abs( a.start - b.start );
double maxInt=1;
if ( ( repIndicators[a.factorIdx ] || repIndicators[ b.factorIdx ] ) && ( repIndicators[a.factorIdx ] && repIndicators[ b.factorIdx ] )   )
maxInt = 1;

else if ( repIndicators[ a.factorIdx ] || repIndicators[ b.factorIdx ]) {
	if (      dist <= Br[0] )                 { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 0 ]; } 
	if (      dist > Br[0] && dist <= Br[1] ) { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 1 ]; }
	if (      dist > Br[1] && dist <= Br[2] ) { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 2 ]; }
	if (      dist > Br[2] && dist <= Br[3] ) { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 3 ]; }
	if (      dist > Br[3] )                 { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 4 ]; }

}
else {
	if (      dist <= B[0] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 0 ]; ExprFunc::Bc[0]++; } 
	if (      dist > B[0] && dist <= B[1] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 1 ]; ExprFunc::Bc[1]++;}
	if (      dist > B[1] && dist <= B[2] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 2 ]; ExprFunc::Bc[2]++;}
	if (      dist > B[2] && dist <= B[3] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 3 ]; ExprFunc::Bc[3]++;}
	if (      dist > B[3] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 4 ]; ExprFunc::Bc[4]++;}

}

     bool orientation = ( a.strand == b.strand ); 
    return  intFunc->compFactorInt( maxInt, dist, orientation );	
}



bool ExprFunc::testRepression( const Site& a, const Site& b ) const
{

    double dist = abs( a.start - b.start );
    return repressionMat( a.factorIdx, b.factorIdx ) && ( dist <= repressionDistThr );
}




ExprPredictor::ExprPredictor(const vector< vector< SiteVec > >& _seqSitesb, const vector< vector< double > >& _bindingData ,  vector< SiteVec >& _seqSitesa,  vector< SiteVec >& _seqSites, const vector< int >& _seqLengths, Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const FactorIntFunc* _intFunc, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr,const vector< int >& _binBord ,const vector< int >& _binBordr, const vector < bool >& _indicator_bool, SeqAnnotator _anny, const string& _file,  vector< SiteVec >& _seqSitesbot,  vector< SiteVec >& _seqSitesm1,  vector< SiteVec >& _seqSitesm2 , vector< SiteVec >& _seqSitesf2 ,vector< SiteVec >& _seqSitesbotf2, vector< SiteVec >& _seqSitesm1f2, vector< SiteVec >& _seqSitesm2f2, vector< SiteVec >& _seqSitesf3 ,  vector< SiteVec >& _seqSitesbotf3,vector< SiteVec >& _seqSitesm1f3 , vector< SiteVec >& _seqSitesm2f3, vector< vector< vector< Site > > >& _d ): seqSitesb( _seqSitesb ), bindingData( _bindingData), seqSitesa(_seqSitesa), seqSites( _seqSites ), seqLengths( _seqLengths ), exprData( _exprData ), motifs( _motifs ), factorExprData( _factorExprData ), intFunc( _intFunc ), coopMat( _coopMat ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), binBord(_binBord), binBordr(_binBordr),  indicator_bool ( _indicator_bool ), anny( _anny) , file( _file ) , seqSitesbot( _seqSitesbot ), seqSitesm1( _seqSitesm1 ), seqSitesm2( _seqSitesm2 ),seqSitesf2( _seqSitesf2 ) ,seqSitesbotf2(_seqSitesbotf2), seqSitesm1f2( _seqSitesm1f2) ,seqSitesm2f2( _seqSitesm2f2), seqSitesf3( _seqSitesf3) , seqSitesbotf3( _seqSitesbotf3),seqSitesm1f3( _seqSitesm1f3 ), seqSitesm2f3( _seqSitesm2f3),  d(_d ), seqSitesm1d1()
{   


    
    ExprPar::modelOption = modelOption;
    ExprFunc::modelOption = modelOption;
    ExprPar::nbins = _binBord.size() +1;  
 
    
    if ( modelOption != LOGISTIC && modelOption != DIRECT && modelOption !=BINS ) {
        ExprPar::min_effect_Thermo = 0.99;
        ExprPar::min_interaction = 0.99;
    }

    
    ExprPar::estBindingOption = estBindingOption;
}

double ExprPredictor::objFunc( const ExprPar& par ) const
{
    if ( objOption == SSE ) return compRMSE( par );	
    if ( objOption == CORR ) return -compAvgCorr( par );
    if ( objOption == CROSS_CORR ) return -compAvgCrossCorr( par ); 
}

double ExprPredictor::objFunc2(  ExprPar& par ) 
{
    return   1 ; 
}

double ExprPredictor::objFuncborder(  ExprPar& par ) 
{
    return compAvgCorrborder( par );
}
double ExprPredictor::objFuncborder2(  ExprPar& par ) 
{
    return compAvgCorrborder2( par );  
}
int ExprPredictor::train( const ExprPar& par_init )
{
   par_model = par_init;
 
 
       if ( nAlternations == 0 ) {obj_model = objFuncborder2( par_model );return 0;}
    
      
    
 
    ExprPar par_result;
    double obj_result;
   for ( int i = 0; i < nAlternations; i++ ) {
        simplex_minimize( par_result, obj_result );
if (!testPar(par_result) ) { par_result.adjust() ;}
        par_model = par_result; 
        
        gradient_minimize( par_result, obj_result );
if (!testPar(par_result) ) { par_result.adjust() ;}
        par_model = par_result;
	
      
    }
	
   
 par_model = par_result; 
    obj_model = obj_result;

printPar( par_model );
    return 0;	
}


int ExprPredictor::train( const ExprPar& par_init, const gsl_rng* rng )
{
   
 train( par_init );
    
	ExprPar par_best = par_model;
	double obj_best = obj_model;
	for ( int i = 0; i < nRandStarts; i++ ) {
        	ExprPar par_curr = par_init; 

	
		
		randSamplePar( rng, par_curr ); 
		
		
		train( par_curr );
       		
		 printPar( par_model );
       
		if ( obj_model < obj_best ) {
			par_best = par_model;
			obj_best = obj_model;	
		}
	}    

    
    if ( nRandStarts ) train( par_best ); 
 
par_model = par_best;
obj_model = obj_best;
printPar( par_model );
    return 0;
}
int ExprPredictor::train3( const ExprPar& par_init, const gsl_rng* rng )
{
   par_model = par_init;

    if ( ExprPar::searchOption == CONSTRAINED ) par_model.adjust();
    obj_model = objFunc( par_model );
       if ( nAlternations == 0 ) return 0;


    return 0;
}
void ExprPredictor::printFile( ostream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    
    
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
 
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( modelOption == LOGISTIC || modelOption == DIRECT ) os << par.txpEffects[i] << "\t";
        else {
            if ( actIndicators[i] ) os << par.txpEffects[i] << "\t";
        }
    }
 
    os << par.basalTxp << "\t"; 
    
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) os << par.factorIntMat( i, j ) << "\t";
        }
    }

   
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( repIndicators[i] ) os << par.repEffects[i] << "\t"; 
        }
    }
    os << jRMSE.getObj() << endl;
   
}
		

void ExprPredictor::printFilePar_KfoldCV( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
    
    
os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
os << endl<<" txpEffects : "<< endl;
 
    for ( int i = 0; i < nFactors(); i++ ) {
       os << par.txpEffects[i] << "\t";
       
    }
os << endl<<  " coopertivity : " << endl;

 for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if ( coopMat(i,j)){
			for(int k =0; k < par.theV[i][j].size(); k++) { 
		
			os <<  "\t" << par.theV[i][j][k] ;
			}
		}	
			
  		
	}
}
os<< endl << " quenching : " << endl;


for(int i =0; i < nFactors(); i++){
	for(int j =0; j <i; j++){
		if(repIndicators[i])  { 
			for(int k =0; k < ExprPar::nbins ; k++){ 
			os <<  "\t" << par.theVr[i][j][k] ;
				
  			}
		}
	}
}

os << endl << "objective function : " << endl;
    os << jRMSE.getObj() << endl;
   
ofstream fo( "Data/format.tex",ios::app );
int j=0;

vector< Site > tsites(0);
for(int m = 0; m < seqSites.size(); m++ ) {
tsites = seqSites[m];
  
	ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
	foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { foo << "b"; j++;  }
			
			else {
				

				if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  
				if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  
				if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	foo << endl;
	foo.close();



 j=0;

	for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
		if( i == 0 ) { fo << "&&" ; continue; }
		fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
		if(i% 90 == 0 ) {fo << "\\\\&&"; }
	}
fo << "\\\\&&\\\\";

fo << ExprPredictor::seqNmes[m] << "&" << "cell"<< cell[m] << "&" ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		if (tsites.size() == 1 ) { break; }
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}







fo << "\\\\&&\\\\&&\\\\" << endl;
}
fo.close();


}


int ExprPredictor::load( const string& fila )
{
   ofstream fout( fila.c_str() );
    vector< string > motifNames;

string dorsal, twist, snail;
dorsal="dl";
twist="tw";
snail="sn";
        motifNames.push_back( dorsal );
	motifNames.push_back( twist );
	motifNames.push_back( snail );
         fout << setprecision(2) ;
    for ( int i = 0; i < nFactors(); i++ ) {
	
        fout << motifNames[i] << '\t'<< par_model.maxBindingWts[i]<< '\t' << par_model.txpEffects[i] << '\n';
   
    }

    map< string, int > factorIdxMap;
    for ( int i = 0; i < nFactors(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }

   
 for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            	if ( coopMat( i, j ) ) {
			fout << motifNames[j]<< '\t' << motifNames[i];
			for(int k=0; k< ExprPar::nbins; k++){
			   fout << '\t' << par_model.theV[i][j][k];               
                	 }
			fout << '\n' ;
            	}  
         }
    }

for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j < i; j++ ) {
            if ( repIndicators[ i] ) {
		fout << motifNames[j]<< '\t' << motifNames[i];
		for(int k=0; k< ExprPar::nbins; k++){  
			fout<< '\t'<< par_model.theVr[i][j][k]  ;             
                 }
		fout << '\n' ;
            }  
        }
    }
fout.close();
return 0;
}	



void ExprPredictor::printFiled( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
    
    
os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t"; 
        }        
    }
os << endl<<" txpEffects : "<< endl;
 
    for ( int i = 0; i < nFactors(); i++ ) {
       os << par.txpEffects[i] << "\t";
       
    }
os << endl<<  " coopertivity : " << endl;

 for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if ( coopMat(i,j)){
			for(int k =0; k < par.theV[i][j].size(); k++) { 
		
			os <<  "\t" << par.theV[i][j][k] ;
			}
		}	
			
  		
	}
}
os<< endl << " quenching : " << endl;


for(int i =0; i < nFactors(); i++){
	for(int j =0; j <i; j++){
		if(repIndicators[i])  {          
			for(int k =0; k < ExprPar::nbins ; k++){ 
			os <<  "\t" << par.theVr[i][j][k] ;
				
  			}
		}
	}
}

os << endl << "objective function : " << endl;
    os << jRMSE.getObj() << endl;
   
ofstream fo( "Data/format.tex",ios::app );
int j=0;

vector< Site > tsites;
vector< Site > tsitesbot;
for(int m = 0; m < seqSites.size(); m++ ) {
tsites = seqSites[m];
tsitesbot = seqSitesbot[m]  ;
 
	ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
	
	foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		
		for (;;) {
			if ( j < tsites[i].start ) { foo << "b"; j++;  }
			
			else {
				

				if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  
				if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  
				if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	
	foo << endl;
	foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
		
		for (;;) {
			if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
			
			else {
				

				if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  
				if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  
				if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	foo << endl;
	foo.close();



 j=0;

	for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
		if( i == 0 ) { fo << "&&" ; continue; }
		fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
		if(i% 90 == 0 ) {fo << "\\\\&&"; }
	}
fo << "\\\\&&\\\\";
ExprFunc* func = this->createExprFunc( par_model );	
   
        vector< double > concs = factorExprData.getCol( cell[2*m] );
        double predicted = func->predictExpr(tsites, 5, concs );
fo << ExprPredictor::seqNmes[m]<< predicted  << "&" << "pcell "<< cell[2*m] << "&" ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
	
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
	
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< 0 << "&" ;
	for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
	
	
		for (;;) {
			if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< 10 << "&" ;
	for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
	
	
		for (;;) {
			if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

for( int e = 0; e < d[m].size() ; e++ ){  
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< 0 << "&" ;
	for( int i = 0; i < d[m][e].size() ; i++ ) {
	
	
		for (;;) {
			if ( j < d[m][e][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( d[m][e][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( d[m][e][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( d[m][e][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

} 
fo << "\\\\&&\\\\&&\\\\" << endl;

}
fo.close();


}
		
void ExprPredictor::printFile2( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
    
 
ofstream fo( "format.tex",ios::app );
int j=0;

vector< Site > tsites;
vector< Site > tsitesbot;
 double bottomofdborder = .3;
 double topofdborder = .8;
 int nrow = factorExprData.nRows();    
 int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  
int tits;
int tibs;
for (int i=0; i<nrow;i++) {
	vector< double > reD;	
	reD = factorExprData.getRow(i);	
	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	
	double m;
	m = gsl_vector_max(rowexprData);      
	
	double exptrace;
	double exppeek;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           
		gsl_vector_set(transition_Indicests, i , tits);  
		gsl_vector_set(transition_Indicesbs, i , tibs);
		
		break;
	}
	else{

	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  
			
			if(exppeek < topofdborder ) { 
				        tits =j;
					
					gsl_vector_set(transition_Indicests, i , tits); 
					
					
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}
					} 
						
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   
			break;
			} 
			
			
			
		}
	}
}

}

 




        	  
	
          

for(int m = 0; m < seqSites.size(); m++ ) {
tsites = seqSites[m];
tsitesbot = seqSitesbot[m]  ;
 
	
	ofstream foo( "/home/jacobc/Desktop/format2.txt",ios::app );
	foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		
		for (;;) {
			if ( j < tsites[i].start ) { foo << "b"; j++;  }
			
			else {
				

				if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  
				if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  
				if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	
	foo << endl;
	foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
		
		for (;;) {
			if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
			
			else {
				

				if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  
				if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  
				if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {
		foo << "b"; j++; 
	}
	foo << endl;
	foo.close();


 j=0;

	for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
		if( i == 0 ) { fo << "&&" ; continue; }
		fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
		if(i% 90 == 0 ) {fo << "\\\\&&"; }
	}
fo << "\\\\&&\\\\";
   
     
       
fo << ExprPredictor::seqNmes[m]<< " "  << "&" << "cell "<< jRMSE.cell[2*m] << "&" ;  
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
	
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;

fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
	for( int i = 0; i < tsitesbot.size() ; i++ ) {
	
	
		for (;;) {  
			if ( j < tsitesbot[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<<"m1" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
	for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
	
	
		for (;;) {
			if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "m2" << "&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
	for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
	
	
		for (;;) {
			if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "d+" << "&" << "cell "<< cell[2*m]<< "&" ;
	for( int i = 0; i < seqSitesf2[m].size() ; i++ ) {
	
	
		for (;;) {
			if ( j < seqSitesf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( seqSitesf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( seqSitesf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "d+"  << "&" << "cell "<< cell[2*m+1]<< "&" ;
	for( int i = 0; i < seqSitesbotf2[m].size() ; i++ ) {
	
	
		for (;;) {
			if ( j < seqSitesbotf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesbotf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( seqSitesbotf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( seqSitesbotf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "m1"<< "d+" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
	for( int i = 0; i < seqSitesm1f2[m].size() ; i++ ) {
	
	
		for (;;) {
			if ( j < seqSitesm1f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm1f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( seqSitesm1f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( seqSitesm1f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] <<"m2"<< "d+" <<"&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
	for( int i = 0; i < seqSitesm2f2[m].size() ; i++ ) {
	
	
		for (;;) {
			if ( j < seqSitesm2f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( seqSitesm2f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( seqSitesm2f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( seqSitesm2f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
				
				
			 }
		}
	}
	while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}


fo << "\\\\&&\\\\&&\\\\" << endl;

}
fo.close();


}

void ExprPredictor::printFile24( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
ofstream tw( "Data/TwistCondit24.txt");
ofstream du( "Data/DorsalUnCondit24.txt" );
ofstream dc( "Data/DorsalCondit24.txt" );
ofstream dcom( "DorsalCombin24.txt" );
ofstream du3( "Data/DUCondit3.txt" );
ofstream dc3( "Data/DCondit3.txt" );
vector< Site > tsites;
for(int m = 0; m < ExprPredictor::seqsya.size(); m++ ) {
	for(int k = 0 ; k < motifs.size() ; k++ ){
	
	tsites = d[m][k];

	for( int i = 0; i < tsites.size() ; i++ ) {			

	
	 			
				if( tsites[i].factorIdx == 0 ) { 
				if(  tsites[i].start-5  <= 0 ){continue;}

				if(  tsites[i].start+motifs[ tsites[i].factorIdx ].length() + 5 > ExprPredictor::seqsya[m].size() ){continue;}
					 int j = 0;
					dc3 << tsites[i].energy << endl;
					dc << ">" << ExprPredictor::seqNmesa[m] << "_" << tsites[i].energy<<endl ;
					int shift = motifs[ tsites[i].factorIdx ].length()+7;
					for( int p =0 ; p < shift; p++ ){
						j = tsites[i].start -3 ;
						if(tsites[i].strand){
						
						dc <<  ALPHABET[ ExprPredictor::seqsya[m][ j+p ]  ] ;
						}
					else{ dc << ALPHABET[complement(symbolToInt(ALPHABET[ ExprPredictor::seqsya[m][j+ shift -1 - p ] ]))];} 				
						} 
					
					dc << endl;
				 }  
				if( tsites[i].factorIdx == 2 ) {
					 int j = 0;
					tw << ">" << ExprPredictor::seqNmesa[m] << "_" << tsites[i].start<< "_" << tsites[i].strand <<endl ;
					for( int p = 0; p < motifs[ tsites[i].factorIdx ].length(); p++ ){
						j = tsites[i].start ;
					if(tsites[i].strand){
					
					tw <<  ALPHABET[ ExprPredictor::seqsya[m][ j+p ]  ] ;
					}
				else{ tw << ALPHABET[complement(symbolToInt(ALPHABET[ ExprPredictor::seqsya[m][ j+ motifs[tsites[i].factorIdx ].length() -1 - p ] ]))];} 				
					} 
					tw << endl;
				 } 
				if( tsites[i].factorIdx == 1 ) { 
					 int j = tsites[i].start ;
					du3 << tsites[i].energy << endl;
					du << ">" << ExprPredictor::seqNmesa[m] << "_" << tsites[i].energy<<endl ;
					for( int p = 0; p < motifs[ tsites[i].factorIdx ].length() ; p++ ){
						
					if(tsites[i].strand){
					
					du <<  ALPHABET[ ExprPredictor::seqsya[m][ j+p ]  ] ;
					}
					else{ du << ALPHABET[complement( symbolToInt( ALPHABET[ ExprPredictor::seqsya[m][j+ motifs[ tsites[i].factorIdx ].length()-1 - p ] ]) )];} 				
					} 
					du << endl;
				 }


				if( tsites[i].factorIdx == 3 ) { 
					 int j = 0;
					dcom << ">" << ExprPredictor::seqNmesa[m] << "_" << tsites[i].start<< "_" << tsites[i].strand <<endl ;
					for( int p = 0; p < motifs[ tsites[i].factorIdx ].length() ; p++ ){
						j = tsites[i].start;
					if(tsites[i].strand){
					
					dcom <<  ALPHABET[ ExprPredictor::seqsya[m][ j+p ]  ] ;
					}
				else{ dcom <<ALPHABET[ complement(symbolToInt(ALPHABET[ ExprPredictor::seqsya[m][j+motifs[ tsites[i].factorIdx ].length() -1 - p ] ]))];} 				
					} 
					dcom << endl;
				 }
	} 
       } 
} 
dc.close();
tw.close();
du.close();
dcom.close();
du3.close();
dc3.close();

}



void ExprPredictor::printFile22( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{

ofstream tw( "TwistCondit.txt", ios::app );
ofstream du( "DorsalUnCondit.txt" , ios::app);
ofstream dc( "DorsalCondit.txt", ios::app );

vector< Site > tsites;


for(int m = 0; m < seqSitesa.size(); m++ ) {
tsites = seqSitesa[m];
	for( int i = 0; i < tsites.size() ; i++ ) {			

				if( tsites[i].factorIdx == 0 ) { 
					int j = 0;
					dc << ">" << ExprPredictor::seqNmesa[m]<< "_"  << tsites[i].start  <<endl ;
					for( int p = 0; p < motifs[ tsites[i].factorIdx ].length() ; p++ ){
						j = tsites[i].start + p;
					dc <<  ALPHABET[ ExprPredictor::seqsya[m][ j ]  ] ;
					} 
					dc << endl;
				 }  
				if( tsites[i].factorIdx == 1 ) {
					
					int j = 0;
					tw << ">" << ExprPredictor::seqNmesa[m]<< "_"  <<  tsites[i].start << endl ;
					for( int p = 0; p < motifs[ tsites[i].factorIdx ].length() ; p++ ){
						j = tsites[i].start + p;
						tw <<  ALPHABET[ ExprPredictor::seqsya[m][ j] ] ;
					}
					tw << endl;
				 } 
				if( tsites[i].factorIdx == 2 ) { 
					 int j = 0;
					du << ">" << ExprPredictor::seqNmesa[m] << "_" << tsites[i].start <<endl ;
					for( int p = 0; p < motifs[ tsites[i].factorIdx ].length() ; p++ ){
						j = tsites[i].start + p;
					du <<  ALPHABET[ ExprPredictor::seqsya[m][ j ]  ] ;
					} 
					du << endl;
				 }
	}
}
dc.close();
tw.close();
du.close();
}
void ExprPredictor::printFile23( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 ); 
    
 
ofstream fo( "format.tex",ios::app );
int j=0;

vector< Site > tsites;
for(int m = 0; m < seqSitesa.size(); m++ ) {
tsites = seqSitesa[m];
 j=0;

	for(int i = 0; i < ExprPredictor::seqsya[m].size(); i++ ) {
		if( i == 0 ) { fo << "&&" ; continue; }
		fo << ALPHABET[ExprPredictor::seqsya[m][i]] ;
		if(i% 90 == 0 ) {fo << "\\\\&&"; }
	}
fo << "\\\\&&\\\\";
   
     
       
fo << ExprPredictor::seqNmesa[m]<< " "  << "&" << "cell neuro " << "&" ;  
	for( int i = 0; i < tsites.size() ; i++ ) {
	
		
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }

				if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  
				if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{D}\\color{black}"; j++ ; break; }
	
			 }
		}
	}
	while( j < ExprPredictor::seqsya[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

fo << "\\\\&&\\\\&&\\\\" << endl;

}
fo.close();


}
void ExprPredictor::printFile25( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
  os.setf( ios::fixed );
    os.precision( 3 ); 
    
 
ofstream fo( "format.tex",ios::app );
int j=0;
vector< Site > tsites;

for(int m = 0; m < ExprPredictor::seqsya.size(); m++ ) {
	for(int k = 0 ; k < 2 ; k++ ){
	
	
	 j=0;
	if( k == 0 ){
		tsites = d[m][k];
		for(int i = 0; i < ExprPredictor::seqsya[m].size(); i++ ) {
			if( i == 0 ) { fo << "&&" ; continue; }
			fo << ALPHABET[ExprPredictor::seqsya[m][i]] ;
			if(i% 90 == 0 ) {fo << "\\\\&&"; }
		}
		fo << "\\\\&&\\\\";
	}
	if( k ==1){
	fo << ExprPredictor::seqNmesa[m]<< " "  << "&" << "   " << "&" ;  
	for( int i = 0; i < d[m][k].size() ; i++ ) {	
		tsites.push_back( d[m][k][i] );
	}
	for( int i = 0; i < d[m][2].size() ; i++ ) {	
		tsites.push_back( d[m][2][i] );
	}
	for( int i = 0; i < d[m][3].size() ; i++ ) {	
		tsites.push_back( d[m][3][i] );
	}
	std::sort(tsites.begin(), tsites.end(), siteSortPredicate); 

	for( int i = 0; i < tsites.size() ; i++ ) {
	
		
		for (;;) {
			if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
			
			else {
				if( j % 90 == 0 ) { fo << "\\\\&&"; }
				if( i > 0 ){
					if ( tsites[i].start == tsites[i-1].start  ) { fo << "!"; j++;}
				}
				if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  
				if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{u}\\color{black}"; j++ ; break; }  
				if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{t}\\color{black}"; j++ ; break; }
				if( tsites[i].factorIdx == 3 ) { fo << "\\color{blue}{c}\\color{black}"; j++ ; break; }
	
			 }
		}
	}
	while( j < ExprPredictor::seqsya[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
		fo << "e"; j++; 
	}

	fo << "\\\\&&\\\\&&\\\\" << endl;
	} 
    } 
}
fo.close();

}

int ExprPredictor::train()
{	
    
	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		
    
    
    ExprPar par_default( nFactors() );
    train( par_default, rng ); 
    
    return 0;	
}

int ExprPredictor::predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs ) 
{
    targetExprs.clear();
    
    
     
    
    ExprFunc* func = this->createExprFunc( par_model );	
    for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( j );
        double predicted = func->predictExpr( targetSites, targetSeqLength, concs );
        targetExprs.push_back( predicted );
    }
    
    return 0; 
}



ModelType ExprPredictor::modelOption = LOGISTIC;
int ExprPredictor::estBindingOption = 1;    
ObjType ExprPredictor::objOption = SSE;

double ExprPredictor::exprSimCrossCorr( const vector< double >& x, const vector< double >& y )
{
    vector< int > shifts; 
    for ( int s = -maxShift; s <= maxShift; s++ ) {
        shifts.push_back( s ); 
    }

    vector< double > cov; 
    vector< double > corr; 
    cross_corr( x, y, shifts, cov, corr ); 
    double result = 0, weightSum = 0; 
    result = *max_element( corr.begin(), corr.end() );

    return result; 
}

int ExprPredictor::maxShift = 5; 
double ExprPredictor::shiftPenalty = 0.8; 
int ExprPredictor::nAlternations = 4;
int ExprPredictor::nRandStarts = 0;
double ExprPredictor::min_delta_f_SSE = 1.0E-8;
double ExprPredictor::min_delta_f_Corr = 1.0E-8;
double ExprPredictor::min_delta_f_CrossCorr = 1.0E-8;
int ExprPredictor::nSimplexIters = 3;
int ExprPredictor::nGradientIters = 6;

int ExprPredictor::randSamplePar( const gsl_rng* rng, ExprPar& par ) const
{
    int counter= 0 ;
    
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
	    if( indicator_bool[counter] ) {
            double rand_weight = exp( gsl_ran_flat( rng, log( ExprPar::min_weight ), log( ExprPar::max_weight ) ) ); 
            par.maxBindingWts[i] = rand_weight;
	    counter++;
	    }
	    else{ counter++; } 
        }        
    }

if(modelOption == BINS){
for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) {
		
		for(int k=0; k< ExprPar::nbins; k++){  
			 
		if( indicator_bool[counter] ) {
		double rand_interaction = exp( gsl_ran_flat( rng, log( ExprPar::min_interaction ), log( ExprPar::max_interaction ) ) );
                par.theV[ i][ j][k] = rand_interaction;	   
		counter++;
		}       
		else{counter++;}
                 } 
		
            }  
           
        }
    }
   
for ( int i = 0; i < nFactors(); i++ ) {
       
	for ( int j = 0; j < i; j++ ) {
            if (repIndicators[i] ) {
		
		for(int k=0; k< ExprPar::nbins; k++){  
		if( indicator_bool[counter] ) {
			 
double rand_interaction = exp( gsl_ran_flat( rng, log( ExprPar::min_interactionr ), log( ExprPar::max_interactionr ) ) );
                par.theVr[ i][ j][k] = rand_interaction;
		counter++;
		}
		else{counter++;}
		}
			              
                
            }  
            
        }
    }
    
}
  
    
    for ( int i = 0; i < nFactors(); i++ ) {
       if( indicator_bool[counter] ) {
            double rand_effect =  gsl_ran_flat( rng,  ExprPar::min_effect_Thermo , ExprPar::max_effect_Thermo  );
            par.txpEffects[i] = rand_effect;
	    counter++;
	    }
	else{counter++;}
       
    }
   
    return 0;
}

bool ExprPredictor::testPar( const ExprPar& par ) const
{
    
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( par.maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) || par.maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) )
                return false; 
        }        
    }

 
   if( modelOption == BINS) {
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
           if(coopMat(i,j)){
             for(int k =0; k < ExprPar::nbins; k++){
                    if ( par.theV[i][j][k] < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) || par.theV[i][j][k] > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) { 
			return false; }
             }
           }
        }
     }
     for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
           if(repIndicators[i]){
             for(int k =0; k < ExprPar::nbins; k++){
                  if ( par.theVr[i][j][k] < ExprPar::min_interactionr * ( 1.0 + ExprPar::delta ) || par.theVr[i][j][k] > ExprPar::max_interactionr * ( 1.0 - ExprPar::delta ) ) { 
			return false; }
             }
           }
        }
     }
    }
 
    for ( int i = 0; i < nFactors(); i++ ) {
       
            if ( par.txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) || par.txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) )
                return false;
        
    }

    return true;    
}
void ExprPredictor::printPar2(  ) 
{ 

}
void ExprPredictor::printPar( const ExprPar& par ) const
{
    
     
    
    
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
		if(i == 0){
            		
		}
		if(i==1){
			
		}
		if(i==2){
			
		} 
        }        
    }
   



  for(int i =0; i < nFactors(); i++){
	for(int j =0; j <= i; j++){
		if(coopMat(i,j))  { 
			for(int k =0; k < ExprPar::nbins ; k++){   
				if( i == 1 && j ==0 ) { 
				}
				
				if( i == 2 && j ==0 ) {  
				}
				
				
				if( i== 2 && j == 1) { 
				}
				
			}
			
  		}
  	}
  }

for(int i =0; i < nFactors(); i++){
	for(int j =0; j <=i; j++){
		if(repIndicators[i])  { 
			for(int k =0; k < ExprPar::nbins ; k++){   
				if( i == 1 && j ==0 ) { 
				}
				
				if( i == 2 && j == 0) {  
				}
				
				
				if( i== 2 && j == 1) { 
				}
				
			}
			
  		}
  	}
  }


   
    for ( int i = 0; i < nFactors(); i++ ) {
      
        
    }


}

ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par ) const
{	
    return new ExprFunc( motifs, intFunc, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, par,  binBord, binBordr );
}
ExprFunc* ExprPredictor::createExprFunc2(  ExprPar& par )  
{	
    return new ExprFunc( motifs, intFunc, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, par, binBord , binBordr);
}

double ExprPredictor::compRMSE( const ExprPar& par ) const
{

    
    ExprFunc* func = createExprFunc( par );
            double rss = 0;
    
    double squaredErr = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            
            
            double predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            
            
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }


        double beta;
        squaredErr += least_square( predictedExprs, observedExprs, beta );
    	
      for ( int i = 0; i < nConds(); i++ ) {  
        rss += ( predictedExprs[i] -  observedExprs[i] ) * ( predictedExprs[i] -  observedExprs[i] ); 
    }   
   }
double rmse =  rss/ ( nSeqs() * nConds() ) ; 
 
    return rmse;
  
}

double ExprPredictor::compAvgCorr( const ExprPar& par ) const
{
	ExprFunc* func = createExprFunc( par );
	vector< double > fOcc;
	vector< double > corrs; 	
	for ( int i = 0; i < seqSitesb.size(); i++ ) {   
		vector< double > predicted; 
		for ( int j = 0; j < seqSitesb[ i ].size(); j++ ) {			
			
			vector< double > concs = factorExprData.getCol( 0 );  

			func->predictOcc( seqSitesb[ i ][ j ], i, concs,fOcc );
 			predicted.push_back( fOcc[i] );
		}
		corrs.push_back( correlation( predicted, bindingData[ i ] ) );
	}	
	
	return mean( corrs ); 
}

void ExprPredictor::predict2(vector< vector< double > >& occs ) 
{
	ExprFunc* func = createExprFunc( par_model );
	vector< double > fOcc;
	vector< double > corrs;
	vector< vector< double > > occstemp(  seqSitesb.size(),vector<double>(seqSitesb[0].size()) ); 
	for ( int i = 0; i < seqSitesb.size(); i++ ) {   
		
		for ( int j = 0; j < seqSitesb[ i ].size(); j++ ) {		
			vector< double > concs = factorExprData.getCol( 0 );  

			func->predictOcc( seqSitesb[ i ][ j ], i, concs,fOcc );
 
			
			occstemp[i][j] =  fOcc[i] ;
					
		}
	}	
	
occs = occstemp;
}
double ExprPredictor::compAvgCorr2(  ExprPar& par, int thres ) 
{
    
    ExprFunc* func = createExprFunc( par );
         vector< double > corrs; 

Matrix f = factorExprData;
Matrix e = exprData;
           
	   anny.annotydorsalZ2( seqsya, d, f, e, *func , thres, ExprPredictor::seqNmesa);
  

return .5;

}

double ExprPredictor::compAvgCorr3(  ExprPar& par, int thres ) 
{
    
    ExprFunc* func = createExprFunc( par );
         vector< double > corrs; 

Matrix f = factorExprData;
Matrix e = exprData;
ofstream he2("Data/hel2.txt");
he2 << d << endl;
he2.close();
           
	   anny.annotydorsalZarend2u( seqsya, d, f, e, *func , thres, ExprPredictor::seqNmesa, seqSitesbot, seqSitesm1);
  

ofstream he("Data/hel.txt");
he << d << endl;
he.close();

return .5;

}









double ExprPredictor::compAvgCorrborder(  ExprPar& par ) 
{
double expmin = .35;
double expmax = .65;
double bottomofdborder = .35;
double topofdborder = .65;


par_model = par;


    ExprFunc* func = createExprFunc( par );
         vector< double > corrs; 
vector< string > rowLabels(0);
 vector< string > colLabels(0);
 vector< vector< double > > data(0) ;

    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
gsl_rng* rng;
gsl_rng_env_setup();
const gsl_rng_type * TT = gsl_rng_default;	
rng = gsl_rng_alloc( TT );
gsl_rng_set( rng, time( 0 ) );		
double  rand_site_index;	
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    
    string line, first, label;
	
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
   
    while ( ss >> label ) {
	
	
	
	
	
	
	
	
	colLabels.push_back( label );
    }
    
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){ 
		
		
		 	

	vals.push_back( val  );		
	
			

	}
	data.push_back( vals );
    }
Matrix data1( data );
Matrix f = factorExprData;
exprData = data;
Matrix e = data;

int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indicesb = gsl_vector_alloc(nrow);  
gsl_vector *transition_Indicest = gsl_vector_alloc(nrow);  
int tit;
int tib;
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	
	double m;
	m = gsl_vector_max(rowexprData);      
	
	double exptrace;
	double exppeek;
	tit = 0;
	tib =0;
	if( m < bottomofdborder ) {                           
		gsl_vector_set(transition_Indicest, i , tit);  
		gsl_vector_set(transition_Indicesb, i , tib);
		
		break;
	}
	else{

	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  
			
			if(exppeek < topofdborder ) { 
				        tit =j;
					 
					gsl_vector_set(transition_Indicest, i , tit); 
					
					
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tit+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tit + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tib=tit + counter;
									
									gsl_vector_set(transition_Indicesb, i , tib);
									break;	
								}
					} 
						
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tit == 0) { 
				  tit = minindex;            
				  tib = minindex;
				gsl_vector_set(transition_Indicest, i , tit); 
				gsl_vector_set(transition_Indicesb, i , tib); 
			} 	   
		 	                   
			break;
			} 
			
			
			
		}
	}
}
}

 vector< double > predictedExprs;
        vector< double > observedExprs;
    
    double totalSim = 0;
double rms = 0;



int c = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {  
		c++;
              
		cell.push_back( gsl_vector_get(transition_Indicest,i) );
		cell.push_back( gsl_vector_get(transition_Indicesb,i) );
		if( gsl_vector_get(transition_Indicest,i) == 0 ) {
			double predicted = 0;   
            		predictedExprs.push_back( predicted );
         
           	
            		double observed = .7;  
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
			double predicted = 0;   
            		predictedExprs.push_back( predicted );
         
           
            		double observed = .7;  
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		else{
	
		   vector< double > concst = factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );
        	   anny.annotydorsal( seqsy[ i], seqSites[ i ], f, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);

  
        	    double predictedt = func->predictExpr( seqSites[ i ], seqLengths[i], concst );
        	    predictedExprs.push_back( predictedt );
  
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - observedt );

		    vector< double > concs2 = factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );
        	    anny.annotydorsal( seqsy[ i], seqSitesbot[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
             
        	    double predicted2 = func->predictExpr( seqSitesbot[ i ], seqLengths[i], concs2 );
        	    predictedExprs.push_back( predicted2 );
         
           
        	    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
           
		     observedExprs.push_back( observed2 );
 		     rms += abs( predicted2 - observed2 );   
		    
       		}
	}  
 bottomofdborder = .5;
 topofdborder = .8;
 nrow = factorExprData.nRows();    
 ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  
int tits;
int tibs;
for (int i=0; i<nrow;i++) {
	vector< double > reD;	
	reD = factorExprData.getRow(i);	
	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	
	double m;
	m = gsl_vector_max(rowexprData);      
	
	double exptrace;
	double exppeek;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           
		gsl_vector_set(transition_Indicests, i , tits);  
		gsl_vector_set(transition_Indicesbs, i , tibs);
		
		break;
	}
	else{

	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  
			
			if(exppeek < topofdborder ) { 
				        tits =j;
					
					gsl_vector_set(transition_Indicests, i , tits); 
					
					
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}
					} 
						
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   
			break;
			} 
			
			
			
		}
	}
}

}




 for ( int i = 0; i < nSeqs(); i++ ) {  
cell.push_back( gsl_vector_get(transition_Indicests,2) );
		
		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
        	  
	
          
		    anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
        	    double predictedt = func->predictExpr( seqSitesm1[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - observedt );
}


for ( int i = 0; i < nSeqs(); i++ ) {

		cell.push_back( gsl_vector_get(transition_Indicesbs,2) );
		 vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicesbs,2) );
        	 
		   anny.annoty3( seqsy[ i], seqSitesm2[ i ], f, e, *func , gsl_vector_get(transition_Indicesbs,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
           
        	    double predictedt = func->predictExpr( seqSitesm2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - observedt );
}



Matrix f2 =factorExprData;
vector< double > temp = factorExprData.getRow(0);
double noise = .1;
for( int i = 0 ; i < temp.size(); i++ ){
temp[i]= temp[i] + noise;
if(temp[i] > 1 ) { temp[i] = 1 ;}
}
f2.setRow(0,temp);

 for ( int i = 0; i < nSeqs(); i++ ) {  
               	
		cell.push_back( gsl_vector_get(transition_Indicest,i) );
		cell.push_back( gsl_vector_get(transition_Indicesb,i) );
		if( gsl_vector_get(transition_Indicest,i) == 0 ) {
			double predicted = 0;   
            		predictedExprs.push_back( predicted );
     
           
            		double observed = .7;  
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
			double predicted = 0;   
            		predictedExprs.push_back( predicted );
         
           
            		double observed = .7;  
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		else{
		

		   vector< double > concst = f2.getCol( gsl_vector_get(transition_Indicest,i) );
        	   anny.annotydorsal( seqsy[ i], seqSitesf2[ i ], f2, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
           
        	    double predictedt = func->predictExpr( seqSitesf2[ i ], seqLengths[i], concst );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - observedt );

		    vector< double > concs2 = f2.getCol( gsl_vector_get(transition_Indicesb,i) );
        	    anny.annotydorsal( seqsy[ i], seqSitesbotf2[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
           
        	    double predicted2 = func->predictExpr( seqSitesbotf2[ i ], seqLengths[i], concs2 );
        	    predictedExprs.push_back( predicted2 );
         
           
        	    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
           
		     observedExprs.push_back( observed2 );
 		     rms += abs( predicted2 - observed2 );   
		    
       		}
	}  
 bottomofdborder = .5;
 topofdborder = .8;
 nrow = factorExprData.nRows();    
 ncol =factorExprData.nCols();
for (int i=0; i<nrow;i++) {
	vector< double > reD;	
	reD = factorExprData.getRow(i);	
	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	
	double m;
	m = gsl_vector_max(rowexprData);      
	
	double exptrace;
	double exppeek;
	tits = 0;
	tibs =0;
	if( m < bottomofdborder ) {                           
		gsl_vector_set(transition_Indicests, i , tits);  
		gsl_vector_set(transition_Indicesbs, i , tibs);
		
		break;
	}
	else{

	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  
			
			if(exppeek < topofdborder ) { 
				        tits =j;
				
					gsl_vector_set(transition_Indicests, i , tits); 
					
					
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tits+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tits + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tibs=tits + counter;
									
									gsl_vector_set(transition_Indicesbs, i , tibs);
									break;	
								}
					} 
						
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tits == 0) { 
				  tits = minindex;            
				  tibs = minindex;
				gsl_vector_set(transition_Indicests, i , tits); 
				gsl_vector_set(transition_Indicesbs, i , tibs); 
			} 	   
		 	                   
			break;
			} 
			
			
			
		}
	}
}

}




 for ( int i = 0; i < nSeqs(); i++ ) {  
cell.push_back( gsl_vector_get(transition_Indicests,2) );
		
		 vector< double > concsm = f2.getCol( gsl_vector_get(transition_Indicests,2) );  
        	  
	
          
		    anny.annoty3( seqsy[ i], seqSitesm1f2[ i ], f2, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSitesf2[i]);
        	    double predictedt = func->predictExpr( seqSitesm1f2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - observedt );
}


for ( int i = 0; i < nSeqs(); i++ ) {

		cell.push_back( gsl_vector_get(transition_Indicesbs,2) );
		 vector< double > concsm = f2.getCol( gsl_vector_get(transition_Indicesbs,2) );
        	 
		   anny.annoty3( seqsy[ i], seqSitesm2f2[ i ], f2, e, *func , gsl_vector_get(transition_Indicesbs,2), ExprPredictor::seqNmes[i],i, seqSitesf2[i]);
           
        	    double predictedt = func->predictExpr( seqSitesm2f2[ i ], seqLengths[i], concsm );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - observedt );
}
return rms/(nSeqs()*8);

}


double ExprPredictor::compAvgCorrborder2(  ExprPar& par ) 
{
double expmin = .35;
double expmax = .65;
double bottomofdborder = .35;
double topofdborder = .65;


par_model = par;


    ExprFunc* func = createExprFunc( par );
         vector< double > corrs; 
vector< string > rowLabels(0);
 vector< string > colLabels(0);
 vector< vector< double > > data(0) ;

    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
gsl_rng* rng;
gsl_rng_env_setup();
const gsl_rng_type * TT = gsl_rng_default;	
rng = gsl_rng_alloc( TT );
gsl_rng_set( rng, time( 0 ) );		
double  rand_site_index;	
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    
    string line, first, label;
	
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
   
    while ( ss >> label ) {
	
	
	
	
	
	
	
	
	colLabels.push_back( label );
    }
    
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){ 
		
		
		 	

	vals.push_back( val  );		
	
			

	}
	data.push_back( vals );
    }
Matrix data1( data );
Matrix f = factorExprData;
exprData = data;
Matrix e = data;

int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indicesb = gsl_vector_alloc(nrow);  
gsl_vector *transition_Indicest = gsl_vector_alloc(nrow);  
int tit;
int tib;
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	
	double m;
	m = gsl_vector_max(rowexprData);      
	
	double exptrace;
	double exppeek;
	tit = 0;
	tib =0;
	if( m < bottomofdborder ) {                           
		gsl_vector_set(transition_Indicest, i , tit);  
		gsl_vector_set(transition_Indicesb, i , tib);
		
		break;
	}
	else{

	 
        for (int j=mi; j< ncol; j++)  {  
		exptrace = gsl_vector_get(rowexprData,j);
		exppeek = gsl_vector_get(rowexprData,j+1);		        
		if ( exppeek > topofdborder ) {
			continue;
	        }
		else {
	        	if( exptrace == exppeek) { continue; }  
			
			if(exppeek < topofdborder ) { 
				        tit =j;
					
					gsl_vector_set(transition_Indicest, i , tit); 
					
					
					int counter=0;	
					for(;;) {
								exptrace = gsl_vector_get(rowexprData,tit+counter);
								counter++;
								exppeek = gsl_vector_get(rowexprData,tit + counter);	
								if ( exppeek > bottomofdborder ) {
									continue;
								}
								else {
									tib=tit + counter;
									
									gsl_vector_set(transition_Indicesb, i , tib);
									break;	
								}
					} 
						
			int minindex;				
			minindex = gsl_vector_min_index(rowexprData);  
		     	if (tit == 0) { 
				  tit = minindex;            
				  tib = minindex;
				gsl_vector_set(transition_Indicest, i , tit); 
				gsl_vector_set(transition_Indicesb, i , tib); 
			} 	   
		 	                   
			break;
			} 
			
			
			
		}
	}
}
}

 vector< double > predictedExprs;
        vector< double > observedExprs;
    
    double totalSim = 0;
double rms = 0;



Matrix f2 =factorExprData;
vector< double > temp = factorExprData.getRow(0);
double noise = .1;
for( int i = 0 ; i < temp.size(); i++ ){
temp[i]= temp[i] + noise;
if(temp[i] > 1 ) { temp[i] = 1 ;}
}

 for ( int i = 0; i < nSeqs(); i++ ) {  
               	
		
		
		if( gsl_vector_get(transition_Indicest,i) == 0 ) {
			double predicted = 0;   
            		predictedExprs.push_back( predicted );
        
           
            		double observed = .7;  
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
			double predicted = 0;   
            		predictedExprs.push_back( predicted );
         
           
            		double observed = .7;  
           
	     		observedExprs.push_back( observed );
 	    		rms += abs( predicted - 0 );
		}
		else{
		

		   vector< double > concst = f2.getCol( gsl_vector_get(transition_Indicest,i) );
        	   anny.annotydorsal( seqsy[ i], seqSitesf2[ i ], f2, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
           
        	    double predictedt = func->predictExpr( seqSitesf2[ i ], seqLengths[i], concst );
        	    predictedExprs.push_back( predictedt );
         
           
        	    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
           
		     observedExprs.push_back( observedt );
 		     rms += abs( predictedt - observedt );

		    vector< double > concs2 = f2.getCol( gsl_vector_get(transition_Indicesb,i) );
        	    anny.annotydorsal( seqsy[ i], seqSitesbotf2[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
           
        	    double predicted2 = func->predictExpr( seqSitesbotf2[ i ], seqLengths[i], concs2 );
        	    predictedExprs.push_back( predicted2 );
         
           
        	    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
           
		     observedExprs.push_back( observed2 );
 		     
		    
       		}
	}  

return rms/nSeqs();

}
int Random(int n)
{
	return rand() % n ;
}

int ExprPredictor::createccdata()
{
double expmin = .35;
double expmax = .65;

int k=3;



    
vector< string > rowLabels(0);
 vector< string > colLabels(0);
 vector< vector< double > > data(0) ;

    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
gsl_rng* rng;
gsl_rng_env_setup();
const gsl_rng_type * TT = gsl_rng_default;	
rng = gsl_rng_alloc( TT );
gsl_rng_set( rng, time( 0 ) );		
double  rand_site_index;	
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) colLabels.push_back( label );
    
    
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){ 
	 rand_site_index = gsl_ran_flat( rng, -.15, .15 );
	if( val + rand_site_index > 1 ) {
	 vals.push_back( 1 );
	}
	
	if( val + rand_site_index < 0 ) {
	 vals.push_back( 0 );
	}
	
	if( val + rand_site_index >= 0  && val + rand_site_index <= 1 ) {
	vals.push_back( val + rand_site_index );
	}
	} 
        data.push_back( vals );
    }
Matrix data1( data );
Matrix f = factorExprData;
exprData = data;
Matrix e = data;

int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indices = gsl_vector_alloc(nrow);  
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	
	double m;
	m = gsl_vector_max(rowexprData);
	int ti;
	ti = 0;
if( m < expmin ) {
	gsl_vector_set(transition_Indices, i , ti);  
	
	break;
	}
else{
	 
        for (int j=mi; j< ncol; j++)  {  
		
		if (expmax < gsl_vector_get(rowexprData,j) ) {
		
			continue;
	        }
		
		else {

			ti = j; break;
	        }
	}
        int minindex;				
	minindex = gsl_vector_min_index(rowexprData);  
	
     	if (ti == 0) { 
	          ti = minindex;            
	} 	   
  
 	gsl_vector_set(transition_Indices, i , ti);                    
}
}

for (int i=0; i<nrow;i++) {

}
ifstream seq_file("efastanotw10txt");  
	ifstream expr_file( "expre10.tab");
assert (seq_file.is_open() &&  expr_file.is_open());
string header;
	
	getline(expr_file, header);
	
	string temp;
	vector <string> seq_name;
	vector <string> seq;
	vector <string> expr;
	
	seq.clear();
	expr.clear();
	
	
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
		
		string name1 (temp, 1, temp.length() - 1);
		seq_name.push_back(name1);
		
		getline(seq_file, temp);
		seq.push_back(temp);
				
		getline(expr_file, temp);
		expr.push_back(temp);
	}
	
	int size = seq.size();
		ofstream meso_seq_file ("mesoseq.txt");
		ofstream meso_expr_file ("mesoexpr.txt");
		
		ofstream neuro_seq_file("neuroseq.txt");
		ofstream neuro_expr_file("neuroexpr.txt");
		
		assert (meso_seq_file.is_open() &&  meso_expr_file.is_open() && neuro_seq_file.is_open() && neuro_expr_file.is_open());
		
		meso_expr_file << header << endl;
		neuro_expr_file << header << endl;
	for(int i = 0; i < size; i++){
		
		
		if( gsl_vector_get(transition_Indices,i) < 10 ) {
				
				meso_seq_file << ">" << seq_name[i] << endl;
				meso_seq_file << seq[i] << endl;
				meso_expr_file << expr[i] << endl;	
		}
		else{
				neuro_seq_file << ">" << seq_name[i] << endl;
				neuro_seq_file << seq[i] << endl;
				neuro_expr_file << expr[i] << endl;
		
		}
	}
		neuro_seq_file.close();
		neuro_expr_file.close();
		
		meso_seq_file.close();
		meso_expr_file.close();
		seq_file.close();
		expr_file.close();


	
	ifstream seq_filem("mesoseq.txt");
	ifstream expr_filem("mesoexpr.txt");
	
	string headerm;
	
	getline(expr_filem, headerm);
	
	string tempm;
	vector <string> seq_namem;
	vector <string> seqm;
	vector <string> exprm;
	
	seqm.clear();
	exprm.clear();
	
	
	while(!seq_filem.eof()){
		tempm = "";
		getline(seq_filem, tempm);
		if(tempm.length() == 0){
			break;
		}
		
		string namem (tempm, 1, tempm.length() - 1);
		seq_namem.push_back(namem);
		
		getline(seq_filem, tempm);
		seqm.push_back(tempm);
				
		getline(expr_filem, tempm);
		exprm.push_back(tempm);
	}
	
	int sizem = seqm.size();
	
	
	for(int i = 0; i < seqm.size(); i++){
		
	}
	
	
	

	vector <int> indicesm(sizem, 0);
	for(int i = 0; i < indicesm.size(); i++)
		indicesm[i] = i;
	std::srand(std::time(0));
	random_shuffle(indicesm.begin(), indicesm.end(), Random);		
	
	
		

	int min_sizem = sizem/k;
	vector <int> partition_sizesm(k, min_sizem);
	
	int residuem = sizem - k * min_sizem;
	for(int i = 0; i < residuem; i++){
		partition_sizesm[i] ++;
	}	
	
	
	int summ = 0;
	for(int i = 0; i < k; i++){
		summ += partition_sizesm[i];
	}
	
	assert (summ == sizem);
	
	
	int startm = 0;
	vector <int> startsm;
	startsm.clear();
	for(int i = 0; i < k; i++){
		startsm.push_back(startm);
		startm += partition_sizesm[i];
	}
	
	for(int i = 0; i < k; i++){
		
	}
	
	for(int i = 0; i < k; i++){
		ostringstream numberm;
		numberm << i;
	
		ofstream train_seq_file (("train_seq_" + numberm.str()).c_str());
		ofstream train_expr_file (("train_expr_" + numberm.str()).c_str());
		
		ofstream test_seq_file(("test_seq_" + numberm.str()).c_str());
		ofstream test_expr_file(("test_expr_" + numberm.str()).c_str());
		
		assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
		
		train_expr_file << headerm << endl;
		test_expr_file << headerm << endl;
		
		for(int j = 0; j < sizem; j++){
			int indexm = indicesm[j];
			if(j >= startsm[i] && j < startsm[i] + partition_sizesm[i]){
				test_seq_file << ">" << seq_namem[indexm] << endl;
				test_seq_file << seqm[indexm] << endl;
				test_expr_file << exprm[indexm] << endl;				
			}
			else{
				train_seq_file << ">" << seq_namem[indexm] << endl;
				train_seq_file << seqm[indexm] << endl;
				train_expr_file << exprm[indexm] << endl;	
			}
		}
		
		
		train_seq_file.close();
		train_expr_file.close();
		
		test_seq_file.close();
		test_expr_file.close();
		
	}


	ifstream seq_filen("neuroseq.txt");
	ifstream expr_filen("neuroexpr.txt");

	string headern;
	
	getline(expr_filen, headern);
	
	string tempn;
	vector <string> seq_namen;
	vector <string> seqn;
	vector <string> exprn;
	
	seqn.clear();
	exprn.clear();
	
	
	while(!seq_filen.eof()){
		tempn = "";
		getline(seq_filen, tempn);
		if(tempn.length() == 0){
			break;
		}
		
		string namen (tempn, 1, tempn.length() - 1);
		seq_namen.push_back(namen);
		
		getline(seq_filen, tempn);
		seqn.push_back(tempn);
				
		getline(expr_filen, tempn);
		exprn.push_back(tempn);
	}
	
	int sizen = seqn.size();
	
	
	for(int i = 0; i < seqn.size(); i++){
		
	}
	
	
	

	vector <int> indicesn(sizen, 0);
	for(int i = 0; i < indicesn.size(); i++)
		indicesn[i] = i;
	std::srand(std::time(0));
	random_shuffle(indicesn.begin(), indicesn.end(), Random);		
	
	
	

	int min_sizen = sizen/k;
	vector <int> partition_sizesn(k, min_sizen);
	
	int residuen = sizen - k * min_sizen;
	for(int i = 0; i < residuen; i++){
		partition_sizesn[i] ++;
	}	
	
	
	int sumn = 0;
	for(int i = 0; i < k; i++){
		sumn += partition_sizesn[i];
	}
	
	assert (sumn == sizen);
	
	
	int startn = 0;
	vector <int> startsn;
	startsn.clear();
	for(int i = 0; i < k; i++){
		startsn.push_back(startn);
		startn += partition_sizesn[i];
	}
	
	for(int i = 0; i < k; i++){
		
	}
	
	for(int i = 0; i < k; i++){
		ostringstream numbern;
		numbern << i;
		
		ofstream train_seq_file (("train_seq_" + numbern.str()).c_str(), ios::app);
		ofstream train_expr_file (("train_expr_" + numbern.str()).c_str(), ios::app);
		
		ofstream test_seq_file(("test_seq_" + numbern.str()).c_str(), ios::app);
		ofstream test_expr_file(("test_expr_" + numbern.str()).c_str(), ios::app);
		
		assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
		
		
		
		
		for(int j = 0; j < sizen; j++){
			int indexn = indicesn[j];
			if(j >= startsn[i] && j < startsn[i] + partition_sizesn[i]){
				test_seq_file << ">" << seq_namen[indexn] << endl;
				test_seq_file << seqn[indexn] << endl;
				test_expr_file << exprn[indexn] << endl;				
			}
			else{
				train_seq_file<< ">" << seq_namen[indexn] << endl;
				train_seq_file << seqn[indexn] << endl;
				train_expr_file << exprn[indexn] << endl;	
			}
		}
		
		
		train_seq_file.close();
		train_expr_file.close();
		
		test_seq_file.close();
		test_expr_file.close();
		
	}
return 1;
}
int ExprPredictor::createccdata2()
{
double expmin = .35;
double expmax = .65;

int k=3;



    
vector< string > rowLabels(0);
 vector< string > colLabels(0);
 vector< vector< double > > data(0) ;

    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
gsl_rng* rng;
gsl_rng_env_setup();
const gsl_rng_type * TT = gsl_rng_default;	
rng = gsl_rng_alloc( TT );
gsl_rng_set( rng, time( 0 ) );		
double  rand_site_index;	
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    
    
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) colLabels.push_back( label );
    
    
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){ 
	 rand_site_index = gsl_ran_flat( rng, -.15, .15 );
	if( val + rand_site_index > 1 ) {
	 vals.push_back( 1 );
	}
	
	if( val + rand_site_index < 0 ) {
	 vals.push_back( 0 );
	}
	
	if( val + rand_site_index >= 0  && val + rand_site_index <= 1 ) {
	vals.push_back( val + rand_site_index );
	}
	} 
        data.push_back( vals );
    }
Matrix data1( data );
Matrix f = factorExprData;
exprData = data;
Matrix e = data;

int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indices = gsl_vector_alloc(nrow);  
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	
	double m;
	m = gsl_vector_max(rowexprData);
	int ti;
	ti = 0;
if( m < expmin ) {
	gsl_vector_set(transition_Indices, i , ti);  
	
	break;
	}
else{
	 
        for (int j=mi; j< ncol; j++)  {  
		
		if (expmax < gsl_vector_get(rowexprData,j) ) {
		
			continue;
	        }
		
		else {

			ti = j; break;
	        }
	}
        int minindex;				
	minindex = gsl_vector_min_index(rowexprData);  
	
     	if (ti == 0) { 
	          ti = minindex;            
	} 	   
  
 	gsl_vector_set(transition_Indices, i , ti);                    
}
}

for (int i=0; i<nrow;i++) {
}
ifstream seq_file("efastanotw10txt");  
	ifstream expr_file( "expre10.tab");
assert (seq_file.is_open() &&  expr_file.is_open());
string header;
	
	getline(expr_file, header);
	
	string temp;
	vector <string> seq_name;
	vector <string> seq;
	vector <string> expr;
	
	seq.clear();
	expr.clear();
	
	
	while(!seq_file.eof()){
		temp = "";
		getline(seq_file, temp);
		if(temp.length() == 0){
			break;
		}
		
		string name1 (temp, 1, temp.length() - 1);
		seq_name.push_back(name1);
		
		getline(seq_file, temp);
		seq.push_back(temp);
				
		getline(expr_file, temp);
		expr.push_back(temp);
	}
	
	int size = seq.size();
		ofstream meso_seq_file ("mesoseq.txt");
		ofstream meso_expr_file ("mesoexpr.txt");
		
		ofstream neuro_seq_file("neuroseq.txt");
		ofstream neuro_expr_file("neuroexpr.txt");
		
		assert (meso_seq_file.is_open() &&  meso_expr_file.is_open() && neuro_seq_file.is_open() && neuro_expr_file.is_open());
		
		meso_expr_file << header << endl;
		neuro_expr_file << header << endl;
	for(int i = 0; i < size; i++){
		
		
		ostringstream numberm;
		numberm << i;
		
		
		if( gsl_vector_get(transition_Indices,i) < 10 ) {
				
				meso_seq_file << ">" << seq_name[i] << endl;
				meso_seq_file << seq[i] << endl;
				meso_expr_file << expr[i] << endl;	
		}
		else{
				neuro_seq_file << ">" << seq_name[i] << endl;
				neuro_seq_file << seq[i] << endl;
				neuro_expr_file << expr[i] << endl;
		
		}
	}
		neuro_seq_file.close();
		neuro_expr_file.close();
		
		meso_seq_file.close();
		meso_expr_file.close();
		seq_file.close();
		expr_file.close();


	
	ifstream seq_filem("mesoseq.txt");
	ifstream expr_filem("mesoexpr.txt");
	
	string headerm;
	
	getline(expr_filem, headerm);
	
	string tempm;
	vector <string> seq_namem;
	vector <string> seqm;
	vector <string> exprm;
	
	seqm.clear();
	exprm.clear();
	
	
	while(!seq_filem.eof()){
		tempm = "";
		getline(seq_filem, tempm);
		if(tempm.length() == 0){
			break;
		}
		
		string namem (tempm, 1, tempm.length() - 1);
		seq_namem.push_back(namem);
		
		getline(seq_filem, tempm);
		seqm.push_back(tempm);
				
		getline(expr_filem, tempm);
		exprm.push_back(tempm);
	}
	
	int sizem = seqm.size();
	
	
	for(int i = 0; i < seqm.size(); i++){
		
	}
	
	
	

	vector <int> indicesm(sizem, 0);
	for(int i = 0; i < indicesm.size(); i++)
		indicesm[i] = i;
	std::srand(std::time(0));
	random_shuffle(indicesm.begin(), indicesm.end(), Random);		
	
	
		

	int min_sizem = sizem/k;
	vector <int> partition_sizesm(k, min_sizem);
	
	int residuem = sizem - k * min_sizem;
	for(int i = 0; i < residuem; i++){
		partition_sizesm[i] ++;
	}	
	
	
	int summ = 0;
	for(int i = 0; i < k; i++){
		summ += partition_sizesm[i];
	}
	
	assert (summ == sizem);
	
	
	int startm = 0;
	vector <int> startsm;
	startsm.clear();
	for(int i = 0; i < k; i++){
		startsm.push_back(startm);
		startm += partition_sizesm[i];
	}
	
	for(int i = 0; i < k; i++){
		
	}
	
	for(int i = 0; i < k; i++){
		ostringstream numberm;
		numberm << i;
		
		ofstream train_seq_file (("train_seq_" + numberm.str()).c_str());
		ofstream train_expr_file (("train_expr_" + numberm.str()).c_str());
		
		ofstream test_seq_file(("test_seq_" + numberm.str()).c_str());
		ofstream test_expr_file(("test_expr_" + numberm.str()).c_str());
		
		assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
		
		train_expr_file << headerm << endl;
		test_expr_file << headerm << endl;
		
		for(int j = 0; j < sizem; j++){
			int indexm = indicesm[j];
			if(j >= startsm[i] && j < startsm[i] + partition_sizesm[i]){
				test_seq_file << ">" << seq_namem[indexm] << endl;
				test_seq_file << seqm[indexm] << endl;
				test_expr_file << exprm[indexm] << endl;				
			}
			else{
				train_seq_file << ">" << seq_namem[indexm] << endl;
				train_seq_file << seqm[indexm] << endl;
				train_expr_file << exprm[indexm] << endl;	
			}
		}
		
		
		train_seq_file.close();
		train_expr_file.close();
		
		test_seq_file.close();
		test_expr_file.close();
		
	}


	ifstream seq_filen("neuroseq.txt");
	ifstream expr_filen("neuroexpr.txt");

	string headern;
	
	getline(expr_filen, headern);
	
	string tempn;
	vector <string> seq_namen;
	vector <string> seqn;
	vector <string> exprn;
	
	seqn.clear();
	exprn.clear();
	
	
	while(!seq_filen.eof()){
		tempn = "";
		getline(seq_filen, tempn);
		if(tempn.length() == 0){
			break;
		}
		
		string namen (tempn, 1, tempn.length() - 1);
		seq_namen.push_back(namen);
		
		getline(seq_filen, tempn);
		seqn.push_back(tempn);
				
		getline(expr_filen, tempn);
		exprn.push_back(tempn);
	}
	
	int sizen = seqn.size();
	
	
	for(int i = 0; i < seqn.size(); i++){
		
	}
	
	
	

	vector <int> indicesn(sizen, 0);
	for(int i = 0; i < indicesn.size(); i++)
		indicesn[i] = i;
	std::srand(std::time(0));
	random_shuffle(indicesn.begin(), indicesn.end(), Random);		
	
	
		

	int min_sizen = sizen/k;
	vector <int> partition_sizesn(k, min_sizen);
	
	int residuen = sizen - k * min_sizen;
	for(int i = 0; i < residuen; i++){
		partition_sizesn[i] ++;
	}	
	
	
	int sumn = 0;
	for(int i = 0; i < k; i++){
		sumn += partition_sizesn[i];
	}
	
	assert (sumn == sizen);
	
	
	int startn = 0;
	vector <int> startsn;
	startsn.clear();
	for(int i = 0; i < k; i++){
		startsn.push_back(startn);
		startn += partition_sizesn[i];
	}
	
	for(int i = 0; i < k; i++){
		
	}
	
	for(int i = 0; i < k; i++){
		ostringstream numbern;
		numbern << i;
		
		ofstream train_seq_file (("train_seq_" + numbern.str()).c_str(), ios::app);
		ofstream train_expr_file (("train_expr_" + numbern.str()).c_str(), ios::app);
		
		ofstream test_seq_file(("test_seq_" + numbern.str()).c_str(), ios::app);
		ofstream test_expr_file(("test_expr_" + numbern.str()).c_str(), ios::app);
		
		assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
		
		
		
		
		for(int j = 0; j < sizen; j++){
			int indexn = indicesn[j];
			if(j >= startsn[i] && j < startsn[i] + partition_sizesn[i]){
				test_seq_file << ">" << seq_namen[indexn] << endl;
				test_seq_file << seqn[indexn] << endl;
				test_expr_file << exprn[indexn] << endl;				
			}
			else{
				train_seq_file << ">" << seq_namen[indexn] << endl;
				train_seq_file << seqn[indexn] << endl;
				train_expr_file << exprn[indexn] << endl;	
			}
		}
		
		
		train_seq_file.close();
		train_expr_file.close();
		
		test_seq_file.close();
		test_expr_file.close();
		
	}
return 1;
}
void ExprPredictor::compvar(vector< double >& vars) 
{
	vector< double > predictedExprs;
        vector< double > stepobservedExprs;

    
    ExprFunc* func = createExprFunc( par_model );
            double rss = 0;
    
double step = .001;
    vector< double > concs = factorExprData.getCol( 0 );

vector< double > concs2 = concs;
concs2[0] = concs[0] - step;
    for ( int i = 0; i < nSeqs(); i++ ) {
        
       
         
            
            
            double predicted = func->predictExpr( seqSites[ i ], 5, concs );
            predictedExprs.push_back( predicted );
           
            
            double stepobserved = func->predictExpr(seqSites[ i ], 5, concs2 );       
            stepobservedExprs.push_back( stepobserved );
        
     }
double partialcu = .5 * ( concs2[0] + concs[0] );
assert( predictedExprs.size() == stepobservedExprs.size());
    for ( int i = 0; i <  predictedExprs.size(); i++ ) {
 
 
     vars.push_back(  partialcu * (predictedExprs[i]  -  stepobservedExprs[i]  ) / step )   ; 
    }
}
double ExprPredictor::compAvgCrossCorr( const ExprPar& par ) const
{
    
    ExprFunc* func = createExprFunc( par );
            
    
    double totalSim = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            
            
            double predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            
            
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        totalSim += exprSimCrossCorr( predictedExprs, observedExprs ); 
    }	

    return totalSim / nSeqs();
}


void ExprPredictor::compOccMat(const gsl_vector* v, void* params  ) 
{
int nrow = exprData.nRows();    
int ncol =exprData.nCols();
gsl_vector *transition_Indices = gsl_vector_alloc(nrow);  
for (int i=0; i<nrow;i++) {
	vector< double > reD;					 
	reD = exprData.getRow(i);	
	int mi;				
	gsl_vector *rowexprData = gsl_vector_alloc(ncol);
	rowexprData = vector2gsl(reD);			 
	mi = gsl_vector_max_index(rowexprData);  
	double m;
	m = gsl_vector_max(rowexprData);
	int ti;
	ti = 0;
        for (int j=mi; j< ncol-mi; j++)  {  
	        double temp;
		temp = m/2;     
	        if (temp < gsl_vector_get(rowexprData,j) ) {
			continue;
	        }
		
		else {
			ti = j; break;
	        }
	}
        int minindex;				
	minindex = gsl_vector_min_index(rowexprData);  
	
     	if (ti == 0) { 
	          ti = minindex;            
	} 	   
  
 	gsl_vector_set(transition_Indices, i , ti);                    
}     

ExprFunc* funcmat = createExprFunc2( par_model );
gsl_matrix * Occupancy = gsl_matrix_alloc(nrow,par_model.nFactors() );  
for ( int i = 0; i < nSeqs(); i++ ) {
                                                                  
        vector< double > concs = factorExprData.getCol( gsl_vector_get(transition_Indices,i) );
        gsl_vector * Occvector =gsl_vector_alloc(motifs.size()); 
        vector <double > fOcc(3);
        funcmat->predictOcc( seqSites[ i ], seqLengths[i], concs, fOcc); 
        Occvector = vector2gsl(fOcc);  
        gsl_matrix_set_row(Occupancy , i , Occvector);           
}


int fixedsize = fix_pars.size();  
int rows = Occupancy -> size1;
int columns = Occupancy -> size2;
int i,j,k;
k=0;
gsl_matrix *X = gsl_matrix_alloc(columns,columns);
gsl_matrix *V = gsl_matrix_alloc(columns,columns);
gsl_vector *S =gsl_vector_alloc(columns);
gsl_vector *xx =gsl_vector_alloc(columns);  
gsl_vector *b =gsl_vector_alloc(rows);     
gsl_vector_set_all( b,0 );
gsl_vector *work=gsl_vector_alloc(columns);
int rsvd;		
int rsvds;               
rsvd=gsl_linalg_SV_decomp(Occupancy,V,S,work);
	if(gsl_vector_get(S,2) < .1 ) {
	
	
			
			
			
		}
	
	
for(int i =0; i < columns;i++){  
			
		}
gsl_matrix_free( Occupancy );
gsl_matrix_free( X );
gsl_matrix_free( V );
gsl_vector_free( xx );
gsl_vector_free( b );
gsl_vector_free( S );
   
}

   

int ExprPredictor::simplex_minimize( ExprPar& par_result, double& obj_result )   
{

 

    vector< double > pars;
  

    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators ); 
       
int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			free_pars.push_back( pars[ index ]);
		}
		else{
			fix_pars.push_back( pars[ index ] );
		}
	}
       
	
	
	
	pars.clear();
	pars = free_pars;
	
for (int i = 0; i < pars.size() ; i++ ) {
}
    gsl_multimin_function my_func;
    my_func.f = &gsl_obj_f;
    my_func.n = pars.size();
    my_func.params = (void*)this;    
    
    gsl_vector* x = vector2gsl( pars );

   const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex;
    gsl_vector* ss = gsl_vector_alloc( my_func.n );
    gsl_vector_set_all( ss, 1.0 );

    
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc( T, my_func.n );

    gsl_multimin_fminimizer_set( s, &my_func, x, ss ); 

    
    size_t iter = 0;
    int status;
    double size;	
    do {
        double f_prev = iter ? s->fval : 1.0E6;     
   
	iter++;

        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_fminimizer_iterate( s );
   
        
        if ( status ) {  break;}
	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}
        ExprPar par_curr = ExprPar ( pars, coopMat, actIndicators, repIndicators );
	
	  
      
   
        
        
         double f_curr = s->fval;            
         double delta_f = abs( f_curr - f_prev ); 
         if ( objOption == SSE && delta_f < min_delta_f_SSE ) {  break;}
         if ( objOption == CORR && delta_f < min_delta_f_Corr ) { break;}
        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_test_size( size,1e-6 ); 
 
 
 		if ( status == GSL_SUCCESS ) {  }

        
     
  
   
    } while (  status == GSL_CONTINUE && iter < nSimplexIters );
 
free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}
 par_result = ExprPar( pars, coopMat, actIndicators, repIndicators );
    obj_result = s->fval;

    return 0;
}

int ExprPredictor::gradient_minimize( ExprPar& par_result, double& obj_result )   
{
	
    
    vector< double > pars;
    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators ); 
       
	int pars_size = pars.size();
	fix_pars.clear();
	free_pars.clear();
	for( int index = 0; index < pars_size; index++ ){
		if( indicator_bool[ index ] ){
			
			free_pars.push_back( pars[ index ]);
		}
		else{
			
			fix_pars.push_back( pars[ index ] );
		}
	}
  
	
	pars.clear();
	pars = free_pars;
	
    
    gsl_multimin_function_fdf my_func;
    my_func.f = &gsl_obj_f;
    my_func.df = &gsl_obj_df;
    my_func.fdf = &gsl_obj_fdf;
    my_func.n = pars.size();
    my_func.params = (void*)this;
    
    
    gsl_vector* x = vector2gsl( pars ); 

	
 	
		
	
 	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_pr;	
		
    
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc( T, my_func.n );
    double init_step = .001, tol = .0001;
    gsl_multimin_fdfminimizer_set( s, &my_func, x, init_step, tol );
    
    
    size_t iter = 0;
    int status;
    do {
        double f_prev = iter ? s->f : 1.0E6;     
       
        iter++;
        status = gsl_multimin_fdfminimizer_iterate( s );
 
	if (status ) {

 break; }
   
	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}
 ExprPar par_curr = ExprPar ( pars, coopMat, actIndicators, repIndicators );
  
        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) {   printPar( par_curr );break;}
        double f_curr = s->f;
        double delta_f = abs( f_curr - f_prev ); 
        if ( objOption == SSE && delta_f < min_delta_f_SSE ) { break;}
        if ( objOption == CORR && delta_f < min_delta_f_Corr ) {  break;}
        if ( objOption == CROSS_CORR && delta_f < min_delta_f_CrossCorr ) break;
  
        status = gsl_multimin_test_gradient( s->gradient, 5e-4 ); 
    
      

    } while ( status == GSL_CONTINUE && iter < nGradientIters );

	free_pars = gsl2vector( s-> x);
	pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < pars_size; index ++ ){
		if( indicator_bool[ index ] ){
			pars.push_back( free_pars[ free_par_counter ++ ]);
		}
		else{
			pars.push_back( fix_pars[ fix_par_counter ++ ]);
		}
	}
        par_result = ExprPar ( pars, coopMat, actIndicators, repIndicators );
	
    obj_result = s->f;


  
    gsl_vector_free( x );    
    gsl_multimin_fdfminimizer_free( s );
    
    return 0;
}
double gsl_obj_f( const gsl_vector* v, void* params )
{ 

    
    ExprPredictor* predictor = (ExprPredictor*)params;
            
    
	vector <double> temp_free_pars = gsl2vector(v);
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
		if( predictor ->  indicator_bool[ index ]  ){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}


	
    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
    
    if( (*predictor).clasvar == 0 ) {
double tempobj2 = predictor->objFuncborder( par ); return tempobj2;
}
else { 
double obj = predictor->objFunc2( par );  return obj;	
}

    
double tempobj2 = predictor->objFuncborder( par ); 

	
 
	
	return tempobj2;
}
double gsl_obj_f2( const gsl_vector* v, void* params )
{ 

    
    ExprPredictor* predictor = (ExprPredictor*)params;
            
    
	vector <double> temp_free_pars = gsl2vector(v);
	vector < double > all_pars;
	all_pars.clear();
	int free_par_counter = 0;
	int fix_par_counter = 0;
	for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
		if( predictor ->  indicator_bool[ index ]  ){
			all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
		}
		else{
			all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
		}
	}


	
    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
    
    
    
if( (*predictor).clasvar == 1 ) {
double tempobj2 = predictor->objFuncborder2( par ); return tempobj2;
}
else { 
double obj = predictor->objFunc( par );  return obj;	
}
	
 
	
	
}
void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad )
{
    double step = 1.0E-3;
    numeric_deriv( grad, gsl_obj_f, v, params, step );	
}

void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
    *result = gsl_obj_f( v, params ); 
    gsl_obj_df( v, params, grad );		
}


