
#include "ExprPredictor.h"


int main( int argc, char* argv[] ) 
{
 
    // command line processing
    string seqFile, annFile,duFile,dcFile, twFile, adamiFile, exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, seqFileb, bIFile;
    string outFile, cbFileo, motifFile2, adamiFileFP;     // output file
    double coopDistThr = 150;
    double factorIntSigma = 50.0;   // sigma parameter for the Gaussian interaction function
    double repressionDistThr = 150;
    double energyThr = 2;
    int maxContact = 1;
	double pva=.05;
int fi =0;
double psed=.1; // psedocount
bool free_fix_indicators[] = {1,1,1,
			    1,1,0,0,0,
///*
0,0,0,0,0,
0,0,0,0,0,
1,0,0 };
int nrandStarts =0;
	int nUpdates =7;
	int wid = 2;
	int l = 1;

	int nExps = 1;	// number of experiments
	vector <bool> indicator_bool ( free_fix_indicators, free_fix_indicators + sizeof( free_fix_indicators )/sizeof( bool ));

    int binwidth; // = 30;
      ExprPredictor::nAlternations = 1;
     for ( int i = 1; i < argc; i++ ) {
	
        if ( !strcmp( "-s", argv[ i ] ) )
            seqFile = argv[ ++i ];
	else if ( !strcmp( "-du", argv[ i ] ) )
            duFile = argv[ ++i ];  
	else if ( !strcmp( "-cbo", argv[ i ] ) )
            cbFileo = argv[ ++i ];  
	else if ( !strcmp( "-dc", argv[ i ] ) )
            dcFile = argv[ ++i ]; 
	else if ( !strcmp( "-tw", argv[ i ] ) )
            twFile = argv[ ++i ];   
	else if ( !strcmp( "-sa", argv[ i ] ) )
            adamiFile = argv[ ++i ];    
	else if ( !strcmp( "-saFP", argv[ i ] ) )
            adamiFileFP = argv[ ++i ];  
        else if ( !strcmp( "-a", argv[ i ] ) )
            annFile = argv[ ++i ];      
        else if ( !strcmp( "-e", argv[ i ] ) )
            exprFile = argv[ ++i ];            
        else if ( !strcmp( "-m", argv[ i ] ) )
            motifFile = argv[ ++i ];
	   else if ( !strcmp( "-m2", argv[ i ] ) )
            motifFile2 = argv[ ++i ];
        else if ( !strcmp( "-f", argv[ i ] ) )
            factorExprFile = argv[ ++i ];    
        else if ( !strcmp( "-o", argv[ i ] ) )
            ExprPredictor::modelOption = getModelOption( argv[++i] );
        else if ( !strcmp( "-c", argv[ i ] ) )
            coopFile = argv[ ++i ];
        else if ( !strcmp( "-i", argv[ i ] ) )
            factorInfoFile = argv[ ++i ];            
        else if ( !strcmp( "-r", argv[ i ] ) )
            repressionFile = argv[ ++i ];  
        else if ( !strcmp( "-oo", argv[ i ] ) )
            ExprPredictor::objOption = getObjOption( argv[++i] );    
        else if ( !strcmp( "-mc", argv[i] ) )
            maxContact = atoi( argv[++i] );
        else if ( !strcmp( "-fo", argv[i] ) )
            outFile = argv[++i];    
        else if ( !strcmp( "-p", argv[i] ) )
            parFile = argv[++i]; 
        else if ( !strcmp( "-ct", argv[i] ) )
            coopDistThr = atof( argv[++i] ); 
else if ( !strcmp( "-binwt", argv[i] ) )
            binwidth = atof( argv[++i] ); 
        else if ( !strcmp( "-sigma", argv[i] ) )
            factorIntSigma = atof( argv[++i] );    
        else if ( !strcmp( "-rt", argv[i] ) )
            repressionDistThr = atof( argv[++i] );    
        else if ( !strcmp( "-et", argv[i] ) )
            energyThr = atof( argv[++i] );
        else if ( !strcmp( "-na", argv[i] ) )
            ExprPredictor::nAlternations = atoi( argv[++i] );  
  	else if ( !strcmp( "-n", argv[ i ] ) )
	    nExps = atoi( argv[ ++i ] );
	else if ( !strcmp( "-bs", argv[i] ) )
            seqFileb = argv[++i];    
        else if ( !strcmp( "-bi", argv[i] ) )
            bIFile = argv[++i]; 
		else if ( !strcmp( "-pva", argv[ i ] ) )
	    pva = atof( argv[ ++i ] );
	else if ( !strcmp( "-nu", argv[i] ) )
            nUpdates = atoi( argv[++i] );
else if ( !strcmp( "-nr", argv[i] ) )
            nrandStarts = atoi( argv[++i] );
	else if ( !strcmp( "-fi", argv[ i ] ) )
	    fi = atoi( argv[ ++i ] );
        else if ( !strcmp( "-w", argv[i] ) )
            wid = atoi( argv[++i] );
    
  	else if ( !strcmp( "-l", argv[ i ] ) )
	    l = atoi( argv[ ++i ] );
	else if ( !strcmp( "-pse", argv[ i ] ) )
	    psed = atof( argv[ ++i ] );
	
    }
    if ( seqFile.empty() || exprFile.empty() || motifFile.empty() || factorExprFile.empty() || outFile.empty() || ( ( ExprPredictor::modelOption == QUENCHING || ExprPredictor::modelOption == CHRMOD_UNLIMITED || ExprPredictor::modelOption == CHRMOD_LIMITED ) &&  factorInfoFile.empty() ) || ( ExprPredictor::modelOption == QUENCHING && repressionFile.empty() ) ) {
        cerr << "Usage: " << argv[ 0 ] << " -s seqFile -e exprFile -m motifFile -f factorExprFile -fo outFile [-a annFile -o modelOption -c coopFile -i factorInfoFile -r repressionFile -oo objOption -mc maxContact -p parFile -rt repressionDistThr -et energyThr -na nAlternations -ct coopDistThr -sigma factorIntSigma]" << endl;
        cerr << "modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited" << endl;
        exit( 1 );
    }  

vector< int > mmm(4);
mmm[0] = 1;   // 5
mmm[1] = 5;   //25
mmm[2] = 40;  //50000
mmm[3] = 50001;

vector< int > mmmr(4);
mmmr[0] = 10;
mmmr[1] = 60;
mmmr[2] = 70;
mmmr[3] = 100;

    double gcContent = .5;
    FactorIntType intOption = BINSF;     // type of interaction function
    ExprPar::searchOption = CONSTRAINED;      // search option: unconstrained; constrained. 
    ExprPar::estBindingOption = 1;
    int nbin = mmm.size() + 1;

    ExprPar::nbins = nbin;    // necessary in for loops for theV
   
    ExprPredictor::nRandStarts = 0;
    ExprPredictor::min_delta_f_SSE = 1.0E-10;
    ExprPredictor::min_delta_f_Corr = 1.0E-10;
    ExprPredictor::min_delta_f_CrossCorr = 1.0E-10;
    ExprPredictor::nSimplexIters = 50;
    ExprPredictor::nGradientIters = 50;

    int rval;
    vector< vector< double > > data;    // buffer for reading matrix data
    vector< string > labels;    // buffer for reading the labels of matrix data
    string factor1, factor2;    // buffer for reading factor pairs

    vector< Sequence > seqs;
    vector< string > seqNames;
    rval = readSequences( seqFile, seqs, seqNames );
    assert( rval != RET_ERROR );
    int nSeqs = seqs.size();
//////////////////////////////////////////////////7292011

rval = readSequences( adamiFileFP, ExprPredictor::seqsyaFP, ExprPredictor::seqNmesaFP ); 
rval = readSequences( seqFile, ExprPredictor::seqsy, ExprPredictor::seqNmes );

  // read the sequences
    vector< Sequence > seqsa;
    vector< string > seqNamesa;
    rval = readSequences( adamiFile, seqsa, seqNamesa );
    assert( rval != RET_ERROR );
    int nSeqsa = seqsa.size();

rval = readSequences( adamiFile, ExprPredictor::seqsya, ExprPredictor::seqNmesa );

    // read the expression data
    vector< string > condNames;  
    rval = readMatrix( exprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nSeqs );
    for ( int i = 0; i < nSeqs; i++ ) assert( labels[i] == seqNames[i] );
    Matrix exprData( data ); 
    int nConds = exprData.nCols();
    ExprPredictor::exprData2 = exprData;

    // read the motifs

    vector< Motif > motifs;
    vector< Motif > motifs2;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );


vector< double > energyThrs2( 4, 0 ); 
string b = "dcondit.txt";
string duc = "dconditdu.txt";
string dcs = "dcondits.txt";
////////////////////////////////////////////////////////////////////////////////////
//    DC
////////////////////////////////////////////////////////////////////////////////////
int ncountsdc =10;

Matrix dcm =countmatrixS( dcFile );

Matrix dcm2 = countmatrix2( dcFile , b );

Motif dcmm( dcm2,psed, background);
double pv =1.5;

string aaadc ="dcondit2.txt";
string aaapdc="dccdfcore.txt";
dcmm.countmatrix39( dcFile, aaadc, aaapdc );

energyThrs2[0]=pv;
motifs2.push_back( dcmm ); 

/////////////////////////////////////////////////////////////
//  DU
//////////////////////////////////////////////////////////////////
string dcu="duconditdc.txt";
string a = "ducondit.txt";
string aa ="ducondit2.txt";
string aap = "ducdfcore.txt";

string ducs = "ducondits.txt";
Matrix dum2 = countmatrix2( duFile, a );

Motif dumm( dum2,psed, background);

dumm.countmatrix39( duFile, aa,  aap );

energyThrs2[1]= pv; //2.6
motifs2.push_back( dumm );     //  this is 2 is seqSitesa
scorematrix( motifs2[1].getLLRMat(), dcFile , motifs2[1].getMaxLLR(), dcu );
scorematrix( motifs2[1].getLLRMat(), duFile , motifs2[1].getMaxLLR(), ducs );
scorematrix( motifs2[0].getLLRMat(), duFile , motifs2[0].getMaxLLR(), duc );

///////////////////////////////////////////////////////////////////////////
//  TW
////////////////////////////////////////////////////////////////////////

string c = "twcondit.txt";

Matrix twm = countmatrixtw( twFile );
Motif twmm( twm, background);
double ethresT = 1;
twmm.setEth( ethresT );

energyThrs2[2]=0;
motifs2.push_back( twmm );   
//////////////////////////////////////////////////////////////////////////////
//    CB
////////////////////////////////////////////////////////////////////////////////////
Matrix cbom = countmatrixS( cbFileo );
Motif cbomf( cbom,psed, background);
string cbfi="cbcondit.txt";
string cbfis="cbcondits.txt";
string aaa = "cbcondit2.txt";
string aaap = "cbcdfcore.txt";

Matrix cbfim = countmatrix2( cbFileo , cbfi );

cbomf.countmatrix39( cbFileo, aaa, aaap );

energyThrs2[3]= pv;

motifs2.push_back( cbomf );     //  this is 2 is seqSitesa
scorematrix( motifs2[3].getLLRMat(), cbFileo , motifs2[3].getMaxLLR(), cbfis );

string cbfisz="cbits.txt";
string deric="deric.txt";
string kelly = "kelly.txt";
double len =0;

scorematrix(  motifs2[0].getLLRMatpos()*(-1), dcFile , len, deric );
scorematrix(  motifs2[1].getLLRMatpos()*(-1), duFile , len, kelly );



string kellyfp2="kellyfptw.txt";
string dericfp2="dericfptw.txt";
scorematrix( motifs2[0].getLLRMatpos()*(-1), duFile , len, dericfp2 );
scorematrix( motifs2[1].getLLRMatpos()*(-1), dcFile , len, kellyfp2 );
rval = readMotifs( motifFile, background, motifs, motifNames ); 

    vector< Motif > motifs3;
    vector< string > motifNames3;
    rval = readMotifs( motifFile2, background, motifs3, motifNames3 ); 
    assert( rval != RET_ERROR );
    int nFactors = motifs.size();
	Motif d;
	d.copy(motifs[0]);
	Matrix dtemp=d.getPwm2();
	vector< double > e(4);
	e[0] = 1;
	e[1] = 1;
	e[2] = 1;
	e[3] = 10;
	vector< vector< double > > d2s;
	for(int i = 0 ; i < dtemp.nRows(); i++ ) {

		d2s.push_back( dtemp.getRow(i) ) ;

		if ( i == 4 ){
		d2s.push_back( e );

		}
	}
	Matrix d2sm( d2s );

	Motif d2( d2sm , background) ;

    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < motifNames.size(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
     
    // read the factor expression data
    labels.clear();
    data.clear();
    rval = readMatrix( factorExprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nFactors && condNames.size() == nConds );
    for ( int i = 0; i < nFactors; i++ ) assert( labels[i] == motifNames[i] );
    Matrix factorExprData( data ); 
    assert( factorExprData.nCols() == nConds ); 
    //////////////////////////////////////////////73011
ExprPredictor::factorExprData2 = factorExprData;
//////////////////////////////////////////////////////////////
    // site representation of the sequences
    vector< double > energyThrs( nFactors, energyThr ); 

    vector< SiteVec > seqSites( nSeqs );
    vector< int > seqLengths( nSeqs );
    SeqAnnotator ann( motifs, energyThrs );
    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.annot( seqs[ i ], seqSites[ i ] );
            seqLengths[i] = seqs[i].size();
        }
    } else {    // read the site representation and compute the energy of sites
        rval = readSites( annFile, factorIdxMap, seqSites, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.compEnergy( seqs[i], seqSites[i] );
            seqLengths[i] = seqs[i].size();
        }
    }

vector< SiteVec > seqSitesa( nSeqsa );
vector< int > seqLengthsa( nSeqsa );
  
SeqAnnotator anndt( motifs2, energyThrs2 );

    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < nSeqsa; i++ ) {

            anndt.annot( seqsa[ i ], seqSitesa[ i ] );
            seqLengthsa[i] = seqsa[i].size();
        }
    } else {    // read the site representation and compute the energy of sites
        rval = readSites( annFile, factorIdxMap, seqSitesa, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqsa; i++ ) {
            anndt.compEnergy( seqsa[i], seqSitesa[i] );
            seqLengthsa[i] = seqsa[i].size();
        }
    }

vector< int > en( motifs2[0].length() + 100,0 );
vector< int > en2( motifs2[1].length() + 100,0 );

SeqAnnotator annbb( motifs2, energyThrs2 );
vector< Motif > motifsexu;
vector< string > motifNamesexu;
countmatrix_unique( duFile, motifNamesexu,  motifsexu,  background);
		map< string, int > factorIdxMapexu;
	    for ( int i = 0; i < motifNamesexu.size(); i++ ) {
		factorIdxMapexu[motifNamesexu[i]] = i;
	    }
	     
	  map< int, string > factorIdxMap2u;
	    for ( int i = 0; i < motifNamesexu.size(); i++ ) {
		factorIdxMap2u[i] = motifNamesexu[i];
	    }

vector< double > energyThrsexu( motifsexu.size(), 0.001 ); 
vector< SiteVec > seqSitesaexu( nSeqsa );
vector< int > seqLengthsaexu( nSeqsa );
SeqAnnotator anndtexu( motifsexu, energyThrsexu );


    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < nSeqsa; i++ ) {
            anndtexu.annot( seqsa[ i ], seqSitesaexu[ i ] );
            seqLengthsaexu[i] = seqsa[i].size();
        }
    } else {    // read the site representation and compute the energy of sites
        rval = readSites( annFile, factorIdxMapexu, seqSitesaexu, false );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqsa; i++ ) {
            anndtexu.compEnergy( seqsa[i], seqSitesaexu[i] );
            seqLengthsaexu[i] = seqsa[i].size();
        }
    }
string  cse= "outsitesdu.txt";
ofstream cs( cse.c_str() );
for ( int i = 0; i < nSeqsa; i++ ) {
	cs << ">" + seqNamesa[i] << endl;
	for ( int j = 0; j < seqSitesaexu[i].size(); j++ ) {

cs << seqSitesaexu[i][j] << '\t' <<  factorIdxMap2u[seqSitesaexu[i][j].factorIdx ] << endl;
	}
}
cs.close();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// load dc
    vector< Motif > motifsex;
    vector< string > motifNamesex;

countmatrix_unique( dcFile, motifNamesex,  motifsex,  background);



    map< string, int > factorIdxMapex;
    for ( int i = 0; i < motifNamesex.size(); i++ ) {
        factorIdxMapex[motifNamesex[i]] = i;
    }
     
  map< int, string > factorIdxMap2;
    for ( int i = 0; i < motifNamesex.size(); i++ ) {
        factorIdxMap2[i] = motifNamesex[i];
    }
vector< double > energyThrsex( motifsex.size(), 0.1 ); 
vector< SiteVec > seqSitesaex( nSeqsa );
    vector< int > seqLengthsaex( nSeqsa );
   
 SeqAnnotator anndtex( motifsex, energyThrsex );

    if ( annFile.empty() ) {        // construct site representation
        for ( int i = 0; i < nSeqsa; i++ ) {

            anndtex.annot( seqsa[ i ], seqSitesaex[ i ] );
            seqLengthsaex[i] = seqsa[i].size();
        }
    } else {    // read the site representation and compute the energy of sites
        rval = readSites( annFile, factorIdxMapex, seqSitesaex, false );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqsa; i++ ) {
            anndtex.compEnergy( seqsa[i], seqSitesaex[i] );
            seqLengthsaex[i] = seqsa[i].size();
        }
    }
string  csed= "outsitesdc.txt";
ofstream csd( csed.c_str() );
for ( int i = 0; i < nSeqsa; i++ ) {
	csd << ">" + seqNamesa[i] << endl;
	for ( int j = 0; j < seqSitesaex[i].size(); j++ ) {

csd << seqSitesaex[i][j].factorIdx  << '\t' <<  factorIdxMap2[seqSitesaex[i][j].factorIdx ] << endl;
	}
}
csd.close();

for ( int i = 0; i < nSeqsa; i++ ) {
	
	for ( int j = 0; j < seqSitesaex[i].size(); j++ ) {

seqSitesaex[i][j].factorIdx = 0;

	}
}

for ( int i = 0; i < nSeqsa; i++ ) {
	
	for ( int j = 0; j < seqSitesaexu[i].size(); j++ ) {

seqSitesaexu[i][j].factorIdx = 1;

	}
}


//////////////////////////////////change sitespp initialization to nothing,

vector< vector< vector< Site > > > sitespp(nSeqsa, vector< vector< Site > >(4));
vector< vector< vector< Site > > > sitesptex(nSeqsa, vector< vector< Site > >(4));
for( int e = 0; e < nSeqsa; e++ ) {
	sitesptex[e][0] = seqSitesaex[e];    // this was filled with core DC
	sitesptex[e][1] = seqSitesaexu[e];      // this was filled with core DU
}
/*
for( int e = 0; e < seq.size(); e++ ) {
sitestoverlap( sitesptex[e][0], sitespp[e][0] );  
sitestoverlap( sitesptex[e][1], sitespp[e][1] );
}*/
for ( int e = 0; e < nSeqsa; e++ ) {  // runs over enhnacers
	  for ( int k = 0; k < 3; k++ ) {  // 
		if( k == 0 ) {
				//for( int i=0; i < seqSitesbot[e].size() ; i++ ){
		
					//if( sitespt[e][k].size() > 0 ) {
				//sitesptex[e][k].push_back(  siteMax(  sitesptex[e][0]  )); 
				//}	
			while( !sitesptex[e][k].empty()){
				if( sitesptex[e][k].size() > 0 ) {
				sitespp[e][k].push_back(  anndt.siteMax(   sitesptex[e][0]   ));  // sitespp[e][k] = sitespt[e][k]; siteMax(  sitespt[e][k] 
				
				anndt.sitestoverlap( sitesptex[e][k], sitespp[e][k] );  
				//cout << " sitespp e k " << sitespp[e][k] << endl;
				}//if site
			}// while
		
			//if( k==1 ){
			//cout << " sitespp e k " << sitespp[e][k] << endl;
			//}
		}//if k
		if( k ==1 ) {
			while( !sitesptex[e][k].empty()){
				if( sitesptex[e][k].size() > 0 ) {
				sitespp[e][k].push_back(  anndt.siteMax(   sitesptex[e][1]   ));  
				
				anndt.sitestoverlap( sitesptex[e][k], sitespp[e][k] );  
				
				}//if site
			}// while
		} // if k==1
	  }//for k
	}// for e
 vector< SiteVec > seqSitesbot;
 vector< SiteVec > seqSitesm1;
for ( int e = 0; e < nSeqsa; e++ ) {  // runs over enhnacers
	  for ( int k = 0; k < 3; k++ ) {  // 
		if( k == 0 ) {
				
				seqSitesbot.push_back( sitespp[e][k] ); 
		}//if k
		if( k ==1 ) {
			seqSitesm1.push_back( sitespp[e][k] );
		} // if k==1
	  }//for k
	}// for e

string  csed2= "outsitesdc2.txt";
ofstream csd2( csed2.c_str() );
for ( int i = 0; i < nSeqsa; i++ ) {
	csd2 << ">" + seqNamesa[i] << endl;
	for ( int j = 0; j < seqSitesbot[i].size(); j++ ) {

csd2 << seqSitesbot[i][j]  << endl;
	}
}
csd2.close();

string  cse2= "outsitesdu2.txt";
ofstream cs2( cse2.c_str() );
for ( int i = 0; i < nSeqsa; i++ ) {
	cs2 << ">" + seqNamesa[i] << endl;
	for ( int j = 0; j < seqSitesm1[i].size(); j++ ) {

cs2 << seqSitesm1[i][j] << endl;
	}
}
cs2.close();
 
 vector< SiteVec > seqSitesm2( nSeqs );
seqSitesm2 = seqSites;

//cout << " asd" <<endl;
vector< SiteVec > seqSitesf2( nSeqs );
seqSitesf2 = seqSites;
vector< SiteVec > seqSitesbotf2(nSeqs);
seqSitesbotf2 = seqSites;
vector< SiteVec > seqSitesm1f2(nSeqs);
seqSitesm1f2 = seqSites; 
vector< SiteVec > seqSitesm2f2( nSeqs);
seqSitesm2f2 = seqSites;

vector< SiteVec >  seqSitesf3( nSeqs);
seqSitesf3 = seqSites;
vector< SiteVec >  seqSitesbotf3( nSeqs);
seqSitesbotf3 = seqSites;
vector< SiteVec > seqSitesm1f3( nSeqs );
seqSitesm1f3 = seqSites;
vector< SiteVec >  seqSitesm2f3( nSeqs);
seqSitesm2f3=seqSites;

vector< vector< Site > > seqSitesm1delete1;
seqSitesm1delete1 = seqSites;
////////////////////////////////////// 7/8/12
 vector< vector< vector< Site > > > dpp(nSeqsa, vector< vector< Site > >( 4 ));
///////////////////////////////////////////////////////////////////////////////////
SeqAnnotator anny( motifs, energyThrs , d2);
   
    // read the cooperativity matrix 
    IntMatrix coopMat( nFactors, nFactors, false );
    if ( !coopFile.empty() ) {
        ifstream fcoop( coopFile.c_str() );
        if ( !fcoop ) {
            cerr << "Cannot open the cooperativity file " << coopFile << endl;
            exit( 1 );
        }  
        while ( fcoop >> factor1 >> factor2 ) {
            assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
            int idx1 = factorIdxMap[factor1];
            int idx2 = factorIdxMap[factor2];
            coopMat( idx1, idx2 ) = true;
            coopMat( idx2, idx1 ) = true;
        }        
    } 

    // read the roles of factors
    vector< bool > actIndicators( nFactors, false );
    vector< bool > repIndicators( nFactors, false );
    if ( !factorInfoFile.empty() ) {
        ifstream finfo( factorInfoFile.c_str() );
        if ( !finfo ) {
            cerr << "Cannot open the factor information file " << factorInfoFile << endl;
            exit( 1 );
        }      
        string name;
        int i = 0, actRole, repRole;
        while ( finfo >> name >> actRole >> repRole ) {
            assert( name == motifNames[i] );
            if ( actRole ) actIndicators[i] = true;
            if ( repRole ) repIndicators[i] = true;
            i++;
        }
    }
    
    // read the repression matrix 
    IntMatrix repressionMat( nFactors, nFactors, false );
    if ( !repressionFile.empty() ) {
        ifstream frepr( repressionFile.c_str() );
        if ( !frepr ) {
            cerr << "Cannot open the repression file " << repressionFile << endl;
            exit( 1 );
        }        
        while ( frepr >> factor1 >> factor2 ) {
            assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
            int idx1 = factorIdxMap[factor1];
            int idx2 = factorIdxMap[factor2];
            repressionMat( idx1, idx2 ) = true;
        }        
    }


	// read the sequences
	vector< Sequence > allSeqs;
	vector< string > allNames;
	rval = readSequences( seqFileb, allSeqs, allNames );
    assert( rval != RET_ERROR );

	int nSeqsb = allSeqs.size() / nExps;
	vector< vector< Sequence > > seqs2( nExps );
	vector< vector< string > > names( nExps );	
	int counter = 0;
	for ( int i = 0; i < nExps; i++ ) {
		for ( int j = 0; j < nSeqsb; j++ ) {
			seqs2[ i ].push_back( allSeqs[ counter ] );	
			names[ i ].push_back( allNames[ counter ] );
			counter++;
		}
	}

	// read the binding data: one row per experiment (binding of all sequences in that experiment)
	ifstream fdata( bIFile.c_str() );
    if ( !fdata ) {
        cerr << "Cannot find the binding data file " << bIFile << endl;
        exit( 1 );
    }
	vector< vector< double > > bindingData( nExps, vector<double > (nSeqsb) );
	for ( int i = 0; i < nExps; i++ ) {
		for ( int j = 0; j < nSeqsb; j++ ) {
			string name;
			fdata >> name;			
			if( name != names[ i ][ j ] ) { 
				cerr << "Error: " << names[ i ][ j ] << "\t" << name << endl;
				exit( 1 );
			}
			fdata >> bindingData[ i ][ j ];	
 
		}
	}	


vector< vector< SiteVec > > seqSitesb( seqs2.size() );  // shouldn't this be seqs2.size(), changed to 2 on 1129
	SeqAnnotator annb( motifs, energyThrs );
	for ( int i = 0; i < nExps; i++ ) {
		for ( int j = 0; j < nSeqsb; j++ ) {
// 			cout << "Exp " << i << "\tSeq " << j << endl;
			SiteVec sites;
			annb.annot( seqs2[ i ][ j ], sites );  // annot( sequence, sitevector)  // annot takes empty sitevector
			seqSitesb[ i ].push_back( sites );
			
		}
	}	

    // create the expression predictor
    FactorIntFunc* intFunc; 
    if ( intOption == BINARY ) intFunc = new FactorIntFuncBinary( coopDistThr ); 
    else if ( intOption == GAUSSIAN ) intFunc = new FactorIntFuncGaussian( coopDistThr, factorIntSigma );
    else if ( intOption == BINSF ) intFunc = new FactorIntFuncBinsf( coopDistThr, nbin );
    else {
        cerr << "Interaction Function is invalid " << endl; exit( 1 ); 
    }
 for ( int i = 0; i < nSeqsa; i++ ) {

	for ( int j = 0; j < seqSitesa[ i ].size(); j++ ) {
            Site site = seqSitesa[i][j];
            char strandChar = site.strand ? '+' : '-';
          
        }
        }

   ExprPredictor* predictor = new ExprPredictor(seqSitesb, bindingData, seqSitesa, seqSites, seqLengths, exprData, motifs2, factorExprData, intFunc, coopMat, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, mmm, mmmr, indicator_bool, anndt, exprFile, seqSitesbot, seqSitesm1,seqSitesm2, seqSitesf2 ,seqSitesbotf2, seqSitesm1f2 ,seqSitesm2f2, seqSitesf3, seqSitesbotf3,seqSitesm1f3, seqSitesm2f3, dpp );  //520

    // read the initial parameter values
    ExprPar par_init( nFactors);  // 519
	if ( !parFile.empty() ) {
        rval = par_init.load( parFile, coopMat, repIndicators );
        if ( rval == RET_ERROR ) {
            cerr << "Cannot read parameters from " << parFile << endl;
            exit( 1 );
        } 
    	}


	gsl_rng* rng;
	gsl_rng_env_setup();
	const gsl_rng_type * T = gsl_rng_default;	// create rng type
	rng = gsl_rng_alloc( T );
	gsl_rng_set( rng, time( 0 ) );		// set the seed equal to simulTime(0)

(*predictor).clasvar = 0;

ofstream ma("muteseqEner.txt");


	ma  << dcmm.getMean()<< '\t' ;
	ma  << dumm.getMean()<< '\t' ;
	ma  << twmm.getMean()<< '\t' ;
	ma  << cbomf.getMean()<< '\t' ;
	ma << endl;
ostringstream numberm;
numberm <<fi;
string daFile= ("da.txt"  + numberm.str());
ofstream da( daFile.c_str());
ofstream da2( "lnodds.txt" );
pv = 0;
	 for( int m = 1 ; m< 500 ; m++){

			
				pv = pv + .1; 
		
				energyThrs2[0]=  pv; 
			
				energyThrs2[1]= pv ; //pv/10;
				
			  	energyThrs2[3]= pv ;  //pv/10;
				
				energyThrs2[2]= 0;//pv;
		
				predictor->modanny( energyThrs2 );

				double ot3 = predictor->compAvgCorr3( par_init , 1 );

		}// for m spacers









    vector< double > energyThrsFP( motifs2.size(), energyThr );
    vector< Sequence > seqsaFP;
    vector< string > seqNamesaFP;
    rval = readSequences( adamiFileFP, seqsaFP, seqNamesaFP );
    assert( rval != RET_ERROR );
    int nSeqsaFP = seqsaFP.size();

vector< SiteVec > seqSitesaFP( nSeqsaFP );
vector< int > seqLengthsaFP( nSeqsaFP );
vector< Motif > motifs4;
motifs4 = motifs2;
   
SeqAnnotator anndtFP( motifs4, energyThrsFP );
        for ( int i = 0; i < nSeqsaFP; i++ ) {
            anndtFP.annot( seqsaFP[ i ], seqSitesaFP[ i ] );
            seqLengthsaFP[i] = seqsaFP[i].size();
        }

 vector< vector< vector< Site > > > dpp2(nSeqsaFP, vector< vector< Site > >( 4 ));

   ExprPredictor* predictor2 = new ExprPredictor(seqSitesb, bindingData, seqSitesaFP, seqSites, seqLengths, exprData, motifs4, factorExprData, intFunc, coopMat, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, mmm, mmmr, indicator_bool, anndtFP, exprFile, seqSitesbot, seqSitesm1,seqSitesm2, seqSitesf2 ,seqSitesbotf2, seqSitesm1f2 ,seqSitesm2f2, seqSitesf3, seqSitesbotf3,seqSitesm1f3, seqSitesm2f3, dpp2 );  //520


pv=0;
 	for( int m = 1 ; m< 500 ; m++){

			
				pv = pv + .1;  
		
				energyThrsFP[0]=  pv;   
			
				energyThrsFP[1]= pv ; 
				
				energyThrsFP[3]= pv ;  

				energyThrsFP[2]= 0;
		
			
				predictor2->modanny( energyThrsFP );

				double ot3 = predictor2->compAvgCorrFP( par_init , 1 );

		}// for m spacers




da2.close();
da.close();

ofstream to2("ot2.txt");
	predictor->printFile25(to2,predictor->getPar(),*predictor);  // rememeber to remove format.tex
	to2.close();

ofstream to22("ot22.txt");
	//predictor->printFile24(to22,predictor->getPar(),*predictor);  // rememeber to remove format.tex
	to22.close();
 
  return 0;	
}

