
#include "ExprPredictor.h"

int main( int argc, char* argv[] ) 
{
 
    // command line processing
    string seqFile, annFile,duFile,dcFile, twFile, adamiFile, exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, seqFileb, bIFile, motifFile2;
    string outFile, cbFileo;     // output file
    double coopDistThr = 150;
    double factorIntSigma = 50.0;   // sigma parameter for the Gaussian interaction function
    double repressionDistThr = 150;
    double energyThr = 2;
    int maxContact = 1;
	double pva=.05;
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
	else if ( !strcmp( "-m2", argv[ i ] ) )
            motifFile2 = argv[ ++i ];
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
        else if ( !strcmp( "-a", argv[ i ] ) )
            annFile = argv[ ++i ];      
        else if ( !strcmp( "-e", argv[ i ] ) )
            exprFile = argv[ ++i ];            
        else if ( !strcmp( "-m", argv[ i ] ) )
            motifFile = argv[ ++i ];
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

        else if ( !strcmp( "-w", argv[i] ) )
            wid = atoi( argv[++i] );
    
  	else if ( !strcmp( "-l", argv[ i ] ) )
	    l = atoi( argv[ ++i ] );
    }
    if ( seqFile.empty() || exprFile.empty() || motifFile.empty() || factorExprFile.empty() || outFile.empty() || ( ( ExprPredictor::modelOption == QUENCHING || ExprPredictor::modelOption == CHRMOD_UNLIMITED || ExprPredictor::modelOption == CHRMOD_LIMITED ) &&  factorInfoFile.empty() ) || ( ExprPredictor::modelOption == QUENCHING && repressionFile.empty() ) ) {
        cerr << "Usage: " << argv[ 0 ] << " -s seqFile -e exprFile -m motifFile -f factorExprFile -fo outFile [-a annFile -o modelOption -c coopFile -i factorInfoFile -r repressionFile -oo objOption -mc maxContact -p parFile -rt repressionDistThr -et energyThr -na nAlternations -ct coopDistThr -sigma factorIntSigma]" << endl;
        cerr << "modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited" << endl;
        exit( 1 );
    }  

// 'mmm' vector contains the spacer window borders, this vector is initialized here, and then modified in the for loop after the ExprPredictor object has been created.
vector< int > mmm(4);
mmm[0] = 1;   // 5
mmm[1] = 1;   //25
mmm[2] = 3;  //50000
mmm[3] = 50001;
// 'mmmr' is quenching function for 'thermodynamic modeling'.
vector< int > mmmr(4);
mmmr[0] = 10;
mmmr[1] = 60;
mmmr[2] = 70;
mmmr[3] = 100;
//mmm[4] = 50;


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

    // read the sequences
    vector< Sequence > seqs;
    vector< string > seqNames;
    rval = readSequences( seqFile, seqs, seqNames );
    assert( rval != RET_ERROR );
    int nSeqs = seqs.size();
//////////////////////////////////////////////////7292011

rval = readSequences( seqFile, ExprPredictor::seqsy, ExprPredictor::seqNmes );
/////////////////////////////////////////////////

  // read the sequences
    vector< Sequence > seqsa;
    vector< string > seqNamesa;
    rval = readSequences( adamiFile, seqsa, seqNamesa );
    assert( rval != RET_ERROR );
    int nSeqsa = seqsa.size();
//////////////////////////////////////////////////7292011
//cout << " adamifile " << endl;
rval = readSequences( adamiFile, ExprPredictor::seqsya, ExprPredictor::seqNmesa );
//cout << " adamifile " << endl;



    // read the expression data
    vector< string > condNames;  
    rval = readMatrix( exprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nSeqs );
    for ( int i = 0; i < nSeqs; i++ ) assert( labels[i] == seqNames[i] );
    Matrix exprData( data ); 
    int nConds = exprData.nCols();
    //////////////////////////////////////////////73011 (this needs to be checked for shallow copy complications.
ExprPredictor::exprData2 = exprData;
//////////////////////////////////////////////////////////////
    // read the motifs

	vector< Motif > motifs;
	vector< Motif > motifs2;
	vector< string > motifNames;
	vector< double > background = createNtDistr( gcContent );
	vector< Motif > motifs3;
	vector< string > motifNames3;
	rval = readMotifs( motifFile2, background, motifs3, motifNames3 ); 
	assert( rval != RET_ERROR );

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
	Motif dcmm( dcm, background);
	double pv =0;
	//pv =dcmm.PvalueSum( pva );
	energyThrs2[0]=pv;
	motifs2.push_back( dcmm ); 
	scorematrix( motifs2[0].getLLRMat(), duFile , motifs2[0].getMaxLLR(), duc );
	scorematrix( motifs2[0].getLLRMat(), dcFile , motifs2[0].getMaxLLR(), dcs );
/////////////////////////////////////////////////////////////
//  DU
//////////////////////////////////////////////////////////////////
	string dcu="duconditdc.txt";
	string a = "ducondit.txt";
	string ducs = "ducondits.txt";
	Matrix dum2 = countmatrix2( duFile, a );

	Matrix dum = countmatrixS( duFile );
	Motif dumm( dum, background);
	//pv =dumm.PvalueSum( pva );
	energyThrs2[1]= pv; //2.6
	motifs2.push_back( dumm );     //  this is 2 is seqSitesa
	scorematrix( motifs2[1].getLLRMat(), dcFile , motifs2[1].getMaxLLR(), dcu );
///////////////////////////////////////////////////////////////////////////
//  TW
////////////////////////////////////////////////////////////////////////

	string c = "twcondit.txt";
	Matrix twm2 = countmatrix2( twFile , c );
	Matrix twm = countmatrix( twFile );
	Motif twmm( twm, background);
	double ethresT = .1;//3;
	twmm.setEth( ethresT );
	energyThrs2[2]=0;//.5;
	motifs2.push_back( motifs3[0]); 
//////////////////////////////////////////////////////////////////////////////
//    CB
////////////////////////////////////////////////////////////////////////////////////
	Matrix cbom = countmatrix( cbFileo );
	Motif cbomf( cbom, background);
	string cbfi="cbcondit.txt";
	string cbfis="cbcondits.txt";
	Matrix cbfim = countmatrix2( cbFileo , cbfi );
	//pv =cbomf.PvalueSum( pva );
	energyThrs2[3]= pv;
	motifs2.push_back( cbomf );     
	scorematrix( motifs2[3].getLLRMat(), cbFileo , motifs2[3].getMaxLLR(), cbfis );

	ifstream tws( duFile.c_str() );
	ifstream dco( dcFile.c_str() );

	string temp;
	int leng=0;
	int len2 =0;
while(!dco.eof()){
		temp = "";
		getline(dco, temp);
		if(temp.length() == 0){
			break;
		}
		getline(dco, temp);
		leng=leng + 1;   
}
while(!tws.eof()){
		temp = "";
		getline(tws, temp);
		if(temp.length() == 0){
			break;
		}
		getline(tws, temp);
		len2=len2 + 1; 
}
tws.close();
dco.close();

int lencb = leng + len2;   //  number of sites in cb, which 


string cbfisz="cbits.txt";
string deric="deric.txt";
string kelly = "kelly.txt";
double len =0;
scorematrix( motifs2[3].getInfoMat()*(-1), cbFileo , len, cbfisz );
scorematrix( motifs2[0].getInfoMat()*(-1), dcFile , len, deric );
scorematrix( motifs2[1].getInfoMat()*(-1), duFile , len, kelly );
    rval = readMotifs( motifFile, background, motifs, motifNames ); 
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

		//d2s.setRow( i,  dtemp.getRow(i) );
		if ( i == 4 ){
		d2s.push_back( e );

		//d2s.setRow( i,  e );
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


vector< vector< vector< Site > > > sitespp(nSeqsa, vector< vector< Site > >(4));
vector< vector< vector< Site > > > sitesptex(nSeqsa, vector< vector< Site > >(4));
for( int e = 0; e < nSeqsa; e++ ) {
	sitesptex[e][0] = seqSitesaex[e];    // this was filled with core DC
	sitesptex[e][1] = seqSitesaexu[e];      // this was filled with core DU
}

for ( int e = 0; e < nSeqsa; e++ ) {  // runs over enhnacers
	  for ( int k = 0; k < 3; k++ ) {  // 
		if( k == 0 ) {
				
			while( !sitesptex[e][k].empty()){
				if( sitesptex[e][k].size() > 0 ) {
				sitespp[e][k].push_back(  anndt.siteMax(   sitesptex[e][0]   ));  
				anndt.sitestoverlap( sitesptex[e][k], sitespp[e][k] );  
				}//if site
			}// while
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
             //dc
 vector< SiteVec > seqSitesm1;
            //du
for ( int e = 0; e < nSeqsa; e++ ) {  // runs over enhnacers
	  for ( int k = 0; k < 3; k++ ) {  // 
		if( k == 0 ) {
				
				seqSitesbot.push_back( sitespp[e][k] ); // sitespp[e][k] = sitespt[e][k]; siteMax(  sitespt[e][k] 
			
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
   
        cout << "Cannot open the cooperativity file " << endl;
////////////////////////////////////////////////////////////////////////////////////
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
	vector< vector< double > > bindingData( nExps, vector< double >(nSeqsb) );
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


vector< vector< SiteVec > > seqSitesb( seqs2.size() ); 
	SeqAnnotator annb( motifs, energyThrs );
	for ( int i = 0; i < nExps; i++ ) {
		for ( int j = 0; j < nSeqsb; j++ ) {
// 			cout << "Exp " << i << "\tSeq " << j << endl;
			SiteVec sites;
			annb.annot( seqs2[ i ][ j ], sites );  // annot( sequence, sitevector)  // annot takes empty sitevector
			seqSitesb[ i ].push_back( sites );
			
		}
	}	

////////////////
  
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
	mmm[2] = 4 ;
	double ot3 = predictor->compAvgCorr3( par_init , 1 );
	ma  << dcmm.getMean()<< '\t' ;
	ma  << dumm.getMean()<< '\t' ;
	ma  << twmm.getMean()<< '\t' ;
	ma  << cbomf.getMean()<< '\t' ;
	ma << endl;
	double score1 = ncountsdc*dcmm.getInformation();
	cout << "ncounts for initial score1 = " << ncountsdc << endl;
	double score2 =0;
	double scoretemp =0;
	double scorebest=score1;
	string bestRC="rcReadb.txt";
	//string dacu= "DorsalUnCondit2.txt";
	//align( dacu ,nrandStarts,  5 ,  wid,  l );
	double base = 10.0;
	double expo = -5.0;
	double pvalsmm =pow(base,expo);
ofstream pvalst("Data/pvalstra.txt");
ofstream pvalst2("Data/pvalstra2.txt");
ofstream pvalst3("Data/pvalstra3.txt");
ofstream mutuali("Data/minfo.txt");
ofstream kldiver("Data/kl.txt");

int nUp = 7;
//wid=10;
cout << " l " << l << endl;
cout << " energyThrs2 just before loop " << energyThrs2 << endl;


ifstream tws2( dcFile.c_str() );
ifstream dco2( duFile.c_str() );
ofstream ffff("combined.txt");
string temp2;
int leng2=0;
int len22 =0;
while(!dco2.eof()){
		temp2 = "";
		getline(dco2, temp2);
		ffff << temp << endl; 
		if(temp2.length() == 0){
			break;
		}
		getline(dco2, temp2);
		ffff << temp << endl;  // index 0 refers to dorsal conditional
		leng2=leng2 + 1;   // number of sites in dc/ giving class prior for 'adjacent to Twist'
}
while(!tws2.eof()){
		temp2 = "";
		getline(tws2, temp2);
		ffff << temp << endl; 		
		if(temp2.length() == 0){
			break;
		}
		getline(tws2, temp2);
		ffff << temp << endl;  // index 0 refers to dorsal conditional
		//ffff << motifs2[0].pwmvalue( temp )*double(len2)/double(lencb) << endl;
		len22=len22 + 1;  // number of sites in du/ giving class prior for 'not adjacent to Twist'
}
tws2.close();
dco2.close();
ffff.close();
int lencb2 = leng2 + len22;   //  number of sites in cb, which 




///////////////////////////////////////////////////////////////////////////////
	vector< Motif > _motifs2;
	_motifs2.push_back(dcmm );
	_motifs2.push_back(dumm );
	_motifs2.push_back(twmm );
	_motifs2.push_back(cbomf );
	//mutuali << mutualinformatio(_motifs2, double(leng2)/double(lencb2) ) << endl;
	_motifs2.clear();


cout << " motifs3.size() " << motifs3.size() << endl;
for( int k = 0; k < motifs3.size() ; k++){
//if( k > 0 ){
ostringstream ssk;
	//istringstream ss1;
	ssk << k;
cout << " dcmm in "<< endl << dcmm << endl;
string da21="lnodds" + ssk.str() + ".txt";
ofstream da2( da21.c_str() );
string mdis1="Data/weblogo-3.3/tmp2/dis" + ssk.str() + ".txt";
string infologodc1="Data/weblogo-3.3/tmp2/infodc" + ssk.str() + ".txt";
string infologodu1="Data/weblogo-3.3/tmp2/infodu" + ssk.str() + ".txt";
string infologocb1="Data/weblogo-3.3/tmp2/infocb" + ssk.str() + ".txt";
string nsitelogodc1="Data/weblogo-3.3/tmp2/nsitesdc" + ssk.str() + ".txt";
string nsitelogodu1="Data/weblogo-3.3/tmp2/nsitesdu" + ssk.str() + ".txt";
string nsitelogocb1="Data/weblogo-3.3/tmp2/nsitescb" + ssk.str() + ".txt";

string infocont="infocon.txt";
ofstream mdis( mdis1.c_str() );
ofstream infologodc(infologodc1.c_str());
ofstream infologodu(infologodu1.c_str());
ofstream infologocb(infologocb1.c_str());
ofstream nsitelogodc(nsitelogodc1.c_str());
ofstream nsitelogodu(nsitelogodu1.c_str());
ofstream nsitelogocb(nsitelogocb1.c_str());
ofstream infoconts(infocont.c_str());
		predictor->anny.killmo();
		predictor->anny.addmo( dcmm );  
		predictor->anny.addmo( dumm );
		cout << " motifs3[k] " << endl;
		cout << " k " << k << endl;
		predictor->anny.addmo( motifs3[ k ] );   
		predictor->anny.addmo( cbomf);
	energyThrs2[2]=0;
	
	mmm[2]=30;
for( int i = 1 ; i < nUpdates; i++){		
	if( i==1){ mmm[3]=300; }
	else { mmm[3]=50001;}
	vector< double > energyThrs3( 4, 0 ); 
	mdis << mmm[2]  << '\t';  
	pv =0;
	
	energyThrs3[0]=pv;
	
	energyThrs3[1]= pv;
	energyThrs3[2]=0;
	
	predictor->modanny( energyThrs3 );  
	double ot = predictor->compAvgCorr3( par_init , i );  // this is only for predictTwist2, which sets up mutualinfo data
	ofstream to2("ot2.txt");

	to2.close();
	ofstream to3("ot3.txt");
	
	to3.close();
	string dccc = "DorsalCondit2.txt" ;
	int lengthyc=0;
	ifstream filestrc;

	filestrc.open("DorsalCondit2.txt", ios::binary); // open your file
	filestrc.seekg(0, ios::end); // put the "cursor" at the end of the file
	lengthyc = filestrc.tellg(); // find the position of the cursor
	filestrc.close(); // close your file
	if ( lengthyc == 1 ){ 
		mutuali  << fixed << setprecision(2)<< mmm[2] << '\t' << 0 << endl;
		mmm[2] =30+30*i; mmm[1]=30*i;
		continue;
	}
	alignNEW( dccc , nrandStarts, nUp , wid,  l, i ,k);
	string dcccrcin = "coreOPT.txt"; //"rcRead.txt";
	string dcccrc="coreOPTout.txt";
	fastaflip( dcccrcin, dcccrc );

	Matrix dcmw = countmatrixS( dcccrc ); 
	Motif dcmmw( dcmw, background);
	
	string dcu= "DorsalUnCondit2.txt";
	
	int lengthy;
	ifstream filestr;

	filestr.open("DorsalUnCondit2.txt", ios::binary); 
	filestr.seekg(0, ios::end); // put the "cursor" at the end of the file
	lengthy = filestr.tellg(); // find the position of the cursor
	filestr.close(); 

	if ( lengthy == 0 ){ continue;}
	//else { }
	alignNEW2( dcu ,nrandStarts,  nUp ,  wid,  l,i,k );
	string dcccrc2in = "coreOPT2.txt";
	string dcccrc2="coreOPT2out.txt";
	
	fastaflip( dcccrc2in, dcccrc2);
	
	Matrix dum2 = countmatrixS( dcccrc2 );

	Motif dumm2( dum2, background);
	
	
	string cbpart ="combi.txt";
	ifstream dus( dcccrc2.c_str() );
	ifstream dco( dcccrc.c_str() );
	ofstream fffff(cbpart.c_str());
	string temp;
	int lengdc=0;
	int len2du =0;
	while(!dco.eof()){
			temp = "";
			getline(dco, temp);
			
			if(temp.length() == 0){
				break;
			}
			fffff << temp << endl;
			getline(dco, temp);
			if(temp.length() == 0){
				break;
			}
			fffff << temp << endl;//ffff << temp << endl;  // index 0 refers to dorsal conditional
			lengdc=lengdc + 1;   // number of sites in dc/ giving class prior for 'adjacent to Twist'
	}
	while(!dus.eof()){
			temp = "";
			getline(dus, temp);
			
			if(temp.length() == 0){
				break;
			}
			fffff << temp << endl;
			getline(dus, temp);
			if(temp.length() == 0){
				break;
			}
			fffff << temp << endl;
			len2du=len2du + 1;  // number of sites in du/ giving class prior for 'not adjacent to Twist'
	}
	
	tws.close();
	dco.close();
	fffff.close();
	int lencb = lengdc + len2du;   //  number of sites in cb, which 
	//alignNEW3( cbpart , nrandStarts, nUp , 5,  l, i , k );


	ostringstream ss2;
	//istringstream ss1;
	ss2 << i;
	
	
	string cbpart3 = "coreOPT" +ssk.str() +"cb" + ss2.str() + ".txt";//"rcRead.txt";
	string dcffa="coreOPT"+ ssk.str() + "dc" +ss2.str() +".txt";  // k i
	string duffa="coreOPT" + ssk.str() + "du" + ss2.str() +".txt";
	//string cbpart3shif = "core" + ss2.str() +"OPT1.txt";//"rcRead.txt";
	
	copyf(cbpart,cbpart3);
	copyf(dcccrc,dcffa);
	copyf(dcccrc2,duffa);

	string tmmot3="coreOPT" +ssk.str() +"tw"+ ss2.str() + ".txt";
	
	ofstream twmot3s( tmmot3.c_str() );
	twmot3s << motifs3[ k] << endl;
		
	////////////////////////////////////////////////////// creating cb shift with dc and du shift
string cbpartshif = "corecb" + ss2.str() +"OPT1.txt";//"rcRead.txt";
string dcFile2in2= "coredc"+ ss2.str()+ "OPT1.txt";  // coreshift of 1
string duFile2in2= "coredu"+ ss2.str()+ "OPT1.txt";  // coreshift of 1

//return 0;
	
	fastaflip( duFile2in2, duFile2in2 );
	fastaflip(dcFile2in2,dcFile2in2);

	ifstream dus2( duFile2in2.c_str() );
	ifstream dco2( dcFile2in2.c_str() );
	ofstream fffff2(cbpartshif.c_str());
	string tempe;
	int lengdcb=0;
	int len2dub =0;
	while(!dco2.eof()){
			tempe = "";
			getline(dco2, tempe);
			
			if(tempe.length() == 0){
				break;
			}
			fffff2 << tempe << endl;
			getline(dco2, tempe);
			if(tempe.length() == 0){
				break;
			}
			fffff2 << tempe << endl;//ffff << temp << endl;  // index 0 refers to dorsal conditional
			lengdcb=lengdcb + 1;   // number of sites in dc/ giving class prior for 'adjacent to Twist'
	}
	while(!dus2.eof()){
			tempe = "";
			getline(dus2, tempe);
			
			if(tempe.length() == 0){
				break;
			}
			fffff2 << tempe << endl;
			getline(dus2, tempe);
			if(tempe.length() == 0){
				break;
			}
			fffff2 << tempe << endl;
			//ffff << temp << endl;  // index 0 refers to dorsal conditional
			//ffff << motifs2[0].pwmvalue( temp )*double(len2)/double(lencb) << endl;
			len2dub=len2dub + 1;  // number of sites in du/ giving class prior for 'not adjacent to Twist'
	}

	dus2.close();
	dco2.close();
	fffff2.close();


	nsitelogodc << lengdc << '\t';
	nsitelogodu << len2du << '\t';
	nsitelogocb << lencb << '\t';
	//dcff.close();
	//duff.close();
	twmot3s.close();
	//coreshifdc.close();
	//coreshifdu.close();
	Matrix dcmw2 = countmatrixS( dcffa ); //dcccrc ); //, ncountsdc );
	Motif dcmmw2( dcmw2, background);
	infologodc << fixed << setprecision(1) << dcmmw2.getInformation()  << '\t'; //- 2.164/ double(lengdc) << '\t';
	Matrix dumw2 = countmatrixS( duffa ); //dcccrc ); //, ncountsdc );
	Motif dummw2( dumw2, background);
	infologodu  << fixed << setprecision(1)<< dummw2.getInformation()  << '\t';  // - 2.164/ double(len2du) << '\t';
	Matrix dumw2cb = countmatrixS( cbpart3 ); //dcccrc ); //, ncountsdc );
	Motif dummw2cb( dumw2cb, background);
	infologocb<< fixed  << setprecision(1)<< dummw2cb.getInformation()   << '\t'; // - 2.164/ double(lencb)<< '\t';
	infoconts<<fixed << setprecision(3)<< "dc du cb " << dcmmw2.getInformation()  << '\t'<< dummw2.getInformation()  << '\t'<< dummw2cb.getInformation() << endl;
///////////////////////////////////////////////////////////////////////////////
	vector< Motif > _motifs;
	_motifs.push_back(dcmmw2 );
	_motifs.push_back(dummw2 );
	_motifs.push_back( motifs3[k] );
	_motifs.push_back(dummw2cb );
	SeqAnnotator aaa(_motifs, energyThrs2 );
	double cmi=.5;
 	vector< vector< vector< Site > > > dpp2(nSeqsa, vector< vector< Site > >( 4 ));
	mutuali  << fixed << setprecision(2)<< mmm[2] << '\t' << mutualinformatio9(_motifs,double(lengdc)/double(lencb) ) << endl;
	

	mmm[2] = 30+30*i; 
	mmm[1] = 30*i; 
}
da2.close();
mdis.close();
infologodc.close(); //("Data/weblogo-3.3/tmp2/infodc.txt");
infologodu.close(); //("Data/weblogo-3.3/tmp2/infodu.txt");
infologocb.close(); //("Data/weblogo-3.3/tmp2/infocb.txt");
nsitelogodc.close(); //("Data/weblogo-3.3/tmp2/nsitesdc.txt");
nsitelogodu.close(); //("Data/weblogo-3.3/tmp2/nsitesdu.txt");
nsitelogocb.close(); //("Data/weblogo-3.3/tmp2/nsitescb.txt");
infoconts.close();
}// for k

mutuali.close();
kldiver.close();


ofstream to2("ot2.txt");
	predictor->printFile25(to2,predictor->getPar(),*predictor);  // rememeber to remove format.tex
	to2.close();

ma.close();
pvalst.close() ;
pvalst2.close();
pvalst3.close();



    return 0;	
}








