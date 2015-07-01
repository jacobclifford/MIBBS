
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

rval = readSequences( adamiFileFP, ExprPredictor::seqsyaFP, ExprPredictor::seqNmesaFP ); 
rval = readSequences( seqFile, ExprPredictor::seqsy, ExprPredictor::seqNmes );

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
    //////////////////////////////////////////////73011 (this needs to be checked for shallow copy complications.
ExprPredictor::exprData2 = exprData;
//////////////////////////////////////////////////////////////
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
Motif dcmm( dcm, background);
double pv =1.5;
string aaadc ="dcondit2.txt";
string aaapdc="dccdfcore.txt";

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
Matrix dum2 =  countmatrix2( duFile, a );   
Motif dumm( dum2, background);

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
motifs2.push_back( twmm );      //  this is 1 in seqSitesa   

//////////////////////////////////////////////////////////////////////////////
//    CB
////////////////////////////////////////////////////////////////////////////////////
Matrix cbom = countmatrixS( cbFileo );
Motif cbomf( cbom, background);
string cbfi="cbcondit.txt";
string cbfis="cbcondits.txt";
string aaa = "cbcondit2.txt";
string aaap = "cbcdfcore.txt";
Matrix cbfim = countmatrix2( cbFileo , cbfi );    ////////aaa


energyThrs2[3]= pv;
vector< Sequence > seqsadc;
    vector< string > seqNamesadc;
    rval = readSequences( dcFile, seqsadc, seqNamesadc );
    assert( rval != RET_ERROR );
    double ndc = seqsadc.size();
vector< Sequence > seqsacb;
    vector< string > seqNamesacb;
    rval = readSequences( cbFileo, seqsacb, seqNamesacb );
    assert( rval != RET_ERROR );
    double ncb = seqsacb.size();
double fdc=ndc/ncb;
double fdu=(ncb-ndc)/ncb;

string cbfisz="cbits.txt";
string deric="deric.txt";
string kelly = "kelly.txt";
double len =0;
rval = readMotifs( motifFile, background, motifs, motifNames ); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
///////////////////////// w Matrix
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Matrix pwmdc =motifs2[0].getPwm();
Matrix pwmdu =motifs2[1].getPwm();
Matrix cbmix= pwmdc*fdc+ pwmdu*fdu;                                     //01182015
Motif cbomf2( cbmix, background);
motifs2.push_back( cbomf2);                                           //////01182015

Matrix pwmcb =motifs2[3].getPwm();
Matrix LLRdc =motifs2[0].getLLRMatpos();
Matrix LLRdu =motifs2[1].getLLRMatpos();
Matrix LLRcb =motifs2[3].getLLRMatpos();

int lenw=pwmdc.nCols();

 Matrix wmatdc(9,4,0);
double mi =0.0;
for( int p0 = 0 ; p0 < pwmdc.nRows() ; p0++){
	for( int j = 0 ; j < 4 ; j++){
	LLRdc(p0,j) = LLRcb(p0,j)+  log2( pwmcb(p0,j) / pwmdc(p0,j) );
	wmatdc(p0,j)=log2( pwmcb(p0,j) / pwmdc(p0,j) );
	}
}
 Matrix wmatdu(9,4,0);

for( int p0 = 0 ; p0 < pwmdc.nRows() ; p0++){
	for( int j = 0 ; j < 4 ; j++){
	LLRdu(p0,j) = LLRcb(p0,j)+  log2( pwmcb(p0,j) / pwmdu(p0,j) );
	wmatdu(p0,j)=log2( pwmcb(p0,j) / pwmdu(p0,j) );
	mi = mi - pwmdc(p0,j)*fdc*wmatdc(p0,j)-  pwmdu(p0,j)*fdu*wmatdu(p0,j) ;//P(S|C)P(C)log( P(S|C)/P(S) )
	}
}

scorematrixscanoverlap( LLRdc*(-1), dcFile, len, deric );
scorematrixscanoverlap( LLRdu*(-1), duFile , len, kelly );
scorematrixscanoverlap( LLRcb*(-1), cbFileo , len, cbfis );

string dericfp="dericfp.txt";
string kellyfp="kellyfp.txt";
string cbfp="cbfp.txt";

scorematrixscanoverlapboth( LLRdc*(-1), adamiFileFP , len, dericfp );
scorematrixscanoverlapboth(  LLRdu*(-1), adamiFileFP , len, kellyfp );
scorematrixscanoverlapboth( LLRcb*(-1), adamiFileFP , len, cbfp );


string kellyfp2="kellyfptw.txt";
string dericfp2="dericfptw.txt";

scorematrix( LLRdc*(-1), duFile , len, dericfp2 );
scorematrix( LLRdu*(-1), dcFile , len, kellyfp2 );


  return 0;	
}

