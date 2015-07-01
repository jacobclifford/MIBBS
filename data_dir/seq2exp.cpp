#include "ExprPredictor.h"

int main( int argc, char* argv[] ) 
{

    string seqFile, annFile,duFile,dcFile, twFile, adamiFile, exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, seqFileb, bIFile;
    string outFile, ipfile, seqFilet;     // output file
    double coopDistThr = 150;
    double factorIntSigma = 50.0;   // sigma parameter for the Gaussian interaction function
    double repressionDistThr = 150;
    double energyThr = 2;
	int initp =10;
    int maxContact = 1;
	int nrandStarts =10;
	 double pseudoCount =1;
	//int nStarts = 10;
	int nUpdates =4;
	int wid = 2;
	int D =0;
	int fi =0;
	int motifclass;
bool free_fix_indicators[] = {1,1,1,
			    1,1,0,0,0,
///*
0,0,0,0,0,
0,0,0,0,0,
1,0,0 };

	int l = 1;	// number of experiments
	vector <bool> indicator_bool ( free_fix_indicators, free_fix_indicators + sizeof( free_fix_indicators )/sizeof( bool ));

    int binwidth; // = 30;

    ExprPredictor::nAlternations = 1;
    for ( int i = 1; i < argc; i++ ) {
	
        if ( !strcmp( "-s", argv[ i ] ) )
            seqFile = argv[ ++i ];
	else if ( !strcmp( "-tw", argv[ i ] ) )
            twFile = argv[ ++i ];  
	else if ( !strcmp( "-dc", argv[ i ] ) )
            dcFile = argv[ ++i ]; 
	else if ( !strcmp( "-du", argv[ i ] ) )
            duFile = argv[ ++i ];   
	else if ( !strcmp( "-nu", argv[i] ) )
            nUpdates = atoi( argv[++i] );
	else if ( !strcmp( "-nr", argv[i] ) )
            nrandStarts = atoi( argv[++i] );
	else if ( !strcmp( "-ipf", argv[i] ) )
            ipfile =  argv[++i] ;
        else if ( !strcmp( "-et", argv[i] ) )
            energyThr = atof( argv[++i] );
        else if ( !strcmp( "-w", argv[i] ) )
            wid = atoi( argv[++i] );
        else if ( !strcmp( "-ip", argv[i] ) )
            initp = atoi( argv[++i] );
  	else if ( !strcmp( "-l", argv[ i ] ) )
	    l = atoi( argv[ ++i ] );
	else if ( !strcmp( "-fi", argv[ i ] ) )
	    fi = atoi( argv[ ++i ] );
	else if ( !strcmp( "-mot", argv[ i ] ) )
	    motifclass = atoi(argv[ ++i ] );
	else if ( !strcmp( "-psed", argv[i] ) )
            pseudoCount  = atof( argv[++i] );
    }
    double opt = 1;
    vector< Sequence > seqs;
    vector< string > seqNames;
    int rval;
    rval = readSequences( seqFile, seqs, seqNames );  
    string file = dcFile;

    ifstream seq_file( file.c_str() );
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
		int na = 0;	
		int naRC = 0;
	int countbias=0;

	for( int j = 0; j < seq.size(); j++ ){
		countbias = countbias + 1;
		Sequence readseq(seq[j]);
		Sequence RCreadseq( readseq.compRevCompl() );  //choose the a rich sequence for counting conservation
		//int na = 0;
		if(length==9){
		for(int i = 0;  i< readseq.size(); i++){
			if(readseq[i] == 0)            //const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' }
				{ na = na + 1;}
		}
		//int naRC = 0;
		for(int i = 0;  i< RCreadseq.size(); i++){
			if(RCreadseq[i] == 0)            //const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' }
				{ naRC = naRC + 1;}
		}
		}//if length==9
		else{
		for(int i = 1;  i< readseq.size()-1; i++){
			if(readseq[i] == 0)            //const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' }
				{ na = na + 1;}
		}
		//int naRC = 0;
		for(int i = 1;  i< RCreadseq.size()-1; i++){
			if(RCreadseq[i] == 0)            //const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' }
				{ naRC = naRC + 1;}
		}
		} //else length 9
		
	} // for j
	for( int j = 0; j < seq.size(); j++ ){

		Sequence readseq(seq[j]);
		Sequence RCreadseq( readseq.compRevCompl() ); 
		if( na > naRC ){
			otes << ">" << seq_name[j] << endl << readseq <<endl;
			for( int i = 0;  i< readseq.size(); i++)  // this for loop is miss-counting by one? maybe overloaded size method?
			{
				if(readseq[i] == 0)            //const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' }
				{ m(0,i) = m(0,i) + 1;}
				if(readseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(readseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(readseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				//if ( readseq.nts.empty() ) continue;
			}
		}  // if na
		else {
			otes << ">" << seq_name[j] << endl << RCreadseq <<endl;
			for( int i = 0;  i< RCreadseq.size(); i++)  // this for loop is miss-counting by one? maybe overloaded size method?
			{
				if(RCreadseq[i] == 0)            //const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' }
				{ m(0,i) = m(0,i) + 1;}
				if(RCreadseq[i] == 1)
				{ m(1,i) = m(1,i) + 1;}
				if(RCreadseq[i] == 2)
				{ m(2,i) = m(2,i) + 1;}
				if(RCreadseq[i] == 3)
				{ m(3,i) = m(3,i) + 1;}
				//if ( readseq.nts.empty() ) continue;
			}
		} //else


	} // for j
	Matrix countMatrix = m.transpose();
	assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
	  
	otes.close();
	Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;

        string dcFile2 = "rcRead.txt";
	double gcContent = .5;
	Matrix dcm4 = countmatrixS( dcFile2 );
        vector< double > background = createNtDistr( gcContent );

	Motif dcmm4( dcm4, pseudoCount, background);

        string dcFile3 = "rcRead2.txt";
        ofstream rcR2(dcFile3.c_str());

	rcR2 << dcmm4 << endl;

	rcR2.close();

        string shel= "info.txt";
	ofstream info(shel.c_str());
	info << setprecision(2) << dcmm4.getInformation()  << endl; //- 2.164/ double(countbias) << endl;
	info.close();
	string shel2= "pwm.txt";
	ofstream info2(shel2.c_str());
	for( int i=0; i< m.nRows(); i++){
		for( int j=0; j< m.nCols(); j++){
			info2 << fixed << setprecision(0) << m(i,j) ;
			if(  j==m.nCols()-1) info2 << endl;
			else info2 << '\t' ;
		}
	}
	for( int j=0; j< m.nCols(); j++){
		int sumer=0;
		for( int i=0; i< m.nRows(); i++){
			sumer=sumer+ m(i,j);
		}

	}

info.close();

    return 0;	
}

