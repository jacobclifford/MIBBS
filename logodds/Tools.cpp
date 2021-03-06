#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include "Tools.h"
void copyf( string source, string dest )
{

 ifstream seq_file( source.c_str() );
 ofstream target( dest.c_str() );
    if ( !seq_file ) {
        //return RET_ERROR;
    }
	string temp;
		while(!seq_file.eof()){
		temp = "";

		getline(seq_file, temp);
		target << temp << endl;
		}
  seq_file.close();
  target.close();

}

const log_add_table table( -10.0, 0, 500 );		// global variable

// create a vector from gsl_vector
vector< double > gsl2vector( const gsl_vector* v )
{
    vector< double > result;
    for ( int i = 0; i < v->size; i++ ) result.push_back( gsl_vector_get( v, i ) );
    return result;	
}

gsl_vector* vector2gsl( const vector< double >& v )
{
    gsl_vector* result = gsl_vector_alloc( v.size() );
    for ( int i = 0; i < v.size(); i++ ) gsl_vector_set( result, i, v[ i ] );
    return result;
}

vector< int > gsl_int2vector( const gsl_vector_int* v )
{
    vector< int > result;
    for ( int i = 0; i < v->size; i++ ) result.push_back( gsl_vector_int_get( v, i ) );
    return result;	
}

gsl_vector_int* vector2gsl_int( const vector< int >& v )
{
    gsl_vector_int* result = gsl_vector_int_alloc( v.size() );
    for ( int i = 0; i < v.size(); i++ ) gsl_vector_int_set( result, i, v[ i ] );
    return result;
}

Matrix::Matrix( int nRows, int nCols ) 
{
    assert( nRows > 0 && nCols > 0 );
    data = gsl_matrix_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "Matrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }
}

Matrix::Matrix( int nRows, int nCols, double init ) 
{
    assert( nRows > 0 && nCols > 0 );
    data = gsl_matrix_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "Matrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }
    gsl_matrix_set_all( data, init );
}

Matrix::Matrix( const gsl_matrix* _data ) 
{
    assert( _data );
    assert( _data->size1 > 0 && _data->size2 > 0 ); 
    int nRows = _data->size1, nCols = _data->size2;
    data = gsl_matrix_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "Matrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }
    gsl_matrix_memcpy( data, _data );
}

// Matrix: constructor
Matrix::Matrix( const double **_data, int nRows, int nCols )
{
    assert( _data && nRows > 0 && nCols > 0 );
    data = gsl_matrix_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "Matrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }    
    for ( int i = 0; i < nRows; ++i ) {
        for ( int j = 0; j < nCols; ++j ) {
                gsl_matrix_set( data, i, j, _data[ i ][ j ] ); 	
        }	
    }
}

Matrix::Matrix( const vector< vector< double > >& _data )
{
    if ( _data.size() == 0 ) {
        data = NULL; 
        return;
    }

    int nRows = _data.size();
    int nCols = _data[0].size();
//     cout << "nRows = " << nRows << " nCols = " << nCols << endl;
    data = gsl_matrix_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "Matrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }        
    for ( int i = 0; i < nRows; i++ ) {
//         cout << "Row " << i << " Size = " << _data[i].size() << endl;
        assert( _data[i].size() == nCols );
        for ( int j = 0; j < nCols; j++ ) {
            gsl_matrix_set( data, i, j, _data[i][j] );
        }
    }
}

// destructor
Matrix::~Matrix()
{
    if ( data ) gsl_matrix_free( data );	
}

// get the row vector
vector< double > Matrix::getRow( int row ) const
{
    gsl_vector* v = gsl_vector_alloc( nCols() );
    gsl_matrix_get_row( v, data, row );
    return gsl2vector( v );
}

// get the column vector
vector< double > Matrix::getCol( int col ) const
{
    gsl_vector* v = gsl_vector_alloc( nRows() );
    gsl_matrix_get_col( v, data, col );
    return gsl2vector( v );	
}

// set the matrix dimensions
void Matrix::setDimensions( int nRows, int nCols)
{
    assert( nRows > 0 && nCols > 0 );
    if ( data && ( data->size1 != nRows || data->size2 != nCols ) ) gsl_matrix_free( data );
    //if ( data ) gsl_matrix_free( data );
    data = gsl_matrix_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "Matrix setDimensions() failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }    
}

// equality test
bool Matrix::operator==( const Matrix& other ) const
{
    // compare dimensions
    if ( nRows() != other.nRows() || nCols() != other.nCols() )
        return false;
    
    // compare elements
    for ( int i = 0; i < nRows(); i++ ) {
        for ( int j = 0; j < nCols(); j++ ) {
            if ( getElement( i, j ) != other.getElement( i, j ) ) return false;	
        }	
    }
    
    return true;
}

// check if square
bool Matrix::isSquare() const
{
    if ( nRows() == nCols() ) 
        return true;
    else 
        return false;	
}

bool Matrix::isSymmetric() const
{
    assert( isSquare() );
    for ( int i = 0; i < nRows(); i++ ) {
        for ( int j = 0; j < i; j++ ) {
            if ( abs( getElement( i, j ) - getElement( j, i ) ) > DBL_EPSILON ) return false;
        }
    }

    return true;
}

// Matrix: set row
void Matrix::setRow( int row, const vector< double >& v )
{
    assert( row >= 0 && row < data->size1 );
    
    for ( int j = 0; j < data->size2; j++ ) 
        setElement( row, j, v[ j ] );	
}

// Matrix: set column
void Matrix::setCol( int col, const vector< double >& v )
{
    assert( col >= 0 && col < data->size2 );
    
    for ( int i = 0; i < data->size1; i++ ) 
        setElement( i, col, v[ i ] );	
            
}

// Matrix: set all rows
void Matrix::setRows( const vector< double >& v )
{
    for ( int i = 0; i < data->size1; i++ )
        setRow( i, v );	
}

// Matrix: set all columns
void Matrix::setCols( const vector< double >& v )
{
    for ( int j = 0; j < data->size2; j++ ) 
        setCol( j, v );	
}

// set the identity matrix of dimension n
void Matrix::setIdentityMatrix( int n )		
{
    assert( n > 0 );
    setDimensions( n, n );	
    setZero();
    for ( int i = 0; i < n; i++ ) setElement( i, i, 1.0 );
}

// set the diagonal matrix
void Matrix::setDiagonalMatrix( const vector< double >& diag )
{
    int n = diag.size();
    
    setDimensions( n, n );
    setZero();
    for ( int i = 0; i < n; i++ ) setElement( i, i, diag[ i ] );
}

// Matrix transpose
Matrix Matrix::transpose() const
{
    Matrix result( nCols(), nRows() );
    
    for ( int i = 0; i < nRows(); i++ ) {
        for ( int j = 0; j < nCols(); j++ ) {
            result( j, i ) = getElement( i, j );
        }
    }
    return( result );
}

// matrix addition
Matrix Matrix::operator+( const Matrix& other ) const
{
    // check dimensions
    assert( nRows() == other.nRows() && nCols() == other.nCols() );
    
    // addition
    Matrix result( *this );
    gsl_matrix_add( result.data, other.data );	

    return result;
}

// matrix addition
Matrix& Matrix::operator+=( const Matrix& other )
{
    // check dimensions
    assert( nRows() == other.nRows() && nCols() == other.nCols() );

    // addition
    gsl_matrix_add( data, other.data );
    
    return ( *this );
}

// matrix multiplication by a constant
Matrix Matrix::operator*( const double c ) const
{
    Matrix result( *this );
    gsl_matrix_scale( result.data, c );	

    return result;
}	

// matrix multiplication by a constant
Matrix& Matrix::operator*=( const double c )
{
    gsl_matrix_scale( data, c );	
    
    return ( *this );	
}

// matrix multiplication
Matrix Matrix::operator*( const Matrix& other ) const
{
    // check dimensions
    assert( nCols() == other.nRows() );
    
    // matrix multiplication
    Matrix result( nRows(), other.nCols() );
    gsl_linalg_matmult( data, other.data, result.data );
    
    return result;
}

// output the Matrix
ostream& operator<<( ostream& os, const Matrix& m )
{
    os.setf( ios::fixed );
    
    for ( int i = 0; i < m.nRows(); i++ ) {
        for ( int j = 0; j < m.nCols(); j++ ) {
            os << m( i, j );
            if ( j != m.nCols() - 1 ) os << "   ";
        }
        os << endl;
    }

    return os;
}

// load Matrix from a file
// If readDims = true, read dimensions (first row of the file); o/w, only read data, assuming the correct dimensions are known
int Matrix::load( const string& fileName, bool readDims )
{
    // open the file
    ifstream fin( fileName.c_str() );
    if ( !fin ){ cerr << "Cannot open " << fileName << endl;	exit( 1 ); } 
    
    // read the matrix dimensions
    int nRows, nColumns;
    if ( readDims ) {
        fin >> nRows >> nColumns;
        assert( nRows > 0 && nColumns > 0 );
        setDimensions( nRows, nColumns );		
    } else {
        nRows = data->size1;
        nColumns = data->size2;
    }

    // read the matrix 
    for ( int i = 0; i < nRows; i++ ) {
        for ( int j = 0; j < nColumns; j++ ) {
            double elem;
            fin >> elem;
            gsl_matrix_set( data, i, j, elem );
        }
    }
    
    return 0;
}

// save Matrix to a file
void Matrix::save( const string& fileName ) const
{
    // open the file
    ofstream fout( fileName.c_str(), ofstream::out | ofstream::app );
    
    // write the matrix
    fout << *this;	
}

IntMatrix::IntMatrix( int nRows, int nCols ) 
{
    assert( nRows > 0 && nCols > 0 );
    data = gsl_matrix_int_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "IntMatrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }
}

IntMatrix::IntMatrix( int nRows, int nCols, int init ) 
{
    assert( nRows > 0 && nCols > 0 );
    data = gsl_matrix_int_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "IntMatrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }
    gsl_matrix_int_set_all( data, init );
}

IntMatrix::IntMatrix( const gsl_matrix_int* _data ) 
{
    assert( _data );
    assert( _data->size1 > 0 && _data->size2 > 0 ); 
    int nRows = _data->size1, nCols = _data->size2;
    data = gsl_matrix_int_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "IntMatrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }
    gsl_matrix_int_memcpy( data, _data );
}

// IntMatrix: constructor
IntMatrix::IntMatrix( const int **_data, int nRows, int nCols )
{
    assert( _data && nRows > 0 && nCols > 0 );
    data = gsl_matrix_int_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "IntMatrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }    
    for ( int i = 0; i < nRows; ++i ) {
        for ( int j = 0; j < nCols; ++j ) {
            gsl_matrix_int_set( data, i, j, _data[ i ][ j ] ); 	
        }	
    }
}

IntMatrix::IntMatrix( const vector< vector< int > >& _data )
{
    if ( _data.size() == 0 ) {
        data = NULL; 
        return;
    }

    int nRows = _data.size();
    int nCols = _data[0].size();
    data = gsl_matrix_int_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "Matrix creation failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }        
    for ( int i = 0; i < nRows; i++ ) {
        assert( _data[i].size() == nCols );
        for ( int j = 0; j < nCols; j++ ) {
            gsl_matrix_int_set( data, i, j, _data[i][j] );
        }
    }
}

// destructor
IntMatrix::~IntMatrix()
{
    if ( data ) gsl_matrix_int_free( data );	
}

// get the row vector
vector< int > IntMatrix::getRow( int row ) const
{
    gsl_vector_int* v = gsl_vector_int_alloc( nCols() );
    gsl_matrix_int_get_row( v, data, row );
    return gsl_int2vector( v );
}

// get the column vector
vector< int > IntMatrix::getCol( int col ) const
{
    gsl_vector_int* v = gsl_vector_int_alloc( nRows() );
    gsl_matrix_int_get_col( v, data, col );
    return gsl_int2vector( v );	
}

// set the matrix dimensions
void IntMatrix::setDimensions( int nRows, int nCols )
{
    assert( nRows > 0 && nCols > 0 );
    if ( data && ( data->size1 != nRows || data->size2 != nCols ) ) gsl_matrix_int_free( data );
    //if ( data ) gsl_matrix_int_free( data );
    data = gsl_matrix_int_alloc( nRows, nCols );
    if ( !data ) {
        cerr << "IntMatrix setDimensions() failed. nRows = " << nRows << " nCols = " << nCols << endl;
        exit( 1 );
    }    
}

// equality test
bool IntMatrix::operator==( const IntMatrix& other ) const
{
    // compare dimensions
    if ( nRows() != other.nRows() || nCols() != other.nCols() )
        return false;
    
    // compare elements
    for ( int i = 0; i < nRows(); i++ ) {
        for ( int j = 0; j < nCols(); j++ ) {
            if ( getElement( i, j ) != other.getElement( i, j ) ) return false;	
        }	
    }
    
    return true;
}

// check if square
bool IntMatrix::isSquare() const
{
    if ( nRows() == nCols() ) 
        return true;
    else 
        return false;	
}

bool IntMatrix::isSymmetric() const
{
    assert( isSquare() );
    for ( int i = 0; i < nRows(); i++ ) {
        for ( int j = 0; j < i; j++ ) {
            if ( getElement( i, j ) != getElement( j, i ) ) return false;
        }
    }

    return true;
}

// IntMatrix: set row
void IntMatrix::setRow( int row, const vector< int >& v )
{
    assert( row >= 0 && row < data->size1 );
    
    for ( int j = 0; j < data->size2; j++ ) 
        setElement( row, j, v[ j ] );	
}

// IntMatrix: set column
void IntMatrix::setCol( int col, const vector< int >& v )
{
    assert( col >= 0 && col < data->size2 );
    
    for ( int i = 0; i < data->size1; i++ ) 
        setElement( i, col, v[ i ] );	
}

// IntMatrix: set all rows
void IntMatrix::setRows( const vector< int >& v )
{
    for ( int i = 0; i < data->size1; i++ )
        setRow( i, v );	
}

// IntMatrix: set all columns
void IntMatrix::setCols( const vector< int >& v )
{
    for ( int j = 0; j < data->size2; j++ ) 
        setCol( j, v );	
}

// set the identity matrix of dimension n
void IntMatrix::setIdentityMatrix( int n )		
{
    assert( n > 0 );
    setDimensions( n, n );	
    setZero();
    for ( int i = 0; i < n; i++ ) setElement( i, i, 1 );
}

// set the diagonal matrix
void IntMatrix::setDiagonalMatrix( const vector< int >& diag )
{
    int n = diag.size();
    
    setDimensions( n, n );
    setZero();
    for ( int i = 0; i < n; i++ ) setElement( i, i, diag[ i ] );
}

// IntMatrix transpose
IntMatrix IntMatrix::transpose() const
{
    IntMatrix result( nCols(), nRows() );
    
    for ( int i = 0; i < nRows(); i++ ) {
        for ( int j = 0; j < nCols(); j++ ) {
            result( j, i ) = getElement( i, j );
        }
    }
    return( result );
}

// matrix addition
IntMatrix IntMatrix::operator+( const IntMatrix& other ) const
{
    // check dimensions
    assert( nRows() == other.nRows() && nCols() == other.nCols() );
    
    // addition
    IntMatrix result( *this );
    gsl_matrix_int_add( result.data, other.data );	

    return result;
}

// matrix addition
IntMatrix& IntMatrix::operator+=( const IntMatrix& other )
{
    // check dimensions
    assert( nRows() == other.nRows() && nCols() == other.nCols() );

    // addition
    gsl_matrix_int_add( data, other.data );
    
    return ( *this );
}

// matrix multiplication by a constant
IntMatrix IntMatrix::operator*( const int c ) const
{
    IntMatrix result( *this );
    gsl_matrix_int_scale( result.data, c );	

    return result;
}	

// matrix multiplication by a constant
IntMatrix& IntMatrix::operator*=( const int c )
{
    gsl_matrix_int_scale( data, c );	
    
    return ( *this );	
}

// output the IntMatrix
ostream& operator<<( ostream& os, const IntMatrix& m )
{
    os.setf( ios::fixed );
    
    for ( int i = 0; i < m.nRows(); i++ ) {
        for ( int j = 0; j < m.nCols(); j++ ) {
            os << m( i, j );
            if ( j != m.nCols() - 1 ) os << " ";
        }
        os << endl;
    }

    return os;
}

// load IntMatrix from a file
// If readDims = true, read dimensions (first row of the file); o/w, only read data, assuming the correct dimensions are known
int IntMatrix::load( const string& fileName, bool readDims )
{
    // open the file
    ifstream fin( fileName.c_str() );
    if ( !fin ){ cerr << "Cannot open " << fileName << endl;	exit( 1 ); } 
    
    // read the matrix dimensions
    int nRows, nColumns;
    if ( readDims ) {
        fin >> nRows >> nColumns;
        assert( nRows > 0 && nColumns > 0 );
        setDimensions( nRows, nColumns );		
    } else {
        nRows = data->size1;
        nColumns = data->size2;
    }

    // read the matrix 
    for ( int i = 0; i < nRows; i++ ) {
        for ( int j = 0; j < nColumns; j++ ) {
            int elem;
            fin >> elem;
            gsl_matrix_int_set( data, i, j, elem );
        }
    }
    
    return 0;
}

// save IntMatrix to a file
void IntMatrix::save( const string& fileName ) const
{
    // open the file
    ofstream fout( fileName.c_str(), ofstream::out | ofstream::app );
    
    // write the matrix
    fout << *this;	
}

// log_add_table constructor: compute the table used for log_add algorithm
log_add_table::log_add_table( double a, double b, int N ) : x_array( N + 1 ), T( N + 1 )
{
    delta = ( b - a ) / N;
    x_array[ 0 ] = a;
    T[ 0 ] = log( 1 + exp( a ) );
    for ( int i = 1; i < N; i++ ) {
        double x = x_array[ i - 1 ] + delta;
        x_array[ i ] = x;
        T[ i ] = log( 1 + exp( x ) );	
    }
    x_array[ N ] = b; 
    T[ N ] = log( 1 + exp( b ) );
}

// log_add algorithm
double log_add( double p, double q )
{
    // check if p or q is -inf
    if ( gsl_isinf( p ) == -1 ) return q;
    if ( gsl_isinf( q ) == -1 ) return p;
    
    // log_add
    int N = table.x_array.size();
    double x, val;
    if ( p > q ) {
        x = q - p; val = p;
    }
    else {
        x = p - q; val = q;
    }
    if ( x < table.x_array[ 0 ] ) {
        val += 0;
    } else {
        int k = (int)floor( ( x - table.x_array[ 0 ] ) / table.delta );
        val += table.T[ k ];
        //if ( k == N ) 
        //	val += table.T[ N ];
        //else 
        //	val += table.T[ k ] + ( x - table.x_array[ k ] ) * ( table.T[ k + 1 ] - table.T[ k ] ) / table.delta;	
    }	
    
    return val;
}

// log_add algorithm
double log_add( const vector< double >& p )
{
    double result = p[ 0 ];
    
    for ( int i = 1; i < p.size(); i++ ) {
        if ( gsl_isinf( p[ i ] ) == -1 ) continue;
        result = log_add( result, p[ i ] );	
    }	
    
    return result;
}

// log_add algorithm
double log_add( double p, double q, double r )
{
    double result;
    result = log_add( p, q );
    result = log_add( result, r );
    
    return result;	
}

vector< double > log( const vector< double >& v )
{
    vector< double > result;
    for ( int i = 0; i < v.size(); i++ ) result.push_back( log( v[ i ] ) );
    return result;	
}

vector< double > exp( const vector< double >& v )
{
    vector< double > result;
    for ( int i = 0; i < v.size(); i++ ) result.push_back( exp( v[ i ] ) );
    return result;			
}

// log. of a matrix
Matrix log( const Matrix& M )
{
    int m = M.nRows(), n = M.nCols();
    
    Matrix result( m, n );	
    for ( int i = 0; i < m; i++ ) {
        for ( int j = 0; j < n; j++ ) {
            result( i, j ) = log( M( i, j ) );
        }
    }
            
    return result;
}	

// maximum of a set of numbers
double max( const vector< double >& v, int &arg )
{
    double v_max = GSL_NEGINF;
    for ( int i = 0; i < v.size(); ++i ) {
        if ( v[ i ] > v_max ) {
            arg = i;
            v_max = v[ i ];	
        }	
    }
    
    return v_max;
}

double mean( const vector< double >& x )
{
    if ( x.size() == 0 ) return 0;
    
    double sum = 0; 
    for ( int i = 0; i < x.size(); i++ ) sum += x[ i ];
    return ( sum / x.size() );
}

double median( const vector< double > &x )
{
    if ( x.size() == 0 ) { return 0; }
    
    vector< double > xCopy( x );
    sort( xCopy.begin(), xCopy.end() );
    int middle = (int)( ( xCopy.size() - 1 ) / 2.0 );
    return xCopy[ middle ];	
}

double std_dev( const vector< double >& x )
{
    double* data = new double[ x.size() ];
    for ( int i = 0; i < x.size(); i++ ) data[ i ] = x[ i ];
    return gsl_stats_sd( data, 1, x.size() );
}

double correlation( const vector< double >& v1, const vector< double >& v2 )
{
	assert( v1.size() == v2.size() );
	int n = v1.size();
	
	double* data1 = new double[ n ]; 
	double* data2 = new double[ n ];
	for ( int i = 0; i < n; i++ ) {
		data1[ i ] = v1[ i ];
		data2[ i ] = v2[ i ];	
	}
	
	return gsl_stats_covariance( data1, 1, data2, 1, n ) / ( gsl_stats_sd( data1, 1, n ) * gsl_stats_sd( data2, 1, n) );
}

double corr( const vector< double >& x, const vector< double >& y )
{
    if ( x.size() != y.size() ) return RET_ERROR; 
    if ( x.size() == 0 ) return RET_ERROR;

    // means of X and Y
    double x_bar = mean( x );
    double y_bar = mean( y );
    
    // pseudo-observation at tau = 0
    vector< double > X( x );
    vector< double > Y( y );
    X.insert( X.begin(), 0.0 );
    Y.insert( Y.begin(), 0.0 );
    int n = X.size() - 1;

    // variance of X and Y
    double sum_x = 0; 
    double sum_y = 0;
    for ( int s = 1; s <= n; s++ ) {
        sum_x += ( X[s] - x_bar ) * ( X[s] - x_bar );
        sum_y += ( Y[s] - y_bar ) * ( Y[s] - y_bar );
    }
    double x_var = sum_x / n;
    double y_var = sum_y / n; 
 //   cout << x_var << "xvar  is f" << endl;
    // covariance and correlation
    double sum = 0;
    for ( int s = 1; s <= n ; s++ ) {
        sum += ( X[s] - x_bar ) * ( Y[s] - y_bar );
    }
    double cov_xy = sum / n;
    double corr_xy = cov_xy / sqrt( x_var * y_var );
    if ( corr_xy != corr_xy ) return 1;
    return corr_xy;
}

int cross_corr( const vector< double >& x, const vector< double >& y, const vector< int >& lag, vector< double >& cov, vector< double >& corr )
{
    if ( x.size() != y.size() ) return RET_ERROR; 
    if ( x.size() == 0 ) return RET_ERROR;
    
    cov.clear();
    corr.clear();

    // means of X and Y
    double x_bar = mean( x );
    double y_bar = mean( y );
    
    // pseudo-observation at tau = 0
    vector< double > X( x );
    vector< double > Y( y );
    X.insert( X.begin(), 0.0 );
    Y.insert( Y.begin(), 0.0 );
    int n = X.size() - 1;

    // variance of X and Y
    double sum_x = 0; 
    double sum_y = 0;
    for ( int s = 1; s <= n; s++ ) {
        sum_x += ( X[s] - x_bar ) * ( X[s] - x_bar );
        sum_y += ( Y[s] - y_bar ) * ( Y[s] - y_bar );
    }
    double x_var = sum_x / n;
    double y_var = sum_y / n; 
    
    // covariance and correlation
    for ( int i = 0; i < lag.size(); i++ ) {
        int t = lag[i]; 
        if ( t > n - 1 || t < -( n - 1 ) ) return RET_ERROR;

        double sum = 0;
        if ( t >= 0 ) {
            for ( int s = 1; s <= n - t; s++ ) {
                sum += ( X[s + t] - x_bar ) * ( Y[s] - y_bar );
            }
        } else {
            for ( int s = 1 - t; s <= n; s++ ) {
                sum += ( X[s + t] - x_bar ) * ( Y[s] - y_bar );
            }
        }
        double cov_xy = sum / n;
        cov.push_back( cov_xy );
        double corr_xy = cov_xy / sqrt( x_var * y_var );
        corr.push_back( corr_xy );
    }
    
    return 0;
}

double least_square( const vector< double >& x, const vector< double >& y, double& beta )
{
    assert( x.size() == y.size() );
    int n = x.size();

    double numerator = 0, denom = 0;
    for ( int i = 0; i < n; i++ ) {
        numerator += x[i] * y[i];
        denom += x[i] * x[i];
    }
    beta = numerator / denom;  // beta is cos(theta) , so if we now want to see how far the projection ( p = xcos(t) = xbeta) is from y, just subtract y -xb gives error.
    double rss = 0;
    for ( int i = 0; i < n; i++ ) {
        rss += ( y[i] - beta * x[i] ) * ( y[i] - beta * x[i] );
    }

    return rss;
}

// check if a vector of real numbers is probablity mass function: return false if not
bool isPmf( const vector< double > &p )
{
    double sum = 0; 
    for ( int i = 0; i < p.size(); ++i ) {
        sum += p[ i ];	
    }
   // if ( abs( sum - 1.0 ) > numeric_limits< double >::epsilon() ) return false;
    //else 
return true;
}

// sample one outcome from a multinomial distribution
int sampleMul( const gsl_rng * r, const vector< double > &p )
{
    int k = p.size();	// # of different outcomes
  
    // sample
    double u = gsl_rng_uniform( r );
    double c = 0;
    for ( int i = 0; i < k; ++i ) {
        c += p[ i ];
        if ( u < c ) { return i; }	
    }			
    
    return k - 1;
}

// sample a positive integer from a truncated geometric series: 1, r, r^2, ... r^(n - 1)
int sampleTruncatedGeometric( const gsl_rng* rng, double r, int n )
{
    double factor = ( 1.0 - r ) / ( 1.0 - pow( r, n ) );
    vector< double > probs;
    double curr = 1.0;
    for ( int i = 0; i < n; i++ ) {
        probs.push_back( curr * factor );
        curr *= r;
    }
    
    return 1 + sampleMul( rng, probs );
}

// mixture of two multinomial distributions, where the weight of the first one is "weight"
vector< double > multinomialMixture( const vector< double >& distr1, const vector< double >& distr2, double weight )
{
    int size = distr1.size();
    assert( distr2.size() == size && isPmf( distr1 ) && isPmf( distr2 ) );
    
    vector< double > mix;
    for ( int i = 0; i < size; i++ ) {
        mix.push_back( weight * distr1[ i ] + ( 1 - weight ) * distr2[ i ] );	
    }	
    
    return mix;
}

// the threshold value at a given percentage in a vector of numbers
// if highest is true, then the value at the highest p; otherwise, the value at the lowest p
double elementAt( const vector< double >& v, double p, bool highest )
{
    assert( p >= 0 && p <= 1 );
    
    // create a copy of data because we don't want to modify the original data
    vector< double > vCopy( v );

    // get the threshold
    sort( vCopy.begin(), vCopy.end() );
    int index;
    if ( highest == true ) {
        index = (int)ceil( vCopy.size() * ( 1 - p ) ) - 1;
    } else {
        index = (int)floor( vCopy.size() * p ) - 1;
    }
    if ( index < 0 ) index = 0;
    
    return vCopy[ index ];		
}

double trunc( double x )
{
    double lower = floor( x );
    double upper = ceil( x );
    if ( lower == upper ) return lower;
    
    if ( x - lower < upper - x ) return lower;
    else if ( x - lower > upper - x ) return upper;
    else return x >= 0 ? lower : upper;
}

double logit( double x )
{
    assert( x >= 0 && x <= 1 );
    
    return log( x / ( 1.0 - x ) );	
}

double inv_logit( double x )
{
    return ( exp( x ) / ( 1.0 + exp( x ) ) );
}

double logistic( double x )
{
    return ( 1.0 / ( 1.0 + exp( -x ) ) );
}

double Eucledian_dist( const vector< double >& x, const vector< double >& y )
{
    assert( x.size() == y.size() );
    
    double sum = 0;
    for ( int i = 0; i < x.size(); i++ ) {
        sum += ( x[ i ] - y[ i ] ) * ( x[ i ] - y[ i ] );	
    }	
    
    return sqrt( sum );
}

double infty_transform( double x, double a, double b )
{
    // assert x \in [a,b]
   // assert( x >= a && x <= b ); 
    
    // transformation
    double y = ( x - a ) / ( b - a ); 
    double z = log( y / ( 1.0 - y ) );

    return z;
}

double inverse_infty_transform( double z, double a, double b )
{
    assert( a <= b );

    double y = exp( z ) / ( 1.0 + exp( z ) );
    double x = a + y * ( b - a );

    return x; 
}

vector< double > weight_transform( const vector< double > w )
{
    double w0 = 1.0 - sum( w );
    vector< double > results;
    for ( int k = 0; k < w.size(); k++ ) results.push_back( log( w0 ) - log( w[ k ] ) );	
    
    return results;
}

vector< double > inverse_weight_transform( const vector< double > u )
{
    vector< double > neg_expo;
    double neg_expo_tot = 0;
    for ( int k = 0; k < u.size(); k++ ) {
        double temp = exp( -u[ k ] );
        neg_expo.push_back( temp );	
        neg_expo_tot += temp;
    }
    
    vector< double > results;
    for ( int k = 0; k < u.size(); k++ ) 
        results.push_back( neg_expo[ k ] / ( 1.0 + neg_expo_tot ) );
            
    return results;
}

void numeric_deriv( gsl_vector* grad, double (*f)( const gsl_vector*, void* ), const gsl_vector* v, void* params, double step )
{
    int n = v->size;
    double f_val = (*f)( v, params );
    
    gsl_vector* dv = gsl_vector_alloc( n );	
    for ( int i = 0; i < n; i++ ) {
        gsl_vector_memcpy( dv, v );          // ( dest, src)
        gsl_vector_set( dv, i, gsl_vector_get( v, i ) + step );
        double partial_deriv = ( (*f)( dv, params ) - f_val ) / step;
        gsl_vector_set( grad, i , partial_deriv );
    }		
}

// read parameter file into a <field, value> table. Return 0 if successful, -1 otherwise
// Format: field = value [# ...] where # indicates the start of comments (ignored)
// Ex. kappa = 2.0 # kappa: the transition/transversion bias
int readParams( const string& file, map< string, string >& params, const char* comment )
{
    params.clear();
    
    // open the parameter file
    ifstream fin( file.c_str() );
    if ( !fin ){ cerr << "Cannot open " << file << endl; exit( 1 );}
    
    // read parameters
    string line;
    while ( getline( fin, line ) ) {
        // skip a line if it starts with a comment or non-alphnumerical charcter 
        int first = line.find_first_not_of( " \t\r" );
        if ( !isalnum( line[ first ] ) || strchr( comment, line[ first ] ) ) continue;
                        
        // check if the line contains '=', if not error
        int posEq = line.find( '=' );
        if ( posEq == string::npos ) return RET_ERROR;
        
        // read the parameter name
        string name;
        int pos1 = line.find_first_not_of( " \t\r" );
        int pos2 = line.find_last_not_of( " \t\r", posEq - 1);
        if ( pos1 <= pos2 ) name = line.substr( pos1, pos2 - pos1 + 1 );
        
        // read the value, skipping the part of sentence starting with a comment symbol if any
        string value;
        pos1 = line.find_first_not_of( " \t\r", posEq + 1 );
        int posComment = line.find_first_of( comment );
        if ( posComment != string::npos ) {
            pos2 = line.find_last_not_of( " \t\r", posComment - 1 );
        } else {
            pos2 = line.find_last_not_of( " \t\r" );
        }
        if ( pos1 <= pos2 ) value = line.substr( pos1, pos2 - pos1 + 1 );
        
        // add to the table of parameters
        if ( name.size() && value.size() ) 
            params[name] = value;
        else 
            return RET_ERROR;	
    }
    
    return 0;
}

void tokenize( const string& str, vector< string >& tokens, const string& delimiters )
{
//     // Skip delimiters at beginning.
//     string::size_type lastPos = str.find_first_not_of(delimiters, 0);
//     // Find first "non-delimiter".
//     string::size_type pos  = str.find_first_of(delimiters, lastPos);

//     while (string::npos != pos || string::npos != lastPos)
//     {
//         // Found a token, add it to the vector.
//         tokens.push_back(str.substr(lastPos, pos - lastPos));
//         // Skip delimiters.  Note the "not_of"
//         lastPos = str.find_first_not_of(delimiters, pos);
//         // Find next "non-delimiter"
//         pos = str.find_first_of(delimiters, lastPos);
//     }

// 	tokens.clear();
	
    int MAX_STRING = 5000;
    char buffer[ MAX_STRING ];
    strcpy( buffer, str.c_str() );
    char* token;
    token = strtok( buffer, delimiters.c_str() );
    while ( token ) {
        tokens.push_back( string(token) );
        token = strtok( NULL, delimiters.c_str() );
    }
}

string toupperStr( const string& str )
{
    string result = str;
    transform( str.begin(), str.end(), result.begin(), ::toupper );
    
    return result;
}

string tolowerStr( const string& str )
{
    string result = str;
    transform( str.begin(), str.end(), result.begin(), ::tolower );
    
    return result;
}

void stringToVector( vector< string >& result, const string& str, const string& leftBoundary, const string& rightBoundary, const string& sep )
{
    // find the substring containing the values
    int start = str.find_first_of( leftBoundary );
    int end = str.find_last_of( rightBoundary );
    string sub = str.substr( start + 1, end - start - 1 );

    // split the substring
    tokenize( sub, result, sep );		
}

void stringToVector( vector< double >& result, const string& str, const string& leftBoundary, const string& rightBoundary, const string& sep )
{
    result.clear();
    vector< string > resultStr;
    stringToVector( resultStr, str, leftBoundary, rightBoundary, sep );
    
    for ( int i = 0; i < resultStr.size(); i++ ) {
        result.push_back( atof( resultStr[ i ].c_str() ) );		
    }
}

string IO_Ctrl::vector_delimiter = "\t";
string IO_Ctrl::map_delimiter = "\n";
string IO_Ctrl::field_separator = "\t";

int indexCompl( const vector< int >& indices, int n, vector< int >& complIndices )
{
    complIndices.clear();
    
    // bit vector of 1 to n
    vector< int > bitVect( n );
    for ( int i = 0; i < n; i++ ) bitVect[ i ] = 0;
    
    // mark all indices that are in the input set
    for ( int i = 0; i < indices.size(); i++ ) {
        if ( indices[ i ] >= n ) return RET_ERROR;
        bitVect[ indices[ i ] ] = 1;
    }
    
    // results
    for ( int i = 0; i < n; i++ ) {
        if ( bitVect[ i ] == 0 ) complIndices.push_back( i ); 		
    }
    
    return 0;
}
