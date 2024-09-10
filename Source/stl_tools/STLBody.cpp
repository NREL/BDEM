/* ============================================================================
    This file is part of BDEM.

    Class STLBody
============================================================================ */
#include "STLBody.H"
#define __deleteArrayIfNULL(arr)  if ( arr == NULL )  delete [] arr;  

using namespace BDEM;

/* ============================================================================
    Constructors and destructor
============================================================================ */

STLBody::STLBody( const amrex::ParmParse& input, std::string stlname )
:
    parameters( input ),
    stl_name( stlname ),
    num_tri( 0 ),
    ndata_per_tri( 9 ),
    ndata_per_normal ( 3 ),
    expfac( 2. ),
    nlines_per_facet ( 7 ),
    tri_pts( NULL ),
    tri_pts_vec( NULL ),
    tri_normals( NULL ),
    tri_normals_vec( NULL ),
    tri_mat( NULL ),
    tri_mat_vec( NULL ),
    domndir( NULL ),
    domndir_vec( NULL ),
    rot_tri_pts( NULL ),
    rot_tri_pts_vec( NULL ),
    sorted_indexarray( NULL ),
    sorted_indexarray_vec( NULL )
{

    readFromFile();
}    

STLBody::~STLBody()
{
    __deleteArrayIfNULL( tri_pts );
    __deleteArrayIfNULL( tri_pts_vec );
    __deleteArrayIfNULL( tri_normals );
    __deleteArrayIfNULL( tri_normals_vec );
    __deleteArrayIfNULL( tri_mat );
    __deleteArrayIfNULL( tri_mat_vec );
    __deleteArrayIfNULL( domndir );
    __deleteArrayIfNULL( domndir_vec );
    __deleteArrayIfNULL( rot_tri_pts );
    __deleteArrayIfNULL( rot_tri_pts_vec );
    __deleteArrayIfNULL( sorted_indexarray );
    __deleteArrayIfNULL( sorted_indexarray_vec );
}

/* ============================================================================
    Member functions
============================================================================ */

void STLBody::readFromFile()
{

    // Define temporaries
    std::string tmpline, tmp1, tmp2;
    int nlines( 0 );

    amrex::Real p[9], mat[9], rot_pts[9];
    int dndir;

    // TODO: allow restart!
    std::ifstream infile( stl_name.c_str() );
    amrex::Print() << "Reading STL file:" << stl_name << "\n";

    // Get number of lines in file
    std::getline( infile, tmpline ); 
    while(!infile.eof())
    {
        std::getline(infile,tmpline);
        if(tmpline.find("endsolid")!=std::string::npos)
        {
            break;
        }
        nlines++;
    }

    if( nlines % nlines_per_facet != 0 )
    {
        amrex::Abort( "There could be blank lines in the STL file\n" );
    }

    // Get number of triangles
    num_tri = nlines / nlines_per_facet;
    amrex::Print() << "Number of triangles:" << num_tri << "\n";

    // Allocate GPU arrays
    tri_pts_vec         = new amrex::Gpu::ManagedVector<amrex::Real>;
    tri_normals_vec     = new amrex::Gpu::ManagedVector<amrex::Real>;
    tri_mat_vec         = new amrex::Gpu::ManagedVector<amrex::Real>;
    domndir_vec         = new amrex::Gpu::ManagedVector<int>;
    rot_tri_pts_vec     = new amrex::Gpu::ManagedVector<amrex::Real>;

    // Adjust array lenght
    tri_pts_vec->resize( num_tri * ndata_per_tri );
    tri_normals_vec->resize( num_tri * ndata_per_normal );
    tri_mat_vec->resize( num_tri * 9 );
    rot_tri_pts_vec->resize( num_tri * 9 );
    domndir_vec->resize( num_tri );

    // Go back to the beginning of the file, let's read again
    infile.seekg( 0 );
    std::getline( infile, tmpline );

    // Read triangle by triangle
    for( int i = 0; i < num_tri; i++ )
    {
        // Normal to tri face
        std::getline( infile, tmpline );  
        std::istringstream fcnormal( tmpline );
        fcnormal >> tmp1 >> tmp2
            >> (*tri_normals_vec)[i*ndata_per_normal+0]
            >> (*tri_normals_vec)[i*ndata_per_normal+1]
            >> (*tri_normals_vec)[i*ndata_per_normal+2];

        std::getline( infile, tmpline );

        // First vertex
        std::getline( infile, tmpline );
        std::istringstream vertex1( tmpline );
        vertex1 >> tmp1
            >> (*tri_pts_vec)[i*ndata_per_tri+0]
            >> (*tri_pts_vec)[i*ndata_per_tri+1]
            >> (*tri_pts_vec)[i*ndata_per_tri+2];

        // Second vertex
        std::getline( infile, tmpline );
        std::istringstream vertex2( tmpline );
        vertex2 >> tmp1 
            >> (*tri_pts_vec)[i*ndata_per_tri+3]
            >> (*tri_pts_vec)[i*ndata_per_tri+4]
            >> (*tri_pts_vec)[i*ndata_per_tri+5];

        // Third vertex
        std::getline( infile, tmpline );
        std::istringstream vertex3( tmpline );
        vertex3 >> tmp1
            >> (*tri_pts_vec)[i*ndata_per_tri+6]
            >> (*tri_pts_vec)[i*ndata_per_tri+7]
            >> (*tri_pts_vec)[i*ndata_per_tri+8];

        std::getline(infile,tmpline);
        std::getline(infile,tmpline);

    }

    // TODO: figure out what's happening here
    for( int i = 0; i < num_tri; i++ )
    {
        for( int n = 0; n < 9; n++ )
        {
            p[n] = (*tri_pts_vec)[i*ndata_per_tri+n];
        }

        set_tri_mat_and_domndir( p, mat, rot_pts, dndir);

        for( int n=0; n<9; n++ )
        {
            (*tri_mat_vec)[i*9+n] = mat[n];
            (*rot_tri_pts_vec)[i*9+n] = rot_pts[n];
        }

        (*domndir_vec)[i]=dndir;
    }

    // Get pointers to data of AoS 
    // !Is the destructor still ok? Probably not 
    tri_pts     = tri_pts_vec->dataPtr();
    tri_normals = tri_normals_vec->dataPtr();

    domndir = domndir_vec->dataPtr();
    tri_mat  = tri_mat_vec->dataPtr();

    rot_tri_pts=rot_tri_pts_vec->dataPtr();

    update_bounding_box();

    // Initialize sorted array with ordered ids
    sorted_indexarray_vec = new amrex::Gpu::ManagedVector<int>;
    sorted_indexarray_vec->resize(num_tri);
    boxmap.clear();

    for( int i=0; i < num_tri; i++ )
    {
        (*sorted_indexarray_vec)[i]=i;
    }
    sorted_indexarray = sorted_indexarray_vec->dataPtr();
    
    orb_of_triangulation( 0, num_tri-1, 0, sorted_indexarray);

}