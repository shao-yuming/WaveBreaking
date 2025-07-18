/**
Various output functions for extracting a surface from field data using the build in VOF PLIC surface or the isosurface reconstruction in bview
Three different formats are supported:
- .ply
- .vtu
- .vtp
all of which can be loaded into Paraview. for the two later alternatives, field data may be added by using the functions with "_w_fielddata".
This is an updated version of the previous fractions_output.h


*/

#include "geometry.h"
#include "fractions.h"

#if dimension == 1
coord mycs (Point point, scalar c) {
    coord n = {1.};
    return n;
}
#elif dimension == 2
# include "myc2d.h"
#else // dimension == 3
# include "myc.h"
#endif

struct OutputFacets_scalar {
    scalar c;
    FILE * fp;     // optional: default is stdout
    //scalar * list;  // List of scalar fields to include when writing vtu surface to file
    vector * vlist; // List of vector fields to include.
    face vector s; // optional: default is none
};


/**
 Outputs ISO surface with fielddata in .vtp format 
 */
void output_vtp (struct OutputFacets_scalar p)
{
    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
    scalar c = p.c;
    //face vector s = p.s;
    if (!p.fp) p.fp = stdout;

    // print header text
    fputs ("<?xml version=\"1.0\"?>\n", p.fp);
    fputs ("<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n", p.fp);
    fputs ("\t<PolyData>\n", p.fp);

    // Start by creating the vertex and smoothed normal field
    vertex scalar v[];
    foreach_vertex()
        v[] = (c[] + c[-1] + c[0,-1] + c[-1,-1] +
        c[0,0,-1] + c[-1,0,-1] + c[0,-1,-1] + c[-1,-1,-1])/8.;


    /** Loop through all surface cells
     *       The point of this first round is to count the number of isosurface triangles. Should be improved...
     */
    int nverts = 0;
    int nfacets = 0;
    //shao change
    //foreach(reduction(+:nfacets) reduction(+:nverts)) {
     foreach(serial) {
        //if (c[] > 1e-7 && c[] < 1. - 1e-7) {
        double val[8] = {
            v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
            v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
        };
        double t[5][3][3];
        int nt = polygonize (val, 0.5, t);
        nfacets += nt;
        nverts += nt*3;
    }

    fprintf (p.fp, "\t\t<Piece NumberOfPoints=\"%i\" NumberOfPolys=\"%i\">\n", nverts, nfacets);

    fputs ("      <CellData>\n", p.fp);
    fputs ("      </CellData>\n", p.fp);

    // Write list of scalar field values to file

    fputs ("\t\t\t <PointData Normals=\"Normals\">\n", p.fp);
/*
    for (scalar s in p.list) {
        fprintf (p.fp,"\t\t\t\t <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", s.name);
        foreach() {
	      	// Rearranging v[]
	      	double val[8] = {
	      		v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
	      		v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
	     	 };
	      	double t[5][3][3];
	      	int nt = polygonize (val, 0.5, t);
	      	for (int i = 0; i < nt; i++) {
		      	for (int j = 0; j < 3; j++) {
	      	  		coord v = {t[i][j][0], t[i][j][1], t[i][j][2]};
	      	  		fprintf (p.fp, "%g\n", interp (point, v, s));
	      		}
	      	}
    	}
        fputs ("\t\t\t\t </DataArray>\n", p.fp);
    }
*/
    // Write list of vector field values to file
    for (vector ve in p.vlist) {
        fprintf (p.fp,"\t\t\t\t <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", ve.x.name);
        foreach(serial) {
	      	// Rearranging v[]
	      	double val[8] = {
	      		v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
	      		v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
	     	 };
	      	double t[5][3][3];
	      	int nt = polygonize (val, 0.5, t);
	      	for (int i = 0; i < nt; i++) {
		      	for (int j = 0; j < 3; j++) {
	      	  		coord v = {t[i][j][0], t[i][j][1], t[i][j][2]};
	      	  		#if dimension == 2
	                fprintf (p.fp, "%g %g 0.\n", interp (point, v, ve.x), interp (point, v, ve.y));
	                #endif
	                #if dimension == 3
	                fprintf (p.fp, "%g %g %g\n", interp (point, v, ve.x), interp (point, v, ve.y), interp (point, v, ve.z));
	                #endif
	      	  		//fprintf (p.fp, "%g\n", interp (point, v, s));
	      		}
	      	}
    	}
        fputs ("\t\t\t\t </DataArray>\n", p.fp);
    }
    fputs ("\t\t\t </PointData>\n", p.fp);


    // Write points to file
    fputs ("      <Points>\n", p.fp);
    fputs ("        <DataArray type=\"Float32\" Name=\"vertices\" NumberOfComponents=\"3\" format=\"ascii\">\n", p.fp);


    foreach(serial){
        //if (c[] > 1e-7 && c[] < 1. - 1e-7) {

        double val[8] = {
            v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
            v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
        };
        double t[5][3][3];
        int nt = polygonize (val, 0.5, t);
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < 3; j++) {
                //coord v = {t[i][j][0], t[i][j][1], t[i][j][2]}, np;
                //shao change
                 coord v = {t[i][j][0], t[i][j][1], t[i][j][2]};
                fprintf (p.fp, "%g %g %g\n",
                         x + v.x*Delta, y + v.y*Delta, z + v.z*Delta);
            }
        }
    }

    fputs ("        </DataArray>\n", p.fp);
    fputs ("      </Points>\n", p.fp);
    fputs ("      <Polys>\n", p.fp);

    fputs ("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n", p.fp);

    // print vert numbers
    for (int ivert = 0; ivert < nverts; ivert++)
        fprintf (p.fp, "%i ", ivert);

    fputs ("        </DataArray>\n", p.fp);
    fputs ("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n", p.fp);

    // print offsets
    for (int ifacet = 0; ifacet < nfacets; ifacet++)
        fprintf (p.fp, "%i ", ifacet*3+3);


    fputs ("        </DataArray>\n", p.fp);
    fputs ("      </Polys>\n", p.fp);
    fputs ("    </Piece>\n", p.fp);
    fputs ("  </PolyData>\n", p.fp);
    fputs ("</VTKFile>\n", p.fp);

    fflush (p.fp);
    #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
    #endif
}

//void output_pvtp(scalar * list, vector * vlist,  FILE * fp, char * subname)
void output_pvtp(vector * vlist,  FILE * fp, char * subname)
{

    fputs("<?xml version=\"1.0\"?>\n", fp);
    fputs ("<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs("\t<PPolyData GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PPointData Normals=\"Normals\">\n", fp);
   /* for(scalar s in list){
    fprintf (fp,"\t\t\t\t <PDataArray type=\"Float32\" Name=\"%s\"/>\n", s.name);
    }*/
    for (vector v in vlist)
    {
      fprintf(fp, "\t\t\t\t <PDataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"%s\"/>\n", v.x.name);
    }
    fputs ("\t\t\t </PPointData>\n", fp);  
    fputs ("\t\t\t <PPoints>\n", fp);
    fputs ("\t\t\t\t <PDataArray type=\"Float32\" Name=\"vertices\" NumberOfComponents=\"3\"/>\n", fp);
    fputs ("\t\t\t </PPoints>\n", fp);
    for (int i = 0; i < npe(); i++)
      fprintf(fp, "\t\t <Piece Source=\"%s_n%3.3d.vtp\"/> \n", subname, i);
    fputs("\t </PPolyData>\n", fp);
    fputs("</VTKFile>\n", fp);

}

