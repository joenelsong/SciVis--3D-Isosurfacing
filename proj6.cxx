/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
//  Function: EvaluateVectorFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//
//     rv (output): the interpolated field value. (0,0) if the location is out of bounds.
//
// ****************************************************************************

class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

unsigned int IdentifyCase(unsigned int bl, unsigned int br, unsigned int tl, unsigned int tr)
{
	br = br<<1; // ^2
	tl = tl<<2; // ^4
	tr = tr<<3; // ^8
	//printf("IdentifyCase(): bl = %d, br =%d, tl = %d, tr = %d, \n", bl, br, tl, tr);

	unsigned int casenum = bl+br+tl+tr;
	return casenum;
}
bool InterpolatePointLocation(int edgeIdx, float *A, float *B, float valA, float valB, float isoVal, float *OutputCoords)
{
	float t;

	if(edgeIdx == 0 || edgeIdx == 2) {
		t = (isoVal - valA) / (valB-valA) ;
		OutputCoords[0] = A[0] + t*(B[0]-A[0]);
		if (A[1] != B[1]) {throw 20; }
		OutputCoords[1] = (A[1]+B[1])/2; // A[1] and B[1] should be the same...
	}
	else if(edgeIdx == 1 || edgeIdx == 3 )
	{
		t = (isoVal - valA) / (valB-valA);
		if (A[0] != B[0]) {throw 20; }
		OutputCoords[0] = (A[0]+B[0])/2; // A[0] and B[0] should be the same...
		OutputCoords[1] = A[1] + t*(B[1]-A[1]);
	}
	else {
		return false; // Invalid Edge Error
	}

	return true;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!
    const int NUMSEGMENTS[16] = {0, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 0};
    // [16] = cases [4] points VALUE = {-1, 0, 1, 2, 3} where -1 means no edge and 0-3 are edges
    int LOOKUP[16][4];
    			LOOKUP[0][0] = -1;  LOOKUP[0][1] = -1; 	LOOKUP[0][2] = -1;	LOOKUP[0][3] = -1; //
    			LOOKUP[1][0] = 0;   LOOKUP[1][1] = 3; 	LOOKUP[1][2] = -1;	LOOKUP[1][3] = -1; //case 01				0,3,-1,-1}  //case 1
    			LOOKUP[2][0] = 0;   LOOKUP[2][1] = 1; 	LOOKUP[2][2] = -1;	LOOKUP[2][3] = -1; //case 02LOOKUP[2] = {				0,1,-1,-1}  //case 2
    			LOOKUP[3][0] = 1;  LOOKUP[3][1] = 3; 	LOOKUP[3][2] = -1;	LOOKUP[3][3] = -1; 
    			LOOKUP[4][0] = 2;  LOOKUP[4][1] = 3; 	LOOKUP[4][2] = -1;	LOOKUP[4][3] = -1; //				2,3,-1,-1}  //case 4
    			LOOKUP[5][0] = 0;  LOOKUP[5][1] = 2; 	LOOKUP[5][2] = -1;	LOOKUP[5][3] = -1; //			0,2,-1,-1}  //case 5
    			LOOKUP[6][0] = 0;  LOOKUP[6][1] = 1; 	LOOKUP[7][2] = 2;	LOOKUP[8][3] = 3; //			0,1,2,3}  //case 6
    			LOOKUP[7][0] = 1;  LOOKUP[7][1] = 2; 	LOOKUP[7][2] = -1;	LOOKUP[7][3] = -1; //				1,2,-1,-1}  //case 7
    			LOOKUP[8][0] = 1;  LOOKUP[8][1] = 2; 	LOOKUP[8][2] = -1;	LOOKUP[8][3] = -1; //1,2,-1,-1}  //case 8
    			LOOKUP[9][0] = 0;  LOOKUP[9][1] = 3; 	LOOKUP[9][2] = 1;	LOOKUP[9][3] = 2; //		0,3,1,2}  //case 9
    			LOOKUP[10][0] = 0;  LOOKUP[10][1] = 2; 	LOOKUP[10][2] = -1;	LOOKUP[10][3] = -1; //			1,2,-1,-1}  //case 10
    			LOOKUP[11][0] = 2;  LOOKUP[11][1] = 3; 	LOOKUP[11][2] = -1;	LOOKUP[11][3] = -1; //			2,3,-1,-1}  //case 11
    			LOOKUP[12][0] = 1;  LOOKUP[12][1] = 3; 	LOOKUP[12][2] = -1;	LOOKUP[12][3] = -1; //			1,3,-1,-1}  //case 12
    			LOOKUP[13][0] = 0;  LOOKUP[13][1] = 1; 	LOOKUP[13][2] = -1;	LOOKUP[13][3] = -1; //			0,1,-1,-1}  //case 13
    			LOOKUP[14][0] = 0;  LOOKUP[14][1] = 3; 	LOOKUP[14][2] = -1;	LOOKUP[14][3] = -1; //		0,3,-1,-1}  //case 14
    			LOOKUP[15][0] = -1;  LOOKUP[15][1] = -1; 	LOOKUP[15][2] = -1;	LOOKUP[15][3] = -1; //		-1,-1,-1,-1};//case 15
    const float ISOVALUE = 3.200000; // Given by instructions
    
    // Iterate through cells starting at 0 and going to numCells
    int x = dims[0];    int y = dims[1];
    int numCells = (x-1) * (y-1);
    float tl, tr, bl, br; // Corners topleft, topright, bottomleft, bottomright
    unsigned int blFlag, brFlag, tlFlag, trFlag; // flags = t true if corner is above ISOVALUE
    blFlag = brFlag = tlFlag= trFlag = -99; //
    
    int caseCount[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    for (j =0; j <(y-1); j++)	{ // j = y index
    	for (i =0; i <(x-1); i++ )	{ // i = x index}
    		// Get Indices for each Cell
    		bl =  F [dims[0] * j + i ];
    		br =  F [dims[0] * j + i + 1];
    		tl =  F [dims[0] * (j+1) + i];
    		tr =  F [dims[0] * (j+1) + i + 1];
    		   		
    		// Set Corner Flags for Case identification
    		blFlag = (bl >= ISOVALUE) ? 1 : 0;
    		brFlag = (br >= ISOVALUE) ? 1 : 0;
    		tlFlag = (tl >= ISOVALUE) ? 1 : 0;
    		trFlag = (tr >= ISOVALUE) ? 1 : 0;
    		
    		// Look up table for precomputed answers to all cases
    		int icase = IdentifyCase(blFlag, brFlag, tlFlag, trFlag);
		//if (brFlag && !blFlag && !tlFlag && !trFlag) {
    			//printf("Vertice Flags: %d %d %d %d\n", blFlag, brFlag, tlFlag, trFlag);
    			//printf("icase = %d\n", icase);
    		//}
    		
    		
	//if (icase != 10)	{continue; } // Test Case by Cases
    		// Working: { 0, 15 }
    		// Probably Working: {1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14 }
    		// Not Working: 
						//   10- Lines OOB mostly right though}
						
    		caseCount[icase]++;

		    int nsegements = NUMSEGMENTS[icase];

		    int edge1 = -99;
		    //int edge2 = -99;
		    float pt1[2][2]; // interpolated position on edge 1
		    //float pt2[2]; // interpolated position on edge 2

		    // if there are no segments then continue to the next iteration (i.e. case 0 & 16)
	    //printf("nsegements %d\n", nsegements);
		    if (nsegements < 0)	{	
		    	continue;	
		    }
		    int s;
		    for (s=0; s<nsegements; s++)
		    {
		    	float A[2];
		    	float B[2];
		    	float valA =-99;
		    	float valB =-99;

		    	int p;
		    	for (p =0; p <2; p++) // Iterate Twice to Interpolate both points so we can draw a line segement
		    	{


			    	edge1 = LOOKUP[icase][2*s+p];
		    	//printf("edge1 = %d\n", edge1);
			    	switch (edge1) {
			    		case 0:
			    			A[0] = X[i];
			    			A[1] = Y[j];
			    			B[0] = X[i+1];
			    			B[1] = Y[j];
			    			valA=bl;
			    			valB=br;
			    			break;
			    		case 1:
			    			A[0] = X[i+1];
			    			A[1] = Y[j];
			    			B[0] = X[i+1];
			    			B[1] = Y[j+1];
			    			valA=br;
			    			valB=tr;
			    			break;
			    		case 2:
			    			A[0] = X[i];
			    			A[1] = Y[j+1];
			    			B[0] = X[i+1];
			    			B[1] = Y[j+1];
			    			valA=tl;
			    			valB=tr;
			    			break;
			    		case 3:
			    			A[0] = X[i];
			    			A[1] = Y[j];
			    			B[0] = X[i];
			    			B[1] = Y[j+1];
			    			valA=bl;
			    			valB=tl;
			    			break;
			    	}

		    		InterpolatePointLocation(edge1, A, B, valA, valB, ISOVALUE, pt1[p]);
		    	}

		    	//printf("Interpolate Point on Edge %d Between (%f,%f) and (%f,%f) = (%f,%f)\n",edge1, A[0],A[1],B[0],B[1],pt1[0],pt1[1]);
		    	//printf("");
		    	//printf("AddSegment((%f, %f) to (%f, %f)\n)",pt1[0][0], pt1[0][1], pt1[1][0], pt1[1][1]);
		    	sl.AddSegment(pt1[0][0], pt1[0][1], pt1[1][0], pt1[1][1]);
		    	//sl.AddSegment(-10, -10, +10, -10);
		    }

    		// ~TestCode- should create diagonal lines from bottom left to top right of every cell
    		//printf("bl=%f4, br=%f4, tl=%f4, tr=%f4",X[i], Y[j], X[i+1], Y[j+1] );
    		//sl.AddSegment(X[i], Y[j], X[i+1], Y[j+1]);

    	}
    }
    /*
    printf("Case Counts\n");
    printf("0 %d\n",caseCount[0]);
    printf("1 %d\n",caseCount[1]); 
	printf("2 %d\n",caseCount[2]); 
	printf("3 %d\n",caseCount[3]);
	printf("4 %d\n",caseCount[4]);
	printf("5 %d\n",caseCount[5]); 
	printf("6 %d\n",caseCount[6]);
	printf("7 %d\n",caseCount[7]);
	printf("8 %d\n",caseCount[8]); 
	printf("9 %d\n",caseCount[9]); 
	printf("10 %d\n",caseCount[10]);
	printf("11 %d\n",caseCount[11]);
	printf("12 %d\n",caseCount[12]);
	printf("13 %d\n",caseCount[13]); 
	printf("14 %d\n",caseCount[14]); 
	printf("15 %d\n",caseCount[15]);
	*/

    



    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
