//      DICMesh.cxx
//      
//      Copyright 2011 Seth Gilchrist <seth@mech.ubc.ca>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.

#ifndef DICMESH_H
#define DICMESH_H

#include <cstring>
#include <ctime>
#include "DIC.cxx"
#include "itkMesh.h"
#include "itkTetrahedronCell.h"
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkQuadraticTetra.h>
#include <vtkQuadraticTriangle.h>
#include <vtkQuadraticEdge.h>
#include <vtkVertex.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkPKdTree.h>
#include <vtkIdList.h>
#include <vtkCellDerivatives.h>
#include <vtkCellDataToPointData.h>
#include <vtkMath.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageData.h>

#include <itkImageFileWriter.h>
#include <itkShrinkImageFilter.h>


template<typename TFixedImage, typename TMovingImage>
class DICMesh : public DIC<TFixedImage, TMovingImage>
{
public:

/** Inherited typedefs from DIC.*/
typedef typename DIC<TFixedImage,TMovingImage>::FixedImageType				FixedImageType;
typedef typename DIC<TFixedImage,TMovingImage>::MovingImageType				MovingImageType;
typedef typename DIC<TFixedImage,TMovingImage>::FixedImagePointer			FixedImagePointer;
typedef typename DIC<TFixedImage,TMovingImage>::MovingImagePointer			MovingImagePointer;
typedef typename DIC<TFixedImage,TMovingImage>::FixedImageType::IndexType 	FixedImageIndexType;
typedef typename DIC<TFixedImage,TMovingImage>::MovingImageType::IndexType	MovingImageIndexType;
typedef	typename DIC<TFixedImage,TMovingImage>::FixedImageRegionType		FixedImageRegionType;
typedef typename DIC<TFixedImage,TMovingImage>::MovingImageRegionType		MovingImageRegionType;
typedef typename DIC<TFixedImage,TMovingImage>::FixedImageRegionListType	FixedImageRegionListType;
typedef	typename DIC<TFixedImage,TMovingImage>::MovingImageRegionListType	MovingImageRegionListType;

typedef typename DIC<TFixedImage,TMovingImage>::ImageRegistrationMethodType	ImageRegistrationType;
typedef typename ImageRegistrationType::ParametersType							RegistrationParametersType;	
typedef typename DIC<TFixedImage,TMovingImage>::TransformInitializerType	TransformInitializerType;
typedef typename DIC<TFixedImage,TMovingImage>::TransformType				TransformType;

typedef		vtkSmartPointer<vtkUnstructuredGrid>	DataImagePointer;
typedef		vtkSmartPointer<vtkPoints>				DataImagePointsPointer;
typedef		vtkSmartPointer<vtkDoubleArray>			DataImagePixelPointer;



/* Methods. **/

/** Constructor **/
DICMesh()
{
	m_DataImage = 0; // must be provided by user or read in using the ReadMeshFromGmshFile method
	m_KDTree = vtkSmartPointer<vtkPKdTree>::New();
	m_errorRadius = 4; // initialize the error search radius to 3 image units
	m_errorTolerance = 2; // difference of a pixel from its neighbours to be considered erronious, in standard deviations from the mean
	m_pointsList = vtkSmartPointer<vtkIdList>::New(); // the points list for analysis
	m_maxMeticValue = -0.00; // TODO: make this setable using a method
	m_GlobalRegDownsampleValue = 3; // This value is the default downsample when preforming the global registration.
}

/** Destructor **/
~DICMesh() {}

/** This function will compile a list of valid fixed image regions.*/
void CalculateInitialMovingImageRegionList()
{
	vtkIdType	numberOfNodes = m_DataImage->GetNumberOfPoints();
	this->m_MovingImageRegionList.clear();
	
	for ( int i = 0; i < numberOfNodes; ++i){ // Iterate through the points in the point set
		double *currentLocation = new double[3];
		currentLocation = this->CalculateFixedRegionLocationFromIndex( i );
		//~ m_DataImage->GetPoint( i, currentPoint ); // get the current node location
		
		MovingImageRegionType *currentRegion = new MovingImageRegionType;
		this->GetMovingImageRegionFromLocation( currentRegion, currentLocation ); // get the current region
		this->PushRegionOntoMovingImageRegionList( currentRegion );
	}
}

double *CalculateMovingImageRegionLocationFromIndex( vtkIdType i)
{
	double *currentPoint	= new double[3];
	double *currentValue	= new double[3];
	double *currentLocation	= new double[3];
	m_DataImage->GetPoint( i, currentPoint ); // get the current point location
	this->m_DataImage->GetPointData()->GetArray("Displacement")->GetTuple( i, currentValue );
	*currentLocation		= *currentPoint		+ *currentValue; // set the current point to the original local to the disp
	*(currentLocation +1)	= *(currentPoint +1)+ *(currentValue +1);
	*(currentLocation +2)	= *(currentPoint +2)+ *(currentValue +2);
	
	return currentLocation;
}

/** A function to calculate the initial fixed image region list.**/
void CalculateInitialFixedImageRegionList()
{
	vtkIdType	numberOfNodes = m_DataImage->GetNumberOfPoints();
	this->m_FixedImageRegionList.clear();
	this->m_pointsList->Reset();
	
	for (int i = 0; i < numberOfNodes; ++i){
		double *currentLocation = new double[3];
		m_DataImage->GetPoint( i, currentLocation );
		//~ currentLocation = this->CalculateFixedRegionLocationFromIndex( i );
		
		FixedImageRegionType	*currentRegion = new FixedImageRegionType;
		this->GetFixedImageRegionFromLocation( currentRegion, currentLocation );
		this->PushRegionOntoFixedImageRegionList( currentRegion );
		this->m_pointsList->InsertNextId( i );
	}
}

/** A function that will calculate the fixed image region initial
 * location from an index. */
double* CalculateFixedRegionLocationFromIndex( vtkIdType i )
{
	double *currentPoint	= new double[3];
	double *currentValue	= new double[3];
	double *currentLocation	= new double[3];
	m_DataImage->GetPoint( i, currentPoint ); // get the current point location
	this->m_DataImage->GetPointData()->GetArray("Displacement")->GetTuple( i, currentValue );
	*currentLocation		= *currentPoint		+ *currentValue; // set the current point to the original local to the disp
	*(currentLocation +1)	= *(currentPoint +1)+ *(currentValue +1);
	*(currentLocation +2)	= *(currentPoint +2)+ *(currentValue +2);
	
	return currentLocation;
}

/** A function to read a gmsh file. */
void ReadMeshFromGmshFile( std::string gmshFileName )
{
	std::stringstream msg("");
	
	std::ifstream gmshFileInput(gmshFileName.c_str());
	//open file for reading
	if (!gmshFileInput){
		msg.str(" ");
		msg << "Cannot open Gmsh file for reading."<<std::endl <<
			"Please check the filename and try again."<<std::endl;
		this->WriteToLogfile( msg.str() );
		std::exit(1);
	}
	char str[255];
	gmshFileInput.getline(str,255);
	
	// read in the nodes /
	while ( strcmp(str,"$Nodes") ){
		gmshFileInput.getline(str,255); // skip all lines until we get to "$Nodes"
	}
	unsigned int numberOfNodes;
	gmshFileInput >> numberOfNodes; // the next line will be the number of nodes
	
	// create the data image
	DataImagePointsPointer		points				= DataImagePointsPointer::New();
	DataImagePointer			meshImage			= DataImagePointer::New();
	DataImagePixelPointer		displacementData	= DataImagePixelPointer::New();
	DataImagePixelPointer		optimizerData		= DataImagePixelPointer::New();
	
	meshImage->SetPoints( points );
	
	displacementData->SetNumberOfComponents(3);
	displacementData->SetName("Displacement");
	meshImage->GetPointData()->AddArray( displacementData );
	
	optimizerData->SetNumberOfComponents(1);
	optimizerData->SetName("Optimizer Value");
	meshImage->GetPointData()->AddArray( optimizerData );
	
	meshImage->GetPoints()->SetNumberOfPoints( (vtkIdType)numberOfNodes );
	for (unsigned int i = 0; i < numberOfNodes; ++i){ // run through the number of nodes and put them in the pointset
		unsigned int pointNumber;
		float x, y, z;
		gmshFileInput >> pointNumber >> x >> y >> z;
		meshImage->GetPoints()->SetPoint(i,x,y,z);
	}
	
	// read in the elements. /
	while ( strcmp(str,"$Elements") ){
		gmshFileInput.getline(str,255); // skip all lines until we get to "$Elements"
	}
	unsigned int numberOfElements;
	gmshFileInput >> numberOfElements; // the next line will be the number of elements
	
	for ( unsigned int i = 0; i < numberOfElements; ++i ){
		unsigned int elNo, elTyp;
		gmshFileInput >> elNo >> elTyp;
		
		/*if (elTyp == 1){ // 2 node line
			unsigned int elNTags;
			gmshFileInput >> elNTags;
			for (unsigned int j = 0; j < elNTags+1; ++j){
				gmshFileInput.ignore(255,' ');
			}
			int nodes[2];
			gmshFileInput >> nodes[0] >> nodes[1];
			vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1 };
			meshImage->InsertNextCell(3, 2, ptIds);
		}*/
		
		/*if (elTyp == 2){ // 3 node triangle
			unsigned int elNTags;
			gmshFileInput >> elNTags;
			for (unsigned int j = 0; j < elNTags+1; ++j){
				gmshFileInput.ignore(255,' ');
			}
			int nodes[3];
			gmshFileInput >> nodes[0] >> nodes[1] >> nodes[2];
			vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1, nodes[2]-1 };
			meshImage->InsertNextCell(5, 3, ptIds);
		}*/
		
		if (elTyp == 4){ // linear tet element
			unsigned int elNTags;
			gmshFileInput >> elNTags;
			for (unsigned int j = 0;j < elNTags+1; ++j){ // ignore tags
				gmshFileInput.ignore(255,' ');
			}
			int nodes[4];
			gmshFileInput >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3]; 
			vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1, nodes[2]-1, nodes[3]-1 };
			meshImage->InsertNextCell(10, 4, ptIds);
		}
		
		/*if (elTyp == 8){ // quadratic line element
			std::stringstream msg("");
			msg << "Warning: Qadratic line element in gmsh file."<<std::endl << "Qadratic lines are not supported" <<std::endl <<
				"Continuing with the file input."<<std::endl; // not supported because they give bad jaboians in the strain calcs
			this->WriteToLogfile( msg.str() );
			gmshFileInput.ignore(255,'\n');
			continue;
		}*/
		
		/*if (elTyp == 9){ //quadratic triangle element
			unsigned int elNTags;
			gmshFileInput >> elNTags;
			for (unsigned int j = 0; j< elNTags+1; ++j){ //ignore tags
				gmshFileInput.ignore(255,' ');
			}
			int nodes[6];
			gmshFileInput >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5];
			vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1, nodes[2]-1, nodes[3]-1, nodes[4]-1, nodes[5]-1};
			meshImage->InsertNextCell(22,6,ptIds);
		}*/
		
		if (elTyp == 11){ // quatratic tet element
			unsigned int elNTags;
			gmshFileInput >> elNTags;
			for (unsigned int j = 0;j < elNTags+1; ++j){ // ignore tags
				gmshFileInput.ignore(255,' ');
			}
			int nodes[10];//, node1, node2, node3, node4, node5, node6, node7, node8, node9;
			gmshFileInput >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5] >> nodes[6] >> nodes[7] >> nodes[9] >> nodes[8]; // nodes[8] and modes[9] intensionally switched due to vtk/gmsh node numbering differences
			vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1, nodes[2]-1, nodes[3]-1, nodes[4]-1, nodes[5]-1, nodes[6]-1, nodes[7]-1, nodes[8]-1, nodes[9]-1 };
			meshImage->InsertNextCell(24, 10, ptIds);
		}
		
		/*if (elTyp == 15){ // 3D vertex
			unsigned int elNTags;
			gmshFileInput >> elNTags;
			for (unsigned int j = 0; j< elNTags+1; ++j){ //ignore tags
				gmshFileInput.ignore(255,' ');
			}
			int nodes[1];
			gmshFileInput >> nodes[0];
			vtkIdType ptIds[] = { nodes[0]-1 };
			meshImage->InsertNextCell(1,1,ptIds);
		}*/
		
		/*if (elTyp != 1 && elTyp != 2 && elTyp != 4 && elTyp != 8 && elTyp != 9 && elTyp != 11 && elTyp != 15){ // only support the above element types
			std::stringstream msg("");
			msg << "Unsupported element types detected!"<<std::endl << "Unsupported type = "<<elTyp<<std::endl <<
				"Continuing with the file input."<<std::endl;
			this->WriteToLogfile( msg.str() );
			gmshFileInput.ignore(255,'\n');
			continue;
		}*/
		if (elTyp != 4 && elTyp != 11){ // only support the above element types
			std::stringstream msg("");
			msg << "Unsupported element types detected!"<<std::endl << "Unsupported type = "<<elTyp<<std::endl <<
				"Continuing with the file input."<<std::endl;
			this->WriteToLogfile( msg.str() );
			gmshFileInput.ignore(255,'\n');
			continue;
		}	
	}
		
	for ( unsigned int i = 0; i < numberOfNodes; ++i ){
		meshImage->GetPointData()->GetArray("Displacement")->InsertNextTuple3( 0, 0, 0 );
		meshImage->GetPointData()->GetArray("Optimizer Value")->InsertNextTuple1( 0 );
	}
	
	this->SetDataImage( meshImage );
}

/** A function to fill a mesh with an single value. */
void SetMeshInitialValue( double initialData[] )
{
	vtkSmartPointer<vtkDataArray> pointData = this->m_DataImage->GetPointData()->GetArray("Displacement"); // get the data array by name
	vtkIdType numberOfNodes = pointData->GetNumberOfTuples();
	
	for ( unsigned int i = 0; i < numberOfNodes; ++i ){
		pointData->SetTuple( i , initialData );
	}

}

/** Get a pixel value from the mesh using an index */
void GetMeshPixelValueFromIndex( vtkIdType index, double *pixel )
{
	this->m_DataImage->GetPointData()->GetArray("Displacement")->GetTuple( index, pixel );
}

/** Get the value of the optimizer at a certain point. **/
void GetMeshPixelOptimizerFromIndex( vtkIdType index, double *opt )
{
	this->m_DataImage->GetPointData()->GetArray("Optimizer Value")->GetTuple( index, opt );
}

/** Set a pixel value in the mesh using an index*/
void SetMeshPixelValueFromIndex(vtkIdType index, double *pixel )
{
	this->m_DataImage->GetPointData()->GetArray("Displacement")->SetTuple( index, pixel );
}

/** Set the value of the optimizer at a certain point. **/
void SetMeshPixelOptimizerFromIndex( vtkIdType index, double *opt )
{
	this->m_DataImage->GetPointData()->GetArray("Optimizer Value")->SetTuple( index, opt );
}

/** Get a point by index from the mesh.**/
void GetMeshPointLocationFromIndex( vtkIdType index, double  *point )
{
	this->m_DataImage->GetPoint(index,point);
}

/** Set the data Image. **/
void SetDataImage( vtkUnstructuredGrid *initialDataImage )
{
	if (this->m_DataImage.GetPointer() != initialDataImage){
		this->m_DataImage = initialDataImage;
	}
	this->KDTreeSetAndBuild();
}

/** get the data image. **/
DataImagePointer GetDataImage()
{
	return this->m_DataImage;
}

/** a function to execute the analysis. */
void ExecuteDIC()
{
	std::stringstream msg("");
	if( !this->m_Registration ){
		msg.str("");
		msg << "Registration not set.  Please define and set the\nregistration using the SetRegistrationMethod() method.\n";
		this->WriteToLogfile( msg.str() );
		std::abort();
	}
	if( !this->m_FixedImage ){
		msg.str("");
		msg << "Fixed image not set.  Please define and set the\n fixed image using the SetFixedImage() method.\n";
		this->WriteToLogfile( msg.str() );
		std::abort();
	}
	if( !this->m_MovingImage ){
		msg.str("");
		msg << "Moving image not set.  Please define and set the\n moving image using the SetMovingImage() method.\n";
		this->WriteToLogfile( msg.str() );
		std::abort();		
	}
	
	if( !this->m_DataImage ){
		msg.str("");
		msg << "The data image not set.  Please define and set the\n initial data image using either the ReadMeshFromGmshFile or\n SetDataImage() methods.\n";
		this->WriteToLogfile( msg.str() );
		std::abort();		
	}
	
	if( this->m_FixedImageRegionList.empty() ){
		this->CalculateInitialFixedImageRegionList();
	}
	if( this->m_MovingImageRegionList.empty() ){
		this->CalculateInitialMovingImageRegionList();
	}
	
	unsigned int nMeshPoints = this->m_pointsList->GetNumberOfIds();
	struct tm * timeValue;
	std::time_t rawTime;
	std::time( &rawTime );
	timeValue = std::localtime( &rawTime );
	msg.str("");
	msg << "Starting DIC at: "<<std::asctime( timeValue )<<std::endl;
	this->WriteToLogfile( msg.str() );
	for( unsigned int i = 0; i<nMeshPoints; ++i){
		vtkIdType pointId = this->m_pointsList->GetId( i );
				
		msg.str("");
		msg << "Starting image registraion for point: "<<i+1<<" of "<<nMeshPoints<<" (mesh index "<<pointId<<")"<<std::endl;
		this->WriteToLogfile( msg.str() );
				
		std::time_t rawTime;
		struct tm * timeValue;
		std::time( &rawTime );
		timeValue = std::localtime( &rawTime );
		msg.str("");
		msg << "Time: "<<std::asctime( timeValue );
		this->WriteToLogfile( msg.str() );
		

		FixedImageRegionType	*fixedRegion = new FixedImageRegionType ;
		FixedImagePointer		fixedImage;
		fixedRegion = this->GetFixedImageRegionFromIndex( i );
		fixedImage = this->GetFixedROIAsImage( fixedRegion );
		this->SetFixedROIImage( fixedImage );
		
		MovingImageRegionType	*movingRegion = new MovingImageRegionType;
		MovingImagePointer		movingImage;
		movingRegion = this->GetMovingImageRegionFromIndex( i );
		movingImage = this->GetMovingROIAsImage( movingRegion );
		this->SetMovingROIImage( movingImage );
		
		this->SetTransformToIdentity();
		
		this->m_TransformInitializer->SetFixedImage( this->m_Registration->GetFixedImage() );
		this->m_TransformInitializer->SetMovingImage( this->m_Registration->GetMovingImage() );
		this->m_TransformInitializer->SetTransform( this->m_Transform );
		this->m_TransformInitializer->InitializeTransform();
		this->m_Registration->SetInitialTransformParameters( this->m_Transform->GetParameters() );
			
		double	*displacementData = new double[3];
		this->GetMeshPixelValueFromIndex( pointId, displacementData );
		this->SetInitialDisplacement( displacementData );
		
		msg.str("");
		msg <<"Current transform: "<<this->m_Registration->GetInitialTransformParameters()<<std::endl;
		this->WriteToLogfile( msg.str() );
		
		if( !strcmp(this->m_Registration->GetOptimizer()->GetNameOfClass(),"LBFGSBOptimizer") ){
			itk::Array< long > boundSelect( this->m_Registration->GetTransform()->GetNumberOfParameters() );
			itk::Array< double > upperBound( this->m_Registration->GetTransform()->GetNumberOfParameters() );
			itk::Array< double > lowerBound( this->m_Registration->GetTransform()->GetNumberOfParameters() );
			
			boundSelect.Fill( 0 );
			boundSelect[6] = 2;
			boundSelect[7] = 2;
			boundSelect[8] = 2;
			upperBound.Fill( 0 );
			typename FixedImageType::SpacingType imageSpacing = this->m_Registration->GetFixedImage()->GetSpacing();
			upperBound[6] = *displacementData + 2*imageSpacing[0];
			upperBound[7] = *(displacementData+1) + 2*imageSpacing[1];
			upperBound[8] = *(displacementData+2) + 2*imageSpacing[2];
			lowerBound.Fill( 0 );
			lowerBound[6] = *displacementData - 2*imageSpacing[0];
			lowerBound[7] = *(displacementData+1) - 2*imageSpacing[1];
			lowerBound[8] = *(displacementData+2) - 2*imageSpacing[2];
			
			reinterpret_cast<itk::LBFGSBOptimizer *>(this->m_Registration->GetOptimizer())->SetBoundSelection( boundSelect );
			reinterpret_cast<itk::LBFGSBOptimizer *>(this->m_Registration->GetOptimizer())->SetUpperBound( upperBound );
			reinterpret_cast<itk::LBFGSBOptimizer *>(this->m_Registration->GetOptimizer())->SetLowerBound( lowerBound );
		}
			
		//~ this->m_Metric->ReinitializeSeed();
		std::cout<<"Number of Fixed Image Samples: "<<this->m_Metric->GetNumberOfPixelsCounted()<<std::endl;

			
		this->UpdateRegionRegistration();
		
		double *lastDisp = new double[3];
		this->GetLastDisplacement( lastDisp );
		this->SetMeshPixelValueFromIndex( pointId, lastDisp );
		double *lastOpt = new double;
		this->GetLastOptimizer( lastOpt );
		this->SetMeshPixelOptimizerFromIndex( pointId, lastOpt );
		msg.str("");
		msg << "Final displacement value: ("<<*lastDisp<<", "<<*(lastDisp+1)<<", "<<*(lastDisp+2)<<")"<<std::endl <<
			"Final optimizer value: "<<*lastOpt<<std::endl<<std::endl;
		this->WriteToLogfile( msg.str() );
		std::string debugFile = this->m_OutputDirectory + "/debug.vtk";
		this->WriteMeshToVTKFile( debugFile );
	}
}

/** A function to write the mesh data to a VTK ASCII file. **/
void WriteMeshToVTKFile(std::string outFile)
{
	vtkSmartPointer<vtkUnstructuredGridWriter>	writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetFileName( outFile.c_str() );
	writer->SetInput( this->m_DataImage );
	writer->Update();
}

/** A function to build the KDTree locator from the data image.  This method
 * is called when the data image is set. */
void KDTreeSetAndBuild()
{
	this->m_KDTree->SetDataSet( this->m_DataImage );
	this->m_KDTree->BuildLocatorFromPoints(this->m_DataImage->GetPoints());
}

/** A function to find the values that are not the same as the neighbours. */
void CreateNewRegionListFromBadPixels()
{	
	this->m_FixedImageRegionList.clear(); // clear the regions for calculation
	this->m_MovingImageRegionList.clear();
	this->m_pointsList->Reset();
	unsigned int numberOfPoints = this->m_DataImage->GetNumberOfPoints();
	for( unsigned int i = 0; i < numberOfPoints; ++i ){
		// find points within radius
		double *currentPoint = new double[3];
		this->GetMeshPointLocationFromIndex( i, currentPoint );
		vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New();
		this->m_KDTree->FindPointsWithinRadius( this->m_errorRadius, currentPoint, pointList );
		// calc the average and stDev of the points in the list
		double *averageValue = new double[3];
		double *averageMag = new double;
		double *stDevValue = new double[3];
		double *stDevMag = new double;
		this->CalculateStats(i, pointList, averageValue, averageMag, stDevValue, stDevMag );
		// get the value of the current point
		double *currentValue = new double[3];
		this->GetMeshPixelValueFromIndex( i, currentValue );
		double *currentOptimizer = new double;
		this->GetMeshPixelOptimizerFromIndex( i, currentOptimizer);
		
		// check if there the pixel is valid - if any component is out of value it is considered invalid
		if( std::fabs(*currentValue - *averageValue) > this->m_errorTolerance*(*stDevValue) ||
			std::fabs(*(currentValue+1) - *(averageValue+1)) > this->m_errorTolerance*(*(stDevValue+1)) ||
			std::fabs(*(currentValue+2) - *(averageValue+2)) > this->m_errorTolerance*(*(stDevValue+2)) ||
			*currentOptimizer > this->m_maxMeticValue ){
			
			MovingImageRegionType *currentMovingRegion = new MovingImageRegionType; // New Moving Region
			double *currentPointLocation = this->m_DataImage->GetPoint( i );
			
			double *movingImageCenterLocation = new double[3];
			*movingImageCenterLocation = *currentPointLocation + *averageValue;
			*(movingImageCenterLocation+1) = *(currentPointLocation+1) + *(averageValue+1);
			*(movingImageCenterLocation+2) = *(currentPointLocation+2) + *(averageValue+2);
			this->GetMovingImageRegionFromLocation( currentMovingRegion, movingImageCenterLocation );
			this->PushRegionOntoMovingImageRegionList( currentMovingRegion );
			
			FixedImageRegionType *currentFixedRegion = new FixedImageRegionType;  // New Fixed Region
			this->GetFixedImageRegionFromLocation( currentFixedRegion, currentPointLocation );
			this->PushRegionOntoFixedImageRegionList( currentFixedRegion );
			this->m_pointsList->InsertNextId( i );
			
			this->SetMeshPixelValueFromIndex( i, averageValue );  // apply pixel value to the mesh
		}
	}
}

/** A function to calcluate the component and mangnitued averages of point
 * values from a list of point IDs. **/ 
void CalculateStats(vtkIdType c_point, vtkSmartPointer<vtkIdList> points, double *vectorAverage, double *magAverage, double *vectorStDev, double *magStDev)
{
	unsigned int nPoints = points->GetNumberOfIds();
	*vectorAverage = 0;
	*(vectorAverage+1) = 0;
	*(vectorAverage+2) = 0;
	*magAverage = 0;
	*vectorStDev = 0;
	*(vectorStDev+1) = 0;
	*(vectorStDev+2) = 0;
	*magStDev = 0;
	for( unsigned int i = 0; i < nPoints; ++i){
		vtkIdType pointId = points->GetId( i );
		if ( pointId == c_point) continue; // don't include the current point in the average.
		double *cPoint = new double[3];
		this->m_DataImage->GetPointData()->GetArray("Displacement")->GetTuple( pointId, cPoint );
		
		*vectorAverage = *vectorAverage + *cPoint;
		*(vectorAverage+1) = *(vectorAverage+1) + *(cPoint+1);
		*(vectorAverage+2) = *(vectorAverage+2) + *(cPoint+2);
		
		*magAverage = *magAverage + std::sqrt( std::pow(*cPoint,2) + std::pow(*(cPoint+1),2) + std::pow(*(cPoint+2),2) );
	}
	
	*vectorAverage = *vectorAverage/nPoints;
	*(vectorAverage+1) = *(vectorAverage+1)/nPoints;
	*(vectorAverage+2) = *(vectorAverage+2)/nPoints;
	
	*magAverage = *magAverage/nPoints-1;
	
	double vectDiff[3] = {0, 0, 0};
	double magDiff = 0;
	for( unsigned int i = 0; i <nPoints; ++i){
		vtkIdType pointId = points->GetId( i );
		double *cPoint = new double[3];
		this->m_DataImage->GetPointData()->GetArray("Displacement")->GetTuple( pointId, cPoint );
		double cMag = std::sqrt( std::pow(*cPoint,2) + std::pow(*(cPoint+1),2) + std::pow(*(cPoint+2),2) );
		
		vectDiff[0] = vectDiff[0] + std::pow( (*cPoint-*vectorAverage), 2);
		vectDiff[1] = vectDiff[1] + std::pow( (*(cPoint+1)-*(vectorAverage+1)) ,2);
		vectDiff[2] = vectDiff[2] + std::pow( (*(cPoint+2)-*(vectorAverage+2)), 2);
		
		magDiff = magDiff + std::pow( (cMag - *magAverage), 2);
	}
	
	*vectorStDev = std::sqrt(vectDiff[0]/nPoints);
	*(vectorStDev+1) = std::sqrt(vectDiff[1]/nPoints);
	*(vectorStDev+2) = std::sqrt(vectDiff[2]/nPoints);
	
	*magStDev = std::sqrt(magDiff/nPoints);
}

/** A function to calculate the strains. */
void GetStrains()
{
	// calculate derivatives
	vtkCellDerivatives *derivativeCalculator = vtkCellDerivatives::New();
	derivativeCalculator->SetInput( this->m_DataImage );
	derivativeCalculator->SetTensorModeToComputeStrain();
	derivativeCalculator->SetVectorModeToComputeGradient();
	
	derivativeCalculator->SetInputArrayToProcess(1,0,0,0,"Displacement");
	derivativeCalculator->Update();
	this->SetDataImage( derivativeCalculator->GetUnstructuredGridOutput() );
		
	// move the data from cells to nodes
	vtkCellDataToPointData	*dataTransformer = vtkCellDataToPointData::New();
	dataTransformer->SetInput( this->m_DataImage );
	dataTransformer->PassCellDataOn();
	dataTransformer->SetInputArrayToProcess(0,0,1,1,"Strain");
	dataTransformer->Update();
	
	this->SetDataImage( dataTransformer->GetUnstructuredGridOutput() );
}

/** A function to calculate the principal strains. */
void GetPrincipalStrains()
{
	vtkSmartPointer<vtkMath>		mathAlgorithm = vtkSmartPointer<vtkMath>::New();
	// First find the principal strains for the point data
	vtkIdType nPoints = this->m_DataImage->GetPointData()->GetArray("Strain")->GetNumberOfTuples();
	
	DataImagePixelPointer V0;
	DataImagePixelPointer V1;
	DataImagePixelPointer V2;
	DataImagePixelPointer val0;
	DataImagePixelPointer val1;
	DataImagePixelPointer val2;
	
	// create the eigenvector containers if they don't exist
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Vector 1") ){
		V0 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray( "Principal Strain Vector 1" ) );

	}
	else{
		V0 = DataImagePixelPointer::New();
		V0->SetNumberOfComponents(3);
		V0->SetName("Principal Strain Vector 1");
		V0->Allocate(nPoints);
		V0->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( V0 );
	}
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Vector 2") ){
		V1 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray( "Principal Strain Vector 2" ) );
	}
	else{
		V1 = DataImagePixelPointer::New();
		V1->SetNumberOfComponents(3);
		V1->SetName("Principal Strain Vector 2");
		V1->Allocate(nPoints);
		V1->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( V1 );
	}
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Vector 3") ){
		V2 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray( "Principal Strain Vector 3" ) );
	}
	else{
		V2 = DataImagePixelPointer::New();
		V2->SetNumberOfComponents(3);
		V2->SetName("Principal Strain Vector 3");
		V2->Allocate(nPoints);
		V2->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( V2 );
	}

	// create the eigenvalue containers if they don't exist
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 1") ){
		val0 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray( "Principal Strain Value 1" ) );
	}
	else{
		val0 = DataImagePixelPointer::New();
		val0->SetNumberOfComponents(1);
		val0->SetName("Principal Strain Value 1");
		val0->Allocate(nPoints);
		val0->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( val0 );
	}
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 2") ){
		val1 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 2") );
	}
	else{
		val1 = DataImagePixelPointer::New();
		val1->SetNumberOfComponents(1);
		val1->SetName("Principal Strain Value 2");
		val1->Allocate(nPoints);
		val1->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( val1 );
	}
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 3") ){
		val2 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 3") );
	}
	else{
		val2 = DataImagePixelPointer::New();
		val2->SetNumberOfComponents(1);
		val2->SetName("Principal Strain Value 3");
		val2->Allocate(nPoints);
		val2->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( val2 );
	}	
	
	for ( unsigned int i = 0; i < nPoints; ++i ){
		// Get the point tensor
		double *cTensorRaw = this->m_DataImage->GetPointData()->GetArray("Strain")->GetTuple( i );
		// reform tensor into an appropreate array
		double *cTensor[3];
		double cT0[3];
		double cT1[3];
		double cT2[3];

		cTensor[0] = cT0; cTensor[1] = cT1; cTensor[2] = cT2;

		cT0[0] = *cTensorRaw;
		cT0[1] = *(cTensorRaw+1);
		cT0[2] = *(cTensorRaw+2);	
		cT1[0] = *(cTensorRaw+3);
		cT1[1] = *(cTensorRaw+4);
		cT1[2] = *(cTensorRaw+5);
		cT2[0] = *(cTensorRaw+6);
		cT2[1] = *(cTensorRaw+7);
		cT2[2] = *(cTensorRaw+8);
				
		// Get the eigen-values and -vectors
		double eigenValues[3]; // The Eigenvalues
		double *eigenVectors[3]; // vector of pointers to eigenvectors
		double v0[3]; // eigen vector 1
		double v1[3]; // eigen vector 2
		double v2[3]; // eigen vector 3
		eigenVectors[0] = v0; eigenVectors[1] = v1; eigenVectors[2] = v2;
		
		// perform the calculation of the principal strains
		mathAlgorithm->Jacobi(cTensor,eigenValues,eigenVectors);
		
		// save the values and vectors into the data image
		V0->SetTuple( i, eigenVectors[0] );
		val0->SetTuple( i, &eigenValues[0] );
		V1->SetTuple( i, eigenVectors[1] );
		val1->SetTuple( i, &eigenValues[1] );
		V2->SetTuple( i, eigenVectors[2] );
		val2->SetTuple( i, &eigenValues[2] );
	}
	
	// next find the principal strains for the cell data
	nPoints = this->m_DataImage->GetCellData()->GetArray("Strain")->GetNumberOfTuples();
	
	// create the eigenvector containers if they don't exist
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 1") ){
		V0 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 1") );
	}
	else{
		V0 = DataImagePixelPointer::New();
		V0->SetNumberOfComponents(3);
		V0->SetName("Principal Strain Vector 1");
		V0->Allocate(nPoints);
		V0->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( V0 );
	}
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 2") ){
		V1 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 2") );
	}
	else{
		V1 = DataImagePixelPointer::New();
		V1->SetNumberOfComponents(3);
		V1->SetName("Principal Strain Vector 2");
		V1->Allocate(nPoints);
		V1->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( V1 );
	}
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 3") ){
		V2 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 3") );
	}
	else{
		V2 = DataImagePixelPointer::New();
		V2->SetNumberOfComponents(3);
		V2->SetName("Principal Strain Vector 3");
		V2->Allocate(nPoints);
		V2->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( V2 );
	}

	// create the eigenvalue containers if they don't exist
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 1") ){
		val0 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 1") );
	}
	else{
		val0 = DataImagePixelPointer::New();
		val0->SetNumberOfComponents(1);
		val0->SetName("Principal Strain Value 1");
		val0->Allocate(nPoints);
		val0->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( val0 );
	}
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 2") ){
		val1 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 2") );
	}
	else{
		val1 = DataImagePixelPointer::New();
		val1->SetNumberOfComponents(1);
		val1->SetName("Principal Strain Value 2");
		val1->Allocate(nPoints);
		val1->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( val1 );
	}
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 3") ){
		val2 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 3") );
	}
	else{
		val2 = DataImagePixelPointer::New();
		val2->SetNumberOfComponents(1);
		val2->SetName("Principal Strain Value 3");
		val2->Allocate(nPoints);
		val2->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( val2 );
	}
	
	for ( unsigned int i = 0; i < nPoints; ++i ){
		// Get the point tensor
		double *cTensorRaw = this->m_DataImage->GetCellData()->GetArray("Strain")->GetTuple( i );

		// reform tensor into an appropreate array
		double *cTensor[3];
		double cT0[3];
		double cT1[3];
		double cT2[3];

		cTensor[0] = cT0; cTensor[1] = cT1; cTensor[2] = cT2;

		cT0[0] = *cTensorRaw;
		cT0[1] = *(cTensorRaw+1);
		cT0[2] = *(cTensorRaw+2);	
		cT1[0] = *(cTensorRaw+3);
		cT1[1] = *(cTensorRaw+4);
		cT1[2] = *(cTensorRaw+5);
		cT2[0] = *(cTensorRaw+6);
		cT2[1] = *(cTensorRaw+7);
		cT2[2] = *(cTensorRaw+8);

		// Get the eigen-values and -vectors
		double eigenValues[3]; // The Eigenvalues
		double *eigenVectors[3]; // vector of pointers to eigenvectors
		double v0[3]; // eigen vector 1
		double v1[3]; // eigen vector 2
		double v2[3]; // eigen vector 3
		eigenVectors[0] = v0; eigenVectors[1] = v1; eigenVectors[2] = v2;

		// perform the calculation of the principal strains
		mathAlgorithm->Jacobi(cTensor,eigenValues,eigenVectors);

		// save the values and vectors into the data image
		V0->SetTuple( i, eigenVectors[0] );
		val0->SetTuple( i, &eigenValues[0] );
		V1->SetTuple( i, eigenVectors[1] );
		val1->SetTuple( i, &eigenValues[1] );
		V2->SetTuple( i, eigenVectors[2] );
		val2->SetTuple( i, &eigenValues[2] );
	}
	
}

/** A function to smooth the mesh image.  This method takes the N closest 
 * points (includeing the point in question) and takes the point with the
 * median magnitude as the value. N must be odd or 1 will be added to it. */
void MedianImageFilter(unsigned int N)
{
	// if N is even, add one
	N % 2 == 1 ? : N += 1;
	// create a duplicate of the data image
	DataImagePointer tempImage = DataImagePointer::New();
	tempImage->DeepCopy(this->m_DataImage);
	// loop through the points of both the images
	unsigned int nPoints = tempImage->GetNumberOfPoints();
	for( unsigned int i = 0; i < nPoints; ++i ){
		double *cPoint = tempImage->GetPoint( i );
		// find the 7 closest points to the position of each point
		vtkIdList *pointList = vtkIdList::New();
		this->m_KDTree->FindClosestNPoints( N, cPoint, pointList );
		// calculate magnitudes
		double neighbourhoodMagnitudes[N];
		std::vector<double> sortingArray; // array for magnitude sorting
		for( unsigned int j = 0; j < N; ++j){
			vtkIdType cId = pointList->GetId( j );
			double *tempPixel = tempImage->GetPointData()->GetArray("Displacement")->GetTuple( cId );
			double cMag = sqrt( pow(*tempPixel,2) + pow(*(tempPixel+1),2) + pow(*(tempPixel+2),2) );
			neighbourhoodMagnitudes[j] = cMag ;
			sortingArray.push_back( cMag );
		}
		// sort the magnitudes
		std::sort( sortingArray.begin(),sortingArray.end() );
		double medianMag = sortingArray[ (N >> 1) ];
		// find the magnitude in the neighbourhoodMagnitudes
		unsigned int j = 0;
		while( neighbourhoodMagnitudes[j] != medianMag ) ++j;
		// set the data image pixel to the value of the temp image
		vtkIdType medianPoint = pointList->GetId( j );
		double *medianPixel = tempImage->GetPointData()->GetArray("Displacement")->GetTuple( medianPoint );
		this->m_DataImage->GetPointData()->GetArray("Displacement")->SetTuple( i, medianPixel );
	}
}

/** This function will use the image registration method from DIC to 
 * align the images in the region bounded by m_DataImage's bounding box.*/
void GlobalRegistration()
{
	// Use an ShrinkImageFilter to blur and downsample fixed and moving images to improve radius of convergance
	typedef itk::ShrinkImageFilter< FixedImageType, FixedImageType > FixedResamplerType;
	typedef itk::ShrinkImageFilter< MovingImageType, MovingImageType > MovingResamplerType;
	
	typename FixedResamplerType::Pointer	fixedResampler = FixedResamplerType::New();
	typename MovingResamplerType::Pointer	movingResampler = MovingResamplerType::New();		
	fixedResampler->SetInput( this->m_FixedImage );
	movingResampler->SetInput( this->m_MovingImage );
	fixedResampler->SetShrinkFactors( this->m_GlobalRegDownsampleValue );
	movingResampler->SetShrinkFactors( this->m_GlobalRegDownsampleValue );
	fixedResampler->SetNumberOfThreads( this->m_Registration->GetNumberOfThreads() );
	movingResampler->SetNumberOfThreads( this->m_Registration->GetNumberOfThreads() );
	std::stringstream msg("");
	msg <<"Resampling for global registration"<<std::endl<<std::endl;
	this->WriteToLogfile(msg.str());
	fixedResampler->Update();
	movingResampler->Update();
	
	// global registration - rotation is centred on the body
	this->m_Registration->SetFixedImage( fixedResampler->GetOutput() );
	this->m_Registration->SetMovingImage( movingResampler->GetOutput() );
	this->m_TransformInitializer->SetFixedImage( fixedResampler->GetOutput() );
	this->m_TransformInitializer->SetMovingImage( movingResampler->GetOutput() );
	this->m_TransformInitializer->SetTransform( this->m_Transform );
	this->m_TransformInitializer->GeometryOn();
	this->m_TransformInitializer->InitializeTransform();
	this->m_Registration->SetInitialTransformParameters( this->m_Transform->GetParameters() );
	
	// restrict the Global registration to just the bounding box of the mesh
	double meshBBox[6];// = new double[6];
	this->m_DataImage->GetBounds( meshBBox );
	typename FixedImageType::PointType meshMinPt; // Get mesh minimum pt - this will be used for the registration region start index 
	meshMinPt[0] = meshBBox[0];
	meshMinPt[1] = meshBBox[2];
	meshMinPt[2] = meshBBox[4];
	double meshSize[3]; 				// Get mesh dimension - this will be used for the registration size
	meshSize[0] = meshBBox[1]-meshBBox[0];
	meshSize[1] = meshBBox[3]-meshBBox[2];
	meshSize[2] = meshBBox[5]-meshBBox[4];
		
	typename FixedImageType::IndexType fixedImageROIStart;
	fixedResampler->GetOutput()->TransformPhysicalPointToIndex(meshMinPt,fixedImageROIStart); // convert min point to start index
	
	typename FixedImageType::SpacingType fixedSpacing = fixedResampler->GetOutput()->GetSpacing(); // convert dimensinos to size in pixels
	typename FixedImageType::SizeType fixedImageROILengths;
	fixedImageROILengths[0] = (int)std::floor(meshSize[0]/fixedSpacing[0]);
	fixedImageROILengths[1] = (int)std::floor(meshSize[1]/fixedSpacing[1]);
	fixedImageROILengths[2] = (int)std::floor(meshSize[2]/fixedSpacing[2]);
	typename FixedImageType::RegionType fixedAnalysisRegion;
	fixedAnalysisRegion.SetIndex( fixedImageROIStart );
	fixedAnalysisRegion.SetSize( fixedImageROILengths );
	this->m_Registration->SetFixedImageRegion( fixedAnalysisRegion ); // set the limited analysis region
	
	this->m_Registration->SetFixedImageRegionDefined( true );
	msg.str("");
	msg << "Global registration in progress"<<std::endl<<std::endl;
	this->WriteToLogfile( msg.str() );
	this->m_Registration->Update();
	this->m_Registration->SetFixedImageRegionDefined( false );
	
	msg.str("");
	msg << "Global Registration complete."<<std::endl;
	this->WriteToLogfile( msg.str() );
	
	
	double globalRegResults[3];// = new double[3];
	RegistrationParametersType	finalParameters = this->m_Registration->GetLastTransformParameters();
	msg.str("");
	msg << "Final Params:"<< finalParameters<<std::endl;
	this->WriteToLogfile( msg.str() );
	globalRegResults[0] = finalParameters[6];
	globalRegResults[1] = finalParameters[7];
	globalRegResults[2] = finalParameters[8];	
	
	msg.str("");
	msg << "Global registration finished.\n Resulting displacement: ("<<globalRegResults[0]<<
		", "<<globalRegResults[1]<<", "<<globalRegResults[2]<<")"<<std::endl<<std::endl;
	this->WriteToLogfile( msg.str() );
	
	this->SetMeshInitialValue( globalRegResults );
}

/** A function to set the downsample value used by the global registration.
 * The default value is 3. */
void SetGlobalRegistrationDownsampleValue( unsigned int value )
{
	if (m_GlobalRegDownsampleValue != value){
		this->m_GlobalRegDownsampleValue = value;
	}
}

/** A function to set the error tolerance in units of standard deviations 
 * from the average of the pixels in the error radius. */
void SetErrorTolerance(double tolerance)
{
	this->m_errorTolerance = tolerance;
}

/** A fucntion to set the error in image dimensions (eg. mm, inches). */
void SetErrorRadius(double radius)
{
	this->m_errorRadius = radius;
}

void SetMaxMetricValue( double maxVal )
{
	if (m_maxMeticValue != maxVal){
		m_maxMeticValue = maxVal;
	}
}

double GetMaxMetricValue()
{
	return m_maxMeticValue;
}

private:

DataImagePointer			m_DataImage;
vtkSmartPointer<vtkPKdTree>	m_KDTree;
double						m_errorRadius;
double						m_errorTolerance;
double						m_maxMeticValue; // this because I'm using the normalized x-correlation coefficient metric
vtkSmartPointer<vtkIdList>	m_pointsList;
RegistrationParametersType	m_GlobalRegistrationParameters;
unsigned int				m_GlobalRegDownsampleValue;
	
}; // end class DICMesh


#endif // DICMESH_H

